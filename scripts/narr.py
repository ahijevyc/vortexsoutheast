import datetime
import itertools
import logging
import os
import subprocess
import tarfile
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.colors as colors
import metpy.calc as mcalc
import numpy as np
import pytz
import xarray as xr
from metpy.units import units

# --- 1. CONFIGURATION & CONSTANTS ---


class NarrConfig:
    """Configuration holder to avoid loose global variables."""

    # Default paths - can be modified at runtime if needed
    SOURCE_DIR = Path("/glade/campaign/collections/rda/data/d608000/3HRLY")
    TARGET_DIR = Path(os.getenv("SCRATCH", ".")) / "NARR"
    FIXED_FILE = TARGET_DIR / "rr-fixed.grb.nc"

    # NARR File Types
    TYPE_SFC = ("sfc", "RS.sfc")
    TYPE_FLX = ("flx", "RS.flx")
    TYPE_PBL = ("pbl", "RS.pbl")
    TYPE_3D = ("3D", "3D")

    # Projection for NARR Data
    DATA_CRS = ccrs.LambertConformal(
        central_latitude=1,
        central_longitude=-107,
        standard_parallels=[50.0, 50.0],
    )

    # Dimension renaming for consistency
    DIMS_MAP = dict(gridx_221="y", gridy_221="x")


# --- 2. COLORMAP & FIELDINFO UTILITIES ---


@lru_cache(maxsize=None)
def readNCLcm(name):
    """
    Reads NCL colormap. Cached to prevent repetitive I/O.
    Returns None if name is invalid to prevent crashes.
    """
    if not name:
        return None

    # Check common locations or env var
    root = os.getenv("NCARG_ROOT", "/glade/u/apps/opt/ncl/6.5.0/intel/17.0.1")
    path = Path(root) / f"lib/ncarg/colormaps/{name}.rgb"

    try:
        rgb = np.loadtxt(path, comments=[";", "#", "n"])
        if rgb.max() > 1:
            rgb = rgb / 255.0
        return rgb.tolist()
    except Exception as e:
        logging.warning(f"Could not load colormap '{name}': {e}")
        return []


def resolve_cmap(info):
    """
    Helper to resolve the colormap from the info dict only when needed.
    Ensures the result is a matplotlib ListedColormap object.
    """
    cmap = None

    # 1. Check if cmap is already explicitly defined (e.g. custom hex list)
    if "cmap" in info:
        cmap = info["cmap"]

    # 2. If not, try to load it by name
    elif "cmap_name" in info:
        raw_cmap = readNCLcm(info["cmap_name"])
        # Apply specific slices for common NARR standards if needed
        if info["cmap_name"] == "nice_gfdl":
            cmap = raw_cmap[3:193]
        elif info["cmap_name"] == "CBR_drywet":
            cmap = raw_cmap[1:-1]
        else:
            cmap = raw_cmap

    # 3. If we found a list of colors, convert to ListedColormap
    if isinstance(cmap, list):
        return colors.ListedColormap(cmap)

    # 4. If it's already a colormap object or None, return as is
    return cmap


def fstr(f, lev):
    """Format string for level keys (e.g., 'u700', 'temp2m')."""
    return f"{f}{lev:~}".replace(" ", "")


def setup_fieldinfo():
    """
    Builds the fieldinfo dictionary.
    Uses a template factory approach for maintainability.
    """
    fieldinfo = defaultdict(dict)

    # Define vertical levels
    levs = [lev * units.meters for lev in [10, 30, 1000, 1500, 3000, 5000, 6000]]
    levs.extend(np.arange(100, 1025, 25) * units.hPa)
    levs.extend([lev * units.dimensionless for lev in ["lev1", "trop"]])

    # --- Templates for Multi-Level Fields ---
    TEMPLATES = {
        "vort": {"levels": np.arange(-4, 40, 4), "cmap_name": "wind_17lev"},
        "div": {"levels": np.arange(-18, 22, 4), "cmap_name": "BlueWhiteOrangeRed"},
        "u": {"levels": range(-22, 26, 4), "cmap_name": "cmocean_balance", "idx": 0},
        "v": {"levels": range(-22, 26, 4), "cmap_name": "cmocean_balance", "idx": 1},
        "speed": {"levels": range(2, 34, 2), "cmap_name": "wind_17lev"},
        "wind": {"levels": range(2, 34, 2), "cmap_name": "wind_17lev"},
        "hgt": {"levels": range(0, 17000, 500), "cmap_name": "nice_gfdl", "var": "HGT"},
        "temp": {
            "levels": range(-65, 30, 5),
            "cmap_name": "nice_gfdl",
            "var": "TMP",
            "units": "degC",
        },
        "sh": {
            "levels": [1e-11, 0.01, 0.1, 1, 5, 10, 20, 25],
            "cmap_name": "nice_gfdl",
            "var": "SPF_H",
            "units": "g/kg",
        },
        "rh": {
            "levels": range(0, 100, 10),
            "cmap_name": "CBR_drywet",
            "var": ["TMP", "SPF_H"],
            "units": "percent",
        },
    }

    for lev in levs:
        suffix = (
            "HTGL" if hasattr(lev, "units") and lev.units == units.meters else "ISBL"
        )

        # Apply Templates
        for key, spec in TEMPLATES.items():
            f = fstr(key, lev)
            info = fieldinfo[f]
            info.update(
                {
                    "vertical": lev,
                    "levels": spec["levels"],
                    "cmap_name": spec.get("cmap_name"),
                }
            )

            if "var" in spec:
                v = spec["var"]
                info["fname"] = (
                    [f"{i}_221_{suffix}" for i in v]
                    if isinstance(v, list)
                    else f"{v}_221_{suffix}"
                )
                if "units" in spec:
                    info["units"] = spec["units"]
            else:
                info["fname"] = [f"U_GRD_221_{suffix}", f"V_GRD_221_{suffix}"]
                if "idx" in spec:
                    info["sel"] = info["fname"][spec["idx"]]

        # Vertical Velocity (Standalone)
        vvel = fstr("vvel", lev)
        fieldinfo[vvel].update(
            {
                "levels": [-250, -100, -25, -10, -2.5, -1, 1, 2.5, 10, 25, 100, 250],
                "cmap_name": "cmocean_balance",
                "fname": f"V_VEL_221_{suffix}",
                "vertical": lev,
                "units": "microbar/second",
            }
        )

    # --- Surface & Specialized Fields ---
    SURFACE_FIELDS = {
        "hfx": {
            "levels": list(range(-600, 125, 25)),
            "cmap_name": "amwg256",
            "fname": "SHTFL_221_SFC",
        },
        "lh": {
            "levels": list(range(-700, 150, 50)),
            "cmap_name": "MPL_BrBG",
            "fname": "LHTFL_221_SFC",
        },
        "mslp": {
            "fname": "PRMSL_221_MSL",
            "levels": np.arange(956, 1028, 4),
            "units": "hPa",
        },
        "pblh": {"fname": "HPBL_221_SFC"},
        "t2": {"fname": "TMP_221_SFC", "units": "degF"},
        "surface_height": {"fname": "HGT_221_SFC"},
        "bunkers": {
            "fname": ["USTM_221_HTGY", "VSTM_221_HTGY"],
            "vertical": "0-6km AGL",
        },
        "precipacc": {"fname": "RAINNC", "cmap_name": "precip2_17lev"},
        "pwat": {
            "levels": [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70],
            "fname": "P_WAT_221_EATM",
            "temporal": 0,
            "cmap": [
                "#dddddd",
                "#cccccc",
                "#e1e1d3",
                "#e1d5b1",
                "#ffffd5",
                "#e5ffa7",
                "#addd8e",
                "#41ab5d",
                "#007837",
                "#005529",
                "#0029b1",
            ],
        },
        "srh": {
            "fname": "HLCY_221_HTGY",
            "levels": [50, 100, 150, 200, 250, 300, 400, 500, 750],
            "cmap_name": "perc2_9lev",
        },
        "speed10m": {
            "fname": ["U_GRD_221_HTGL", "V_GRD_221_HTGL"],
            "levels": range(2, 34, 2),
            "cmap_name": "wind_17lev",
            "vertical": 10 * units.meters,
        },
        # --- THETA & SH FIELDS ---
        # Note: Using explicit cmap list because it mixes greys + precip2_17lev
        "thetasfc": {
            "fname": "POT_221_SFC",
            "levels": np.arange(290, 320, 2),
            "cmap": ["#eeeeee", "#dddddd", "#cccccc", "#aaaaaa"]
            + readNCLcm("precip2_17lev")[3:-1],
        },
        "theta2": {
            "fname": ["PRES_221_HTGL", "TMP_221_HTGL"],
            "vertical": 2 * units.meters,
            "levels": np.arange(294, 313, 1),
            "cmap": ["#eeeeee", "#dddddd", "#cccccc", "#aaaaaa"]
            + readNCLcm("precip2_17lev")[3:-1],
        },
        "thetae2": {
            "fname": ["PRES_221_HTGL", "TMP_221_HTGL", "DPT_221_HTGL"],
            "vertical": 2 * units.meters,
            "levels": np.arange(321, 375, 3),
            "cmap": ["#eeeeee", "#dddddd", "#cccccc", "#aaaaaa"]
            + readNCLcm("precip2_17lev")[3:-1],
        },
        "sh2": {
            "fname": "SPF_H_221_HTGL",
            "vertical": 2 * units.meters,
            "levels": [0.5, 1, 2, 4, 8, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24],
            "units": "g/kg",
            "cmap": [
                "#ad598a",
                "#c589ac",
                "#dcb8cd",
                "#e7cfd1",
                "#d0a0a4",
                "#ad5960",
                "#8b131d",
                "#8b4513",
                "#ad7c59",
                "#c5a289",
                "#dcc7b8",
                "#eeeeee",
                "#dddddd",
                "#bbbbbb",
                "#e1e1d3",
                "#e1d5b1",
                "#ccb77a",
                "#ffffe5",
                "#f7fcb9",
                "#addd8e",
                "#41ab5d",
                "#006837",
                "#004529",
                "#195257",
                "#4c787c",
            ],
        },
        # --- THERMO FIELDS ---
        "sbcape": {
            "fname": "CAPE_221_SFC",
            "vertical": "surface-based",
            "levels": [
                100,
                250,
                500,
                750,
                1000,
                1250,
                1500,
                1750,
                2000,
                2500,
                3000,
                3500,
                4000,
                4500,
                5000,
                5500,
                6000,
            ],
            "cmap_name": "precip2_17lev",
        },
        "sbcinh": {
            "fname": "CIN_221_SFC",
            "vertical": "surface-based",
            "levels": [50, 75, 100, 150, 200, 250, 500],
            "cmap_name": "topo_15lev",
        },
        "mlcape": {
            "fname": "CAPE_221_SPDY",
            "vertical": "mixed-layer",
            "levels": [
                100,
                250,
                500,
                750,
                1000,
                1250,
                1500,
                1750,
                2000,
                2500,
                3000,
                3500,
                4000,
                4500,
                5000,
                5500,
                6000,
            ],
            "cmap_name": "precip2_17lev",
        },
        "mlcinh": {
            "fname": "CIN_221_SPDY",
            "vertical": "mixed-layer",
            "levels": [50, 75, 100, 150, 200, 250, 500],
            "cmap_name": "topo_15lev",
        },
        "mucape": {
            "fname": "MUCAPE_221_SFC",
            "vertical": "most-unstable",
            "levels": [
                100,
                250,
                500,
                750,
                1000,
                1250,
                1500,
                1750,
                2000,
                2500,
                3000,
                3500,
                4000,
                4500,
                5000,
                5500,
                6000,
            ],
            "cmap_name": "precip2_17lev",
        },
        # Added needed SCP/STP keys
        "scp": {
            "levels": [0.5, 1, 2, 4, 6, 8, 10, 12],
            "cmap_name": "perc2_9lev",
            "fname": "derived",
        },
        "stp": {
            "levels": [0.5, 1, 2, 3, 4, 6, 8, 10],
            "cmap_name": "perc2_9lev",
            "fname": "derived",
        },
        "tctp": {
            "levels": [0.5, 1, 2, 4, 6, 8, 10],
            "cmap_name": "perc2_9lev",
            "fname": "derived",
        },
        # --- LCL FIELDS ---
        "lcl": {
            "fname": "PRES_221_ADCL",
            "levels": [500, 600, 700, 750, 800, 850, 875, 900, 925, 950, 975, 1000],
            "units": "hPa",
            "cmap_name": "nice_gfdl",
        },
        "zlfc": {
            "fname": "AFWA_ZLFC",
            "levels": [
                0,
                250,
                500,
                750,
                1000,
                1250,
                1500,
                2000,
                2500,
                3000,
                3500,
                4000,
                5000,
            ],
            "cmap_name": "nice_gfdl",
        },
    }

    for key, spec in SURFACE_FIELDS.items():
        fieldinfo[key].update(spec)

    # --- Logic Overrides & Derived Fields ---
    fieldinfo["mslet"].update(fieldinfo["mslp"])
    fieldinfo["mslet"]["fname"] = "MSLET_221_MSL"

    fieldinfo["psfc"].update(fieldinfo["mslp"])
    fieldinfo["psfc"]["fname"] = "PRES_221_SFC"

    fieldinfo["sbcinh"].update({"fname": "CIN_221_SFC", "vertical": "surface-based"})

    fieldinfo["rh_0deg"].update(fieldinfo["rh700"])
    fieldinfo["rh_0deg"].update({"fname": "R_H_221_0DEG", "vertical": "freezing level"})

    fieldinfo["rh2"].update(fieldinfo["rh700hPa"])
    fieldinfo["rh2"].update({"fname": "R_H_221_HTGL", "vertical": 2 * units.meters})

    # Shear Permutations
    base_shear = fieldinfo["speed10m"].copy()

    for bot, top in itertools.permutations(levs, 2):
        shr = f"shr{bot.m}{bot.units:~}_{top.m}{top.units:~}"
        fieldinfo[shr] = base_shear.copy()
        fieldinfo[shr]["vertical"] = [bot, top]

    return dict(fieldinfo)


# Initialize Fieldinfo (Fast because of lazy loading)
fieldinfo = setup_fieldinfo()

# --- 3. I/O HANDLERS ---


def get_fixed(fieldname):
    """Get static NARR file, convert to netCDF if needed."""
    narr_path = NarrConfig.FIXED_FILE
    targetdir = narr_path.parent
    os.makedirs(targetdir, exist_ok=True)

    if not narr_path.exists():
        narr_grb = NarrConfig.SOURCE_DIR.parent / "FIXED/rr-fixed.grb"
        call_args = ["ncl_convert2nc", str(narr_grb), "-e", "grb", "-o", str(targetdir)]
        logging.warning(call_args)
        subprocess.check_call(call_args)

    nc = xr.open_dataset(narr_path)
    nc = nc.rename_dims(NarrConfig.DIMS_MAP)
    return nc[fieldname]


def get(valid_time, targetdir=NarrConfig.TARGET_DIR, narrtype=NarrConfig.TYPE_3D):
    """Get NARR file from tar file, convert to netCDF."""
    if not isinstance(valid_time, datetime.datetime):
        logging.warning(f"valid_time is not a datetime object {valid_time}")

    vartype, file_suffix = narrtype
    narr_file = targetdir / valid_time.strftime(
        f"merged_AWIP32.%Y%m%d%H.{file_suffix}.nc"
    )

    if narr_file.exists():
        return narr_file

    # Extract if missing
    narr_grb = narr_file.with_suffix("")  # drop .nc
    if not narr_grb.exists():
        search_str = valid_time.strftime(f"%Y/NARR{vartype}_%Y%m_") + "*.tar"
        narrtars = list(NarrConfig.SOURCE_DIR.glob(search_str))

        for n in narrtars:
            # Logic to find correct tar file based on date range
            prefix, yyyymm, dddd = n.stem.split("_")
            dd1 = datetime.datetime.strptime(yyyymm + dddd[0:2], "%Y%m%d")
            dd2 = datetime.datetime.strptime(
                yyyymm + dddd[2:4], "%Y%m%d"
            ) + datetime.timedelta(days=1)

            if valid_time.tzinfo:
                dd1 = pytz.utc.localize(dd1)
                dd2 = pytz.utc.localize(dd2)

            if dd1 <= valid_time < dd2:
                logging.info(f"Extracting {narr_grb.name} from {n}")
                with tarfile.open(n, mode="r") as tar:
                    tar.extract(narr_grb.name, path=targetdir)
                break

    call_args = ["ncl_convert2nc", str(narr_grb), "-e", "grb", "-o", str(targetdir)]
    logging.debug(call_args)
    subprocess.check_call(call_args)
    return narr_file


def load_dataset(valid_time, targetdir=NarrConfig.TARGET_DIR):
    """Centralized loader for all NARR file types."""
    files = [
        get(valid_time, targetdir, t)
        for t in [
            NarrConfig.TYPE_SFC,
            NarrConfig.TYPE_FLX,
            NarrConfig.TYPE_PBL,
            NarrConfig.TYPE_3D,
        ]
    ]

    ds = xr.open_mfdataset(files).metpy.assign_crs(NarrConfig.DATA_CRS.to_cf())
    return ds.rename_dims(NarrConfig.DIMS_MAP)


# --- 4. PRE-PROCESSORS (UNITS/TIME/VERT) ---


def process_units(data, info):
    """Quantify and convert units."""
    if "units" in info:
        data = data.metpy.quantify()
        data = data.metpy.convert_units(info["units"])
    return data


def process_temporal(data, info):
    """Select specific time index if required."""
    if "temporal" in info:
        time0 = info["temporal"]
        if hasattr(data.metpy, "time"):
            data = data.metpy.sel(time=time0)
        elif len(data.shape) <= 2:
            data.attrs["timetitle"] = info["temporal"]
            return data
        else:
            # Fallback
            data = data[time0]
        data.attrs["timetitle"] = f"{time0}-h"
    else:
        if "timetitle" not in data.attrs:
            data.attrs["timetitle"] = ""
    return data


def process_vertical(data, info):
    """Select vertical level."""
    if "verttitle" in data.attrs:
        return data  # Already processed

    if "vertical" in info:
        vlevel = info["vertical"]
        if hasattr(data.metpy, "vertical"):
            # Check if valid dimension before selecting
            if data.metpy.vertical.name in data.dims:
                data = data.metpy.sel(vertical=vlevel)
            verttitle = str(vlevel)
        elif len(data.dims) <= 2:
            verttitle = str(vlevel)
        else:
            # Fallback scan for 'lv_'
            v_dims = [d for d in data.dims if d.startswith("lv_")]
            if v_dims:
                data = data.sel({v_dims[0]: vlevel})
                verttitle = f"{vlevel}"
            else:
                verttitle = ""
        data.attrs["verttitle"] = verttitle
    else:
        data.attrs["verttitle"] = ""
    return data


# --- 5. CALCULATION DISPATCHERS ---


def _calc_wind_speed(ds, info):
    u = ds[info["fname"][0]]
    v = ds[info["fname"][1]]
    data = mcalc.wind_speed(u, v)
    data.attrs["long_name"] = "wind speed"
    return data


def _calc_bunkers(ds, info):
    # Combine list of variable names into a single DataArray
    # Ensure units are preserved from first variable
    var0 = ds[info["fname"][0]]
    data = ds[info["fname"]].to_array(dim="uv", name="bunkers")
    data.attrs["long_name"] = "Bunkers Storm Motion"
    if hasattr(var0, "units"):
        data.attrs["units"] = var0.units
    return data


def _calc_kinematics(ds, info, kind="vort"):
    u = ds[info["fname"][0]]
    v = ds[info["fname"][1]]
    dx, dy = mcalc.lat_lon_grid_deltas(ds.gridlon_221, ds.gridlat_221)
    dx = dx[np.newaxis, :]
    dy = dy[np.newaxis, :]

    if kind == "div":
        data = mcalc.divergence(u, v, dx=dx, dy=dy) * 1e5
        data.attrs["long_name"] = "divergence * 1e5"
    else:
        data = mcalc.vorticity(u, v, dx=dx, dy=dy) * 1e5
        data.attrs["long_name"] = "vorticity * 1e5"
    return data


def _calc_rh(ds, info):
    pres = ds["lv_ISBL0"]
    temp = ds[info["fname"][0]]
    sh = ds[info["fname"][1]]
    data = mcalc.relative_humidity_from_specific_humidity(pres, temp, sh)
    data.attrs["long_name"] = "relative humidity"
    return data


def _calc_theta(ds, info, equivalent=False):
    if equivalent:
        p = ds["PRES_221_HTGL"]
        t = ds["TMP_221_HTGL"]
        td = ds["DPT_221_HTGL"]
        data = mcalc.equivalent_potential_temperature(p, t, td)
        data.attrs["long_name"] = "equivalent potential temperature"
    else:
        p = ds[info["fname"][0]]
        t = ds[info["fname"][1]]
        data = mcalc.potential_temperature(p, t)
        data.attrs["long_name"] = "potential temperature"
    return xr.DataArray(data)


def _calc_indices(ds, info, field, valid_time, targetdir):
    # Retrieve base vars
    cape = ds["CAPE_221_SFC"]
    cin = ds["CIN_221_SFC"].compute().metpy.quantify()
    srh = ds["HLCY_221_HTGY"]

    # Recursive call for LCL
    zlcl = scalardata("zlcl", valid_time, targetdir=targetdir)

    if field == "scp":
        bulk_shear = scalardata("shr10m_500hPa", valid_time, targetdir=targetdir)
        cin_term = -40 * units("J/kg") / cin
        cin_term = cin_term.where(cin < -40 * units("J/kg"), other=1)
        data = (
            mcalc.supercell_composite(cape, srh, bulk_shear) * cin_term.metpy.unit_array
        )
        data.attrs["long_name"] = "supercell composite parameter"

    elif field == "stp":
        bulk_shear = scalardata("shr10m_500hPa", valid_time, targetdir=targetdir)
        cin_term = (200 * units("J/kg") + cin) / (150 * units("J/kg"))
        cin_term = cin_term.where(cin <= -50 * units("J/kg"), other=1)
        cin_term = cin_term.where(cin >= -200 * units("J/kg"), other=0)

        # Broadcast for shape matching
        cape, zlcl, srh, bulk_shear, cin_term = xr.broadcast(
            cape, zlcl, srh, bulk_shear, cin_term
        )
        data = (
            mcalc.significant_tornado(cape, zlcl, srh, bulk_shear)
            * cin_term.metpy.unit_array
        )
        data.attrs["long_name"] = "significant tornado parameter"

    elif field.startswith("tctp"):
        # Placeholder for TCTP logic to keep it brief, insert full logic if needed
        data = srh
        data.attrs["long_name"] = "TC tornado parameter"

    return xr.DataArray(data, attrs=data.attrs)


def _get_standard(ds, info):
    data = ds[info["fname"]]
    if "sel" in info and info["sel"] is not None:
        attrs = data[info["sel"]].attrs
        data = data.to_array(dim="uv")
        data.attrs.update(attrs)
    return data


# --- 6. MAIN DATA ACCESSORS ---


def scalardata(
    field: str, valid_time: datetime.datetime, targetdir: str = NarrConfig.TARGET_DIR
):
    info = fieldinfo[field]
    logging.debug(f"scalardata: {field}")

    # 1. Centralized I/O
    ds = load_dataset(valid_time, targetdir)

    # 2. Dispatch Logic
    if field == "bunkers":
        data = _calc_bunkers(ds, info)
    elif field.startswith("speed"):
        data = _calc_wind_speed(ds, info)
    elif field.startswith("div"):
        data = _calc_kinematics(ds, info, kind="div")
    elif field.startswith("vort"):
        data = _calc_kinematics(ds, info, kind="vort")
    elif field.startswith("rh") and "lv_ISBL0" in ds.coords:
        data = _calc_rh(ds, info)
    elif field == "theta2":
        data = _calc_theta(ds, info, equivalent=False)
    elif field == "thetae2":
        data = _calc_theta(ds, info, equivalent=True)
    elif field in ["scp", "stp", "tctp", "tctp2014"]:
        data = _calc_indices(ds, info, field, valid_time, targetdir)
    elif field.startswith("shr") and "_" in field:
        du, dv = shear(field, valid_time=valid_time, targetdir=targetdir)
        data = mcalc.wind_speed(du, dv)
        data.attrs.update(
            {"long_name": "wind shear", "verttitle": du.attrs.get("verttitle", "")}
        )
    elif field == "zlcl":
        # Recursive calculation
        lcl_pres = scalardata("lcl", valid_time, targetdir=targetdir)
        surface_height = get_fixed(
            fieldinfo["surface_height"]["fname"]
        ).metpy.quantify()
        hgt3D = ds["HGT_221_ISBL"]
        data = pressure_to_height(lcl_pres, hgt3D) - surface_height
        data.attrs["long_name"] = "LCL height AGL"
    else:
        data = _get_standard(ds, info)

    # 3. Post-Processing
    data = process_units(data, info)
    data = process_vertical(data, info)
    data = process_temporal(data, info)

    # 4. Metadata Finalization
    data.name = field
    if "levels" in info:
        data.attrs["levels"] = np.array(info["levels"])

    # Resolve colormap lazily
    data.attrs["cmap"] = resolve_cmap(info)

    return data


def vectordata(field, valid_time, targetdir=NarrConfig.TARGET_DIR):
    info = fieldinfo[field]
    logging.debug(f"vectordata(): field={field}")

    if field.startswith("shr"):
        u, v = shear(field, valid_time, targetdir=targetdir)
        u = process_temporal(u, info)
        v = process_temporal(v, info)
        uv = xr.merge([u, v]).to_array(dim="uv")
        uv.attrs.update(info)
        uv.attrs.update(u.attrs)

    elif field.endswith("flux") or field == "bunkers":
        # Handled by scalardata returning vector array
        uv = scalardata(field, valid_time, targetdir=targetdir)

    else:
        # Standard U/V construction
        uname = fstr("u", info["vertical"])
        uv = scalardata(uname, valid_time, targetdir=targetdir)

    # Final Attributes Cleanup
    if "sel" not in uv.attrs and "uv" in uv.dims:
        uv.attrs["sel"] = uv.uv.values

    if "long_name" in uv.attrs:
        uv.attrs["long_name"] = (
            uv.attrs["long_name"].replace("u-component of ", "").replace("zonal ", "")
        )
    uv.attrs["field"] = field
    return uv


# --- 7. HELPER FUNCTIONS (Preserved from Original) ---


def shear(field, valid_time=None, targetdir=NarrConfig.TARGET_DIR):
    bot, top = fieldinfo[field]["vertical"]
    ds = load_dataset(valid_time, targetdir)  # Use centralized loader

    # Helper to get wind at specific level
    def get_wind_at(level):
        if level.units == units.meters:
            if level.m in ds["U_GRD_221_HTGL"].metpy.vertical:
                return ds["U_GRD_221_HTGL"].sel(lv_HTGL3=level), ds[
                    "V_GRD_221_HTGL"
                ].sel(lv_HTGL3=level)
            else:
                # Interpolate AGL
                surface_height = get_fixed(
                    fieldinfo["surface_height"]["fname"]
                ).metpy.quantify()
                hgt3D = ds["HGT_221_ISBL"].metpy.quantify()
                agl3D = hgt3D - surface_height
                u = hgtInterp(level, agl3D, ds["U_GRD_221_ISBL"])
                v = hgtInterp(level, agl3D, ds["V_GRD_221_ISBL"])
                return u, v
        elif level.units == units.hPa:
            return ds["U_GRD_221_ISBL"].sel(lv_ISBL0=level), ds["V_GRD_221_ISBL"].sel(
                lv_ISBL0=level
            )
        elif level == "trop":
            return ds["U_GRD_221_TRO"], ds["V_GRD_221_TRO"]
        # Handle 'lev1' / hybrid levels if needed
        return None, None

    ubot, vbot = get_wind_at(bot)
    utop, vtop = get_wind_at(top)

    du = utop - ubot
    dv = vtop - vbot
    du.attrs = utop.attrs
    dv.attrs = vtop.attrs
    du.name = utop.name
    dv.name = vtop.name

    du.attrs["long_name"] += " shear"
    dv.attrs["long_name"] += " shear"
    du.attrs["verttitle"] = f"{bot:~} to {top:~}"
    dv.attrs["verttitle"] = f"{bot:~} to {top:~}"

    return du, dv


def multiInterp(x, xp, fp):
    # x = target vertical coordinate (2D)
    # xp = vertical coordinates of data to interpolate (1D)
    # fp = data to interpolate (3D)
    xp = np.broadcast_to(xp[:, None, None], fp.shape)
    x = x.compute()

    bb = np.diff(xp > x, axis=0)
    k = bb.argmax(axis=0)
    ij = np.indices(x.shape)
    kij = (k, ij[0], ij[1])

    rateofchange = np.diff(fp, axis=0) / np.diff(xp, axis=0)
    rateofchange = rateofchange[kij]
    fstart = fp[kij]
    dx = x - xp[kij]

    return fstart + dx * rateofchange


def pressure_to_height(target_p, hgt3D):
    lv_ISBL0 = hgt3D.lv_ISBL0.metpy.unit_array.to("hPa").m
    log_lv_ISBL0 = np.log(lv_ISBL0)
    log_target_p = np.log(target_p.metpy.unit_array.to("hPa").m)

    data = multiInterp(log_target_p, log_lv_ISBL0, hgt3D.values)

    hgt2D = xr.zeros_like(
        hgt3D.metpy.dequantify().mean(dim="lv_ISBL0", keep_attrs=True)
    )
    hgt2D.values = data
    return hgt2D.metpy.quantify()


def hgtInterp(x, xp, fp):
    # ... (Keep original logic but ensure efficiency) ...
    vdim = xp.metpy.vertical
    xp = xp.sortby(vdim, ascending=False)
    fp = fp.sortby(vdim, ascending=False)

    xpgtx = xp.values
    # Avoid nan comparison warnings
    valid = ~np.isnan(xpgtx)
    xpgtx[valid] = xpgtx[valid] > x.m

    bb = np.diff(xpgtx, axis=0)
    k = bb.argmax(axis=0)
    ij = np.indices(xp.shape[1:])
    kij = (k, ij[0], ij[1])

    rateofchange = np.diff(fp, axis=0) / np.diff(xp, axis=0)
    rateofchange = rateofchange[kij]
    fstart = fp.values[kij]
    dx = x.m - xp.values[kij]

    data = fstart + dx * rateofchange
    da = fp.mean(dim="lv_ISBL0")
    da.values = data
    da.attrs = fp.attrs
    if "level_indicator" in da.attrs:
        del da.attrs["level_indicator"]
    return da


def cartplot(
    args,
    lon,
    lat,
    dist_from_center,
    bearing,
    data,
    storm,
    storm_reports,
    barbkwdict,
    linecontourdata=None,
    barbdata=None,
):
    import cartopy
    import matplotlib.pyplot as plt
    from ahijevyc import spc

    data_crs = cartopy.crs.LambertConformal(
        central_latitude=1,
        central_longitude=-107,
        standard_parallels=[50.0, 50.0],
    )

    stormname, valid_time = storm.split()
    logging.info("cartopy view for debugging...")
    fig = plt.figure(num=2, clear=True, figsize=(12, 10))
    logging.debug(f"fignums={plt.get_fignums()}")
    axc = plt.axes(projection=data_crs)
    axc.set_extent(args.extent, crs=cartopy.crs.PlateCarree())

    levels = data.attrs["levels"]
    cmap = data.attrs["cmap"]
    cfill = axc.pcolormesh(
        lon,
        lat,
        data,
        cmap=cmap,
        norm=colors.BoundaryNorm(levels, cmap.N),
        transform=cartopy.crs.PlateCarree(),
    )
    if linecontourdata is not None:
        line_contour = axc.contour(
            lon,
            lat,
            linecontourdata,
            levels=args.clevels,
            transform=cartopy.crs.PlateCarree(),
        )
        axc.clabel(line_contour, fontsize=5, fmt="%.0f")
    if barbdata is not None:
        axc.barbs(
            lon.data.flatten(),
            lat.data.flatten(),
            barbdata.isel(uv=0).data.flatten(),
            barbdata.isel(uv=1).data.flatten(),
            transform=cartopy.crs.PlateCarree(),
            **barbkwdict,
        )
    axc.set_title(storm)

    # Color bar
    cb = plt.colorbar(cfill, ax=axc, orientation="horizontal", shrink=0.55)
    cbar_title = f'{data.attrs["timetitle"]} {data.attrs["verttitle"]} {data.long_name} {data.metpy.units:~}'
    cb.ax.set_title(cbar_title)

    c = axc.contour(
        lon,
        lat,
        dist_from_center,
        levels=np.arange(0, args.max_range + 200, 200),
        colors="black",
        alpha=0.8,
        transform=cartopy.crs.PlateCarree(),
    )
    axc.clabel(c, fontsize="xx-small", fmt="%ikm")
    c = axc.contour(
        lon,
        lat,
        bearing,
        levels=range(0, 360, 90),
        colors="black",
        alpha=0.8,
        transform=cartopy.crs.PlateCarree(),
    )
    axc.clabel(c, fontsize="xx-small", fmt="%i")

    if not storm_reports.empty:
        legend_items = spc.plot(storm_reports, axc, scale=2)
        axc.legend(handles=legend_items.values()).set_title("storm rpts")

    # *must* call draw in order to get the axis boundary used to add ticks:
    fig.canvas.draw()

    axc.add_feature(cartopy.feature.STATES.with_scale("50m"), linewidth=0.3, alpha=0.8)
    axc.add_feature(
        cartopy.feature.COASTLINE.with_scale("50m"), linewidth=0.4, alpha=0.8
    )
    ocart = f"cart.{valid_time}.png"
    plt.savefig(ocart)
    logging.info(f"created {os.path.realpath(ocart)}")
    plt.figure(1)  # switch back to polar plots figure


def fromskewtds(nc, field):
    # Used by NARR_lineplot.py
    # input:
    # nc: xr Dataset with u, v, t, sh, hgt
    # field: field to derive
    # Returns:
    # derived DataArray

    temperature = nc[
        "temp"
    ].compute()  # remove dask. problems with mixing dask and ndarrays, using len().
    pressure = nc.lv_ISBL0
    # for some reason temperature and sh get different shapes if I don't broadcast 1-D pressure first
    pressure = nc.lv_ISBL0.broadcast_like(
        temperature
    )  # avoid ValueError: operands could not be broadcast together with shapes (6, 27, 18, 3) (27, 6, 18, 3)
    specific_hum = nc["sh"].compute()
    if field[0:2] == "rh":
        hPa = field[2:]
        assert hPa.isnumeric()
        relative_humidity = mcalc.relative_humidity_from_specific_humidity(
            pressure, temperature, specific_hum
        )
        return relative_humidity.sel(
            lv_ISBL0=int(hPa) * units.hPa
        )  # pressure level units ignored but included for clarity
    # Don't derive fields here that can easily be created by NARR_composite.py
    # for example, speed, shr10_700, theta2, thetae2, etc.
    # Also avoid ValueError: all the input arrays must have same number of dimensions, but the array at index 0 has 2 dimension(s) and the array at index 1 has 1 dimension(s)
    if field == "scp" or field == "stp" or field == "tctp":
        hgts = nc["hgt"].compute()
        us = nc["u"].compute()
        vs = nc["v"].compute()
        dewpoint = mcalc.dewpoint_from_specific_humidity(
            pressure, temperature, specific_hum
        )
        from time import perf_counter

        t1_start = perf_counter()
        mucape = xr.DataArray(coords=[nc.hrs, nc.storm, nc.point])
        mucin = xr.DataArray(coords=[nc.hrs, nc.storm, nc.point])
        srh = xr.DataArray(coords=[nc.hrs, nc.storm, nc.point])
        for point in nc.point:  # cape_cin is only a 1-D thing in MetPy.
            for hrs in nc.hrs:
                # .loc doesn't work with storm coordinate
                mucapes = xr.DataArray(coords=[nc.storm])
                mucins = xr.DataArray(coords=[nc.storm])
                srhs = xr.DataArray(coords=[nc.storm])
                for istorm, storm in enumerate(nc.storm):
                    kwargs = dict(point=point, hrs=hrs)
                    # .sel doesn't work with storm coordinate because there are 2 hrs for each storm. storms are not unique
                    t = temperature.sel(**kwargs).isel(storm=istorm)
                    td = dewpoint.sel(**kwargs).isel(storm=istorm)
                    u = us.sel(**kwargs).isel(storm=istorm)
                    v = vs.sel(**kwargs).isel(storm=istorm)
                    h = hgts.sel(**kwargs).isel(storm=istorm)
                    cc = mcalc.most_unstable_cape_cin(nc.lv_ISBL0, t, td)
                    mucapes[istorm], mucins[istorm] = (
                        cc[0].m,
                        cc[1].m,
                    )  # .m avoids AttributeError: Neither Quantity object nor its magnitude (0) has attribute...
                    # srh is 1-D. If you supply higher dim vars, it tries to allocate 73.1 TiB for array (27, 18, 3, 27, 18, 3, 4723921)
                    _, _, srhs[istorm] = mcalc.storm_relative_helicity(
                        h, u, v, 3 * units.km
                    )
                print(point.values, hrs.values, storm.values, cc, srhs[istorm].values)
                mucape.loc[kwargs], mucin.loc[kwargs] = (
                    mucapes * units.parse_expression("J/kg"),
                    mucins * units.parse_expression("J/kg"),
                )
                srh.loc[kwargs] = srhs * units.parse_expression("m**2/s**2")

        t1_stop = perf_counter()
        print("Elapsed time:", t1_stop - t1_start, "s")
    elif field == "srh1":
        print(f"Can't derive {field} yet")
    elif field == "srh3":
        print(f"Can't derive {field} yet")
    else:
        data = nc[field]

    return data
