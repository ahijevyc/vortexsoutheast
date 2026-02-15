import datetime
import itertools
import logging
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import cartopy
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import metpy.calc as mcalc
import numpy as np
import pytz
import xarray
from ahijevyc import spc
from fieldinfo import fieldinfo, readNCLcm
from metpy.units import units

# (vartype, file_suffix)
narrSfc = ("sfc", "RS.sfc")
narrFlx = ("flx", "RS.flx")
narrPBL = ("pbl", "RS.pbl")
narr3D = ("3D", "3D")
targetdir = Path(os.getenv("SCRATCH")) / "NARR"
narrFixed = targetdir / "rr-fixed.grb.nc"

data_crs = cartopy.crs.LambertConformal(
    central_latitude=1,
    central_longitude=-107,
    standard_parallels=[50.0, 50.0],
)

# Modify fieldinfo dictionary for NARR.
# fieldinfo keys are "human-readable" string nicknames like "u700" and "speed10m".
# Their values are also dictionaries with information about how to find the data, read the data, extract vertical level(s), and plot the data.
#   For example the following key : value pairs
#   key      : value
#   cmap     : color map
#   fname    : field variable name
#   levels   : contour levels
#   vertical : vertical level(s) to extract. These are either pint quantities with units. For layers: [bot, top], or string descriptions (e.g. "surface", "freezing level").


def fstr(f, lev):
    return f"{f}{lev:~}".replace(" ", "")


# wind, hgt, rh, sh, temp at all possible vertical levels
levs = [lev * units.meters for lev in [10, 30, 1000, 1500, 3000, 5000, 6000]]
levs.extend(
    np.arange(100, 1025, 25) * units.hPa
)  # large range useful for Chris Rozoff's CM1 model. Use wide range to prevent "out of contour range" error in NARR_composite.py.
levs.extend([lev * units.dimensionless for lev in ["lev1", "trop"]])

# --- 1. CONFIGURATION TEMPLATES ---
TEMPLATES = {
    "vort": {"levels": np.arange(-4, 40, 4), "cmap": "wind_17lev"},
    "div": {"levels": np.arange(-18, 22, 4), "cmap": "BlueWhiteOrangeRed"},
    "u": {"levels": range(-22, 26, 4), "cmap": "cmocean_balance", "idx": 0},
    "v": {"levels": range(-22, 26, 4), "cmap": "cmocean_balance", "idx": 1},
    "speed": {"levels": range(2, 34, 2), "cmap": "wind_17lev"},
    "wind": {"levels": range(2, 34, 2), "cmap": "wind_17lev"},
    "hgt": {"levels": range(0, 17000, 500), "cmap": "nice_gfdl", "var": "HGT"},
    "temp": {
        "levels": range(-65, 30, 5),
        "cmap": "nice_gfdl",
        "var": "TMP",
        "units": "degC",
    },
    "sh": {
        "levels": [1e-11, 0.01, 0.1, 1, 5, 10, 15, 20, 25],
        "cmap": "nice_gfdl",
        "var": "SPF_H",
        "units": "g/kg",
    },
    "rh": {
        "levels": range(0, 100, 10),
        "cmap": "CBR_drywet",
        "var": ["TMP", "SPF_H"],
        "units": "percent",
    },
}

# --- 2. DYNAMIC LEVEL GENERATION ---
for lev in levs:
    suffix = "HTGL" if hasattr(lev, "units") and lev.units == units.meters else "ISBL"

    for key, spec in TEMPLATES.items():
        f = fstr(key, lev)
        info = fieldinfo[f]
        info["vertical"] = lev
        info["levels"] = spec["levels"]

        raw_cmap = readNCLcm(spec["cmap"])
        if spec["cmap"] == "nice_gfdl":
            info["cmap"] = raw_cmap[3:193]
        elif spec["cmap"] == "CBR_drywet":
            info["cmap"] = raw_cmap[1:-1]
        else:
            info["cmap"] = raw_cmap

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

    vvel = fstr("vvel", lev)
    fieldinfo[vvel] = {
        "levels": [-250, -100, -25, -10, -2.5, -1, 1, 2.5, 10, 25, 100, 250],
        "cmap": readNCLcm("cmocean_balance")[::-1],
        "fname": f"V_VEL_221_{suffix}",
        "vertical": lev,
        "units": "microbar/second",
    }
    fieldinfo[vvel]["cmap"][127] = "white"

# --- 3. SURFACE & SPECIALIZED FIELDS ---
SURFACE_FIELDS = {
    "hfx": {
        "levels": list(range(-600, 125, 25)),
        "cmap": readNCLcm("amwg256")[::-1],
        "fname": "SHTFL_221_SFC",
    },
    "lh": {
        "levels": list(range(-700, 150, 50)),
        "cmap": readNCLcm("MPL_BrBG")[127:30:-1],
        "fname": "LHTFL_221_SFC",
    },
    "mslp": {
        "fname": "PRMSL_221_MSL",
        "levels": np.arange(956, 1028, 4),
        "units": "hPa",
    },
    "pblh": {"fname": "HPBL_221_SFC"},
    "t2": {"fname": "TMP_221_SFC", "units": "degF"},
    "sbcape": {"fname": "CAPE_221_SFC", "vertical": "surface-based"},
    "pwat": {
        "levels": [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70],
        "fname": "P_WAT_221_EATM",
        "temporal": 0,
    },
    "surface_height": {"fname": "HGT_221_SFC"},
    "bunkers": {
        "fname": ["USTM_221_HTGY", "VSTM_221_HTGY"],
        "levels": range(0, 60, 10),
        "vertical": "0-6km AGL",
    },
    "lcl": {
        "fname": "PRES_221_ADCL",
        "levels": [500, 600, 700, 750, 800, 850, 875, 900, 925, 950, 975, 1000],
        "units": "hPa",
    },
    "precipacc": {"fname": "RAINNC"},
    "shrtrop": {
        "levels": np.array([2, 5, 10, 15, 20, 30, 50]) * 1e-3,
        "fname": "VWSH_221_TRO",
        "vertical": "tropopause",
    },
}

for key, spec in SURFACE_FIELDS.items():
    fieldinfo[key].update(spec)

# Logic Overrides
fieldinfo["mslet"].update(fieldinfo["mslp"])
fieldinfo["mslet"]["fname"] = "MSLET_221_MSL"

fieldinfo["psfc"].update(fieldinfo["mslp"])
fieldinfo["psfc"]["fname"] = "PRES_221_SFC"

fieldinfo["zlcl"]["fname"] = ["PRES_221_ADCL", "HGT_221_ISBL"]
fieldinfo["lcl"]["cmap"] = [
    readNCLcm("nice_gfdl")[i]
    for i in [3, 20, 37, 54, 72, 89, 106, 123, 141, 158, 175, 193]
][::-1]

fieldinfo["sbcinh"].update({"fname": "CIN_221_SFC", "vertical": "surface-based"})
fieldinfo["sbcinh"]["cmap"].reverse()
fieldinfo["sbcinh"]["levels"] = [
    -np.round(x / 2) for x in reversed(fieldinfo["sbcinh"]["levels"])
]

# CAPE/CIN Mixes
fieldinfo["mlcape"].update({"fname": "CAPE_221_SPDY", "vertical": "mixed-layer"})
fieldinfo["mlcinh"].update(fieldinfo["sbcinh"])
fieldinfo["mlcinh"].update({"fname": "CIN_221_SPDY", "vertical": "mixed-layer"})
fieldinfo["mucape"]["vertical"] = "most unstable"

# Theta / Potential Temp
thcmap = ["#eeeeee", "#dddddd", "#cccccc", "#aaaaaa"] + readNCLcm("precip2_17lev")[3:-1]
fieldinfo["thetasfc"].update(
    {"levels": np.arange(290, 320, 2), "cmap": thcmap, "fname": "POT_221_SFC"}
)
fieldinfo["theta2"].update(
    {
        "levels": np.arange(294, 313, 1),
        "cmap": thcmap,
        "fname": ["PRES_221_HTGL", "TMP_221_HTGL"],
        "vertical": 2 * units.meters,
    }
)
fieldinfo["thetae2"].update(
    {
        "levels": np.arange(321, 375, 3),
        "cmap": thcmap,
        "fname": ["PRES_221_HTGL", "TMP_221_HTGL", "DPT_221_HTGL"],
        "vertical": 2 * units.meters,
    }
)

# Copy-Link Fields (Level Specific)
fieldinfo["rh_0deg"].update(fieldinfo["rh700"])  # defaultdict handles existence
fieldinfo["rh_0deg"].update({"fname": "R_H_221_0DEG", "vertical": "freezing level"})

fieldinfo["rh2"].update(fieldinfo["rh700hPa"])
fieldinfo["rh2"].update({"fname": "R_H_221_HTGL", "vertical": 2 * units.meters})

fieldinfo["rhlev1"].update(fieldinfo["rh700hPa"])
fieldinfo["rhlev1"].update({"fname": "R_H_221_HYBL", "vertical": "lowest model level"})

fieldinfo["vvellev1"].update(fieldinfo["vvel700hPa"])
fieldinfo["vvellev1"]["fname"] = "V_VEL_221_HYBL"

# Specific Humidity 2m / Level 1
sh_lvl_cfg = {
    "levels": [0.5, 1, 2, 4, 8, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24],
    "cmap": fieldinfo["td2"]["cmap"],
    "units": "g/kg",
}
fieldinfo["sh2"].update(sh_lvl_cfg)
fieldinfo["sh2"].update({"fname": "SPF_H_221_HTGL", "vertical": 2 * units.meters})
fieldinfo["shlev1"].update(sh_lvl_cfg)
fieldinfo["shlev1"].update(
    {"fname": "SPF_H_221_HYBL", "vertical": "lowest model level"}
)

# Convective Indices
fieldinfo["scp"].update(fieldinfo["stp"])
fieldinfo["tctp"].update(fieldinfo["stp"])
fieldinfo["tctp"]["levels"] = [2 * x for x in fieldinfo["tctp"].get("levels", [])]
fieldinfo["tctp2014"].update(fieldinfo["tctp"])

# --- 4. SHEAR, SRH & FLUXES ---
base_shear = fieldinfo.get("speed10m", {}).copy()
for bot, top in itertools.permutations(levs, 2):
    shr = f"shr{bot.m}{bot.units:~}_{top.m}{top.units:~}"
    fieldinfo[shr] = base_shear.copy()
    fieldinfo[shr]["vertical"] = [bot, top]

fieldinfo["shrtrop"]["cmap"] = fieldinfo.get("speed700", {}).get("cmap")

fieldinfo["srh"].update(fieldinfo.get("srh1", {}))
fieldinfo["srh"].update(
    {
        "fname": "HLCY_221_HTGY",
        "levels": list(fieldinfo["srh"].get("levels", [])) + [750],
        "cmap": list(fieldinfo["srh"].get("cmap", [])) + readNCLcm("wind_17lev")[-6:-4],
    }
)

# Fluxes
for flux in ["wvflux", "wcflux"]:
    prefix = "WV" if flux == "wvflux" else "WC"
    fieldinfo[flux].update(
        {
            "fname": [f"{prefix}UFLX_221_ISBY_acc3h", f"{prefix}VFLX_221_ISBY_acc3h"],
            "cmap": [],
            "levels": [],
        }
    )
    fieldinfo[flux]["sel"] = fieldinfo[flux]["fname"]

fieldinfo["wvfluxconv"].update(
    {
        "fname": "WVCONV_221_ISBY_acc3h",
        "levels": np.array(fieldinfo["precip"]["levels"]) * 10,
        "cmap": readNCLcm("prcp_1"),
    }
)
fieldinfo["wcfluxconv"].update(
    {
        "fname": "WCCONV_221_ISBY_acc3h",
        "levels": np.array(fieldinfo["precip"]["levels"]) * 2,
        "cmap": readNCLcm("prcp_1"),
    }
)


#######################################################################
idir = Path("/glade/campaign/collections/rda/data/d608000/3HRLY")  # path to NARR
#######################################################################

# apply this to rename_dims to avoid UserWarning: Horizonal dim numbers not found...Defaulting to (..., Y, X) order.
# Careful -- ncl_convert2nc names the x dimension gridy_221, and the y dimension gridx_221.
dims_dict = dict(
    gridx_221="y",  # Ny = 277
    gridy_221="x",  # Nx = 349
)


# Get static NARR file, convert to netCDF.
# Return netCDF filename.
def get_fixed(fieldname):
    narr = narrFixed
    targetdir = os.path.dirname(narr)
    os.makedirs(targetdir, exist_ok=True)
    # Convert to netCDF if netCDF file doesn't exist.
    if not narr.exists():
        narr_grb = idir.parent / "FIXED/rr-fixed.grb"
        call_args = ["ncl_convert2nc", narr_grb, "-e", "grb", "-o", targetdir]
        logging.warning(call_args)
        subprocess.check_call(call_args)

    nc = xarray.open_dataset(narr)
    logging.debug(f"rename {narr} {dims_dict}")
    nc = nc.rename_dims(dims_dict=dims_dict)
    return nc[fieldname]


surface_height = get_fixed(fieldinfo["surface_height"]["fname"]).metpy.quantify()


# Get NARR file from tar file, convert to netCDF.
# Return netCDF filename.
def get(valid_time, targetdir=targetdir, narrtype=narr3D):
    if narrtype == narrFixed:
        return narrFixed
    if not isinstance(valid_time, datetime.datetime):
        logging.warning(f"valid_time is not a datetime object {valid_time}")
    logging.debug(
        f"narr.get(): valid_time={valid_time} targetdir={targetdir} narrtype={narrtype} idir={idir}"
    )
    vartype, file_suffix = narrtype  # 3D, clm, flx, pbl, or sfc
    narr = os.path.join(
        targetdir, valid_time.strftime("merged_AWIP32.%Y%m%d%H." + file_suffix + ".nc")
    )
    # Convert to netCDF if netCDF file doesn't exist.
    if os.path.exists(narr):
        # TODO: Make sure file is complete. Another instance of narr.get() may still be writing it.
        logging.debug(f"narr.get(): found previously existing narr file {narr}")
        pass
    else:
        logging.warning(f"narr.get(): no previously existing narr file {narr}")
        narr_grb, ext = os.path.splitext(narr)  # drop '.nc' suffix
        # Extract NARR grib file from tar file if it doesn't exist.
        if not os.path.exists(narr_grb):
            search_str = valid_time.strftime("%Y/NARR" + vartype + "_%Y%m_") + "*.tar"
            logging.debug(f"narr.get(): search_str={search_str}")
            narrtars = idir.glob(search_str)
            for n in narrtars:
                prefix, yyyymm, dddd = str(n).split("_")
                dd1 = datetime.datetime.strptime(yyyymm + dddd[0:2], "%Y%m%d")
                # NARR3D_201704_0103.tar starts on 1st at 0Z and ends just before 4th at 0Z. so add one to eday
                dd2 = datetime.datetime.strptime(
                    yyyymm + dddd[2:4], "%Y%m%d"
                ) + datetime.timedelta(days=1)
                # make dd1 and dd2 timezone-aware if valid_time is timezone-aware
                if (
                    valid_time.tzinfo is not None
                    and valid_time.tzinfo.utcoffset(valid_time) is not None
                ):
                    logging.debug(
                        "making dd1 and dd2 timezone-aware cause valid_time is"
                    )
                    dd1 = pytz.utc.localize(dd1)
                    dd2 = pytz.utc.localize(dd2)
                if valid_time >= dd1 and valid_time < dd2:
                    break
            narrtar = n
            tar = tarfile.open(narrtar, mode="r")
            logging.info(
                f"extracting {os.path.basename(narr_grb)} from tar file {narrtar}"
            )
            tar.extract(os.path.basename(narr_grb), path=targetdir)
        call_args = ["ncl_convert2nc", narr_grb, "-e", "grb", "-o", targetdir]
        logging.debug(call_args)
        # Tried using pygrib or xarray.open_dataset(engine="cfgrib"), but neither was as good as ncl_convert2nc, even though
        # ncl_convert2nc calls the x dimension gridy_221 and y dimension gridx_221
        subprocess.check_call(call_args)
    return narr


def myunits(data, info):
    # Quantify xarray.
    data = data.metpy.quantify()  # quantify in MetPy>=1
    # Convert to info["units"].
    if "units" in info:
        logging.debug(f'converting to {info["units"]}')
        data = data.metpy.convert_units(
            info["units"]
        )  # Don't .sel(lv_HTGL3=10) before .metpy.convert_units('kt'). You get None in return.
    return data


def temporal(data, info):
    logging.debug(f"narr.temporal(): info={info}")

    # Define 'timetitle' attribute of data, no matter what.
    if "timetitle" not in data.attrs:
        data.attrs["timetitle"] = ""

    if "temporal" in info:
        time0 = info["temporal"]
        if hasattr(data.metpy, "time"):
            (temporal,) = data.metpy.coordinates("time")
            data = data.metpy.sel(time=time0)
        elif (
            len(data.shape) <= 2
        ):  # If data has only 2 dimensions assume it has no temporal dimension (like tropopause-level, or max-wind-level)
            data.attrs["timetitle"] = info["temporal"]
            logging.debug(
                f"narr.temporal(): assume no temporal dimension. setting timetitle={info['temporal']}"
            )
            return data
        else:
            temporal_str = data.dims[0]
            assert "time" in temporal_str
            logging.debug(
                f"narr.temporal(): metpy does not identify temporal coordinate. assume 1st dim is ({temporal_str})"
            )
            data = data[time0]
        data.attrs["timetitle"] = (
            str(time0) + "-h"
        )  # TODO: only works with zero. convert [ns] to hours.

    return data


def vertical(data, info):
    logging.debug(f"narr.vertical(): data {data} info {info}")

    if "verttitle" in data.attrs:
        return data

    # Define 'verttitle' attribute of data, no matter what.
    data.attrs["verttitle"] = ""

    # If data has a vertical coordinate dimension, select the appropriate level.

    if "vertical" in info:
        vlevel = info["vertical"]
        if hasattr(data.metpy, "vertical"):
            (vertical,) = data.metpy.coordinates("vertical")
            if vertical.name in data.dims:
                # If vertical has already been run on this DataArray, the vertical dimension may be length 1.
                # This zero-dimensional dimension was dropped and can't be selected again.
                # Tried DataArray.expand_dims() but the degenerate vertical dimension confused the ax.contour command.
                # contour expects a 2D array not 3D.
                data = data.metpy.sel(vertical=vlevel)
            else:
                logging.debug(
                    f"narr.vertical(): {data} {vertical.name} not in {data.dims}. assuming vertical has already had its way"
                )
            verttitle = str(vlevel)
        elif (
            len(data.dims) <= 2
        ):  # If data has only 2 dimensions assume it has no vertical dimension (like tropopause-level, or max-wind-level)
            logging.debug(f"narr.vertical(): {data} is 2D already.")
            logging.debug(f"narr.vertical(): setting verttitle={vlevel}")
            data.attrs["verttitle"] = vlevel
            return data
        else:
            # using list comprehension  # to get element with substring
            res = [i for i in data.dims if i.startswith("lv_")]
            vertical = res[0]
            logging.debug(
                f"narr.vertical(): metpy does not identify vertical coordinate. assume it is ({vertical})"
            )
            data = data.sel(vertical=vlevel)
            vertical = data.coords[vertical]
            verttitle = str(vlevel) + vertical.units
        logging.debug(f"narr.vertical(): setting verttitle {verttitle}")
        data.attrs["verttitle"] = verttitle

    return data


def shear(field, valid_time=None, targetdir=targetdir):
    # bottom and top vertical level are in fieldinfo[field][vertical]
    bot, top = fieldinfo[field]["vertical"]

    logging.debug(f"narr.shear(): bot={bot} top={top}")

    # winds are found in the flx or 3D file. Open both.
    ifiles = [
        get(valid_time, targetdir=targetdir, narrtype=narrtype)
        for narrtype in [narrFlx, narr3D]
    ]
    ds = xarray.open_mfdataset(ifiles)

    logging.debug(f"rename {dims_dict}")
    ds = ds.rename_dims(dims_dict=dims_dict)

    # ubot and vbot

    if bot.units == units.meters:
        ubot = ds["U_GRD_221_HTGL"].sel(lv_HTGL3=bot)
        vbot = ds["V_GRD_221_HTGL"].sel(lv_HTGL3=bot)
    elif bot.units == units.hPa:
        ubot = ds["U_GRD_221_ISBL"].sel(lv_ISBL0=bot)
        vbot = ds["V_GRD_221_ISBL"].sel(lv_ISBL0=bot)
    elif bot == "lev1":  # lowest model level
        ubot = ds["U_GRD_221_HYBL"]
        vbot = ds["V_GRD_221_HYBL"]
        ubot.attrs["verttitle"] = "bottom model level"
        vbot.attrs["verttitle"] = "bottom model level"
    else:
        logging.error(f"narr.shear(): unexpected bot {bot}")
        sys.exit(1)

    # utop and vtop
    if top.units == units.meters:
        if top.m in ds["U_GRD_221_HTGL"].metpy.vertical:
            utop = ds["U_GRD_221_HTGL"].sel(lv_HTGL3=top)
            vtop = ds["V_GRD_221_HTGL"].sel(lv_HTGL3=top)
        else:
            # height in meters AGL
            hgt3D = ds["HGT_221_ISBL"].metpy.quantify()
            agl3D = hgt3D - surface_height
            utop = hgtInterp(top, agl3D, ds["U_GRD_221_ISBL"])
            vtop = hgtInterp(top, agl3D, ds["V_GRD_221_ISBL"])
    elif top.units == units.hPa:
        utop = ds["U_GRD_221_ISBL"].sel(lv_ISBL0=top)
        vtop = ds["V_GRD_221_ISBL"].sel(lv_ISBL0=top)
    elif top == "trop":  #  tropopause
        utop = ds["U_GRD_221_TRO"]
        vtop = ds["V_GRD_221_TRO"]
        utop.attrs["verttitle"] = "tropopause"
        vtop.attrs["verttitle"] = "tropopause"
    else:
        print("narr.shear(): unexpected top {top}")
        sys.exit(1)

    du = utop - ubot
    dv = vtop - vbot
    du.attrs = utop.attrs
    dv.attrs = vtop.attrs
    if du.name is None:  # xarray.DataArray.name used in vectordata()
        du.name = utop.name
    if dv.name is None:
        dv.name = vtop.name
    du.attrs["long_name"] += " shear"
    dv.attrs["long_name"] += " shear"
    du.attrs["verttitle"] = f"{bot:~} to {top:~}"
    dv.attrs["verttitle"] = f"{bot:~} to {top:~}"
    return du, dv


def hgtInterp(x, xp, fp):
    assert x.units == xp.metpy.units, f"{x} and {xp} should be same units"
    # x = target height (scalar)
    # xp = input heights (3D)
    # fp = data to interpolate to target (3D)

    # If vertical dimension of variable goes from high altitude to low altitude (pressure dimension sorted in increasing order),
    # the np.diff() function produces a column of zeros and a -1, not a bunch of zeros and a 1, as expected. Then argmax() function misses the target layer.
    vdim = xp.metpy.vertical
    xp = xp.sortby(vdim, ascending=False)
    fp = fp.sortby(vdim, ascending=False)
    assert (
        xp.mean(dim=xp.dims[1:])[0] < xp.mean(dim=xp.dims[1:])[1]
    ), "mean xp heights should go low to high"
    # xp>x is False below the vertical layer that encompasses target
    # Once xp>x turns True, the np.diff function keeps the first occurrence of True in the vertical.
    xpgtx = xp.values
    xpgtx[~np.isnan(xpgtx)] = (
        xpgtx[~np.isnan(xpgtx)] > x.m
    )  # avoid RuntimeWarning: invalid value encountered in greater. Tried .fillna, .ffill, lots of stuff.
    bb = np.diff(
        xpgtx, axis=0
    )  # 3d boolean array. True at start of vertical layer that encompasses target
    k = bb.argmax(
        axis=0
    )  # vertical index of start of vertical layer that encompasses target
    ij = np.indices(xp.shape[1:])  # 2-d indices of spatial coordinates
    kij = (k, ij[0], ij[1])  # 3d indices of start of vertical layer
    rateofchange = np.diff(fp, axis=0) / np.diff(
        xp, axis=0
    )  # rate of change of f with respect to x
    rateofchange = rateofchange[kij]  # just for the layer that encompasses target
    fstart = fp.values[kij]  # f at start of layer
    dx = x.m - xp.values[kij]  # difference betweeen target and start of layer
    data = fstart + dx * rateofchange
    da = fp.mean(dim="lv_ISBL0")  # eliminate vertical dimension
    da.values = data
    # Tried being clever by expanding dims, making 1-element vertical dimension, but you have to use .sel later to eliminate it.
    da.attrs = fp.attrs
    del da.attrs["level_indicator"]

    return da


def multiInterp(x, xp, fp):
    # x = target vertical coordinate (2D)
    # xp = vertical coordinates of data to interpolate (1D)
    # fp = data to interpolate (3D)
    xp = np.broadcast_to(
        xp[:, None, None], fp.shape
    )  # broadcast 1D xp array across 2 new spatial dimensions of fp
    assert (
        xp.shape == fp.shape
    ), f"narr.multiInterp(): shapes of xp and fp differ {xp.shape} {fp.shape}"
    # xp>x is False below the vertical layer that encompasses target
    # Once xp>x turns True, the np.diff function keeps the first occurrence of True in the vertical.
    x = x.compute()  # compute dask array to avoid TypeError: unsupported operand type(s) for -: 'Array' and 'Array' in np.diff
    bb = np.diff(
        xp > x, axis=0
    )  # 3d boolean array. True at start of vertical layer that encompasses target
    k = bb.argmax(
        axis=0
    )  # vertical index of start of vertical layer that encompasses target
    ij = np.indices(x.shape)
    kij = (k, ij[0], ij[1])  # 3d indices of start of vertical layer
    rateofchange = np.diff(fp, axis=0) / np.diff(
        xp, axis=0
    )  # rate of change of f with respect to x
    rateofchange = rateofchange[kij]  # just for the layer that encompasses target
    fstart = fp[kij]  # f at start of layer
    dx = x - xp[kij]  # difference betweeen target and start of layer
    data = fstart + dx * rateofchange

    return data


def pressure_to_height(target_p, hgt3D):
    # target_p and vertical coordinate of hgt3D are both in hPa before taking natural log.
    lv_ISBL0 = hgt3D.lv_ISBL0.metpy.unit_array.to("hPa").m
    log_lv_ISBL0 = np.log(lv_ISBL0)
    log_target_p = np.log(target_p.metpy.unit_array.to("hPa").m)

    data = multiInterp(log_target_p, log_lv_ISBL0, hgt3D.values)

    # numpy array to xarray
    hgt2D = xarray.zeros_like(
        hgt3D.metpy.dequantify().mean(dim="lv_ISBL0", keep_attrs=True)
    )  # dequantify moves units to attributes
    hgt2D.values = data
    # with quantified data
    return hgt2D.metpy.quantify()


def scalardata(field: str, valid_time: datetime.datetime, targetdir: str = targetdir):
    # Get color map, levels, and netCDF variable name appropriate for requested variable (from fieldinfo dictionary).
    info = fieldinfo[field]

    # Make cmap a colors.ListedColormap, if it is not already.
    if "cmap" in info and not isinstance(info["cmap"], (colors.ListedColormap)):
        info["cmap"] = colors.ListedColormap(info["cmap"])
    logging.debug(f"scalardata: {field} info={info}")

    # Get narr file and filename.
    ifiles = [
        get(valid_time, targetdir=targetdir, narrtype=narrtype)
        for narrtype in [narrSfc, narrFlx, narrPBL, narr3D]
    ]

    logging.debug(f"about to open {ifiles}")
    # Would like to use parse_cf() but file is not CF compliant
    # not sure if this helps anything. x and y still are not projection coordinates; they are index values.
    nc = xarray.open_mfdataset(ifiles).metpy.assign_crs(data_crs.to_cf())

    # Avoid UserWarning: Horizontal dimension numbers not found.
    logging.debug(f"rename {dims_dict}")
    nc = nc.rename_dims(dims_dict=dims_dict)

    data = nc[info["fname"]]
    # Define data array. Speed and shear derived differently.
    # Define 'long_name' attribute

    if field.startswith("speed"):
        u = data[info["fname"][0]]
        v = data[info["fname"][1]]
        data = mcalc.wind_speed(u, v)
        data.attrs["long_name"] = "wind speed"
    elif field == "bunkers":
        # get units from 1st DataArray
        units = data[list(data.data_vars)[0]].units
        data = data.to_array(dim="uv", name=field)
        data.attrs["verttitle"] = info["vertical"]
        data.attrs["long_name"] = "Bunkers Storm Motion"
        data.attrs.update(units=units)
    elif field.startswith("div") or field.startswith("vort"):
        # TODO: why don't mcalc.div and vort work with vertical dimension? Have to sel a single level before calling them.
        u = data[info["fname"][0]]
        v = data[info["fname"][1]]
        lats = nc.gridlat_221
        lons = nc.gridlon_221
        dx, dy = mcalc.lat_lon_grid_deltas(lons, lats)
        dx = dx[np.newaxis, :]  # allow for vertical dimension
        dy = dy[
            np.newaxis, :
        ]  # avoid IndexError: too many indices for array: array is 2-dimensional, but 3 were indexed
        if field.startswith("div"):
            data = mcalc.divergence(u, v, dx=dx, dy=dy) * 1e5
            data.attrs["long_name"] = "divergence * 1e5"
        if field.startswith("vort"):
            data = mcalc.vorticity(u, v, dx=dx, dy=dy) * 1e5
            data.attrs["long_name"] = "vorticity * 1e5"
    elif field.startswith("shr") and "_" in field:
        du, dv = shear(field, valid_time=valid_time, targetdir=targetdir)
        data = mcalc.wind_speed(du, dv)
        data.attrs.update(
            {"long_name": "wind shear", "verttitle": du.attrs["verttitle"]}
        )
    elif (
        field.startswith("rh") and "lv_ISBL0" in data.coords
    ):  # could be 2-m RH or rh_0deg
        pres = data["lv_ISBL0"]
        temp = data[info["fname"][0]]
        sh = data[info["fname"][1]]
        data = mcalc.relative_humidity_from_specific_humidity(pres, temp, sh)
        data.attrs["long_name"] = "relative humidity"
    elif field == "theta2":
        prs = data[info["fname"][0]]
        tmp = data[info["fname"][1]]
        data = mcalc.potential_temperature(
            prs, tmp
        )  # Tried being clever and using *data, but complains about no units
        data = xarray.DataArray(data=data)
        data.attrs["long_name"] = "potential temperature"
    elif field == "thetae2":
        prs = data["PRES_221_HTGL"]
        tmp = data["TMP_221_HTGL"]
        dpt = data["DPT_221_HTGL"]
        data = mcalc.equivalent_potential_temperature(prs, tmp, dpt)
        data = xarray.DataArray(data=data)
        data.attrs["long_name"] = "equivalent potential temperature"
    elif field == "scp" or field == "stp" or field.startswith("tctp"):
        cape, cin, srh = data.data_vars.values()
        # compute to avoid ValueError: Cannot compare PlainQuantity and <class 'dask.array.core.Array'> in cin_term
        # quantify to avoid pint.errors.DimensionalityError: Cannot convert from 'dimensionless' to 'joule / kilogram' in cin_term
        cin = cin.compute().metpy.quantify()
        lifted_condensation_level_height = scalardata(
            "zlcl", valid_time, targetdir=targetdir
        )
        if field == "scp":
            bulk_shear = scalardata("shr10m_500hPa", valid_time, targetdir=targetdir)
            # In SPC help, cin is positive in SCP formulation.
            cin_term = -40 * units.parse_expression("J/kg") / cin
            cin_term = cin_term.where(
                cin < -40 * units.parse_expression("J/kg"), other=1
            )  # muCIN term based on Gropp and Davenport (2018), Aug issue of WaF. Set to 1.0 when muCIN > -40 kg-1.
            scp = (
                mcalc.supercell_composite(cape, srh, bulk_shear)
                * cin_term.metpy.unit_array
            )
            attrs = {"long_name": "supercell composite parameter"}
            data = xarray.DataArray(data=scp, attrs=attrs)
        elif field == "stp":
            bulk_shear = scalardata("shr10m_500hPa", valid_time, targetdir=targetdir)
            cin_term = (200 * units.parse_expression("J/kg") + cin) / (
                150 * units.parse_expression("J/kg")
            )
            cin_term = cin_term.where(
                cin <= -50 * units.parse_expression("J/kg"), other=1
            )
            cin_term = cin_term.where(
                cin >= -200 * units.parse_expression("J/kg"), other=0
            )
            # CAPE, srh, bulk_shear, cin may be one vertical level, but LCL may be multiple heights.
            # xarray.broadcast() makes them all multiple heights with same shape, so significant_tornado doesn't
            # complain about expecting lat/lon 2 dimensions and getting 3 dimensions..
            (cape, lifted_condensation_level_height, srh, bulk_shear, cin_term) = (
                xarray.broadcast(
                    cape, lifted_condensation_level_height, srh, bulk_shear, cin_term
                )
            )
            # Caveat, NARR storm relative helicity (srh) is 0-3 km AGL, while STP expects 0-1 km AGL.
            # Tried to ignore non-finite elements to avoid RuntimeWarning: invalid value encountered in greater/less but couldn't use 2-d boolean indexing with cape
            # cape and bulk_shear have different nans
            # metpy.calc.significant_tornado is the same as the official SPC mesoanalysis definition, but without the cin term.
            stp = (
                mcalc.significant_tornado(
                    cape, lifted_condensation_level_height, srh, bulk_shear
                )
                * cin_term.metpy.unit_array
            )
            attrs = {
                "long_name": "significant tornado parameter"
            }  # , 'verttitle':lifted_condensation_level_height.attrs['verttitle']} # don't want "2 meter" verttitle
            data = xarray.DataArray(data=stp, attrs=attrs)
        elif field.startswith("tctp"):
            bulk_shear = scalardata("shr10m_3000m", valid_time, targetdir=targetdir)
            if field == "tctp":
                tctp = (
                    srh
                    / (40 * units.parse_expression("m**2/s**2"))
                    * bulk_shear
                    / (12 * units.parse_expression("m/s"))
                    * (2000 * units.meters - lifted_condensation_level_height)
                    / (1400 * units.meters)
                )
                # NARR storm relative helicity (srh) is 0-3 km AGL, while TCTP expects 0-1 km AGL.
                attrs = {"long_name": "TC tornado parameter"}
            elif field == "tctp2014":
                RH24 = (
                    scalardata("rh800hPa", valid_time)
                    + scalardata("rh750hPa", valid_time)
                    + scalardata("rh700hPa", valid_time)
                    + scalardata("rh650hPa", valid_time)
                    + scalardata("rh600hPa", valid_time)
                ) / 5
                tctp = (
                    bulk_shear
                    / (11 * units.parse_expression("m/s"))
                    * srh
                    / (40 * units.parse_expression("m**2/s**2"))
                    * RH24
                    / (80 * units.percent)
                )
                attrs = {"long_name": "TC tornado parameter Eastin et al. 2014"}
            data = xarray.DataArray(data=tctp, attrs=attrs)
    elif field == "lcl":
        data.attrs["long_name"] = (
            "adiabatic lifting condensation level pressure"  # zlcl long_name based on this
        )
    elif field == "zlcl":
        LCL_pressure = scalardata("lcl", valid_time, targetdir=targetdir)
        hgt3D = data["HGT_221_ISBL"]
        data = pressure_to_height(LCL_pressure, hgt3D)
        data = data - surface_height
        data.attrs["long_name"] = LCL_pressure.attrs["long_name"].replace(
            "pressure", "height AGL"
        )
    else:
        if "sel" in info and info["sel"] is not None:  # this is a component of a vector
            attrs = (
                data[info["sel"]].attrs
            )  # remember attributes of sel component before .to_array removes them
            data = data.to_array(
                dim="uv"
            )  # convert from Dataset with 2 DataArrays to single DataArrray
            data.attrs.update(attrs)  # for long_name

    data = myunits(data, info)
    data = vertical(data, info)
    data = temporal(data, info)

    data.name = field
    # data.attrs['field'] = field  # maybe delete. stored in DataArray.name
    # use np.array to allow for levels to be a range
    levels = np.array(info["levels"])
    data.attrs["levels"] = levels
    data.attrs.update(info)

    return data


def vectordata(field, valid_time, targetdir=targetdir):
    # Get color map, levels, and netCDF variable name appropriate for requested variable (from fieldinfo module).
    info = fieldinfo[field]
    logging.debug(f"vectordata(): field={field} info={info}")
    if field.startswith("shr"):
        u, v = shear(field, valid_time, targetdir=targetdir)
        u = temporal(u, info)  # shear doesn't apply temporal like scalardata does.
        v = temporal(v, info)
        uv = xarray.merge(
            [u, v]
        ).to_array(
            dim="uv"
        )  # Tried concat, but didn't preserve the dataarray names or uv coordinate values (needed for uvsel).
        uv.attrs.update(
            info
        )  # shear() doesn't copy over attributes like scalardata does
        uv.attrs.update(u.attrs)
    elif field.endswith("flux"):
        uv = scalardata(field, valid_time, targetdir=targetdir)
    elif "fname" in info and isinstance(info["fname"], (list, tuple)):
        # handles 'bunkers'
        uv = scalardata(field, valid_time, targetdir=targetdir)
    else:
        uname = fstr("u", info["vertical"])
        uv = scalardata(uname, valid_time, targetdir=targetdir)

    # Fix sel attribute
    uv.attrs["sel"] = uv.uv.values  # select all (both) dimensions
    # Fix long_name, which was copied from u-component DataArray when you requested scalardata(uname).
    uv.attrs["long_name"] = (
        uv.attrs["long_name"].replace("u-component of ", "").replace("zonal ", "")
    )
    uv.attrs["field"] = (
        field  # 'field' attribute should have been added to u and v separately in scalardata().
    )
    return uv


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
    # nc: xarray Dataset with u, v, t, sh, hgt
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
        mucape = xarray.DataArray(coords=[nc.hrs, nc.storm, nc.point])
        mucin = xarray.DataArray(coords=[nc.hrs, nc.storm, nc.point])
        srh = xarray.DataArray(coords=[nc.hrs, nc.storm, nc.point])
        for point in nc.point:  # cape_cin is only a 1-D thing in MetPy.
            for hrs in nc.hrs:
                # .loc doesn't work with storm coordinate
                mucapes = xarray.DataArray(coords=[nc.storm])
                mucins = xarray.DataArray(coords=[nc.storm])
                srhs = xarray.DataArray(coords=[nc.storm])
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
