import itertools

import numpy as np
from ahijevyc.fieldinfo import readNCLcm
from metpy.units import units


def fstr(f, lev):
    """Format string for level keys (e.g., 'u700', 'temp2m')."""
    return f"{f}{lev:~}".replace(" ", "")


def setup_fieldinfo(fieldinfo):
    """
    Builds the fieldinfo dictionary.
    Uses a template factory approach for maintainability.
    """

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
