import argparse
import datetime
import logging
import os

import ibtracs
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import VSE
from ahijevyc import spc
from pandas.api.types import CategoricalDtype

sns.set_theme()


def countblock(x):
    s = x.copy()
    s[x == 0] = "0"
    s[(x >= 1) & (x <= 4)] = "1-4"
    s[(x >= 5) & (x <= 9)] = "5-9"
    s[(x >= 10) & (x <= 19)] = "10-19"
    s[(x >= 20) & (x <= 49)] = "20-49"
    s[x > 50] = "50+"
    return s


# =============Arguments===================
parser = argparse.ArgumentParser(
    description="Plot NARR", formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-d", "--debug", action="store_true")
parser.add_argument(
    "--spctd",
    type=float,
    default=1.5,
    help="hours on both sides of valid time to show tornadoes",
)
parser.add_argument(
    "ifiles",
    nargs="+",
    help="input text files. each line starts with a storm name, followed by a yyyymmddhh time",
)
parser.add_argument("--ofile", help="output file name")

# Assign arguments to simple-named variables
args = parser.parse_args()
debug = args.debug
ifiles = args.ifiles
ofile = args.ofile
spc_td = pd.to_timedelta(args.spctd, unit="hours")

level = logging.DEBUG if debug else logging.INFO
logging.basicConfig(format="%(asctime)s %(message)s", level=level)

logging.debug(args)

logging.info("load tornado reports (subsample later).")
all_storm_reports = spc.get_storm_reports(
    start=pd.to_datetime("19900101", utc=True),
    end=pd.to_datetime("20220101", utc=True),
    event_types=["torn"],
)

# round storm report times to nearest narr
epoch = pd.to_datetime("19700101", utc=True)
f = ((all_storm_reports["time"] + spc_td - epoch) / (spc_td * 2)).astype(int)
all_storm_reports["time3h"] = epoch + f * (spc_td * 2)

TCTOR = spc.getTCTOR()

logging.info("Load ibtracs (subsample later).")
all_tracks, best_track_file = ibtracs.get_df(
    basin="NA"
)  # North Atlantic basin guaranteed for NARR

logging.info(f"got {len(all_tracks)} track data from {best_track_file}")

fmt = "%Y%m%d %H%Mz"
notorn = []
stormrpts = []
for ifile in ifiles:
    logging.info(ifile)
    df = pd.read_csv(ifile, delim_whitespace=True, names=["stormname", "narrtime"])
    df["narrtime"] = pd.to_datetime(df["narrtime"], format="%Y%m%d%H", utc=True)
    df["year"] = df["narrtime"].dt.year
    category, ext = os.path.splitext(os.path.basename(ifile))

    for (year, stormname), group in df.groupby(["year", "stormname"]):
        assert len(group) == 1, (
            "stormname/year combo should be unique in each input file."
            "Do not use input file with all times of diurnal cycle. just 1st one"
        )
        stormyear = f"{stormname} {year}"
        # first narrtime and 3-h intervals covering a diurnal cycle
        narrtime0 = pd.to_datetime(group.narrtime.values[0])
        narrtimes = pd.date_range(
            start=narrtime0,
            end=narrtime0 + pd.Timedelta(21, unit="hour"),
            freq="3H",
            tz="UTC",
        )
        imatch = (
            (all_tracks["stormname"] == stormname.upper())
            & (all_tracks["valid_time"].dt.year == year)
            & (all_tracks["rad"] == 34)
        )
        assert imatch.sum(), f"no matching IBTRACS for {stormname.upper()} {year}"
        trackdf = all_tracks[imatch]

        in_narr = trackdf.valid_time.isin(narrtimes)
        if in_narr.sum() < 8:
            logging.warning(
                f"only {in_narr.sum()} track times in {narrtime0.strftime(fmt)} diurnal cycle. checking extensions"
            )
            trackdf = ibtracs.getExt(stormname, year, trackdf, narrtimes)

        logging.debug(
            f"{stormname} {year} {len(trackdf)} times "
            f"[{trackdf.valid_time.min().strftime(fmt)},{trackdf.valid_time.max().strftime(fmt)}]"
        )

        start = trackdf.valid_time.min() - spc_td
        end = trackdf.valid_time.max() + spc_td
        logging.debug(f"storm reports [{start.strftime(fmt)},{end.strftime(fmt)})")
        twin = (all_storm_reports.time >= start) & (all_storm_reports.time < end)
        if not any(twin):
            logging.info("no storm reports")
            notorn.append((category, stormyear))
            continue
        stormrpts_twin = all_storm_reports[twin]

        # restrict reports to within 800 km of TC
        s = stormrpts_twin.groupby("time3h", group_keys=False).apply(
            VSE.tc_coords, trackdf
        )
        s = s[s.dist_from_origin < 800]
        # restrict distance for some cases to match TCTOR
        if stormname == "Erin" and year == 2007:
            s = s[s.dist_from_origin < 700]
        if stormname == "Olga" and year == 2007:
            s = s[s.dist_from_origin < 688]
        if stormname == "Hermine" and year == 2004:
            s = s[s.dist_from_origin < 475]
        s["category"] = category
        s["stormyear"] = stormyear

        logging.info(f"{len(s)} storm reports near {stormname} {year}")

        TCTOR_twin = TCTOR[
            (TCTOR.time >= start)
            & (TCTOR.time < end)
            & (TCTOR["Tor-Center Dist km (from Worksheet)"] < 800)
            & (TCTOR["TC Name"] == f"{stormname}-{year%100:02d}")
        ]
        if len(s) != len(TCTOR_twin):
            logging.warning(f"{len(s)} SPC, but {len(TCTOR_twin)} TCTOR reports")
            if abs(len(s) - len(TCTOR_twin)) > 1:
                print(s[["time", "slat", "slon", "dist_from_origin"]])
                print(
                    TCTOR_twin[
                        ["time", "slat", "slon", "Tor-Center Dist km (from Worksheet)"]
                    ]
                )
        stormrpts.append(s)

stormrpts = pd.concat(stormrpts)

stormrpts["day"] = (stormrpts.time.dt.time < datetime.time(hour=1, minute=30)) | (
    stormrpts.time.dt.time >= datetime.time(hour=13, minute=30)
)
stormrpts.loc[stormrpts.day, "time"] = "day"
stormrpts.loc[~stormrpts.day, "time"] = "night"
hh = stormrpts.groupby(["category", "stormyear", "time"]).size()
for x in notorn:
    hh.loc[(*x, "n/a")] = 0  # Assign zero-count storms to time="n/a"
hh.name = "tornadoes"

hh = countblock(hh)
cat_type = CategoricalDtype(
    categories=["0", "1-4", "5-9", "10-19", "20-49", "50+"], ordered=True
)

hh = hh.astype(cat_type)
sns.displot(
    data=hh.reset_index(),
    x="tornadoes",
    kind="hist",
    multiple="stack",
    hue="time",
    hue_order=["night", "day", "n/a"],
    col="category",
    col_wrap=3,
)
plt.yticks(np.arange(0, 20, 4))

if ofile is None:
    ofile = os.path.realpath("torn_hist.ps")
plt.savefig(ofile)
logging.info(f"created {ofile}")
