import argparse
import os

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(
    description="given stormtime file, return storm list for NARR_composite.py with filtered times of diurnal cycle",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("ifile", help="input file")
parser.add_argument("hr0", type=int, help="anchor hour")
parser.add_argument("twin", type=int, default=6, help="time window in hours")
parser.add_argument("i", type=int, default=0, help="group number")
args = parser.parse_args()

ifile = args.ifile
hr0 = args.hr0
twin = args.twin
i = args.i

fmt = "%Y%m%d%H"
df = []
with open(ifile, "r") as f:
    for line in f:
        storm, start = line.split()
        start = pd.to_datetime(start, format=fmt)
        narrtimes = pd.date_range(
            start=start,
            end=start + pd.Timedelta(1, unit="day"),
            freq="3H",
            inclusive="left",
        )
        for narrtime in narrtimes:
            df.append([storm, narrtime])

df = pd.DataFrame(df, columns=["storm", "narrtime"])


def hh(hr0, twin, dt=3):
    x = np.array(range(0, 24, dt))
    x = x + hr0
    x = x % 24
    n = int(twin / dt)
    for i in range(0, len(x), n):
        yield x[i : i + n]
    return x


groups = hh(hr0, twin)
hh = list(groups)[i]
zero_padded_hours = [f"{x:02d}" for x in hh]
df = df[df.narrtime.dt.hour.isin(hh)]
base = os.path.basename(ifile)
base, ext = os.path.splitext(base)
ofile = os.path.join(
    os.getenv("TMPDIR"), base + "." + "".join(zero_padded_hours) + "z.txt"
)
df = df.sort_values(["storm", "narrtime"])
df["narrtime"] = df["narrtime"].dt.strftime(fmt)
df.to_csv(ofile, index=False, header=False, sep=" ")
