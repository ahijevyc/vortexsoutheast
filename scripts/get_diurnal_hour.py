import datetime
import pandas as pd
import sys

"""
Given a (list of) storm and first time of 24-h time window,
Output storm name and times for requested hour[s]-of-day (0, 3, 6, 9, 12, 15, 18, or 21).
"""

ifile = sys.argv[1]
diurnal_hours = sys.argv[2:]
diurnal_hours = [int(x) for x in diurnal_hours]

df0 = pd.read_csv(
    ifile,
    delim_whitespace=True,
    names=["stormname", "start"],
    parse_dates=["start"],
    date_parser=lambda x: datetime.datetime.strptime(x, "%Y%m%d%H"),
)
df = df0.copy()

# append start time +3 hour, +6 hour, ... , +21 hour for each stormname.
for hours in range(3, 24, 3):
    df_plus = df0.copy()
    df_plus["start"] = df0["start"] + datetime.timedelta(hours=hours)
    df = pd.concat([df, df_plus], ignore_index=True)
df["hr"] = df.start.dt.hour

# grab requested hours
out = df[df["hr"].isin(diurnal_hours)]

# sort values
out = out.sort_values(["stormname", "start"], axis="index", ignore_index=True)

# output stormname and yyyymmddhh for requested hours
for index, row in out.iterrows():
    print(f'{row["stormname"]} {row["start"].strftime("%Y%m%d%H")}')
