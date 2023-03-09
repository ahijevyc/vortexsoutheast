import argparse
import atcf
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import gridlines
import ibtracs
import itertools
import logging
import matplotlib.pyplot as plt
from metpy.units import units
import os
import pandas as pd 
import pdb
import spc

def tc_coords(df, trackdf):
    """
    get distance from TC center for each storm report in df
    """
    originlon = trackdf.loc[trackdf.valid_time == df.name,"lon"].values[0] * units.degrees_E
    originlat = trackdf.loc[trackdf.valid_time == df.name,"lat"].values[0] * units.degrees_N
    dist_from_origin, heading = spc.gdist_bearing(originlon, originlat, df["slon"].values * units.degrees_E, df["slat"].values * units.degrees_N)
    df["dist_from_origin"] = dist_from_origin
    return df

# =============Arguments===================
parser = argparse.ArgumentParser(description = "Plot NARR", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--debug", action='store_true')
parser.add_argument("-e", "--extent", type=float, nargs=4, help="debug plot extent lonmin, lonmax, latmin, latmax", 
        default=[-100.4, -75.6, 24.7, 39.5])
parser.add_argument("--spctd", type=float, default=1.5, help="hours on either side of valid time to show SPC storm reports")
parser.add_argument("--onetrackcolor", action="store_true", help="instead of coloring by intensity, color by track")
parser.add_argument("ifile", help="text file. each line starts with a storm name, followed by a yyyymmddhh time")

# Assign arguments to simple-named variables
args = parser.parse_args()
debug = args.debug
extent = args.extent
ifile = args.ifile
onetrackcolor = args.onetrackcolor
spc_td = pd.to_timedelta(args.spctd, unit='hours')

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
if debug:
    logger = logging.getLogger().setLevel(logging.DEBUG)

logging.debug(args)

logging.info(f"ifile {ifile}")
df = pd.read_csv(ifile, delim_whitespace=True, names=["stormname", "narrtime"])
df["narrtime"] = pd.to_datetime(df["narrtime"], format="%Y%m%d%H", utc=True)

logging.info("Load all storm reports and track data (subsample later).")
all_storm_reports = spc.get_storm_reports(start=pd.to_datetime("20000101",utc=True), end=pd.to_datetime("20200101",utc=True), event_types=["torn"])

# round storm report times to nearest narr
epoch = pd.to_datetime("19700101", utc=True)
f = ((all_storm_reports["time"] + spc_td - epoch) / (spc_td*2) ).astype(int)
all_storm_reports["time"] = epoch + f * (spc_td*2)

all_tracks, best_track_file = ibtracs.get_df(basin="NA") # North Atlantic basin guaranteed for NARR

logging.info(f"got {len(all_tracks)} track data from {best_track_file}")
fig = plt.figure(figsize=(11,8.5))
axc = plt.axes(projection=cartopy.crs.LambertConformal())
base, ext = os.path.splitext(ifile)
axc.set_title(os.path.basename(base))
stormrpts = []

colors = [None]
if onetrackcolor:
    colors = plt.cm.tab20.colors
color = itertools.cycle(colors)

for stormname, group in df.groupby("stormname"):
    onecolor = next(color)
    narrtime=group.narrtime
    year = narrtime.dt.year.min()
    logging.info(f"{stormname} {year}")
    imatch = (all_tracks['stormname'] == stormname.upper()) & (all_tracks['valid_time'].dt.year == year) & (all_tracks["rad"] == 34)
    trackdf = all_tracks[imatch]
    extension = ibtracs.extension(stormname, year)
    trackdf = pd.concat([trackdf, extension], ignore_index=True)
    trackdf["valid_time"] = pd.to_datetime(trackdf["valid_time"], utc=True)
    trackdf = trackdf[trackdf.valid_time.isin(narrtime)] # drop track valid_times not in narrtime array
    mintime = group.narrtime.min()
    maxtime = group.narrtime.max()
    logging.info(f"plot_track {mintime}-{maxtime}")
    trackdf = trackdf[(trackdf.valid_time >= mintime) & ( trackdf.valid_time <= maxtime)]
    atcf.plot_track(axc, "o", trackdf, stormname, label_interval_hours=0, scale=1.2, onecolor=onecolor)
    start = narrtime.min()-spc_td
    end = narrtime.max()+spc_td
    logging.info(f"get storm reports during narr {start} - {end}")
    twin = (all_storm_reports.time >= start ) & (all_storm_reports.time < end)
    if not any(twin):
        continue
    stormrpts_twin = all_storm_reports[twin]
    # restrict reports to within 800 km of TC
    s = stormrpts_twin.groupby("time", group_keys=False).apply(tc_coords,trackdf)
    s = s[s.dist_from_origin < 800]
    stormrpts.append(s)
    legend_items = spc.plot(s, axc, scale=5, alpha=1, onecolor=onecolor)

if not onetrackcolor:
    TCleg = atcf.TClegend(axc)
stormrpts = pd.concat(stormrpts)
legend_items = spc.plot(stormrpts, axc, scale=0)
axc.legend(handles=legend_items.values())

axc.set_extent(extent, crs=cartopy.crs.PlateCarree()) 
# *must* call draw in order to get the axis boundary used to add ticks:
fig.canvas.draw()

# Define gridline locations and draw the lines using cartopy's built-in gridliner:
xticks = list(range(-160,-50,10))
yticks = list(range(0,65,5))
axc.gridlines(xlocs=xticks, ylocs=yticks, linewidth=0.4, alpha=0.8, linestyle='--')
axc.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
axc.yaxis.set_major_formatter(LATITUDE_FORMATTER)
gridlines.lambert_xticks(axc, xticks)
gridlines.lambert_yticks(axc, yticks)
axc.set_xlabel('')
axc.set_ylabel('')
axc.add_feature(cartopy.feature.STATES.with_scale('50m'), linewidth=0.35, alpha=0.6)
axc.add_feature(cartopy.feature.COASTLINE.with_scale('50m'), linewidth=0.5, alpha=0.5)

ofile = os.path.realpath(os.path.basename(base)+".png")
if onetrackcolor:
    ofile = os.path.realpath(os.path.basename(base)+".colorbytrack.png")
plt.savefig(ofile)
logging.info(f'created {ofile}')
