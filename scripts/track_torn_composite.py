import argparse
import atcf
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import gridlines
import ibtracs
import logging
import matplotlib.pyplot as plt
import os
import pandas as pd 
import pdb
import spc

# =============Arguments===================
parser = argparse.ArgumentParser(description = "Plot NARR", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--debug", action='store_true')
parser.add_argument("-e", "--extent", type=float, nargs=4, help="debug plot extent lonmin, lonmax, latmin, latmax", 
        default=[-99, -76, 24, 40])
parser.add_argument("--spctd", type=float, default=1.5, help="hours on either side of valid time to show SPC storm reports")
parser.add_argument("ifile", help="text file. each line starts with a storm name, followed by a yyyymmddhh time")

# Assign arguments to simple-named variables
args = parser.parse_args()
debug = args.debug
extent = args.extent
spc_td = pd.to_timedelta(args.spctd, unit='hours')
ifile = args.ifile

logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
if debug:
    logger = logging.getLogger().setLevel(logging.DEBUG)

logging.debug(args)

logging.info(f"ifile {ifile}")
df = pd.read_csv(ifile, delim_whitespace=True, names=["stormname", "narrtime"])
df["narrtime"] = pd.to_datetime(df["narrtime"], format="%Y%m%d%H", utc=True)

logging.info("Load all storm reports and track data (subsample later).")
all_storm_reports = spc.get_storm_reports(start=pd.to_datetime("20000101",utc=True), end=pd.to_datetime("20200101",utc=True))
all_tracks, best_track_file = ibtracs.get_df(basin="NA") # North Atlantic basin guaranteed for NARR

logging.info(f"got {len(all_tracks)} track data from {best_track_file}")
fig = plt.figure(figsize=(10,8))
axc = plt.axes(projection=cartopy.crs.LambertConformal())
base, ext = os.path.splitext(ifile)
axc.set_title(os.path.basename(base))
stormrpts = []
for stormname, group in df.groupby("stormname"):
    narrtime=group.narrtime
    year = narrtime.dt.year.min()
    logging.info(f"{stormname} {year}")
    imatch = (all_tracks['stormname'] == stormname.upper()) & (all_tracks['valid_time'].dt.year == year)
    trackdf = all_tracks[imatch]
    extension = ibtracs.extension(stormname, year)
    trackdf = pd.concat([trackdf, extension], ignore_index=True)
    logging.info(f"plot_track")
    atcf.plot_track(axc, "o", trackdf, stormname, label_interval_hours=None)
    start = narrtime.min()-spc_td
    end = narrtime.max()+spc_td
    logging.info(f"get storm reports {start} - {end}")
    twin = (all_storm_reports.time >= start ) & (all_storm_reports.time < end)
    # TODO: restrict reports to within 800 km of TC
    stormrpts.append(all_storm_reports[twin])

stormrpts = pd.concat(stormrpts)
stormrpts = stormrpts[stormrpts["event_type"].str.contains("torn")]
legend_items = spc.plot(stormrpts, axc, scale=3)
axc.legend(handles=legend_items.values())
TCleg = atcf.TClegend(axc)

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
plt.savefig(ofile)
logging.info(f'created {ofile}')
