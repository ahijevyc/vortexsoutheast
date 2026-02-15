import argparse
import datetime
import logging
import os
from pathlib import Path

import matplotlib.pyplot as plt
import seaborn as sns
import VSE
import xarray
from metpy.calc import wind_speed

parser = argparse.ArgumentParser(
    description="line plot of output from NARR_composite.py",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument("-d", "--debug", action="store_true", help="print debug messages")
parser.add_argument(
    "--desc",
    default="tornadoes_near_coast",
    choices=VSE.composites(),
    help="description of composite category",
)
parser.add_argument("field", default="mlcape", help="NARR field")
parser.add_argument("--clobber", action="store_true", help="overwrite old file")
parser.add_argument(
    "--coord",
    choices=["north", "storm motion", "wind shear"],
    default="north",
    help="coordinate system. fields rotated so that this vector points up",
)
parser.add_argument(
    "--dpi", type=int, default=150, help="output resolution in dots per inch"
)
parser.add_argument(
    "--idir",
    default=Path(__file__).parent.parent.absolute() / "output/composites/nc",
    help="input directory",
)
parser.add_argument(
    "--no-fineprint",
    action="store_true",
    help="Don't write additional information at bottom of figure",
)
locations = parser.add_mutually_exclusive_group(required=False)
locations.add_argument(
    "--centroid",
    default="shr10_700 max",
    choices=VSE.centroids(),
    help="use predetermined list of locations defined by coord and centroid",
)
locations.add_argument(
    "--location",
    nargs="+",
    type=str,
    default=None,
    help="label/azimuth/range. azimuth in deg, starting at north and increasing cw. range in km",
)
parser.add_argument(
    "--simplegend",
    action="store_true",
    help="no azimuth and range in legend location labels",
)
parser.add_argument("--twin", type=int, default=3, help="time window in hours")
args = parser.parse_args()

centroid = args.centroid
clobber = args.clobber
coord = args.coord
desc = args.desc
dpi = args.dpi
field = args.field
idir = args.idir
no_fineprint = args.no_fineprint
location = args.location
simplegend = args.simplegend
twin = args.twin

level = logging.DEBUG if args.debug else logging.INFO
logging.basicConfig(format="%(asctime)s - %(message)s", level=level, force=True)
logging.debug(args)

if centroid:
    location = VSE.pointlist[coord][centroid]
azimuths = [ar.split("/")[1].replace("deg", "") for ar in location]
ranges_km = [ar.split("/")[2].replace("km", "") for ar in location]
if simplegend:
    location_labels = [ar.split("/")[0] for ar in location]


hour = xarray.DataArray(
    data=["09z", "12z", "15z", "18z", "21z", "00z", "03z", "06z", "09z", "12z"],
    coords={"hour": range(-1, 9)},
    dims="hour",
    name="hour",
    attrs={"format": "%Hz"},
)
if twin == 6:
    hour = xarray.DataArray(
        data=["0609z", "1215z", "1821z", "0003z", "0609z", "1215z"],
        coords={"hour": [-1, 0, 1, 2, 3, 4]},
        dims="hour",
        name="hour",
        attrs={"format": "H1H2z"},
    )
    # hour = xarray.DataArray(data=["1215z","1821z","0003z","0609z"], coords={"hour":[0,1,2,3]}, dims="hour", name="hour", attrs={"format":"H1H2z"})

print(hour)

dss = []

files = [f"{idir}/{desc}.{field}.{h}.nc" for h in hour.values]
logging.info(f"opening {files}")
for f in files:
    dss.append(xarray.open_dataset(f))

ds = xarray.concat(dss, dim=hour)
# TODO: extend times after filling DataSet
hour.values[0] = f"pre {hour.values[0]}"
hour.values[-1] = f"post {hour.values[-1]}"
if twin == 6:
    # Extend line before 1215z and after 0609z.
    hour.values = ["pre 12-15z", "12-15z", "18-21z", "00-03z", "06-09z", "post 06-09z"]
ds = ds.assign_coords({"hour": hour})

target_azimuths = xarray.DataArray(azimuths, dims=["location"])
target_ranges = xarray.DataArray(ranges_km, dims=["location"])

# Calculate magnitude from u- and v- components.
if "component" in ds[field].coords:
    ds[field] = wind_speed(
        ds[field].sel(component="u"), ds[field].sel(component="v")
    ).metpy.dequantify()  # keep units as attribute
narrds = ds.sel(range=target_ranges, method="nearest", tolerance=20.0)
narrds = narrds.sel(azimuth=target_azimuths, method="nearest", tolerance=0.5)
narrds = narrds.sel(coord=coord)

logging.info(narrds.coords)

narrds = narrds.assign_coords(location=location_labels)

# narrds['storm'] = narrds.storm.astype(str) # convert object to str (does it help?)


fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.18)
# convert xarray to pandas dataframe to avoid ValueError: Could not interpret value `hour` for parameter `x`
# 95% CIs change due to bootstrapping. (unless seed is set)
ci, n_boot = 95, 10000
# tried data=narrds, but ValueError: arrays must all be same length
ax = sns.lineplot(
    data=narrds.to_dataframe(),
    x="hour",
    y=field,
    hue="location",
    errorbar=("ci", ci),
    n_boot=n_boot,
    style="location",
    dashes=False,
    markers=["+", "+"],
    seed=14,
    ax=ax,
)
ax.set(xlim=(0.5, hour.size - 1.5))
ax.set(ylabel=f"{field} [{narrds[field].units}]")
fontsize = "small" if simplegend else "x-small"
plt.setp(ax.get_legend().get_texts(), fontsize=fontsize)
plt.setp(ax.get_legend().get_title(), fontsize=fontsize)
ax.grid(alpha=0.5, lw=0.5)

plt.suptitle(f"{desc}  {coord} points up", wrap=True)
fineprint_text = f"{ci}% confidence interval  n_boot={n_boot}"
fineprint_text += f"\n{narrds.storm.size} storms"
fineprint_text += f"\ncreated {datetime.datetime.now()}"
fineprint = plt.annotate(
    text=fineprint_text,
    xy=(4, 1),
    xycoords=("figure pixels", "figure pixels"),
    va="bottom",
    fontsize=4.9,
)
if no_fineprint:
    fineprint.set_visible(False)

ofile = f"{idir}/{desc}.{coord.replace(' ','_')}.{field}.pdf"
plt.savefig(ofile, dpi=dpi)
logging.info(f"made {os.path.realpath(ofile)}")

narrds.to_dataframe().sort_values(field)
