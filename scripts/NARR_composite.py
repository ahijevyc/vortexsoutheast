"""

Fill contour, line contour, wind barbs, and quiver.

Tried import pygrib and working with grib file directly but the variables were not translated correctly. 
Instead, convert to netCDF with ncl_convert2nc.

"""
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import argparse
import atcf
import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import datetime
import gridlines
import ibtracs
import logging
from metpy.units import units
import numpy as np
import narr # fieldinfo levels and color tables and narr.get()
import ncl 
import os
import pandas as pd # for forward fill of NaNs
import pdb
import pytz
import re
from scipy import spatial
import spc
import sys
import xarray

def isvector(data):
    # Is 'uv' one of the coordinates?
    return 'uv' in data.coords


def uvsel(data):
    # If 'uv' is one of the coordinates, return component(s) listed in "sel" attribute.
    # "sel" attribute should have been set in narr.py
    # Since "sel" may be a one-element list, squeeze the dimension.
    if isvector(data):
        data = data.sel(uv=data.attrs["sel"]).squeeze()
    return data


def rotate_vector(vector,az):
    if not isvector(vector):
        return vector
    az = np.radians(az)
    u = vector.isel(uv=0).copy() # copy() is VERY important or else u remains a reference to uv and is changed by first .loc command below
    v = vector.isel(uv=1).copy()
    vector.values[0] = np.cos(az)*u - np.sin(az)*v
    vector.values[1] = np.sin(az)*u + np.cos(az)*v
    return vector

def roll_rotate_sel(values, heading, debug=False):
    # Get daz from azimuth coordinate of values DataArray.
    daz = values["azimuth"].diff(dim="azimuth").values
    assert daz.min() == daz.max(), 'roll_rotate_sel(): not all daz are the same'
    # assume units are degrees. Remove if we eventually figure out how to put units on coordinates without breaking matplotlib
    daz = daz[0] * units.deg
    nrotate = int(round(heading / daz))
    shifts = {"azimuth": -nrotate}
    values = values.roll(shifts=shifts, roll_coords=False) # rotate azimuthal dimension 
    values = rotate_vector(values, heading) # rotate vector
    values = uvsel(values)
    return values


def histogram2d_weighted(bearing, dist_from_center, azbins, rbins, data):
    bins = [azbins, rbins]
    bearing = bearing.values.flatten()
    dist_from_center = dist_from_center.values.flatten()
    data = data.metpy.dequantify() # move units to attrs dict so we may assign them to new DataArray.
    n, _, _ = np.histogram2d(bearing, dist_from_center, bins=bins)
    drs = np.diff(rbins) # bin widths
    dazs = np.diff(azbins)
    assert np.all(drs == drs[0]), "range bins not equal width"
    assert np.all(dazs == dazs[0]), "azimuth bins not equal width"
    range_center1D, theta_center1D = drs[0]/2+rbins[:-1], dazs[0]/2+azbins[:-1]
    if isvector(data):
        weightedh2 = np.zeros((2, len(bins[0])-1, len(bins[1])-1))
        weightedh2[0], _, _ = np.histogram2d(bearing, dist_from_center, bins=bins, weights=data.isel(uv=0).values.flatten())
        weightedh2[1], _, _ = np.histogram2d(bearing, dist_from_center, bins=bins, weights=data.isel(uv=1).values.flatten())
        coords={"uv": data.coords["uv"], "azimuth":theta_center1D, "range":range_center1D}
    else:
        weightedh2, _, _ = np.histogram2d(bearing, dist_from_center, bins=bins, weights=data.values.flatten())
        coords={"azimuth":theta_center1D, "range":range_center1D}
    n = np.ma.masked_where(n == 0,n)# avoid RuntimeWarning: invalid value encountered in true_divide
    da = xarray.DataArray(data = weightedh2/n, coords=coords, dims=coords.keys(), attrs=data.attrs, name=data.name) 
    # Fill NaNs with neighbors
    da = fill_NaNs_with_neighbors(da)
    da = da.metpy.quantify() # move units from attrs dict back to data array (make it a pint quantity)
    return da

def fill_NaNs_with_neighbors2d(array2d):
    # Fill NaNs with neighbors
    # Input DataArray
    # convert to ndarray
    values = array2d.values
    values = pd.DataFrame(values)
    # operate on DataFrame
    values = values.interpolate(limit_direction='both') # goofy around NARR edges
    if values.isnull().all().any(): # if any columns are all null
        print("***Found nan in filled value array. trying to interpolate across columns***") # Think this is fine but not sure.
        # Happens when you zoom in a lot.
        values = values.interpolate(axis='columns')
        if values[0].isnull().all(): # if column 0 is still all null, backfill from neighboring column
            values.fillna(method="backfill", axis="columns", inplace=True)
            print("***Backfilled missing data in column 0***")
    # Output DataArray
    array2d.values = values
    return array2d


def add_rgrid(ax, ring_interval=200):
    ax.grid(alpha=0.4, linewidth=0.4, color='k') # tried moving outside loop right after figure was created but grid was not visible (even with zorder=3)
    ax.tick_params(labelsize='xx-small')
    #ax.yaxis.set_units('km') # Not sure if this helps anything. Certainly not when axis is normalized.
    start, stop = ax.get_ylim()
    radii = np.arange(start,stop+ring_interval,ring_interval)
    rlabels = [f"{x:.0f}km" for x in radii]
    for i in range(len(rlabels)-1):
        rlabels[i] = ""
    logging.debug(f"add_rgrid: ring_interval={ring_interval} start={start}, stop={stop} radii={radii} {rlabels}")
    lines, labels = ax.set_rgrids(radii, labels=rlabels, angle=0, ha='center',va="bottom", fontsize=4.5) # va="baseline" puts label too close to ring
    return lines, labels
    
def fill_NaNs_with_neighbors(values):        
    if len(values.shape) == 2:
        return fill_NaNs_with_neighbors2d(values)
    if len(values.shape) == 3:
        for i, aa in enumerate(values):
            aa = fill_NaNs_with_neighbors2d(aa)
            values[i] = aa.values
    return values

def xr_list_mean(x):
    # mean of list of xarray DataArrays
    return xarray.concat(x, dim="storm").mean(dim="storm")


# =============Arguments===================
parser = argparse.ArgumentParser(description = "Plot NARR", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-b", "--barb", type=str, default= None, help='variable name for wind barb field')
parser.add_argument("--cart", action='store_true', help="plot cartopy version (not TC-centric)")
parser.add_argument("--clevels", type=float, nargs='+', default= [500, 1000, 2000, 3000, 4000, 5000, 10000], help='line contour levels')
parser.add_argument("--clobber", action='store_true', help="overwrite any old outfile, if it exists")
parser.add_argument("-d", "--debug", action='store_true')
parser.add_argument("-e", "--extent", type=float, nargs=4, help="debug plot extent lonmin, lonmax, latmin, latmax", default=[-98, -74, 20, 38])
parser.add_argument("-f", "--fill", type=str, default= 'shr10m_700hPa', help='variable name for contour fill field')
parser.add_argument("--hail", action='store_true', help="overlay hail reports")
parser.add_argument("-l", "--line", type=str, default=None, help='variable name for line contour field')
parser.add_argument("--max_range", type=float, default=800., help="maximum range in km, or, if normalize_by is set, units of normalize_by.")
parser.add_argument("--netcdf", type=str, default=None, help="netCDF file to write composite to")
parser.add_argument("--no-fineprint", action='store_true', help="Don't write details at bottom of image")
parser.add_argument("--normalize_by", type=str, choices=["rmw","r34","Vt500km"], default=None, help="normalize range from TC center by this scale")
parser.add_argument("--no-torn", action='store_false', help="don't overlay tornado reports")
parser.add_argument("-o", "--ofile", type=str, help="name of final composite image")
parser.add_argument("-q", "--quiver", type=str, default= None, help='variable name for quiver field')
parser.add_argument("--spctd", type=float, default=1.5, help="hours on either side of valid time to show SPC storm reports")
parser.add_argument("ifiles", nargs="+", type=argparse.FileType('r'), help="text file(s). each line starts with an atcf filename, followed by a list of yyyymmddhh times")
parser.add_argument("--torn", action='store_true', help="overlay tornado reports")
parser.set_defaults(torn=True)
parser.add_argument("--wind", action='store_true', help="overlay wind reports")
parser.add_argument('-w', "--workdir", type=str, default="/glade/scratch/"+os.getenv("USER")+"/NARR", help="directory to untar NARR into")


# Assign arguments to simple-named variables
args = parser.parse_args()
barb           = args.barb
cart           = args.cart
clobber        = args.clobber
contour_levels = args.clevels
debug          = args.debug
extent         = args.extent
fill           = args.fill
hail           = args.hail
line           = args.line
max_range      = args.max_range
netcdf         = args.netcdf
no_fineprint   = args.no_fineprint
normalize_by   = args.normalize_by
quiver         = args.quiver
spc_td         = datetime.timedelta(hours=args.spctd)
ifiles         = args.ifiles
torn           = args.torn
wind           = args.wind
workdir        = args.workdir

logger = logging.getLogger()
logging.basicConfig(format='%(asctime)s %(message)s')
if debug:
    logger.setLevel(logging.DEBUG)

logging.debug(args)


if args.ofile:
    ofile = os.path.realpath(args.ofile)
else:
    ofile = "./"+".".join([x for x in [fill, line, barb, quiver] if x]) + '.composite.png'
if netcdf:
    netcdf = os.path.realpath(netcdf)

if normalize_by:
    bname, ext = os.path.splitext(ofile)
    ofile  =  bname + ".normalize_by_" + normalize_by + ext
    if netcdf:
        bname, ext = os.path.splitext(netcdf)
        netcdf = bname + ".normalize_by_" + normalize_by + ext
# If png already exists skip this file
if not clobber and os.path.isfile(ofile) and netcdf and os.path.isfile(netcdf):
    print(ofile,"and",netcdf,"exist. Skipping. Use --clobber option to override.")
    sys.exit(0)

def desc(data):
    assert 'timetitle' in data.attrs, "timetitle not attribute of data"+str(data)
    assert 'verttitle' in data.attrs
    assert 'long_name' in data.attrs
    timetitle = data.attrs['timetitle']
    verttitle = str(data.attrs['verttitle'])
    desc = f"{timetitle} {verttitle} {data.long_name} [{data.metpy.units:~}]" # abbreviate units.
    return desc


# Tried adding units to daz and dr but units really messed up polar pcolor and contour functions downstream.
daz = 1# azimuth bin width (delta azimuth)
dr  = 40# range bin width (delta range)
dr_units = "km"
if normalize_by:
    dr_units = normalize_by
    if normalize_by == "r34" and max_range > 100:
        dr  = 0.25 # normalized range bin width (delta range)
        max_range = 4.
    if normalize_by == "Vt500km":
        dr_units = "normalized km"
    if normalize_by == "rmw" and max_range > 100:
        print("normalize_by", normalize_by, "max_range", max_range)
        print("setting max_range to 10")
        print("setting dr to 0.5")
        dr = 0.5
        max_range = 10.
pbot = 850*units("hPa") # for wind shear coordinate
ptop = 200*units("hPa")
fineprint_string = f"daz: {daz}$\degree$   dr: {dr}{dr_units}   layer for wind shear coordinate:{ptop:~}-{pbot:~}"
fineprint_string += "\nstorm and time(s) text file(s): " + " ".join([os.path.realpath(x.name) for x in ifiles])
azbins = np.arange(0,360+daz,daz)
theta_lines = range(0,360,45)
storm_rpt_legend_kw = dict(fontsize=4.8, bbox_to_anchor=(0.8, -0.01), loc='upper left', borderaxespad=0, frameon=False, title_fontsize='xx-small')
linecontour_alpha = .6
linecontour_fontsize = 4
barblength = 3.5
barbunits = 'm/s'
barb_increments = {'half':2.5, 'full':5, 'flag':25} # for matplotlib.quiver.barb
barbkwdict = dict(length=barblength, linewidth=0.2, alpha=0.6, barb_increments=barb_increments)

rbins  = np.arange(0,max_range+dr,dr)
r2d, theta2d = np.meshgrid(rbins, azbins)
range_center1D, theta_center1D = dr/2+rbins[:-1], daz/2+azbins[:-1]
polar_sz = (len(azbins)-1, len(rbins)-1)

#---Locate gridded wind barbs on polar plot
dx = 100 # wind barb grid spacing
mg = np.arange(-max_range, max_range+dx, dx)
x0, y0 = np.meshgrid(mg, mg)
r_center2D, theta_center2D = np.meshgrid(range_center1D, theta_center1D)
xp = r_center2D * np.cos(np.radians(90-theta_center2D))
yp = r_center2D * np.sin(np.radians(90-theta_center2D))
tree = spatial.KDTree(np.c_[xp.ravel(),yp.ravel()])
dd, ii = tree.query(np.c_[x0.ravel(),y0.ravel()], k=1, distance_upper_bound=dx)
ii = ii[dd < dx/2] # remove neighbors over dx/2 away from nearest point


best_track_files, narr_files = [], [] # used for fineprint
lons, lats, times = [], [], [] # used for netCDF file
stormname_years, narr_valid_times = [], [] # for composite title

figplr, axes = plt.subplots(ncols=2,nrows=2,subplot_kw=dict(projection='polar'))
# default left=0.125 right=0.9 bottom=0.1 top=0.9 wspace=0.2 hspace=0.2
plt.subplots_adjust(left=-0.1, right=1.1, bottom=0.05, top=0.87, hspace=0.3, wspace=0)

# These dicts will have one entry for each coordinate system (axes).
filldict   = {}
linedict   = {}
barbdict   = {}
quiverdict = {}
storm_rpts_overlays = {}
# any reason to keep as 2x2 ndarray?
axes = axes.flatten()
for ax in axes:
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    storm_rpts_overlays[ax] = []# for centroid of reports. list of lists of multiple event type
    filldict[ax] = []
    linedict[ax] = []
    barbdict[ax] = []
    quiverdict[ax] = []
(northax, stormmotionax, blankax, windshearax) = axes
blankax.set_visible(False)
cbar_ax = figplr.add_axes([0.05, 0.32, 0.49, 0.015])
# Empty fineprint_string placeholder for fine print in lower left corner of image.
fineprint = plt.annotate(text="", xy=(4,1), xycoords=('figure pixels','figure pixels'), va="bottom", fontsize=3.2, wrap=True)


logging.info(f"ifiles {ifiles}")
stormlist = []
for ifile in ifiles:
    stormlist.extend(ifile.readlines())
for storm in stormlist:
    words = storm.split()
    stormname = words[0]
    year =  words[1]
    stormname_year = stormname + " " + year
    timelist = [pytz.utc.localize(datetime.datetime.strptime(x, '%Y%m%d%H')) for x in words[2:]]
    # Read TC tracks
    useIBTrACS = True
    useNHCbesttracks = not useIBTrACS
    if useNHCbesttracks:
        atcfname = ibtracs.get_atcfname_from_stormname_year(stormname, year, debug=debug) 
        best_track_file = atcf.archive_path(atcfname) # provide path to atcf archive
        df = atcf.read(best_track_file, debug=debug)
        df = atcf.interpolate(df, '3H', debug=debug) # NARR is every 3 hours; NHC best track is usually every 6.
    if useIBTrACS:
        # Use IBTrACS
        df, best_track_file = ibtracs.get_atcf(stormname, year, basin="NA") # North Atlantic basin guaranteed for NARR
        extension = ibtracs.extension(stormname, year)
        df = df.append(extension, sort=False, ignore_index=True)
        # Considered interpolating to 3H but IBTrACS already has data every 3H for most storms.
        # IBTRaCS is not guaranteed to have every 3H, only every 6H.
    df.valid_time = df.valid_time.apply(pytz.utc.localize) # make timezone-aware even without spc reports
    # ignore 50 and 64 knot rad lines. Keep 0 and 34-knot lines.Yes, there are 0-knot lines. rad is a string 
    df = df[df.rad.astype("float") <= 35] # Gabrielle 2001 atcf bal082001.dat has "35" knot winds, not 34. As a quick-fix, just use 35 as the comparison number.
    stormname_years.append(stormname_year) # used for composite title
    logging.info(f"storm {stormname_year} {df.valid_time.min()} to {df.valid_time.max()}")
    

    # Make sure all times in time list are within the TC track time window.
    for t in timelist:
        if t < df.valid_time.min() or t > df.valid_time.max():
            print("requested time",t,"is outside TC track time window. Exiting.")
            sys.exit(1)



    all_storm_reports = pd.DataFrame()
    if torn or wind or hail:
        all_storm_reports = spc.get_storm_reports(start=min(timelist)-spc_td, end=max(timelist)+spc_td)

    if not hail:
        all_storm_reports = all_storm_reports[all_storm_reports["event_type"] != "hail"]
        all_storm_reports = all_storm_reports[all_storm_reports["event_type"] != "large hail"]
    if not torn:
        all_storm_reports = all_storm_reports[all_storm_reports["event_type"] != "torn"]
        all_storm_reports = all_storm_reports[all_storm_reports["event_type"] != "sigtorn"]
    if not wind:
        all_storm_reports = all_storm_reports[all_storm_reports["event_type"] != "wind"]
        all_storm_reports = all_storm_reports[all_storm_reports["event_type"] != "high wind"]

    for index, (valid_time, lon1, lat1, storm_heading, storm_speed) in df[['valid_time','lon','lat','heading','speed']].iterrows():
        lat1 *= units.degrees_N
        lon1 *= units.degrees_E
        storm_heading *= units.degrees
        storm_speed =  storm_speed*units.parse_expression("knots").to("m/s")
        origin_place_time = f'({lat1.m:.1f}$\degree$N, {lon1.m:.1f}$\degree$E) {valid_time.strftime("%Y-%m-%d %H UTC")}'
        logging.debug(f'({lat1:.1f}, {lon1:.1f}) {valid_time.strftime("%Y-%m-%d %H UTC")}')

        # skip non-3hrly valid times
        if valid_time.hour % 3 != 0:
            logging.info(f"skipping non-3-hrly time {valid_time}")
            continue


        if valid_time not in timelist:
            logging.debug(f'TC valid time not requested. Skipping.')
            continue


        snapshot = os.path.dirname(ofile) + "/" + fill + "." + valid_time.strftime("%Y%m%d%H") + f".{lat1.m:04.1f}N{lon1.m:06.1f}E.png"
        supt = figplr.suptitle(stormname_year + "\n" + origin_place_time, wrap=True, fontsize="small")

        normalize_range_by = None
        if normalize_by:
            normalize_range_by = atcf.get_normalize_range_by(df, index, normalize_by)
            bname, ext = os.path.splitext(snapshot)
            snapshot = bname + ".normalize_by_" + normalize_by + ext
        # Used to skip this file If png already exists but we don't want to have incomplete composite.
        # Finish it if you start it. 
        if not clobber and os.path.isfile(snapshot):
            logging.info(f"overwriting existing file {snapshot}")

        # Define storm report time window. Use time range for legend title.
        # Extract storm reports within time window from all_storm_reports DataFrame.
        storm_report_start = valid_time - spc_td
        storm_report_end   = valid_time + spc_td
        storm_report_time_window_str = "storm reports\n"+storm_report_start.strftime('%m/%d %H%M') +' - '+ storm_report_end.strftime('%m/%d %H%M %Z')
        storm_report_window = (all_storm_reports.time >= storm_report_start) & (all_storm_reports.time < storm_report_end)
        storm_reports = all_storm_reports[storm_report_window]

        if storm_reports.empty:
            storm_report_time_window_str = "no "+storm_report_time_window_str

        data = narr.scalardata(fill, valid_time, targetdir=workdir)
        levels = data.attrs['levels']
        best_track_files.append(best_track_file)
        narr_files.append(data.attrs['ifile'])
        narr_valid_times.append(valid_time)
        lons.append(lon1)
        lats.append(lat1)
        cbar_title = "fill: "+ desc(data)
        cmap = data.attrs['cmap']
        cmap.set_under(color='white')
        if line:
            linecontourdata = narr.scalardata(line, valid_time, targetdir=workdir)
            # Add line contour description above color bar title
            cbar_title += "\nline: " + desc(linecontourdata)
        if barb:
            barbdata = narr.vectordata(barb, valid_time, targetdir=workdir)
            barbdata = barbdata.metpy.convert_units(barbunits)
            cbar_title += f"\nbarbs: {barb_increments} " + desc(barbdata) 
        if quiver:
            quiverdata = narr.vectordata(quiver, valid_time, targetdir=workdir)
            cbar_title += f"\nquiver: " + desc(quiverdata) 


        logging.debug(f"plotting filled contour {fill}...")

        # .load() to avoid Userwarning about passing obj to dask.array that is already a Dask collection.
        lon, lat = data.metpy.longitude.load(), data.metpy.latitude.load()
        dist_from_center, bearing = atcf.dist_bearing(lon1, lat1, lon, lat)
        if dist_from_center.min() > 32*units.km:
            logging.error(f"at {valid_time} {stormname_year} more than {dist_from_center.min()} from nearest NARR grid point")
            sys.exit(1)

        if cart or debug:
            logging.info("cartopy view for debugging...")
            fig = plt.figure(num=2, clear=True)
            logging.debug(f"fignums={plt.get_fignums()}")
            axc = plt.axes(projection=cartopy.crs.LambertConformal())
            axc.set_extent(extent, crs=cartopy.crs.PlateCarree()) 

            cfill = axc.contourf(lon,lat,uvsel(data),levels=levels,cmap=cmap,norm=colors.BoundaryNorm(levels,cmap.N),transform=cartopy.crs.PlateCarree())
            if line:
                line_contour = axc.contour(lon,lat,uvsel(linecontourdata),levels=contour_levels,cmap=linecontourdata.attrs["cmap"],norm=colors.BoundaryNorm(contour_levels,linecontourdata.attrs["cmap"].N),transform=cartopy.crs.PlateCarree())
                line_contour_labels = axc.clabel(line_contour, fontsize=linecontour_fontsize, fmt='%.0f')
            axc.set_title(stormname_year)

            # Color bar
            cb = plt.colorbar(cfill, ax=axc, orientation='horizontal', shrink=0.55)
            cb.ax.set_title(cbar_title, fontsize='xx-small')
            cb.ax.tick_params(labelsize='xx-small')
            cb.outline.set_linewidth(0.5) 

            c = axc.contour(lon, lat, dist_from_center, levels=np.arange(0,max(rbins)+200,200), colors='black', alpha=0.8, transform=cartopy.crs.PlateCarree())
            axc.clabel(c, fontsize='xx-small', fmt='%ikm')
            c = axc.contour(lon, lat, bearing, levels=range(0,360,45), colors='black', alpha=0.8, transform=cartopy.crs.PlateCarree())
            axc.clabel(c, fontsize='xx-small', fmt='%i')

            axc.set_title(axc.get_title() + "  " + valid_time.strftime('%Y-%m-%d %H UTC'), fontsize='x-small')

            if not storm_reports.empty:
                legend_items = spc.plot(storm_reports, axc, drawrange=0, debug=debug)
                axc.legend(handles=legend_items.values(), fontsize='xx-small').set_title(storm_report_time_window_str, prop={'size':'4.5'})

            # *must* call draw in order to get the axis boundary used to add ticks:
            fig.canvas.draw()

            # Define gridline locations and draw the lines using cartopy's built-in gridliner:
            xticks = list(range(-160,-50,10))
            yticks = list(range(0,65,5))
            axc.gridlines(xlocs=xticks, ylocs=yticks, linewidth=0.4, alpha=0.8, linestyle='--')
            axc.xaxis.set_major_formatter(LONGITUDE_FORMATTER) 
            axc.yaxis.set_major_formatter(LATITUDE_FORMATTER)
            gridlines.lambert_xticks(axc, xticks, fontsize='x-small')
            gridlines.lambert_yticks(axc, yticks, fontsize='x-small')
            axc.set_xlabel('')
            axc.set_ylabel('')
            axc.add_feature(cartopy.feature.STATES.with_scale('50m'), linewidth=0.35, alpha=0.55)
            axc.add_feature(cartopy.feature.COASTLINE.with_scale('50m'), linewidth=0.5, alpha=0.55)
            ocart = valid_time.strftime("cart.%Y%m%d%H.png")
            plt.savefig(ocart)
            logging.info(f'created {os.path.realpath(ocart)}')
            plt.figure(1) # switch back to polar plots figure

        if normalize_by:
            print(f"normalizing range from TC center by {normalize_by}. {normalize_range_by:.3f}")
            dist_from_center = dist_from_center / normalize_range_by
            supt = supt.set_text(supt.get_text() + f"  1 km here = {1*units.km * normalize_range_by:.3f}")
        
        logging.debug("creating weighted and unweighted 2d histograms of theta and range...")

        for ax in axes:
            ax.grid(False) #Auto-removal of grids by pcolor() and pcolormesh() is deprecated since 3.5 and will be removed two minor releases later; please call grid(False)

        values = histogram2d_weighted(bearing, dist_from_center, azbins, rbins, data) # Don't apply uvsel() here; we need values for storm motion and wind shear
        polarf = northax.pcolor(np.radians(theta2d), r2d, uvsel(values), cmap=cmap,norm=colors.BoundaryNorm(levels,cmap.N)) # tried DataArray.plot.imshow and pcolormesh, tried 1D coordinates.
        filldict[northax].append(uvsel(values))
        if line: 
            linecontour_values = histogram2d_weighted(bearing, dist_from_center, azbins, rbins, linecontourdata) # Don't apply uvsel() here; we need values for storm motion and wind shear
            polarc = northax.contour(np.radians(theta_center1D), range_center1D, uvsel(linecontour_values).T, contour_levels, colors="black", alpha=linecontour_alpha)
            northax.clabel(polarc, fontsize=linecontour_fontsize, fmt='%.0f')
            linedict[northax].append(uvsel(linecontour_values))
        if barb:
            barb_values = histogram2d_weighted(bearing, dist_from_center, azbins, rbins, barbdata)
            polarb = northax.barbs(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *barb_values.stack(gate=["azimuth","range"]).isel(gate=ii), **barbkwdict)
            barbdict[northax].append(barb_values) # TODO: do we need uvsel() for barbs?
        if quiver:
            quiver_values = histogram2d_weighted(bearing, dist_from_center, azbins, rbins, quiverdata)
            polarq = northax.quiver(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *quiver_values.stack(gate=["azimuth","range"]).isel(gate=ii))
            quiverdict[northax].append(quiver_values)
            powerof10 = 10**int(np.log10(np.sqrt(quiverdata[0]**2 + quiverdata[1]**2).max().values)) # order of magnitude  
            northax.quiverkey(polarq, 0.1, -0.05, powerof10, f"{powerof10} {quiver_values.metpy.units:~}", labelpos='S', fontproperties=dict(size='xx-small'))


        # Rotate so Storm motion vector points up
        storm_motion_values = roll_rotate_sel(values, storm_heading)

        polarf2 = stormmotionax.pcolor(np.radians(theta2d), r2d, storm_motion_values, cmap=cmap,norm=colors.BoundaryNorm(levels,cmap.N))
        filldict[stormmotionax].append(storm_motion_values)
        if line:
            linecontour_storm_motion_values = roll_rotate_sel(linecontour_values, storm_heading)
            polarc = stormmotionax.contour(np.radians(theta_center1D), range_center1D, linecontour_storm_motion_values.T, contour_levels, colors="black", alpha=linecontour_alpha)
            stormmotionax.clabel(polarc, fontsize=linecontour_fontsize, fmt='%.0f')
            linedict[stormmotionax].append(linecontour_storm_motion_values)
        if barb:
            storm_motion_barb_values = roll_rotate_sel(barb_values, storm_heading) 
            polarb = stormmotionax.barbs(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *storm_motion_barb_values.stack(gate=["azimuth","range"]).isel(gate=ii), **barbkwdict)
            barbdict[stormmotionax].append(storm_motion_barb_values)
        if quiver:
            storm_motion_quiver_values = roll_rotate_sel(quiver_values, storm_heading)
            polarq = stormmotionax.quiver(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *storm_motion_quiver_values.stack(gate=["azimuth","range"]).isel(gate=ii))
            quiverdict[stormmotionax].append(storm_motion_quiver_values)

           

        # Rotate so Wind shear vector points up
        TCFlow = ncl.steering_flow(lat0=lat1,lon0=lon1,storm_heading=storm_heading,storm_speed=storm_speed,
                stormname=stormname, rx=4.5, ensmember=stormname, ptop=ptop, pbot=pbot, 
                file_ncl=narr.get(valid_time, narrtype=narr.narr3D, targetdir=workdir), debug=debug)
        wind_shear_heading = TCFlow["wind_shear_heading"] * units.degrees
        wind_shear_mag   = TCFlow["wind_shear_speed"] * units.parse_expression("m/s")
        wind_shear_values = roll_rotate_sel(values, wind_shear_heading) 
        polarf3 = windshearax.pcolor(np.radians(theta2d), r2d, wind_shear_values, cmap=cmap,norm=colors.BoundaryNorm(levels,cmap.N))
        filldict[windshearax].append(wind_shear_values)
        if line:
            linecontour_wind_shear_values = roll_rotate_sel(linecontour_values, wind_shear_heading)
            polarc = windshearax.contour(np.radians(theta_center1D), range_center1D, linecontour_wind_shear_values.T, contour_levels, colors="black", alpha=linecontour_alpha)
            windshearax.clabel(polarc, fontsize=linecontour_fontsize, fmt='%.0f')
            linedict[windshearax].append(linecontour_wind_shear_values)
        if barb:
            wind_shear_barb_values = roll_rotate_sel(barb_values, wind_shear_heading)
            polarb = windshearax.barbs(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *wind_shear_barb_values.stack(gate=["azimuth","range"]).isel(gate=ii), **barbkwdict)
            barbdict[windshearax].append(wind_shear_barb_values)
        if quiver:
            wind_shear_quiver_values = roll_rotate_sel(quiver_values, wind_shear_heading)
            polarb = windshearax.quiver(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *wind_shear_quiver_values.stack(gate=["azimuth","range"]).isel(gate=ii))
            quiverdict[windshearax].append(wind_shear_quiver_values)
            
            

        for ax in axes:
            lines, labels = add_rgrid(ax)

        # Add storm reports. save reports from stormmitionax and windshearax in addition to northax so we can clear axes but replot at the end.
        if not storm_reports.empty:
            logging.info(f"found {len(storm_reports)} storm reports. {storm_reports['significant'].sum()} significant.")
            rptkw = dict(normalize_range_by=normalize_range_by, debug=debug)
            storm_rpts_overlays[northax].append(      spc.polarplot(lon1, lat1, storm_reports, northax,       zero_azimuth=0*units.deg,        **rptkw) )
            storm_rpts_overlays[stormmotionax].append(spc.polarplot(lon1, lat1, storm_reports, stormmotionax, zero_azimuth=storm_heading,      **rptkw) )
            storm_rpts_overlays[windshearax].append(  spc.polarplot(lon1, lat1, storm_reports, windshearax,   zero_azimuth=wind_shear_heading, **rptkw) )
            # legend below northax 
            storm_report_legend = northax.legend(handles=storm_rpts_overlays[northax][-1].values(), **storm_rpt_legend_kw, title=storm_report_time_window_str)
            nhandles = len(storm_report_legend.legendHandles)
            logging.info(f"found {nhandles} storm report legend handles")

        _, _ = northax.set_thetagrids(theta_lines, labels=(f"N","","E","","S","","W",""), fontweight="demi", fontsize="large")
        _, _ = stormmotionax.set_thetagrids(theta_lines, labels=(f"storm heading {storm_heading.m:.0f}$\degree$ at {storm_speed:~.1f}\nfront","","R","","rear","","L",""), fontweight="demi")
        _, _ = windshearax.set_thetagrids(theta_lines, labels=(f"wind shear heading {wind_shear_heading.m:.0f}$\degree$ at {wind_shear_mag:~.1f}\ndownshear","","R","","upshear","","L",""), fontweight="demi")

        # Shared colorbar
        if polarf.cmap != polarf2.cmap or polarf.cmap != polarf3.cmap: # Make sure colorbar applies to all subplots
            print("colormap changed between subplots")
            sys.exit(1)
        cb2 = figplr.colorbar(polarf, cax=cbar_ax, orientation='horizontal')
        cb2.set_label(data.metpy.units, fontsize='xx-small')
        cb2.ax.set_title(cbar_title, fontsize='xx-small')
        cb2.ax.tick_params(labelsize='xx-small')
        cb2.set_ticks(levels)
        cb2.outline.set_linewidth(0.35) # TODO: Does this change anything?

        fineprint.set_text(fineprint_string + f"\n{data.attrs['ifile']}\ncreated {datetime.datetime.now()}")
        if no_fineprint: # hide fineprint
            fineprint.set_visible(False)

        # Save image. 
        plt.savefig(snapshot, dpi=220)
        os.system("mogrify -colors 128 "+snapshot) # severe drop in quality with 64 colors
        print('created ' + os.path.realpath(snapshot))

        # Clear all axes. We have lists of storm report collections for later.
        for ax in axes:
            ax.clear()
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)

# Composite

if len(narr_files) == 0:
    print("no NARR files to composite")
    sys.exit(0)

best_track_and_narr_files = [" ".join(x) for x in zip(best_track_files, narr_files)] # create list of best track and narr files paired with " " space
# Don't list all narr ifiles if there are a lot. show count.
if len(best_track_and_narr_files) > 15:
    best_track_and_narr_files_list = str(len(best_track_and_narr_files))+ " narr files"
else:
    best_track_and_narr_files_list = "\n".join(best_track_and_narr_files)
if debug:
    print("best_track_and_narr_files",best_track_and_narr_files)
fineprint.set_text(fineprint_string + "\n" + best_track_and_narr_files_list + f"\ncreated {datetime.datetime.now()}")

# Create hours UTC string for figure title - used to show min-max, but with a plethora of years that didn't make sense. 
fmt = '%H %Z'
hours_str = (x.strftime(fmt) for x in narr_valid_times)
# First remove duplicates with set(). Then convert to list, and sort.
hours_str = " & ".join(sorted(list(set(hours_str))))
# As we move to 6-hour blocks we have more than one time from the same storm.
# Use set of stormname_years to eliminate duplicates. This destroys the order. 
uniq_stormname_years = list(set(stormname_years))
# sort by year first, then name 
sorted_storms = sorted(uniq_stormname_years, key=lambda x: tuple(x.split()[::-1]))
fontsize = "x-small" if len(stormname_years) < 20 else "xx-small"
figplr.suptitle(", ".join(sorted_storms) + "\n" + hours_str, fontsize=fontsize) # Clean title w/o origin place and time

# avoid spiral artifact when pcolor ignores entire columns of nan. pcolor.get_array() does not include nans in PolyCollection.
# shape of array ends up changing after pcolor is called if one column of values array is all nan.
assert polarf.get_array().size == np.mean(filldict[northax], axis=0).size, "filled contour field changed size. Check for nans"
for ax in [northax, stormmotionax, windshearax]:
    ax.grid(False)
    ax.pcolor(np.radians(theta2d), r2d, xr_list_mean(filldict[ax]), cmap=cmap,norm=colors.BoundaryNorm(levels,cmap.N))
    lines, labels = add_rgrid(ax)
    x = xr_list_mean(filldict[ax])
    print(f"heading and range of maximum {fill}:", end="")
    print(x.isel(x.argmax(...)))

    if line:
        pc = ax.contour(np.radians(theta_center1D), range_center1D, (np.mean(linedict[ax], axis=0)).T, contour_levels, colors="black", alpha=linecontour_alpha)
        ax.clabel(pc, fontsize=linecontour_fontsize, fmt='%.0f')
    if barb:
        ax.barbs(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *xr_list_mean(barbdict[ax]).stack(gate=["azimuth","range"]).isel(gate=ii), **barbkwdict)
    if quiver:
        polarq = ax.quiver(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *xr_list_mean(quiverdict[ax]).stack(gate=["azimuth","range"]).isel(gate=ii)) 

    # Composite storm reports
    # Concatenate values of same keys in a list of dicts.
    merged_dict = {k:[d.get(k) for d in storm_rpts_overlays[ax] if k in d] for k in set().union(*storm_rpts_overlays[ax])}
    # get centroid of each type of storm report
    handles = [] # for legend
    labels = [] # for legend
    for event_type in merged_dict:
        pcs = [ax.add_artist(pc) for pc in merged_dict[event_type]]
        lsr = [x.get_offsets().data for x in merged_dict[event_type]]
        lsr = np.concatenate(lsr)
        theta_deg = np.degrees(lsr[:,0]) * units.deg # maybe assign units back in .get_offsets().data ?
        range_km  = lsr[:,1] * units.km
        handles.append(pcs[0]) # only use first PathCollection element for legend
        labels.append(f"{event_type} ({len(lsr)})") # but append total number of events to legend label
        print(f"centroid of {ax} coord {event_type} reports {spc.centroid_polar(theta_deg, range_km)}")
    if ax == northax:
        storm_report_legend = ax.legend(handles=handles, labels=labels, title='storm reports', **storm_rpt_legend_kw)
        if quiver:
            powerof10 =10**int(np.log10(np.sqrt(quiverdata[0]**2 + quiverdata[1]**2).max().values))
            northax.quiverkey(polarq, 0.1, -0.03, powerof10, f"{powerof10} {quiver_values.metpy.units:~}", labelpos='S', fontproperties=dict(size='xx-small'))

# Tried popping top Text major ticklabel object from thetalabels, altering it, and reinserting into thetalabels, but had no effect.
_, _ = northax.set_thetagrids(theta_lines, labels=(f"N","","E","","S","","W",""), fontweight="demi", fontsize="large")
_, _ = stormmotionax.set_thetagrids(theta_lines, labels=("front","","R","","rear","","L",""), fontweight="demi", fontsize=11)
_, _ = windshearax.set_thetagrids(theta_lines, labels=("downshear","","R","","upshear","","L",""), fontweight="demi", fontsize=11)


# Save image. 
plt.savefig(ofile, dpi=220) # use dpi=200 for Stan
os.system("mogrify -colors 128 "+ofile)
print('created ' + os.path.realpath(ofile))

if netcdf:
    coord_stack = np.stack([filldict[northax], filldict[stormmotionax], filldict[windshearax]])
    # Tried to Quantify coord_stack with data.metpy.units so netCDF will have units attribute, but .to_netcdf can't handle Quantity object.
    ds = xarray.Dataset(
            {
                fill               : (('coord', 'storm', 'azimuth', 'range'),  coord_stack), 
                "best_track_files" : (('storm'), best_track_files),
                "narr_files"       : (('storm'), narr_files),
                "lon"              : (('storm'), lons),
                "lat"              : (('storm'), lats),
                "time"             : (('storm'), narr_valid_times)
            },
            coords={'coord' : ["north","storm motion","wind shear"], 'range' : range_center1D, 'azimuth' : theta_center1D, 'storm':stormname_years}
        )
    # TODO: add line contour and wind barb components, if they exist.
    attrs = data.metpy.dequantify().attrs # dequantify returns variables cast to magnitude and units on attribute. Or else units are lost here. 
    # "cmap" doesn't translate well to netCDF, "ifile" and "initial_time" populated by last NARR file, not list of all NARR files, as it should be
    # "arraylevel" should work but sometimes it is data type O.
    to_del = ["levels", "cmap", "ifile", "initial_time", "arraylevel"]
    for d in to_del:
        if d in attrs:
            del attrs[d]
    # 'vertical' may be a list of pint Quantities. pdb.set_trace()
    if 'vertical' in attrs:
        try:
            if not isinstance(attrs["vertical"], str):
                attrs["vertical"] = [str(x) for x in attrs["vertical"]] # convert list of pint Quantities to list of str or get ValueError: setting an array element with a sequence.
        except:
            attrs["vertical"] = str(attrs["vertical"]) # convert pint Quantity to str or get TypeError: Invalid value for attr ...
    ds[fill].attrs = attrs
    ds['range'].attrs = {"units": "km", "long_name":"range from origin"}
    ds['azimuth'].attrs = {"units":"degrees","long_name":"degrees clockwise from north"}
    ds.attrs['daz'] = daz
    ds.attrs['dr'] = dr
    ds.attrs['ifiles'] = [os.path.realpath(ifile.name) for ifile in ifiles]
    ds["time"].values = ds.time.astype(np.datetime64)
    ds["time"].encoding['units'] = "hours since 1970-01-01"
    ds.to_netcdf(netcdf)
    print("wrote",netcdf)

