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
from scipy.stats import gaussian_kde
import spc
import sys
import xarray

def m(x, w):
    """Weighted Mean"""
    return np.sum(x * w) / np.sum(w)

def cov(x, y, w):
    """Weighted Covariance"""
    return np.sum(w * (x - m(x, w)) * (y - m(y, w))) / np.sum(w)

def corr(x, y, w):
    """Weighted Correlation"""
    x = x.values
    y = y.values
    return cov(x, y, w) / np.sqrt(cov(x, x, w) * cov(y, y, w))

def desc(data):
    assert 'timetitle' in data.attrs, "timetitle not attribute of data"+str(data)
    assert 'verttitle' in data.attrs
    assert 'long_name' in data.attrs
    timetitle = data.attrs['timetitle']
    verttitle = str(data.attrs['verttitle'])
    desc = f"{timetitle} {verttitle} {data.long_name} [{data.metpy.units:~}]" # abbreviate units.
    return desc

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

def roll_rotate_sel(values, heading):
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
parser.add_argument("--dpi", type=int, default=250, help="dpi of image output")
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
parser.add_argument("ifiles", nargs="+", type=argparse.FileType('r'), help="text file(s). each line starts with a storm name, followed by a yyyymmddhh time")
parser.add_argument("--torn", action='store_true', help="overlay tornado reports")
parser.set_defaults(torn=True)
parser.add_argument("--wind", action='store_true', help="overlay wind reports")
parser.add_argument('-w', "--workdir", type=str, default="/glade/scratch/"+os.getenv("USER")+"/NARR", help="directory to untar NARR into")


# Assign arguments to simple-named variables
args = parser.parse_args()
barb            = args.barb
cart            = args.cart
clobber         = args.clobber
lcontour_levels = args.clevels
debug           = args.debug
dpi             = args.dpi
extent          = args.extent
fill            = args.fill
hail            = args.hail
line            = args.line
max_range       = args.max_range
netcdf          = args.netcdf
no_fineprint    = args.no_fineprint
normalize_by    = args.normalize_by
quiver          = args.quiver
spc_td          = datetime.timedelta(hours=args.spctd)
ifiles          = args.ifiles
torn            = args.torn
wind            = args.wind
workdir         = args.workdir

logger = logging.getLogger()
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
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
azbins = np.arange(0,360+daz,daz)
theta_lines = range(0,360,45)
linecontour_alpha = .6
linecontour_fontsize = 4
barblength = 3.5
barbunits = 'm/s'
barb_increments = {'half':2.5, 'full':5, 'flag':25} # for matplotlib.quiver.barb
barbkwdict = dict(length=barblength, linewidth=0.2, alpha=0.6, barb_increments=barb_increments)

rbins  = np.arange(0,max_range+dr,dr)
r2d, theta2d = np.meshgrid(rbins, azbins)
binarea = np.pi * (2*r2d + dr) * dr * daz/360.
range_center1D, theta_center1D = dr/2+rbins[:-1], daz/2+azbins[:-1]
r_center2D, theta_center2D = np.meshgrid(range_center1D, theta_center1D) # one less element in range and azimuth than range_center1D
polar_sz = (len(azbins)-1, len(rbins)-1)

#---Locate gridded wind barbs on polar plot
dx = 100 # wind barb grid spacing
mg = np.arange(-max_range, max_range+dx, dx)
x0, y0 = np.meshgrid(mg, mg)
xp = r_center2D * np.cos(np.radians(90-theta_center2D))
yp = r_center2D * np.sin(np.radians(90-theta_center2D))
tree = spatial.KDTree(np.c_[xp.ravel(),yp.ravel()])
dd, ii = tree.query(np.c_[x0.ravel(),y0.ravel()], k=1, distance_upper_bound=dx)
ii = ii[dd < dx/2] # remove neighbors over dx/2 away from nearest point


best_track_files, narr_files = [], [] # used for netCDF file
lons, lats, times = [], [], [] # used for netCDF file
stormname_years, narr_valid_times = [], [] # for composite title

figplr, axes = plt.subplots(ncols=2,nrows=2,subplot_kw=dict(projection='polar'))
# default left=0.125 right=0.9 bottom=0.1 top=0.9 wspace=0.2 hspace=0.2
plt.subplots_adjust(left=0, right=0.99, bottom=0.06, top=0.82, hspace=0.43, wspace=0)

# These dicts will have one entry for each coordinate system (axes).
filldict   = {}
linedict   = {}
barbdict   = {}
quiverdict = {}
storm_rpt_dict = {}
# any reason to keep as 2x2 ndarray?
axes = axes.flatten()
for ax in axes:
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    storm_rpt_dict[ax] = pd.DataFrame(columns=["event_type"])# for centroid of reports. 
    filldict[ax] = []
    linedict[ax] = []
    barbdict[ax] = []
    quiverdict[ax] = []
(northax, stormmotionax, blankax, windshearax) = axes
blankax.set_visible(False)
cbar_ax = figplr.add_axes([0.05, 0.32, 0.45, 0.015])
# Empty fineprint_string placeholder for fine print in lower left corner of image.
fineprint = plt.annotate(text="", xy=(4,1), xycoords=('figure pixels','figure pixels'), va="bottom", fontsize=3.2, wrap=True)


logging.info(f"ifiles {ifiles}")
stormlist = []
for ifile in ifiles:
    stormlist.extend(ifile.readlines())
for storm in stormlist:
    stormname, narrtime = storm.split()
    narrtime = pytz.utc.localize(datetime.datetime.strptime(narrtime, '%Y%m%d%H'))
    year =  str(narrtime.year)
    stormname_year = stormname + " " + year
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
    if narrtime < df.valid_time.min() or narrtime > df.valid_time.max():
        logging.error(f"requested time {narrtime} is outside TC track time window. Exiting.")
        sys.exit(1)



    all_storm_reports = pd.DataFrame()
    if torn or wind or hail:
        all_storm_reports = spc.get_storm_reports(start=narrtime-spc_td, end=narrtime+spc_td)
        # Until I figure out how to combine significant and non-significant for centroid and kde while leaving separate legend handles, combine everything now.
        all_storm_reports = spc.combine_significant(all_storm_reports)

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


        if valid_time !=  narrtime:
            logging.debug(f'TC valid time {valid_time} not requested. Skipping.')
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
        storm_reports = all_storm_reports.loc[storm_report_window].copy() # avoid SettingWithCopyWarning in spc.polarkde by applying .copy()

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
                line_contour = axc.contour(lon,lat,uvsel(linecontourdata),levels=lcontour_levels,cmap=linecontourdata.attrs["cmap"],norm=colors.BoundaryNorm(lcontour_levels,linecontourdata.attrs["cmap"].N),transform=cartopy.crs.PlateCarree())
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
                legend_items = spc.plot(storm_reports, axc, drawrange=0)
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

        values = spc.histogram2d_weighted(bearing, dist_from_center, azbins, rbins, data) # Don't apply uvsel() here; we need values for storm motion and wind shear

        TCFlow = ncl.steering_flow(lat0=lat1,lon0=lon1,storm_heading=storm_heading,storm_speed=storm_speed,
                stormname=stormname, rx=4.5, ensmember=stormname, ptop=ptop, pbot=pbot, 
                file_ncl=narr.get(valid_time, narrtype=narr.narr3D, targetdir=workdir), debug=debug)
        wind_shear_heading = TCFlow["wind_shear_heading"] * units.degrees
        wind_shear_mag   = TCFlow["wind_shear_speed"] * units.parse_expression("m/s")
        wind_shear_values = roll_rotate_sel(values, wind_shear_heading) 

        ax_headings = {
                northax      :0 * units.deg, # north points up
                stormmotionax: storm_heading, # Rotate so Storm motion vector points up
                windshearax  : wind_shear_heading # Rotate so Wind shear vector points up
                }

        ax_labels = {
                northax : ("N","","E","","S","","W",""),
                stormmotionax: (f"storm heading {storm_heading.m:.0f}$\degree$ at {storm_speed:~.1f}\nfront","","R","","rear","","L",""),
                windshearax : (f"wind shear heading {wind_shear_heading.m:.0f}$\degree$ at {wind_shear_mag:~.1f}\ndownshear","","R","","upshear","","L","")
                }

        for ax,axheading in ax_headings.items():
            ax.grid(False) #Auto-removal of grids by pcolor() and pcolormesh() is deprecated since 3.5 and will be removed two minor releases later; please call grid(False)

            rotated_values = roll_rotate_sel(values, axheading)
            polarf = ax.pcolor(np.radians(theta2d), r2d, uvsel(rotated_values), cmap=cmap,norm=colors.BoundaryNorm(levels,cmap.N)) # tried DataArray.plot.imshow and pcolormesh, tried 1D coordinates.
            filldict[ax].append(uvsel(rotated_values))
            if line: 
                linecontour_values = spc.histogram2d_weighted(bearing, dist_from_center, azbins, rbins, linecontourdata) # Don't apply uvsel() here; we need values for storm motion and wind shear
                linecontour_values = roll_rotate_sel(linecontour_values, axheading)
                polarc = ax.contour(np.radians(theta_center1D), range_center1D, uvsel(linecontour_values).T, lcontour_levels, colors="black", alpha=linecontour_alpha)
                ax.clabel(polarc, fontsize=linecontour_fontsize, fmt='%.0f')
                linedict[ax].append(uvsel(linecontour_values))
            if barb:
                barb_values = spc.histogram2d_weighted(bearing, dist_from_center, azbins, rbins, barbdata)
                barb_values = roll_rotate_sel(barb_values, axheading) 
                polarb = ax.barbs(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *barb_values.stack(gate=["azimuth","range"]).isel(gate=ii), **barbkwdict)
                barbdict[ax].append(barb_values) 
            if quiver:
                quiver_values = spc.histogram2d_weighted(bearing, dist_from_center, azbins, rbins, quiverdata)
                quiver_values = roll_rotate_sel(quiver_values, axheading)
                polarq = ax.quiver(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *quiver_values.stack(gate=["azimuth","range"]).isel(gate=ii))
                quiverdict[ax].append(quiver_values)
                powerof10 = 10**int(np.log10(np.sqrt(quiverdata[0]**2 + quiverdata[1]**2).max().values)) # order of magnitude  
                ax.quiverkey(polarq, 0.1, -0.05, powerof10, f"{powerof10} {quiver_values.metpy.units:~}", labelpos='S', fontproperties=dict(size='xx-small'))
            lines, labels = add_rgrid(ax)

            # Plot storm reports. Append reports to DataFrame for stormmitionax and windshearax in addition to northax so we can clear axes but replot at the end.
            if not storm_reports.empty:
                rptkw = dict(normalize_range_by=normalize_range_by)
                # Return DataFrame filtered for max_range with additional "range" and "heading" columns.
                storm_reports = spc.polarplot(lon1, lat1, storm_reports, ax, zero_azimuth=axheading, add_legend = ax == northax, legend_title=storm_report_time_window_str, **rptkw)
                logging.info(f"found {len(storm_reports)} storm reports within {max_range}. {storm_reports['significant'].sum()} significant.")
                storm_rpt_dict[ax] = pd.concat([storm_rpt_dict[ax], storm_reports])
                logging.info(f"{len(storm_rpt_dict[ax])} total")
                showkde = True
                if showkde:
                    w = r_center2D**2 # weight by area (range squared)
                    rptkde = spc.polarkde(lon1, lat1, storm_reports, ax, azbins, rbins, spc_td, zero_azimuth=axheading, ds=10*units.km, add_colorbar=False, **rptkw) #TODO: avoid additional colorbars pushing axes inward. Can't figure out how to remove them.
                    for event_type in rptkde:
                        rho = corr(uvsel(rotated_values), rptkde[event_type], w)
                        logging.info(f"{fill} {event_type} r={rho}")
                        ax.set_title(f"{desc(data)} {event_type} r={rho:.3f}", fontsize="xx-small")

            fontsize = "large" if ax == northax else "x-small"
            _, _ = ax.set_thetagrids(theta_lines, labels=ax_labels[ax], fontweight="demi", fontsize=fontsize)


        # Shared colorbar
        cb = figplr.colorbar(polarf, cax=cbar_ax, orientation='horizontal')
        cb.set_label(data.metpy.units, fontsize='xx-small')
        cb.ax.set_title(cbar_title, fontsize='xx-small')
        cb.ax.tick_params(labelsize='xx-small')
        cb.set_ticks(levels)
        cb.outline.set_linewidth(0.35) # TODO: Does this change anything?

        fineprint.set_text(fineprint_string + f"\n{best_track_file} {data.attrs['ifile']}\ncreated {datetime.datetime.now()}")
        if no_fineprint: # hide fineprint
            fineprint.set_visible(False)

        # Save image. 
        plt.savefig(snapshot, dpi=dpi)
        os.system("mogrify -colors 128 "+snapshot) # severe drop in quality with 64 colors
        logging.info(f'created {os.path.realpath(snapshot)}')

        # Clear all axes. We have lists of storm report collections for later.
        for ax in axes:
            # Kludge to remove kde colorbars
            #[x.remove() for x in figplr.axes if x.get_label() == '<colorbar>']
            ax.clear()
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)

# Composite

if len(narr_files) == 0:
    print("no NARR files to composite")
    sys.exit(0)

logging.debug(f"{len(stormlist)} storm and times {stormlist}")
ncols = 5 
stormlist_fineprint = f"{len(stormlist)} times\n"
if len(stormlist) < 120:
    for i,s in enumerate(stormlist):
        stormlist_fineprint += f" {s.rstrip().center(25)}"
        if i % ncols == ncols-1:
            stormlist_fineprint += "\n"
fineprint.set_text(fineprint_string + "\n" + stormlist_fineprint.rstrip() + f"\ncreated {datetime.datetime.now()}")

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
fontsize = "xx-small" if len(stormname_years) < 20 else 5
figplr.suptitle(", ".join(sorted_storms) + "\n" + hours_str, fontsize=fontsize) # Clean title w/o origin place and time

# avoid spiral artifact when pcolor ignores entire columns of nan. pcolor.get_array() does not include nans in PolyCollection.
# shape of array ends up changing after pcolor is called if one column of values array is all nan.
assert polarf.get_array().size == np.mean(filldict[northax], axis=0).size, "filled contour field changed size. Check for nans"
for ax in [northax, stormmotionax, windshearax]:
    ax.grid(False)
    fillmean = xr_list_mean(filldict[ax]) # used for pcolor and correlation
    ax.pcolor(np.radians(theta2d), r2d, fillmean, cmap=cmap,norm=colors.BoundaryNorm(levels,cmap.N))
    lines, labels = add_rgrid(ax)
    # don't use logging here. don't want timestamp
    print(f"heading and range of maximum {fill}:")
    print(fillmean.isel(fillmean.argmax(...)))

    if line:
        pc = ax.contour(np.radians(theta_center1D), range_center1D, (np.mean(linedict[ax], axis=0)).T, lcontour_levels, colors="black", alpha=linecontour_alpha)
        ax.clabel(pc, fontsize=linecontour_fontsize, fmt='%.0f')
    if barb:
        ax.barbs(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *xr_list_mean(barbdict[ax]).stack(gate=["azimuth","range"]).isel(gate=ii), **barbkwdict)
    if quiver:
        # save polarq for quiverkey
        polarq = ax.quiver(np.radians(theta_center2D.ravel()[ii]), r_center2D.ravel()[ii], *xr_list_mean(quiverdict[ax]).stack(gate=["azimuth","range"]).isel(gate=ii)) 

    # Composite storm reports
    # get centroid of each type of storm report
    handles = [] # for legend
    labels = [] # for legend
    kwdict = spc.symbol_dict()
    for event_type, xrpts in storm_rpt_dict[ax].groupby("event_type"):
        show_storm_rpt_scatter = False
        if show_storm_rpt_scatter:
            pc = ax.scatter(np.radians(xrpts["heading"]), xrpts["range"], alpha=0.9, **kwdict[event_type])
            handles.append(pc) # PathCollection element for legend
        lsr_heading = xrpts["heading"].values * units.deg 
        lsr_range  = xrpts["range"].values * units.km


        nrpts = len(xrpts)
        labels.append(f"{event_type} ({nrpts})") # append total number of events to legend label
        logging.info(f"centroid of {ax} coord {nrpts} {event_type} reports\n{spc.centroid_polar(lsr_heading, lsr_range)}")
        if nrpts < 3:
            logging.info(f"{nrpts} {event_type} reports not enough for kde (require at least 3)")
            continue
        # Try to plot kde
        rptx = lsr_range * np.cos(lsr_heading)
        rpty = lsr_range * np.sin(lsr_heading)
        ds = 2
        X, Y = np.mgrid[-max_range:max_range+ds:ds, -max_range:max_range+ds:ds] 
        az = pd.DataFrame(np.degrees(np.arctan2(Y,X))) # Make DataFrame because histogram2d_weighted expects a .values attribute
        az[az < 0] = az[az < 0] + 360 # azbins is 0-360 but az was -180 to +180
        dist_from_center = pd.DataFrame(np.sqrt(X**2 + Y**2))
        positions = np.vstack([X.ravel(), Y.ravel()])
        weights=None
        if event_type == "torn":
            weights = xrpts["mag"].fillna(value=0) + 1 # F-scale plus one. (F-sum from McCaul, 1991) # unknown tornado F-scales exist.
            event_type = f"torn F-sum"
        kernel = gaussian_kde(np.vstack([rptx,rpty]), weights=weights)
        Z = np.reshape(kernel(positions).T, X.shape)
        time_window_sum = len(stormlist) * 2 * spc_td.total_seconds()*units.seconds
        time_window_sum = time_window_sum.to("days")
        Z = xarray.DataArray(nrpts * Z  / lsr_range.units**2 / time_window_sum) # histogram2d_weighted expects a DataArray with units.
        # convert Cartesian grid Z to polar coordinates
        pkde = spc.histogram2d_weighted(az, dist_from_center, azbins, rbins, Z)
        pkde["azimuth"] = np.radians(pkde.azimuth) # Polar Axes are in radians not degrees.
        kdelevels = [0.00002, 0.00004, 0.00008, 0.00016, 0.00032, 0.00064]
        if Z.max() > kdelevels[0] * units.parse_expression("1/km**2/day"): # avoid ZeroDivisionError with add_colorbar
            logging.info(f"max density {Z.max()} plotting kde")
            polarc = pkde.plot.contour(x="azimuth",y="range", ax=ax, levels=kdelevels, colors="black", 
                    linewidths=0.5, add_colorbar=True, cbar_kwargs={"shrink":0.75,"pad":0.09})

            cb = polarc.colorbar
            cb.formatter.set_powerlimits((-2,2))
            cb.update_ticks()
            cb.ax.yaxis.offsetText.set(size='xx-small')
            cb.ax.yaxis.offsetText.set_horizontalalignment("left")
            cb.set_label(event_type + " " + cb.ax.yaxis.get_label().get_text(),fontsize="xx-small")
            cb.ax.tick_params(labelsize='xx-small')
        else:
            logging.info(f"max density {Z.max()} <= first kde level {kdelevels[0]}. skipping kdeplot")
    
        ax.set_xlabel('')
        ax.set_ylabel('')
        w = r_center2D**2 # weight by area (range squared)
        rho = corr(fillmean, pkde, w)
        rho_limit = 600*units.km
        w[r_center2D*units.km >= rho_limit] = 0.
        rho_limited = corr(fillmean, pkde, w) 
        logging.info(f"{fill} {nrpts} {event_type} r={rho:.3f} r{rho_limit:~}={rho_limited:.3f}")
        ax.set_title(f"{desc(data)} {nrpts} {event_type} r={rho:.3f} r{rho_limit:~}={rho_limited:.3f}", fontsize=4, wrap=True)
    if ax == northax: # add storm report legend and quiver scale to north axes
        if quiver:
            powerof10 =10**int(np.log10(np.sqrt(quiverdata[0]**2 + quiverdata[1]**2).max().values))
            ax.quiverkey(polarq, 0.1, -0.03, powerof10, f"{powerof10} {quiver_values.metpy.units:~}", labelpos='S', fontproperties=dict(size='xx-small'))

# Tried popping top Text major ticklabel object from thetalabels, altering it, and reinserting into thetalabels, but had no effect.
_, _ = northax.set_thetagrids(theta_lines, labels=(f"N","","E","","S","","W",""), fontweight="demi", fontsize="large")
_, _ = stormmotionax.set_thetagrids(theta_lines, labels=("front","","R","","rear","","L",""), fontweight="demi", fontsize=11)
_, _ = windshearax.set_thetagrids(theta_lines, labels=("downshear","","R","","upshear","","L",""), fontweight="demi", fontsize=11)


# Save image. 
plt.savefig(ofile, dpi=dpi) 
os.system("mogrify -colors 128 "+ofile)
logging.info(f'created {os.path.realpath(ofile)}')

def coord_stack(d):
    # Tried to Quantify coord_stack with data.metpy.units so netCDF will have units attribute, but .to_netcdf can't handle Quantity object.
    return np.stack([d[northax], d[stormmotionax], d[windshearax]])

if netcdf:

    ds_dict = {
                fill               : (('coord', 'storm', 'azimuth', 'range'),  coord_stack(filldict)),
                "best_track_files" : (('storm'), best_track_files),
                "narr_files"       : (('storm'), narr_files),
                "lon"              : (('storm'), lons),
                "lat"              : (('storm'), lats),
                "time"             : (('storm'), narr_valid_times)
              }
    
    # add line contour, wind barb components and quiver components, if they exist. If line, barb, or quiver is same string as fill, fill will be replaced. no diff for line, but barb and quiver have u and v components too.
    if line: ds_dict[line] = (('coord', 'storm', 'azimuth', 'range'), coord_stack(linedict))
    if barb: ds_dict[barb] = (('coord', 'storm', 'component', 'azimuth', 'range'), coord_stack(barbdict))
    if quiver: ds_dict[quiver] = (('coord', 'storm', 'component', 'azimuth', 'range'), coord_stack(quiverdict))
    ds = xarray.Dataset(ds_dict, coords={'coord' : ["north","storm motion","wind shear"], 'range' : range_center1D, 'azimuth' : theta_center1D, 'storm':stormname_years, 'component':["u","v"]})

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
    logging.info(f"wrote {netcdf}")

