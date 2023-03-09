import argparse
import datetime
import logging
import matplotlib.pyplot as plt
from metpy.interpolate import interpolate_1d
from metpy.plots import Hodograph, SkewT # SkewT segfaults when outside of base conda env on casper
from metpy.units import units
import metpy.calc as mpcalc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import os
import pdb
import sys
import VSE
import xarray


def main():
    parser = argparse.ArgumentParser(description='skew t of output from NARR_composite.py', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--coord', type=str, choices=["north", "storm motion", "wind shear"], default = "north", help='coordinate system. fields rotated so that this vector points up')
    parser.add_argument('-d','--debug', action="store_true", help='print debug messages')
    parser.add_argument('--desc', type=str, default="all_LTC", choices=VSE.composites(), help='description')
    parser.add_argument('-f','--force_new', action="store_true", help='overwrite old file')
    parser.add_argument('--hrs', type=str, default="0003z", help='hrs')
    parser.add_argument('--idir', type=str, default = "/glade/scratch/ahijevyc/trier/VSE/nc", help='input directory')
    parser.add_argument('--centroids', type=str, choices = VSE.centroids(), default="torn max", help="predetermined list of points by coord and centroid")
    args = parser.parse_args()


    centroids   = args.centroids
    coord       = args.coord
    debug       = args.debug
    desc        = args.desc
    force_new   = args.force_new
    hrs         = args.hrs
    idir        = args.idir
    points      = VSE.pointlist[coord][centroids]

    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(format='%(asctime)s %(message)s', level=level)

    logging.debug(args)

    skewtfields = ["hgt", "sh", "u", "v", "temp"]
    files = [ f"{idir}/{desc}.{x}.{hrs}.nc" for x in skewtfields]
    logging.debug(f"opening {files}")
    for f in files:
        assert os.path.exists(f), f"{os.path.realpath(f)} does not exist. Did u remember to run ~/bin/stack_NARR.csh?"
            
    #ds = xarray.open_mfdataset(files, combine="by_coords", concat_dim=None)
    ds = xarray.open_mfdataset(files)
    capeds = xarray.open_dataset(f"{idir}/{desc}.mlcape.{hrs}.nc")
    # Load xarray into memory instead of dask-delaying until you need it
    ds.load() # or else .shade_cin() throws NotImplementedError: The out parameter is not fully supported. Received type ndarray, expected Dask Array

    ds = ds.metpy.quantify().mean(dim="storm") # if you don't quantify() you lose units.

    fig = plt.figure(figsize=(8,8))
    for point in points:
        point_label = point.split("/")[0]
        azimuth     = point.split("/")[1].replace("deg","")
        range_km    = point.split("/")[2].replace("km","")

        narr = ds.sel(range=range_km, azimuth=azimuth, method="nearest").sel(coord=coord) # coord is a string. "nearest" method undefined.
        cape = capeds.sel(range=range_km, azimuth=azimuth, method="nearest").sel(coord=coord)
        srh = xarray.open_dataset(f"{idir}/{desc}.srh.{hrs}.nc").sel(range=range_km, azimuth=azimuth, method="nearest").sel(coord=coord)
        p = narr.lv_ISBL0.metpy.unit_array # .metpy.unit_array or shade_cin will break
        T = narr.temp.metpy.unit_array # .metpy.unit_array or shade_cin will break
        height = narr.hgt.metpy.unit_array
        Td = mpcalc.dewpoint_from_specific_humidity(p, T, narr.sh)
        barb_increments = {"flag":25,"full":5,"half":2.5}
        plot_barbs_units = "m/s"
        u = narr.u.metpy.convert_units(plot_barbs_units)
        v = narr.v.metpy.convert_units(plot_barbs_units)

        skew = SkewT(fig, rotation=45)
        skew.plot(p, T, 'r')
        skew.plot(p, Td, 'g')

        skew.ax.set_ylim(1000, 100)
        skew.ax.set_xlim(-25, 55)
        # Calculate LCL height and plot as black dot.
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
        logging.debug(f"lcl_pressure {lcl_pressure} lcl_temperature {lcl_temperature}")

        # Calculate full parcel profile and add to plot as black line


        p_without_lcl = p # remember so we can interpolate u,v,height too later
        # Get parcel temperature at all pressures, including lcl. Return new p, T, and Td with lcl included.
        p, T, Td, prof = mpcalc.parcel_profile_with_lcl(p, T, Td) # not virtual temperature yet
        # Now p, T, Td, and prof are defined at lcl.  Add u, v, height at lcl
        u, v, height = interpolate_1d(p, p_without_lcl, u, v, height)
        
        # environmental relative humidity
        relative_humidity = mpcalc.relative_humidity_from_dewpoint(T, Td)

        # parcel virtual temperature
        # parcel mixing ratio (use environmental relative humidity at the surface)
        mixing_ratio_parcel = mpcalc.mixing_ratio_from_relative_humidity(p[0], T[0], relative_humidity[0])
        prof_mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(p, prof, 1) # saturated mixing ratio
        prof_mixing_ratio[p >= lcl_pressure] = mixing_ratio_parcel # unsaturated mixing ratio (constant up to LCL)
        profTv = mpcalc.virtual_temperature(prof, prof_mixing_ratio)

        # environment virtual temperature
        mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(p, T, relative_humidity)
        Tv = mpcalc.virtual_temperature(T, mixing_ratio)
        skew.plot(p, Tv, 'r', linewidth=0.5, linestyle="dashed")
        import matplotlib.transforms as transforms
        trans = transforms.blended_transform_factory(skew.ax.transAxes, skew.ax.transData)
        skew.ax.plot([0.86,0.88], 2*[lcl_pressure.m], transform=trans)
        skew.ax.text(0.885,lcl_pressure,"sfc LCL", transform=trans, horizontalalignment="left", verticalalignment="center")

        sfcape, sfcin = mpcalc.cape_cin(p,Tv,Td,profTv)
        storm_u = 0. * units("m/s")
        storm_v = 0. * units("m/s")
        right_mover, left_mover, wind_mean = mpcalc.bunkers_storm_motion(p, u, v, height)
        storm_u, storm_v = wind_mean
        srh03_pos, srh03_neg, srh03_tot = mpcalc.storm_relative_helicity(height, u, v, 3 * units.km, storm_u=storm_u, storm_v=storm_v)
        srh01_pos, srh01_neg, srh01_tot = mpcalc.storm_relative_helicity(height, u, v, 1 * units.km, storm_u=storm_u, storm_v=storm_v)

        skew.plot(p, profTv, 'k', linewidth=1.5, linestyle="dashed")

        # Shade areas of CAPE and CIN
        skew.shade_cin(p, T, profTv, Td)
        skew.shade_cape(p, Tv, profTv)

        # An example of a slanted line at constant T -- in this case the 0
        # isotherm
        skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

        # Add the relevant special lines
        skew.plot_dry_adiabats()
        skew.plot_moist_adiabats(lw=0.5)
        skew.plot_mixing_lines()
        suptitle  = f"{desc}    {hrs}    {coord} coord    {point_label}   az={azimuth}{ds.azimuth.units}  r={range_km}{ds.range.units}"
        suptitle += f"\nwind barbs and hodograph in {plot_barbs_units} {barb_increments}"
        suptitle += f"\nNARR mlcape={cape.mlcape.metpy.quantify().mean().data:~.0f}    srh={srh.srh.metpy.quantify().mean().data:~.0f}"
        suptitle += f"\nmetpy sfcape={sfcape:~.0f}   sfcin={sfcin:~.0f}   storm_u={storm_u:~.1f}   storm_v={storm_v:~.1f}"
        suptitle += f"\nmetpy 0-3km srh+={srh03_pos:~.0f}   srh-={srh03_neg:~.0f}   srh(tot)={srh03_tot:.~0f}"
        suptitle += f"\nmetpy 0-1km srh+={srh01_pos:~.0f}   srh-={srh01_neg:~.0f}   srh(tot)={srh01_tot:.~0f}"
        plt.suptitle(suptitle, fontsize="x-small")


        ax_hod = inset_axes(skew.ax, '35%', '35%', loc=1)
        h = Hodograph(ax_hod, component_range=30.)
        h.add_grid(increment=5, linewidth=0.75)
        ax_hod.plot(storm_u,storm_v,"X")
        ax_hod.set_xlabel("")
        ax_hod.set_ylabel("")

        label_hgts = [1, 3, 6, 10] * units("km")
        h.plot_colormapped(u,v,height)


        for label_hgt in label_hgts:
            ax_hod.text(np.interp(label_hgt, height, u), np.interp(label_hgt, height, v), label_hgt.to("km").m, fontsize=7)
            ax_hod.plot(np.interp(label_hgt, height, u), np.interp(label_hgt, height, v), '.', color="white", markersize=2.5)

        bbz = skew.plot_barbs(p, u, v, length=6, plot_units=plot_barbs_units, linewidth=0.6, xloc=1.06, barb_increments=barb_increments)

        text = f"created {datetime.datetime.now()}"
        fineprint = plt.annotate(text=text, xy=(2,1), xycoords=('figure pixels','figure pixels'), va="bottom", fontsize=6)

        ofile = f"{desc}.{hrs}.{coord.replace(' ','_')}.{azimuth}{ds.azimuth.units}.{range_km}{ds.range.units}.skewt.png"
        plt.savefig(ofile)
        print("made", os.path.realpath(ofile))
        fig.clear()
if __name__ == '__main__':
    main()
