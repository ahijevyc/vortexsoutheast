from atcf import dist_bearing
import argparse
import cartopy
import datetime
import ibtracs
import logging
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from metpy.interpolate import interpolate_1d
# SkewT segfaults when outside of base conda env on casper
# errors with different units registries if you use base conda env. use `vortexsoutheast` conda env
from metpy.plots import Hodograph, SkewT
from metpy.units import units, pandas_dataframe_to_unit_arrays
import metpy.calc as mpcalc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import narr
import numpy as np
import os
import pandas as pd
import pdb
from siphon.simplewebservice.wyoming import WyomingUpperAir
import sys
import xarray


logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO, force=True)

def getparser():
    parser = argparse.ArgumentParser(description='skew t of output from NARR', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--debug', action="store_true", help='print debug messages')
    parser.add_argument('-f','--force_new', action="store_true", help='overwrite old file')
    parser.add_argument('ifile', type=argparse.FileType("r"), help='input csv file with initialization time,station')
    parser.add_argument("--no-fineprint", action='store_true', help="Don't write details at bottom of image")
    parser.add_argument('--pmin', type=int, default=100, help='mask values below this pressure')
    return parser

def main():
    args = getparser().parse_args()

    debug       = args.debug
    force_new   = args.force_new
    ifile       = args.ifile.name
    no_fineprint=args.no_fineprint
    pmin        = args.pmin

    level = logging.DEBUG if debug else logging.INFO
    logging.getLogger().setLevel(level)

    logging.debug(args)

    df = pd.read_csv(ifile, names=["storm","itime","place","station"], parse_dates=["itime"])
    assert set(df["place"].unique()) == set(["upshr","dnshr"]), "expected 2 places: upshr and dnshr"
    df["itime"] = pd.to_datetime(df["itime"], utc=True)

    odf = pd.DataFrame()

    fig = plt.figure(figsize=(8,8))

    ibtracs_df, ibtracs_file = ibtracs.get_df(basin="NA")
    # filter out "duplicate" 50 and 64-kt lines or it will break ibtracs.getExt()
    ibtracs_df = ibtracs_df[ibtracs_df.rad == 34]

    for index, row in df.iterrows():
        storm = row.storm
        year = row.itime.year
        trackdf = ibtracs.this_storm(ibtracs_df, storm, year)
        logging.debug(f"got {storm} {year} ({len(trackdf)} lines) from {ibtracs_df}")
        itime = row.itime
        if itime > trackdf.valid_time.max():
            logging.warning(f"requested time {itime} is after {storm} {year} ends")
            trackdf = ibtracs.getExt(storm, year, trackdf, [itime])
        
        TC = trackdf.set_index("valid_time")[["lon","lat"]]
        TC = TC.resample('H').interpolate(method="linear").loc[itime]
        station = row.station
        cache = f"/glade/scratch/ahijevyc/wyocache/{itime.strftime('%Y%m%dT%H%M')}{station}.csv"
        if os.path.exists(cache):
            obs = pd.read_csv(cache, index_col=0)
            obs = obs.drop_duplicates(subset="pressure", ignore_index=True) # Avoid ValueError: Values in `t_eval` are not properly sorted.
            udict = {'pressure': 'hPa', 'height': 'meter', 'temperature': 'degC', 'dewpoint': 'degC', 'direction': 'degrees', 'speed': 'knot', 'u_wind': 'knot', 'v_wind': 'knot', 'station': None, 'station_number': None, 'time': None, 'latitude': 'degrees_north', 'longitude': 'degrees_east', 'elevation': 'meter'}
        else:
            obs = WyomingUpperAir.request_data(itime, station)
            obs.to_csv(cache)
            udict = obs.units
        obs = obs[obs["pressure"] >= pmin]
        obs = pandas_dataframe_to_unit_arrays(obs, column_units=udict)
        slon = obs["longitude"][0]
        slat = obs["latitude"][0]
        # put slon and slat in lists because dist_bearing expects iterables in args 3,4
        distance_from_TCcenter, bearing_from_TCcenter = dist_bearing(TC.lon*units.deg, TC.lat*units.deg, [slon], [slat])
        sigps = [1000, 925, 850, 700, 500, 400, 300, 200, 100] * units.hPa
        od = {}
        for stype in ["NARR", "obs"]:
            plot_barbs_units = "m/s"
            skew = SkewT(fig, rotation=45)
            trans = transforms.blended_transform_factory(skew.ax.transAxes, skew.ax.transData)
            if stype == "obs":
                p = obs["pressure"]
                T = obs["temperature"]
                height = obs["height"]
                Td = obs["dewpoint"]
                u = obs["u_wind"].to(plot_barbs_units)
                v = obs["v_wind"].to(plot_barbs_units)
                od["psfc"] = p.max()
                od["pwat"] = mpcalc.precipitable_water(p,Td)
                # Get rh and rh at 0degC level.
                relative_humidity = mpcalc.relative_humidity_from_dewpoint(T,Td).to("percent")
                od["0degC rh"], = interpolate_1d(0*units.parse_expression("degC"), T, relative_humidity)
                lon, lat  = slon, slat
                zs, T0s, Td0s = interpolate_1d(sigps, p, height, T, Td) # ignore UserWarning: Interp pt out of data bounds... they're set to nan.
                for z,T0,Td0,sigp in zip(zs,T0s,Td0s,sigps):
                    if np.isnan(z): continue
                    skew.ax.text(0.01,sigp,f"{z:~.0f}", fontsize='small', transform=trans, ha="left", va="center")
                    od[f"{sigp:.0f} temp"] = T0
                    od[f"{sigp:.0f} thetae"] = mpcalc.equivalent_potential_temperature(sigp, T0, Td0).to("degC")
            else:
                logging.debug(f"open {itime} NARR")
                ds = xarray.open_dataset(narr.get(itime))
                ds = ds.rename_dims(dims_dict=narr.dims_dict)
                ds = ds.metpy.assign_crs(narr.data_crs.to_cf()) 
                ds = ds.metpy.assign_y_x(tolerance=2*units.m)

                logging.debug("reverse pressure dimension")
                ds = ds.reindex(lv_ISBL0=ds.lv_ISBL0[::-1]) # 1000 to 100mb
                logging.debug("find closest gridpoint to station")

                # transform slat/slon to x/y with cartopy
                x, y = narr.data_crs.transform_point(slon.m, slat.m, src_crs=cartopy.crs.PlateCarree())

                imin = ((ds.x-x)**2 + (ds.y-y)**2).argmin(dim=["x","y"])
                prof = ds.isel(imin)
                nearest = dict(x=x,y=y, method="nearest", tolerance=17000*units.m)
                if not prof.equals(ds.metpy.sel(**nearest)):
                    logging.warning(f"sel method pulled different profile than imin")
                    pdb.set_trace()


                lon = prof.gridlon_221.metpy.quantify().data
                lat = prof.gridlat_221.metpy.quantify().data
                p = prof.lv_ISBL0
                T = prof["TMP_221_ISBL"].metpy.convert_units("degC") # skew.plot() assumes T is in degC
                height = prof["HGT_221_ISBL"]
                Td = mpcalc.dewpoint_from_specific_humidity(p, T, prof["SPF_H_221_ISBL"])
                barb_increments = dict(flag=25,full=5,half=2.5)
                u = prof["U_GRD_221_ISBL"]
                v = prof["V_GRD_221_ISBL"]
                u = u.metpy.convert_units(plot_barbs_units)
                v = v.metpy.convert_units(plot_barbs_units)
                logging.debug("NARR temp and thetae")
                for sigp in sigps:
                    z = prof["HGT_221_ISBL"].sel(dict(lv_ISBL0=sigp))
                    skew.ax.text(0.01,sigp,f"{z.astype(int).metpy.quantify().data:~}", fontsize='small', transform=trans, ha="left", va="center")
                    T0  =  T.sel(dict(lv_ISBL0=sigp))
                    Td0 = Td.sel(dict(lv_ISBL0=sigp))
                    od[f"{sigp:.0f} temp"] = T0.data
                    od[f"{sigp:.0f} thetae"] = mpcalc.equivalent_potential_temperature(sigp, T0, Td0).data.to("degC")
                logging.debug("NARR sbcape")
                od["narr sbcape"]    = narr.scalardata('sbcape', itime).isel(imin).data.compute()
                od["narr sbcinh"]    = narr.scalardata('sbcinh', itime).isel(imin).data.compute()
                od["narr mlcape"]    = narr.scalardata('mlcape', itime).isel(imin).data.compute()
                od["narr mlcinh"]    = narr.scalardata('mlcinh', itime).isel(imin).data.compute()
                od["narr 0degC rh"]  = narr.scalardata('rh_0deg', itime).isel(imin).data.compute()
                # compute() to avoid TypeError: len() of unsized object when printing LCL text
                # seems to be too high up in NARR. derive it from 1000 hPa parcel instead
                #od["narr lcl_pressure"]   = narr.scalardata('lcl', itime).isel(imin).data.compute()
                od["psfc"]           = narr.scalardata('psfc', itime).isel(imin).data.compute()
                od["pwat"]           = narr.scalardata('pwat', itime).isel(imin).data.compute() * units.m**3 / (1000*units.kg)
                od["pwat"]           = od["pwat"].to("mm")
                logging.debug("NARR shears")
                od["narr shr10m_900hPa"] = narr.scalardata('shr10m_900hPa', itime).isel(imin).data.compute()
                od["narr shr10m_700hPa"] = narr.scalardata('shr10m_700hPa', itime).isel(imin).data.compute()
                od["narr shr10m_1000m"] = narr.scalardata('shr10m_1000m', itime).isel(imin).data.compute()
                od["narr shr10m_3000m"] = narr.scalardata('shr10m_3000m', itime).isel(imin).data.compute()

            skew.plot(p, T, 'r')
            skew.plot(p, Td, 'g')

            skew.ax.set_ylim(1050, pmin)
            skew.ax.set_xlim(-25, 55)

            logging.debug("Calculate full parcel profile and plot as black line")

            p_without_lcl = p # remember so we can interpolate u,v,height too later
            logging.debug("Get parcel temperature at all pressures, including lcl. Return new p, T, and Td with lcl included.")
            p, T, Td, profT = mpcalc.parcel_profile_with_lcl(p, T, Td) # not virtual temperature yet
            if (p[0], profT[0]) == (p[1], profT[1]):
                assert (T[0], Td[0]) == (T[1], Td[1]), "bottom pressure repeated but (T,Td) differ"
                # Avoid problem when lcl is at bottom and data are repeated
                #   File "/glade/work/ahijevyc/20201220_casper/lib/python3.7/site-packages/metpy/calc/tools.py", line 425, in _get_bound_pressure_height
                #    if height is not None and not (_less_or_close(bound_height, np.nanmax(height))
                # ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
                p = p[1:]
                T = T[1:]
                Td = Td[1:]
                profT = profT[1:]
            logging.debug("Now p, T, Td, and prof are defined at lcl.  Add u, v, height at lcl")
            u, v, height = interpolate_1d(p, p_without_lcl, u, v, height)
            
            logging.debug("environmental relative humidity")
            relative_humidity = mpcalc.relative_humidity_from_dewpoint(T, Td)

            logging.debug("parcel virtual temperature")
            logging.debug("parcel mixing ratio (use environmental relative humidity at the surface)")
            parcel_level = 0 if stype == "obs" else 0
            mixing_ratio_parcel = mpcalc.mixing_ratio_from_relative_humidity(p[parcel_level], T[parcel_level], relative_humidity[parcel_level])
            prof_mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(p, profT, 1) # saturated mixing ratio
            od["lcl_pressure"], _ = mpcalc.lcl(p[parcel_level], T[parcel_level], Td[parcel_level])
            prof_mixing_ratio[p >= od["lcl_pressure"]] = mixing_ratio_parcel # unsaturated mixing ratio (constant up to LCL)
            profTv = mpcalc.virtual_temperature(profT, prof_mixing_ratio)

            logging.debug("environment virtual temperature")
            mixing_ratio = mpcalc.mixing_ratio_from_relative_humidity(p, T, relative_humidity)
            cold_and_probably_dry = np.isnan(relative_humidity) & (T < -40 * units.degC)
            if np.any(cold_and_probably_dry):
                logging.warning(f"replace {cold_and_probably_dry.sum()} missing mixing ratios with zero where T < -40C")
                mixing_ratio[cold_and_probably_dry] = 0
            Tv = mpcalc.virtual_temperature(T, mixing_ratio)
            skew.plot(p, Tv, 'r', linewidth=0.5, linestyle="dashed")
            skew.ax.plot([0.83,0.85], 2*[od["lcl_pressure"].m], transform=trans)
            skew.ax.text(0.855,od["lcl_pressure"],f'sfc LCL ({od["lcl_pressure"]:~.0f})', transform=trans, ha="left", va="center", fontsize='x-small')

            sfcape, sfcin = mpcalc.cape_cin(p,Tv,Td,profTv)
            storm_u = 0. * units.m / units.s
            storm_v = 0. * units.m / units.s
            right_mover, left_mover, wind_mean = mpcalc.bunkers_storm_motion(p, u, v, height)
            storm_u, storm_v = wind_mean
            srh03_pos, srh03_neg, srh03_tot = mpcalc.storm_relative_helicity(height, u, v, 3 * units.km, storm_u=storm_u, storm_v=storm_v)
            srh01_pos, srh01_neg, srh01_tot = mpcalc.storm_relative_helicity(height, u, v, 1 * units.km, storm_u=storm_u, storm_v=storm_v)
         

            od["time"] = itime
            od["place"] = row.place
            od["station"] = station
            od["sounding type"] = stype
            od["station lon"] = slon
            od["station lat"] = slat
            od["lon"] = lon
            od["lat"] = lat
            od["lowest level temp"] = T[0]
            od["lowest level dwpt"] = Td[0]
            od["lowest level rh"] = mpcalc.relative_humidity_from_dewpoint(T[0], Td[0])
            od["lowest level theta"] = mpcalc.potential_temperature(p[0], T[0])
            od["lowest level thetae"] = mpcalc.equivalent_potential_temperature(p[0], T[0], Td[0])
            od["sfc parcel mixing ratio"] = mixing_ratio_parcel
            od["sfc parcel lcl pressure AGL"] = od["psfc"] - od["lcl_pressure"]
            od["metpy sfcape"] = sfcape
            od["metpy sfcin"] = sfcin
            od["storm_u"] = storm_u
            od["storm_v"] = storm_v
            od["sfc-5km shr mag"] = mpcalc.wind_speed(*mpcalc.bulk_shear(p, u, v, height, bottom=height.min(), depth=5*units.km))
            od["sfc-3km shr mag"] = mpcalc.wind_speed(*mpcalc.bulk_shear(p, u, v, height, bottom=height.min(), depth=3*units.km))
            od["sfc-1km shr mag"] = mpcalc.wind_speed(*mpcalc.bulk_shear(p, u, v, height, bottom=height.min(), depth=1*units.km))
            od["metpy srh03+"] = srh03_pos
            od["metpy srh03-"] = srh03_neg
            od["metpy srh03"] = srh03_tot
            od["metpy srh01+"] = srh01_pos
            od["metpy srh01-"] = srh01_neg
            od["metpy srh01"] = srh01_tot
            #mlcape,mlcin = mpcalc.mixed_layer_cape_cin(p,Tv,Td,height=height,depth=30*units.hPa) # not sure how metpy deals with Tv. 
            #od["metpy mlcape"] = mlcape
            #od["metpy mlcin"] = mlcin
            od["sounding list"] = ifile
            od["distance_from_TCcenter"] = distance_from_TCcenter
            od["bearing_from_TCcenter"] = bearing_from_TCcenter

            skew.plot(p, profTv, 'k', linewidth=1.5, linestyle="dashed")

            logging.debug("Shade areas of CIN")
            skew.shade_cin(p, Tv, profTv)
            logging.debug("Shade areas of CAPE")
            skew.shade_cape(p, Tv, profTv)

            logging.debug("0degC isotherm")
            skew.ax.plot([0,0], [skew.ax.get_ylim()[0],300], color='c', linestyle='--', linewidth=2, alpha=0.6)

            logging.debug("dry and moist adiabats")
            skew.plot_dry_adiabats()
            skew.plot_moist_adiabats(lw=0.5)
            logging.debug("mixing ratio lines")
            lines = skew.plot_mixing_lines(alpha=0.5)
            for q, path in zip([0.4, 1, 2, 4, 7, 10, 16, 24, "32g/kg"], lines.get_paths()):
                skew.ax.text(*path.vertices[-3], q, fontsize='x-small', ha="center", color=lines.get_color()[0], alpha=1, clip_on=True)
            suptitle  = f"{stype} {station} ({lon:~.2f},{lat:~.2f}) {itime}"
            suptitle += f"\nwind barbs and hodograph in {plot_barbs_units} {barb_increments}"
            suptitle += f"\nmetpy sfcape={sfcape:~.0f}   sfcin={sfcin:~.0f}   storm_u={storm_u:~.1f}   storm_v={storm_v:~.1f}"
            suptitle += f"\nmetpy 0-3km srh+={srh03_pos:~.0f}   srh-={srh03_neg:~.0f}   srh(tot)={srh03_tot:.~0f}"
            suptitle += f"\nmetpy 0-1km srh+={srh01_pos:~.0f}   srh-={srh01_neg:~.0f}   srh(tot)={srh01_tot:.~0f}"
            plt.suptitle(suptitle, fontsize="x-small")

            logging.debug("hodograph")
            ax_hod = inset_axes(skew.ax, '35%', '35%', loc=1)
            h = Hodograph(ax_hod, component_range=30.)
            h.add_grid(increment=5, linewidth=0.75)
            ax_hod.plot(storm_u,storm_v,"X")
            ax_hod.set_xlabel("")
            ax_hod.set_ylabel("")

            label_hgts = [1, 3, 6, 10] * units.km
            h.plot_colormapped(u,v,height)


            for label_hgt in label_hgts:
                ax_hod.text(np.interp(label_hgt, height, u), np.interp(label_hgt, height, v), label_hgt.to("km").m, fontsize=7)
                ax_hod.plot(np.interp(label_hgt, height, u), np.interp(label_hgt, height, v), '.', color="white", markersize=2.5)

            bbz = skew.plot_barbs(p, u, v, length=6, plot_units=plot_barbs_units, linewidth=0.6, xloc=1.06, barb_increments=barb_increments)
            text = f"created {datetime.datetime.now()}\nsounding list {ifile}"
            fineprint = plt.annotate(text=text, xy=(2,1), xycoords=('figure pixels','figure pixels'), va="bottom", fontsize=6)
            if no_fineprint: fineprint.set_visible(False)

            odir = "/glade/scratch/ahijevyc/vortexsoutheast/output/skewT"
            ofile = os.path.join(odir, f"{itime.strftime('%Y%m%dT%H%M')}.{station}.{stype}.skewt.png")
            plt.savefig(ofile)
            logging.info(f"made {os.path.realpath(ofile)}")
            fig.clear()
            odf = pd.concat([odf,pd.DataFrame([od])], ignore_index=True) # Avoid FutureWarning about append. Keep od as dictionary (don't convert to DataFrame) 

    logging.info(f"move units from values to column names")
    odf = move_units_from_values_to_column_names(odf)
    odir = os.path.join(os.path.dirname(ifile), "../output")  # output directory 
    if not os.path.exists(odir):
        os.mkdir(odir)
    path, ext = os.path.splitext(os.path.basename(ifile))
    ocsv = os.path.join(odir, path+".data.csv")
    odf.set_index(["time","station","sounding type"]).to_csv(ocsv)
    logging.info(f"made {ocsv}")


def move_units_from_values_to_column_names(df):
    rdict = {}
    for col in df.columns:
        igood = ~df[col].isna()
        values = df.loc[igood,col].values
        # If first good element has units...
        if hasattr(values[0], 'units'):
            # Make sure units are all the same
            us = [x.units for x in values]
            assert len(set(us)) == 1, f"units of {col} not all the same {set(us)}"
            logging.debug(f"move {us[0]} of {col} from values to column name")
            # Take magnitude of column values (remove units)
            df.loc[igood,col] = [x.m for x in values]
            # append units to column name
            rdict[col] = f"{col} [{us[0]:~}]"
    df = df.rename(columns=rdict)
    return df

if __name__ == '__main__':
    main()
