import argparse
import ibtracs
import logging
import metpy
from metpy.units import units
import narr
import ncl
import numpy as np
import os
import pandas as pd
from pathlib import Path
import pdb
import pytz


def main():
    # =============Arguments===================
    parser = argparse.ArgumentParser(description = "NARR env shear", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--debug", action='store_true', help="more debug messages and pdb.set_trace before ending")
    parser.add_argument("-n", action='store_true', help='list missing storm and times but do not fill them')
    parser.add_argument("ifile", help="text file. each line starts with a storm name, followed by a yyyymmddhh time"
            "like $TMPDIR/strong_LTC_many_tornadoes.0003060912151821z.txt "
            "Don't provide the category files because they only have the first hour of the diurnal cycle")

    # Assign arguments to simple-named variables
    args = parser.parse_args()
    debug = args.debug
    ifile = args.ifile

    logging.basicConfig(format='%(asctime)s %(message)s')
    if debug:
        logger = logging.getLogger().setLevel(logging.DEBUG)

    logging.debug(args)

    logging.info(f"ifile {ifile}")
    df = pd.read_csv(ifile, delim_whitespace=True, names=["stormname", "narrtime"])
    df["narrtime"] = pd.to_datetime(df["narrtime"], format="%Y%m%d%H", utc=True)

    all_tracks, best_track_file = ibtracs.get_df(basin="NA") # North Atlantic basin guaranteed for NARR
    all_tracks = all_tracks[all_tracks.rad == 34]

    logging.info(f"got {len(all_tracks)} tracks from {best_track_file}")

    out = df.groupby(["stormname","narrtime"], group_keys=False).apply(getTCFlow, all_tracks, args)
    wind_shear_vector  = metpy.calc.wind_components(out["wind_shear_speed"].values*units.m/units.s, (out["wind_shear_heading"].values +180.)*units.deg)
    storm_motion_vector = metpy.calc.wind_components(out["storm_motion_speed"].values*units.m/units.s, (out["storm_motion_heading"].values +180.)*units.deg)
    # shear_x_motion is positive when wind shear points to right of storm motion.
    shear_x_motion = np.cross(np.array(wind_shear_vector).T,np.array(storm_motion_vector).T)
    out["shear_x_motion"] = shear_x_motion
    out["shear_rt_of_motion"] = out["shear_x_motion"] > 0
    out["shear_x_motion_norm"] = shear_x_motion / out.wind_shear_speed / out.storm_motion_speed
    assert all(df.groupby("stormname").narrtime.count() > 1), (  
            "You should have multiple narrtimes for each storm. Use file like $TMPDIR/strong_LTC_many_tornadoes.0003060912151821z.txt\n"
            "you don't want the category files because they only have the first hour of the diurnal cycle")
    if debug:
        print(out)
    avg = out.mean(numeric_only=True)
    assert len(df) == len(out), f"{df} {out} not equal length"
    # Not sure why I calculated du and dv. they are just the components of the imaginary vector with magnitude equal to average wind shear magnitude
    # at each time and heading equal to the mean of the wind shear headings (not a valid, circular average)
    u, v = metpy.calc.wind_components(avg.wind_shear_speed*units.meter/units.second, avg.wind_shear_heading * units.deg + 180*units.deg)
    out["year"] = pd.to_datetime(out.init,format="%Y-%m-%d_%H:%M:%S").dt.year
    print(f"{os.path.basename(ifile).replace('.0003060912151821z.txt','').center(53)}  "
           # f"n={len(df)} avg shear mag={avg.wind_shear_speed:.2f}  "
           # f"avg cross={avg.shear_x_motion:.1f} avg cross_norm={avg.shear_x_motion_norm:.1}  "
           # f"shear right of motion={avg.shear_rt_of_motion:.1%}  "
            f"{out.groupby(['ensmember','year']).ngroups} cases  "
            f"{np.sqrt(wind_shear_vector[0].mean()**2+wind_shear_vector[1].mean()**2):~5.2f}  "
            f"{metpy.calc.wind_direction(wind_shear_vector[0].mean(), wind_shear_vector[1].mean()):~.0f}")
    base, ext = os.path.splitext(ifile) # separate extension and ignore
    ofile = Path(__file__).parent.parent / "output" / f"{base}.env_windshear.csv"
    # Output numeric columns
    out.select_dtypes(include=[np.number]).to_csv(ofile)
    logging.info(f"write {ofile}")

    if debug:
        pdb.set_trace()

fmt = '%Y%m%d%H'
def getTCFlow(df, all_tracks, args):
    stormname, narrtime = df.name
    year = narrtime.year
    logging.debug(f"{stormname} {narrtime} {year}")
    file_ncl = narr.get(narrtime, narrtype=narr.narr3D, targetdir="/glade/scratch/ahijevyc/NARR")
    imatch = (all_tracks['stormname'] == stormname.upper()) & (all_tracks['valid_time'].dt.year == year)
    trackdf = all_tracks[imatch]
    if narrtime > trackdf.valid_time.max():
        trackdf = ibtracs.getExt(stormname, year, trackdf, [narrtime])
    trackdf = trackdf[trackdf.valid_time == narrtime]
    pbot = 850*units.hPa
    ptop = 200*units.hPa
    rx = 4.5 # no units (like NARR_composite.py)
    TCFlow_csv = os.path.join(os.path.dirname(file_ncl),
            f"{stormname}.{narrtime.strftime(fmt)}.{pbot.m:.0f}-{ptop.m:.0f}hPa.{rx:.1f}degrees.000.csv") 
    if os.path.exists(TCFlow_csv):
        logging.debug(f"use existing TCFlow csv file {TCFlow_csv}")
        TCFlow = pd.read_csv(TCFlow_csv)
    else:
        logging.warning(f"no {stormname} {narrtime.strftime(fmt)} TCFlow csv file {TCFlow_csv} making it.")
        if args.n:
            return
        # make dataframe a series (just like in NARR_composite.py)
        track = trackdf.squeeze()
        lat0=track.lat*units.degrees_N
        lon0=track.lon*units.degrees_E
        storm_heading=track.heading*units.deg
        storm_speed = (track.speed*units.kts).to("m/s")
        TCFlow = ncl.steering_flow(lat0=lat0,lon0=lon0,storm_heading=storm_heading,storm_speed=storm_speed,
                stormname=stormname, rx=rx, ensmember=stormname, ptop=ptop, pbot=pbot, file_ncl=file_ncl)
    return TCFlow
if __name__ == "__main__":
    main()
