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
import pdb
import pytz


def main():
    # =============Arguments===================
    parser = argparse.ArgumentParser(description = "NARR env shear", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--debug", action='store_true')
    parser.add_argument("-n", action='store_true', help='list missing storm and times but do not fill them')
    parser.add_argument("ifile", help="text file. each line starts with a storm name, followed by a yyyymmddhh time")

    # Assign arguments to simple-named variables
    args = parser.parse_args()
    debug = args.debug
    ifile = args.ifile

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.WARNING, force=True)
    if debug:
        logger = logging.getLogger().setLevel(logging.DEBUG)

    logging.debug(args)

    logging.info(f"ifile {ifile}")
    df = pd.read_csv(ifile, delim_whitespace=True, names=["stormname", "narrtime"])
    df["narrtime"] = pd.to_datetime(df["narrtime"], format="%Y%m%d%H", utc=True)

    all_tracks, best_track_file = ibtracs.get_df(basin="NA") # North Atlantic basin guaranteed for NARR
    all_tracks = all_tracks[all_tracks.rad == 34]

    logging.info(f"got {len(all_tracks)} track data from {best_track_file}")

    out = df.groupby(["stormname","narrtime"], group_keys=False).apply(getTCFlow, all_tracks, args)
    avg = out.mean(numeric_only=True)
    assert len(df) == len(out), f"{df} {out} not equal length"
    u, v = metpy.calc.wind_components(avg.wind_shear_speed*units.meter/units.second, avg.wind_shear_heading * units.deg + 180*units.deg)
    print(f"n={len(df)} avg u={u:~.2f} avg v={v:~.2f} avg shear mag={avg.wind_shear_speed:.2f}")

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
