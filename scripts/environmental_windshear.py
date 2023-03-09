import argparse
import atcf
import ibtracs
import logging
import metpy
from metpy.units import units
import narr
import numpy as np
import os
import pandas as pd 
import pdb
import pytz


def main():
    # =============Arguments===================
    parser = argparse.ArgumentParser(description = "NARR env shear", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-d", "--debug", action='store_true')
    parser.add_argument("ifile", help="text file. each line starts with a storm name, followed by a yyyymmddhh time")

    # Assign arguments to simple-named variables
    args = parser.parse_args()
    debug = args.debug
    ifile = args.ifile

    logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
    if debug:
        logger = logging.getLogger().setLevel(logging.DEBUG)

    logging.debug(args)

    logging.info(f"ifile {ifile}")
    df = pd.read_csv(ifile, delim_whitespace=True, names=["stormname", "narrtime"])
    df["narrtime"] = pd.to_datetime(df["narrtime"], format="%Y%m%d%H", utc=True)

    all_tracks, best_track_file = ibtracs.get_df(basin="NA") # North Atlantic basin guaranteed for NARR

    logging.info(f"got {len(all_tracks)} track data from {best_track_file}")

    out = df.groupby(["stormname","narrtime"], group_keys=False).apply(getTCFlow, all_tracks)
    avg = out.mean(numeric_only=True)
    u, v = metpy.calc.wind_components(avg.wind_shear_speed*units.meter/units.second, avg.wind_shear_heading * units.deg + 180*units.deg)
    print(f"avg u={u:~.2f} avg v={v:~.2f} avg shear mag={avg.wind_shear_speed:.2f}")

def getTCFlow(df, all_tracks):
    stormname, narrtime = df.name
    year = narrtime.year
    logging.debug(f"{stormname} {narrtime} {year}")
    file_ncl = narr.get(narrtime, narrtype=narr.narr3D, targetdir=os.path.join("/glade/scratch",os.getenv('USER'),"NARR"))
    imatch = (all_tracks['stormname'] == stormname.upper()) & (all_tracks['valid_time'].dt.year == year)
    trackdf = all_tracks[imatch]
    extension = ibtracs.extension(stormname, year)
    trackdf = pd.concat([trackdf, extension], ignore_index=True)
    trackdf["valid_time"] = trackdf.valid_time.apply(pytz.utc.localize) # make timezone-aware even without spc reports
    trackdf = trackdf[trackdf.rad.astype("float") <= 35]
    trackdf = trackdf[trackdf.valid_time == narrtime]
    pbot = 850
    ptop = 200
    rx = 4.5
    TCFlow_csv = os.path.join(os.path.dirname(file_ncl), 
            f"{stormname}.{narrtime.strftime('%Y%m%d%H')}.{pbot}-{ptop}hPa.{rx:.1f}degrees.000.csv") 
    if os.path.exists(TCFlow_csv):
        logging.debug(f"use existing TCFlow csv file {TCFlow_csv}")
        TCFlow = pd.read_csv(TCFlow_csv)
    else:
        logging.info(f"no TCFlow csv file {TCFlow_csv} making it.")
        TCFlow = ncl.steering_flow(lat0=trackdf.lat,lon0=trackdf.lon,storm_heading=trackdf.heading,storm_speed=track_df.speed,
                stormname=stormname, rx=rx, ensmember=stormname, ptop=ptop, pbot=pbot, file_ncl=file_ncl)
    return TCFlow
if __name__ == "__main__":
    main()
