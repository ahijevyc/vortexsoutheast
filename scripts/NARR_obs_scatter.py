import argparse
import datetime
import logging
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from metpy.interpolate import interpolate_1d
from metpy.units import units, pandas_dataframe_to_unit_arrays
import metpy.calc as mpcalc
import narr
import os
import pandas as pd
import pdb
import seaborn as sns
import sys

sns.set_theme(style="darkgrid")

def main():
    parser = argparse.ArgumentParser(description='scatterplot of NARR vs observed variables', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d','--debug', action="store_true", help='print debug messages')
    parser.add_argument('-f','--force_new', action="store_true", help='overwrite old file')
    parser.add_argument('ifiles', nargs="+", type=argparse.FileType("r"), help='input csv file(s) with data to plot')
    parser.add_argument("-l", "--label", action='store_true', help="label points with time and station")
    parser.add_argument("--no-fineprint", action='store_true', help="Don't write details at bottom of image")
    parser.add_argument('-x', type=str, default="metpy sfcape [J / kg]", help='x-axis variable name')
    parser.add_argument('-y', type=str, default="narr sbcape [J / kg]", help='y-axis variable name')
    args = parser.parse_args()


    debug       = args.debug
    force_new   = args.force_new
    ifiles      = [x.name for x in args.ifiles]
    label       = args.label
    no_fineprint=args.no_fineprint
    varx        = args.x
    vary        = args.y

    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(format='%(asctime)s %(message)s', level=level)

    logging.debug(args)

    logging.info(f"Reading {len(ifiles)} csv files") 
    df = pd.concat( [pd.read_csv(ifile, header=0, parse_dates=["time"]) for ifile in ifiles] )


    if varx not in df.columns or vary not in df.columns:
        if varx not in df.columns: logging.error(f"{varx} not available")
        if vary not in df.columns: logging.error(f"{vary} not available")
        logging.error(f"choices {df.columns}")
        sys.exit(1)

    # remove dirname and extension from sounding list column.
    df["sounding list"] = [os.path.basename(os.path.splitext(slist)[0]).replace(".","\n") for slist in df["sounding list"]]

    # This is kind of kludgy, but better than before
    obs  = df[df["sounding type"] == "obs"].reset_index()
    narr = df[df["sounding type"] == "NARR"].reset_index()
    obs[vary] = narr[vary] # Replace vary in obs with vary in narr
    # Then plot obs
    ax = sns.scatterplot(data=obs, x=varx, y=vary, hue="place", style="sounding list")
    plt.setp(ax.get_legend().get_texts(), fontsize='6') # for legend text
    plt.setp(ax.get_legend().get_title(), fontsize='8') # for legend title
   
    ax.plot(ax.get_xlim(),ax.get_xlim(), color="red", linewidth=0.5)
    ax.set_aspect('equal')

    if label:
        labels = obs[["time","station"]]
        labels = [f"{time.strftime('%Y-%m-%d')}\n{station}" for time,station in labels.values] 
        for x,y,s in zip(obs[varx],obs[vary],labels):
            ax.text(x,y,s,fontsize="xx-small", ha="center", va="center_baseline")

    text = f"created {datetime.datetime.now()}"
    fineprint = plt.annotate(text=text, xy=(2,1), xycoords=('figure pixels','figure pixels'), va="bottom", fontsize=6)
    if no_fineprint: fineprint.set_visible(False)

    odir = "/glade/scratch/ahijevyc/trier/VSE"
    ofile = os.path.join(odir, f"t.scatter.png")
    plt.tight_layout()
    plt.savefig(ofile,dpi=175)
    logging.info(f"made {os.path.realpath(ofile)}")


if __name__ == '__main__':
    main()
