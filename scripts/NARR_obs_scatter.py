import argparse
import datetime
import logging
import matplotlib.pyplot as plt
import os
import pandas as pd
from pathlib import Path
import seaborn as sns
import sys

sns.set_theme(style="darkgrid")


def main():
    parser = argparse.ArgumentParser(
        description="scatterplot of NARR vs observed variables",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-d", "--debug", action="store_true", help="print debug messages"
    )
    parser.add_argument(
        "-f", "--force_new", action="store_true", help="overwrite old file"
    )
    parser.add_argument(
        "ifiles",
        nargs="+",
        type=argparse.FileType("r"),
        help="input csv file(s) with data to plot",
    )
    parser.add_argument(
        "-l", "--label", action="store_true", help="label points with time and station"
    )
    parser.add_argument(
        "--no-fineprint",
        action="store_true",
        help="Don't write details at bottom of image",
    )
    parser.add_argument(
        "-x", type=str, default="metpy sfcape [J / kg]", help="x-axis variable name"
    )
    parser.add_argument(
        "-y", type=str, default="narr sbcape [J / kg]", help="y-axis variable name"
    )
    args = parser.parse_args()

    debug = args.debug
    ifiles = [x.name for x in args.ifiles]
    label = args.label
    no_fineprint = args.no_fineprint
    varx = args.x
    vary = args.y

    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(format="%(asctime)s %(message)s", level=level, force=True)

    logging.debug(args)

    logging.info(f"Reading {len(ifiles)} csv files")
    df = pd.concat(
        [pd.read_csv(ifile, header=0, parse_dates=["time"]) for ifile in ifiles]
    )

    if varx not in df.columns or vary not in df.columns:
        if varx not in df.columns:
            logging.error(f"{varx} not available")
        if vary not in df.columns:
            logging.error(f"{vary} not available")
        logging.error(f"choices {df.columns}")
        sys.exit(1)

    df = df.rename(columns={"sounding list": "LTC category"})
    # remove dirname and extension from category column.
    df["LTC category"] = [
        os.path.basename(os.path.splitext(slist)[0]).replace(".", "\n")
        for slist in df["LTC category"]
    ]

    # This is kind of kludgy, but better than before
    obs = df[df["sounding type"] == "obs"].reset_index()
    narr = df[df["sounding type"] == "NARR"].reset_index()

    narrvalues = narr[varx]
    if varx == vary:
        # differentiate column names
        varx = "narr " + varx
    obs[varx] = narrvalues  # Replace varx in obs with narr values
    # Then plot obs
    blue, orange, green, red, purple, brown, pink, gray, yellow, teal = (
        sns.color_palette()
    )
    colors = dict(upshr=teal, dnshr=red)
    ax = sns.scatterplot(
        data=obs,
        x=varx,
        y=vary,
        hue="place",
        palette=colors,
        style="LTC category",
        markers=["o", "^", "P"],
    )  # Tried "+" instead of "s" (square) but ValueError: Filled and line art markers cannot be mixed. Tried "*" instead of "o" but star is too tiny.
    plt.setp(ax.get_legend().get_texts(), fontsize=5)  # for legend text

    # one-to-one line
    ax.plot(ax.get_xlim(), ax.get_xlim(), color=gray, zorder=ax.get_zorder() - 1)
    ax.set_aspect("equal")

    if label:
        labels = obs[["time", "station"]]
        labels = [
            f"{time.strftime('%Y%m%d')}\n{station}" for time, station in labels.values
        ]
        for x, y, s in zip(obs[varx], obs[vary], labels):
            ax.text(x, y, s, fontsize="xx-small", ha="center", va="center_baseline")

    assert len(obs) == len(narr)
    text = f"n = {len(obs)}\n"
    text += f"Radiosonde Mean = {obs[vary].mean():.2f}\n"
    text += f"NARR Mean Absolute Error = {abs(obs[varx] - obs[vary]).mean():.2f}\n"
    text += f"NARR Bias = {obs[varx].mean() - obs[vary].mean():+.2f}\n"
    text += f"$R^2$ (Radiosonde, NARR) = {obs[varx].corr(obs[vary])**2:.2f}\n"
    text += f"created {datetime.datetime.now()}"
    fineprint = plt.annotate(
        text=text,
        xy=(2, 1),
        xycoords=("figure pixels", "figure pixels"),
        va="bottom",
        fontsize=6,
    )
    if no_fineprint:
        fineprint.set_visible(False)

    odir = Path(__file__).parent.parent.absolute() / "output"
    svarx = shortstr(varx)
    svary = shortstr(vary)
    ofile = os.path.join(odir, f"{svarx}.{svary}.scatter.png")
    plt.tight_layout()
    if label:
        plt.show()
    plt.savefig(ofile, dpi=175)
    logging.info(f"made {os.path.realpath(ofile)}")


def shortstr(s):
    # Turn string variable names with units into short strings appropriate for file name
    i = s.index("[")
    s = s[: i - 1]
    s = s.replace(" ", "_")
    return s


if __name__ == "__main__":
    main()
