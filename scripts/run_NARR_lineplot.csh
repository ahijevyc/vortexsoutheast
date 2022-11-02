#!/bin/csh

conda activate
python NARR_lineplot.py mlcape --twin 6 --desc all_LTC --coord north --centroid "torn max"
python NARR_lineplot.py mlcinh --twin 6 --desc all_LTC --coord north --centroid "torn max"
python NARR_lineplot.py shr10m_500hPa --twin 6 --desc all_LTC --coord north --centroid "torn max"
python NARR_lineplot.py shr10m_5000m --twin 6 --desc all_LTC --coord north --centroid "torn max"
python NARR_lineplot.py rh500hPa --twin 6 --desc all_LTC --coord north --centroid "torn max"
python NARR_lineplot.py rh_0deg --twin 6 --desc all_LTC --coord north --centroid "torn max"
