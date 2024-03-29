#!/bin/csh

conda activate
foreach coord (north "storm motion")
    python NARR_lineplot.py sbcape --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py sbcinh --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py mlcape --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py mlcinh --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py shr10m_500hPa --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py shr10m_5000m --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py shr10m_700hPa --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py shr10m_3000m --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py rh500hPa --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
    python NARR_lineplot.py rh_0deg --twin 6 --desc all_LTC --coord "$coord" --centroid "torn max" --simplegend
end
