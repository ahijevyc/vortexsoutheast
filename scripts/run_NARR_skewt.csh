#!/bin/csh
python NARR_skewt.py ../closest_observed_soundings/coastal.torn.csv
python NARR_skewt.py ../closest_observed_soundings/inland.torn.csv
python NARR_skewt.py ../closest_observed_soundings/strongLTC.manytorn.csv
python NARR_skewt.py ../closest_observed_soundings/weakLTC.fewornotorn.csv
python NARR_skewt.py ../closest_observed_soundings/weakLTC.manytorn.csv
