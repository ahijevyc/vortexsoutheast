#!/bin/csh

foreach f (\
    ../closest_observed_soundings/sometorn.txt \
    ../closest_observed_soundings/strongLTC.manytorn.txt\
    ../closest_observed_soundings/weakLTC.fewornotorn.txt\
    ../closest_observed_soundings/weakLTC.manytorn.txt\
    )
    set csv=../output/`basename $f txt`csv
    python ../closest_observed_soundings/read_origin_email.py $f > $csv
    python NARR_skewt.py $csv
end
