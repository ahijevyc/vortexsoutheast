#!/bin/csh
module load conda
conda activate vse
foreach f (\
    ../closest_observed_soundings/sometorn.txt \
    ../closest_observed_soundings/strongLTC.manytorn.txt\
    ../closest_observed_soundings/weakLTC.fewornotorn.txt\
    ../closest_observed_soundings/weakLTC.manytorn.txt\
    )
    set csv=../output/`basename $f txt`csv
    python read_origin_close_soundings_email.py $f > $csv
    python NARR_skewt.py $csv
end
