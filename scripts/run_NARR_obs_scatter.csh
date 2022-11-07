#!/bin/csh
conda activate
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.* -x "lowest level temp [°C]" -y "lowest level temp [°C]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.* -x "lowest level dwpt [°C]" -y "lowest level dwpt [°C]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.* -x "pwat [mm]" -y "pwat [mm]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.* -x "narr 0degC rh [%]" -y "0degC rh [%]"
