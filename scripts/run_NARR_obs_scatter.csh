#!/bin/csh
conda activate
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "lowest level temp [°C]" -y "lowest level temp [°C]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "lowest level dwpt [°C]" -y "lowest level dwpt [°C]"

python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "pwat [mm]" -y "pwat [mm]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "narr 0degC rh [%]" -y "0degC rh [%]"

python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "metpy sfcape [J / kg]" -y "metpy sfcape [J / kg]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "narr sbcape [J / kg]" -y "metpy sfcape [J / kg]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "metpy srh03 [m ** 2 / s ** 2]" -y "metpy srh03 [m ** 2 / s ** 2]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "narr srh [m ** 2 / s ** 2]" -y "metpy srh03 [m ** 2 / s ** 2]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "metpy sfc-1km shr mag [m / s]" -y "metpy sfc-1km shr mag [m / s]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "narr shr10m_1000m [m / s]"     -y "metpy sfc-1km shr mag [m / s]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "metpy sfc-3km shr mag [m / s]" -y "metpy sfc-3km shr mag [m / s]"
python NARR_obs_scatter.py ../output/strongLTC.manytorn.data.csv ../output/weakLTC.*.data.csv -x "narr shr10m_3000m [m / s]"     -y "metpy sfc-3km shr mag [m / s]"
