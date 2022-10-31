#!/bin/csh
module load nco

cd /glade/scratch/ahijevyc/trier/VSE/nc

# stack netCDF output from NARR_composites.py in the vertical
# add isobaric level dimension
# stack with ncrcat 

# TODO: Make NARR_composites.py add isobaric level dimension and shortened field name without pressure level attached to it.


set ps=(100 125 150 175 200 225 250 275 300 350 400 450 500 550 600 650 700 750 800 825 850 875 900 925 950 975 1000)

#foreach desc (all_LTC)
foreach desc (all_LTC strong_LTC_many_tornadoes_prelandfall weak-to-intermediate_LTC_many_tornadoes_prelandfall no_or_few_tornadoes_prelandfall strong_LTC_many_tornadoes weak-to-intermediate_LTC_many_tornadoes no_or_few_tornadoes tornadoes_well_inland tornadoes_near_coast)
    foreach f (u v sh hgt temp)
    #foreach f (u)

        #foreach t (0003z 0609z 1215z 1821z)
        foreach t (00z 03z 06z 09z 12z 15z 18z 21z)
        #foreach t (0003z )

            set ofile=$desc.${f}.$t.nc
            if (-s $ofile) then
                echo "Found $ofile. Skipping."
                continue
            endif
            foreach p ($ps)
                # define new vertical dimension
                ncap2 -s 'defdim("lv_ISBL0",1);' -O $desc.$f$p.$t.nc $desc.$f$p.$t.x.nc
                # rename field name. Remove pressure part
                ncrename -v $f$p,$f -O $desc.$f$p.$t.x.nc $desc.$f$p.$t.x.nc
                # Define vertical coordinate 
                ncap2 -s "lv_ISBL0[lv_ISBL0]=$p;" -O $desc.$f$p.$t.x.nc $desc.$f$p.$t.x.nc
                # metadata
                ncap2 -s 'lv_ISBL0@units="hPa";lv_ISBL0@long_name="isobaric level";' -O $desc.$f$p.$t.x.nc $desc.$f$p.$t.x.nc
                # add vertical dimension to field. TODO: don't make a record dimension
                ncecat -u lv_ISBL0 -O $desc.$f$p.$t.x.nc $desc.$f$p.$t.x.nc
            end

            ncrcat -O\
                   $desc.${f}1000.$t.x.nc\
                   $desc.${f}975.$t.x.nc\
                   $desc.${f}950.$t.x.nc\
                   $desc.${f}925.$t.x.nc\
                   $desc.${f}900.$t.x.nc\
                   $desc.${f}875.$t.x.nc\
                   $desc.${f}850.$t.x.nc\
                   $desc.${f}825.$t.x.nc\
                   $desc.${f}800.$t.x.nc\
                   $desc.${f}750.$t.x.nc\
                   $desc.${f}700.$t.x.nc\
                   $desc.${f}650.$t.x.nc\
                   $desc.${f}600.$t.x.nc\
                   $desc.${f}550.$t.x.nc\
                   $desc.${f}500.$t.x.nc\
                   $desc.${f}450.$t.x.nc\
                   $desc.${f}400.$t.x.nc\
                   $desc.${f}350.$t.x.nc\
                   $desc.${f}300.$t.x.nc\
                   $desc.${f}275.$t.x.nc\
                   $desc.${f}250.$t.x.nc\
                   $desc.${f}225.$t.x.nc\
                   $desc.${f}200.$t.x.nc\
                   $desc.${f}175.$t.x.nc\
                   $desc.${f}150.$t.x.nc\
                   $desc.${f}125.$t.x.nc\
                   $desc.${f}100.$t.x.nc\
                   $ofile
            echo $ofile
        end
    end
end
