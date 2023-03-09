#!/bin/csh
module load nco

cd /glade/scratch/ahijevyc/trier/VSE/nc

# stack netCDF output from NARR_composites.py in the vertical
# add isobaric level dimension
# stack with ncrcat 

# TODO: Make NARR_composites.py add isobaric level dimension and shortened field name without pressure level attached to it.


set ps=(100 125 150 175 200 225 250 275 300 350 400 450 500 550 600 650 700 750 800 825 850 875 900 925 950 975 1000)

#foreach desc (all_LTC)
foreach desc (strong_LTC_many_tornadoes_prelandfall \
              weak-to-intermediate_LTC_many_tornadoes_prelandfall \
              no_or_few_tornadoes_prelandfall \
              strong_LTC_many_tornadoes \
              weak-to-intermediate_LTC_many_tornadoes \
              no_or_few_tornadoes \
              tornadoes_well_inland \
              tornadoes_near_coast \
              all_LTC \
              all_LTC_prelandfall )
    foreach f (u v sh hgt temp)
    #foreach f (u)

        foreach t (0003z 0609z 1215z 1821z)
        #foreach t (00z 03z 06z 09z 12z 15z 18z 21z)
        #foreach t (0003z )

            set ofile=$desc.${f}.$t.nc
            if (-s $ofile) then
                echo "Found $ofile. Skipping."
                continue
            endif
            foreach p ($ps)
                # define new vertical dimension
                ncap2 -s 'defdim("lv_ISBL0",1);' -O $desc.$f${p}hPa.$t.nc $desc.$f${p}hPa.$t.x.nc
                # rename field name. Remove pressure part
                ncrename -v $f${p}hPa,$f -v wind${p}hPa,wind -O $desc.$f${p}hPa.$t.x.nc $desc.$f${p}hPa.$t.x.nc
                # Define vertical coordinate 
                ncap2 -s "lv_ISBL0[lv_ISBL0]=$p;" -O $desc.$f${p}hPa.$t.x.nc $desc.$f${p}hPa.$t.x.nc
                # metadata
                ncap2 -s 'lv_ISBL0@units="hPa";lv_ISBL0@long_name="isobaric level";' -O $desc.$f${p}hPa.$t.x.nc $desc.$f${p}hPa.$t.x.nc
                # add vertical dimension to field. TODO: don't make a record dimension
                ncecat -u lv_ISBL0 -O $desc.$f${p}hPa.$t.x.nc $desc.$f${p}hPa.$t.x.nc
            end

            ncrcat -x -v narr_files -O\
                   $desc.${f}1000hPa.$t.x.nc\
                   $desc.${f}975hPa.$t.x.nc\
                   $desc.${f}950hPa.$t.x.nc\
                   $desc.${f}925hPa.$t.x.nc\
                   $desc.${f}900hPa.$t.x.nc\
                   $desc.${f}875hPa.$t.x.nc\
                   $desc.${f}850hPa.$t.x.nc\
                   $desc.${f}825hPa.$t.x.nc\
                   $desc.${f}800hPa.$t.x.nc\
                   $desc.${f}750hPa.$t.x.nc\
                   $desc.${f}700hPa.$t.x.nc\
                   $desc.${f}650hPa.$t.x.nc\
                   $desc.${f}600hPa.$t.x.nc\
                   $desc.${f}550hPa.$t.x.nc\
                   $desc.${f}500hPa.$t.x.nc\
                   $desc.${f}450hPa.$t.x.nc\
                   $desc.${f}400hPa.$t.x.nc\
                   $desc.${f}350hPa.$t.x.nc\
                   $desc.${f}300hPa.$t.x.nc\
                   $desc.${f}275hPa.$t.x.nc\
                   $desc.${f}250hPa.$t.x.nc\
                   $desc.${f}225hPa.$t.x.nc\
                   $desc.${f}200hPa.$t.x.nc\
                   $desc.${f}175hPa.$t.x.nc\
                   $desc.${f}150hPa.$t.x.nc\
                   $desc.${f}125hPa.$t.x.nc\
                   $desc.${f}100hPa.$t.x.nc\
                   $ofile
            if ($status == 0) then
                echo $ofile
                rm $desc.${f}[0-9]*[0-9]hPa.$t.x.nc
            else
                echo "error no $ofile made"
            endif
        end
    end
end
