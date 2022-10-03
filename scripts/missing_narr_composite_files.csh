#!/bin/csh

# look for missing PNG files
# list qsub jobs needed to recreate them


cd /glade/scratch/ahijevyc/trier/VSE

set b=b.csh
if (-e $b) rm $b
foreach desc (strong_LTC_many_tornadoes \
              weak-to-intermediate_LTC_many_tornadoes \
              no_or_few_tornadoes \
              strong_LTC_many_tornadoes_prelandfall \
              weak-to-intermediate_LTC_many_tornadoes_prelandfall \
              no_or_few_tornadoes_prelandfall \
              tornadoes_well_inland \
              tornadoes_near_coast \
              all_LTC \
              all_LTC_prelandfall \
              )

    foreach hh (03060912 15182100 0003060912151821)
        foreach filllinebarb (`cat /glade/scratch/ahijevyc/vortexsoutheast/scripts/filllinebarb.txt`)
            set split=($filllinebarb:as,/, ,)
            set fill=$split[1]
            set line=$split[2]
            set barb=$split[3]

            set line="$line."
            if ($line == _.) set line=""
                    
            set barb="$barb."
            if ($barb == _.) set barb=""
                    
            ls $desc.$fill.$line$barb${hh}z.png > /dev/null
            if ($status != 0) echo qsub /glade/scratch/ahijevyc/temp/$fill.$line$barb$hh.pbs >> $b
        end
        foreach f (u v sh hgt temp)
            continue # skip these
            foreach p (100 125 150 175 200 225 250 275 300 350 400 450 500 550 600 650 700 750 800 825 850 875 900 925 950 975 1000)
                set f=$f$p
                ls $desc.$f.${hh}z.nc > /dev/null
                if ($status != 0) echo qsub /glade/scratch/ahijevyc/temp/$f.$hh.sbatch >> $b
            end
        end
    end
end

sort $b | uniq

