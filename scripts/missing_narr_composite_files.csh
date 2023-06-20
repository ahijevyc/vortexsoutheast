#!/bin/csh

set repo=/glade/scratch/ahijevyc/vortexsoutheast

# look for missing PNG files
# list qsub jobs needed to recreate them


cd $repo/output/composites

set b=b.csh
if (-e $b) rm $b
foreach desc ($repo/categories/*.txt)
    set desc=`basename $desc .txt`
    #if ( $desc !~ strong* ) continue
    foreach hh (0003 0609 1215 1821 03060912 15182100 0003060912151821)
        foreach filllinebarb (`cat $repo/scripts/filllinebarb.txt` `cat $repo/CM1_input_fields.txt`)
            set split=($filllinebarb:as,/, ,)
            set fill=$split[1]
            set line=$split[2]
            set barb=$split[3]

            set line="$line."
            if ($line == _.) set line=""
                    
            set barb="$barb."
            if ($barb == _.) set barb=""
                   
            # Look for PNG files 
            ls $desc.$fill.$line$barb${hh}z.png > /dev/null
            if ($status != 0) echo qsub /glade/scratch/ahijevyc/temp/$fill.$line$barb$hh.pbs >> $b

            # look for netCDF files too.
            ls nc/$desc.$fill.${hh}z.nc > /dev/null
            if ($status != 0) echo qsub /glade/scratch/ahijevyc/temp/$fill.$line$barb$hh.pbs >> $b
        end
    end
end

sort $b | uniq

