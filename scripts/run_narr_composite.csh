#!/bin/csh

# Composite NARR for many fields, times, and storms

set repo=/glade/scratch/ahijevyc/vortexsoutheast

# Storm and time lists written to this directory
set stormlistdir=$repo/stormtimelist

# plots written to this directory
set odir=$repo/output/composites

# option to normalize radial dimension
set normalize_str=""
#set normalize_str="--normalize_by Vt500km"

# Hours in time window
set time_window_hours=24
set anchor_hour=0

set CM1=../CM1_input_fields.txt
if (-e $CM1) rm $CM1
foreach f (u v sh hgt temp)
    foreach p (100 125 150 175 200 225 250 275 300 350 400 450 500 550 600 650 700 750 800 825 850 875 900 925 950 975 1000)
        echo $f${p}hPa/_/wind${p}hPa >> $CM1
    end
end

# If you normalize the radial dimension (--normalize option), NARR_composite renames output. 
# You don't need to change --ofile and --netcdf arguments yourself. 


foreach filllinebarb (`cat $repo/scripts/filllinebarb.txt` `cat $CM1`)

    
    # substitution operator. Turn slashes to spaces so array can be set.
    set split=($filllinebarb:as,/, ,)
    set fill=$split[1]
    set line=$split[2]
    set barb=$split[3]

    if ($line == _) then
        set line=""
        set lineargs=""
    else
        # line contour argumements, and add line to output file name
        if ($line == sbcape) set lineargs="--line $line --clev 500 1000 2000" 
        if ($line =~ shr10m_*) set lineargs="--line $line --clev 8 12 16 20" 
        if ($line == srh) set lineargs="--line $line --clev 150 200 300" 
        set line="$line."
    endif

    if ($barb == _) then
        set barb=""
        set barbargs=""
    else
        # wind barb/quiver argumements, and to output file name
        set barbargs="--barb $barb" # Can't make one-line because 2nd element not defined.
        if ($barb =~ shr*) set barbargs="--quiver $barb" # If shear, use quiver, not barb
        set barb="$barb."
    endif

    # Loop through diurnal cycle
    @ i_time_window = 0
    foreach h (`seq 0 $time_window_hours 21`)
        # add anchor_hour and take modulo 24
        set h=`echo "($h+$anchor_hour)%24" |bc`
        # hours UTC after applying offset from 00z
        # used in output filename
        set hh="" # string
        set diurnal_hours=""
        set h_end=`echo $h + $time_window_hours -3|bc`
        foreach n (`seq --format='%02.0f' $h 3 $h_end`)
            set n=`echo "$n % 24" |bc`
            if ($n<10) set n=0$n
            set hh=$hh$n
            set diurnal_hours="$diurnal_hours $n"
        end
        set batch=$TMPDIR/$fill.$line$barb$hh.pbs

        cat<<END>$batch
#!/bin/csh
#PBS -A NMMM0021
#PBS -N $fill$hh$line
#PBS -e $TMPDIR/$fill.$line$hh.err
#PBS -o $TMPDIR/$fill.$line$hh.out
#PBS -q casper
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=4GB

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR

module reset
module load ncl
conda activate
END
        foreach category (../categories/*.txt)
            set desc=`basename $category txt`
            set a=$TMPDIR/$desc${hh}z.txt # Time-saving measure, but If you update stormlists, remove old versions of this.
            # Not sure if diurnal_hours.py is any better or faster than get_diurnal_hour.py
            if (! -s $a) python $repo/scripts/diurnal_hours.py $category $anchor_hour $time_window_hours $i_time_window
            echo python $repo/scripts/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc$fill.${hh}z.nc >> $batch
        end
        @ i_time_window++


        echo qsub $batch
    end
    
end
