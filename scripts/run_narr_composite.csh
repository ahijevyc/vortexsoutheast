#!/bin/csh

# Composite NARR for many fields, times, and storms


# Storm and time lists written to this directory
set stormlistdir=/glade/work/ahijevyc/share/VSE

# plots written to this directory
set odir=/glade/scratch/ahijevyc/trier/VSE

# Best tracks - atcf format
set atcf=/glade/work/ahijevyc/atcf/archive

# option to normalize radial dimension
set normalize_str=""
#set normalize_str="--normalize_by Vt500km"

# Hours in time window
set time_window_hours=24
set anchor_hour=0

module load conda
module load ncl
conda activate npl

set CM1=CM1_input_fields.txt
if (-e $CM1) rm $CM1
foreach f (u v sh hgt temp)
    foreach p (100 125 150 175 200 225 250 275 300 350 400 450 500 550 600 650 700 750 800 825 850 875 900 925 950 975 1000)
        echo $f$p/wind$p >> $CM1
    end
end

# If you normalize the radial dimension (--normalize option), NARR_composite renames output. 
# You don't need to change --ofile and --netcdf arguments yourself. 


foreach fillbarb (div250hPa/wind250hPa div925hPa/wind925hPa mslp/wind10m mlcape/wind10m mlcinh/wind10m pwat/wind10m rh925hPa/wind925hPa rh850hPa/wind850hPa rh700hPa/wind700hPa rh500hPa/wind500hPa rh_0deg/wind500hPa shr10m_500hPa/shr10m_500hPa shr10m_700hPa/shr10m_700hPa shr10m_900hPa/shr10m_900hPa speed10m/wind10m speed700hPa/wind700hPa speed500hPa/wind500hPa sh2/wind10m srh stp theta2/wind10m thetae2/wind10m thetasfc/wind10m vort250hPa/wind250hPa vort925hPa/wind925hPa vvel700hPa/wind700hPa vvel850hPa/wind850hPa speed10m/shr10m_700hPa tctp2014/shr10m_3000m tctp/shr10m_3000m tctp/shr10m_1000m scp/shr10m_700hPa scp/shr10m_900hPa shr10m_900hPa/shr10m_900hPa shr10m_700hPa/shr10m_700hPa) # torn symbols in magenta?
#foreach fillbarb (tctp2014/shr10m_3000m tctp/shr10m_3000m tctp/shr10m_1000m)
#foreach fillbarb (vvel700hPa/wind700hPa) # torn symbols in yellow
#foreach fillbarb (speed10m/shr10m_700hPa tctp/shr10m_700hPa tctp/shr10m_900hPa scp/shr10m_700hPa scp/shr10m_900hPa shr10m_900hPa/shr10m_900hPa shr10m_700hPa/shr10m_700hPa) # torn symbols in magenta?
#foreach fillbarb (`cat $CM1`)
#foreach fillbarb (mlcinh/wind10m mlcape/wind10m srh/shr10m_700hPa) # need mlcape for NARR_composite_skewt.py

    
    # substitution operator. Turn slash to space so array can be set.
    set split=($fillbarb:s,/, ,)
    set fill=$split[1]

    set lineargs=""
    set line=""
    # line contour argumements, and add line to output file name
    if ($fill =~ shr10m_* | $fill =~ vvel700*) then
        # manually switch comments below and rerun to get both sbcape and shear
        set line="sbcape."
        set lineargs="--line sbcape --clev 500 1000 2000" 
        #set line="shr10m_700hPa."
        #set lineargs="--line shr10m_700hPa --clev 8 12 16 20" 
    endif
    if ($fill == "scp" | $fill == "tctp") then
        set line="srh."
        set lineargs="--line srh --clev 150 200 300" 
    endif

    set barbargs=""
    set barb=""
    # wind barb/quiver argumements, and to output file name
    if ($#split >= 2) then
        set barb="$split[2]."
        set barbargs="--barb $split[2]" # Can't make one-line because 2nd element not defined.
        if ($split[2] =~ shr*) set barbargs="--quiver $split[2]" # If shear, use quiver, not barb
    endif

    # Loop through diurnal cycle
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
module load conda
module load ncl
conda activate npl
END

        foreach desc (all_LTC \
                      strong_LTC_many_tornadoes_prelandfall \
                      weak-to-intermediate_LTC_many_tornadoes_prelandfall \
                      no_or_few_tornadoes_prelandfall \
                      strong_LTC_many_tornadoes \
                      weak-to-intermediate_LTC_many_tornadoes \
                      no_or_few_tornadoes \
                      tornadoes_well_inland \
                      tornadoes_near_coast \
                          )

            set a=$stormlistdir/$desc.${hh}z.txt
            if (! -s $a) then
                python /glade/scratch/ahijevyc/vortexsoutheast/scripts/get_diurnal_hour.py /glade/scratch/ahijevyc/vortexsoutheast/categories/$desc.txt $diurnal_hours > $a 
            endif
            echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch
        end



        echo qsub $batch \; sleep 3
    end

end
