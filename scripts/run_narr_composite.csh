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
set time_window_hours=6
set extra_time_window_hours=`expr $time_window_hours - 3`


set CM1=CM1_input_fields.txt
if (-e $CM1) rm $CM1
foreach f (u v sh hgt temp)
    foreach p (100 125 150 175 200 225 250 275 300 350 400 450 500 550 600 650 700 750 800 825 850 875 900 925 950 975 1000)
        echo $f$p/wind$p >> $CM1
    end
end

# If you normalize the radial dimension (--normalize option), NARR_composite renames output. 
# You don't need to change --ofile and --netcdf arguments yourself. 


#foreach fillbarb (div250/wind250 div925/wind925 mslp/wind10m mlcape/wind10m mlcinh/wind10m pwat/wind10m rh925/wind925 rh850/wind850 rh700/wind700 rh500/wind500 rh_0deg/wind500 shr10m_500hPa/shr10m_500hPa shr10m_700hPa/shr10m_700hPa shr10m_900hPa/shr10m_900hPa speed10m/wind10m speed700/wind700 speed500/wind500 sh2/wind10m srh stp theta2/wind10m thetae2/wind10m thetasfc/wind10m vort250/wind250 vort925/wind925 vvel700/wind700 vvel850/wind850 speed10m/shr10m_700hPa tctp/shr10m_3000m tctp/shr10m_1000m scp/shr10m_700hPa scp/shr10m_900hPa shr10m_900hPa/shr10m_900hPa shr10m_700hPa/shr10m_700hPa) # torn symbols in magenta?
foreach fillbarb (tctp2014/shr10m_3000m tctp/shr10m_3000m tctp/shr10m_1000m)
#foreach fillbarb (vvel700/wind700) # torn symbols in yellow
#foreach fillbarb (speed10m/shr10m_700hPa tctp/shr10m_700hPa tctp/shr10m_900hPa scp/shr10m_700hPa scp/shr10m_900hPa shr10m_900hPa/shr10m_900hPa shr10m_700hPa/shr10m_700hPa) # torn symbols in magenta?
#foreach fillbarb (`cat $CM1`)
#foreach fillbarb (mlcinh/wind10m mlcape/wind10m srh/shr10m_700hPa) # need mlcape for NARR_composite_skewt.py

    
    # substitution operator. Turn slash to space so array can be set.
    set split=($fillbarb:s,/, ,)
    set fill=$split[1]

    set lineargs=""
    set line=""
    # line contour argumements, and add line to output file name
    if ($fill =~ shr10m_* | $fill =~ vvel700) then
        set line="sbcape."
        set lineargs="--line sbcape --clev 500 1000 2000" 
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

    # hours offset from 00z
    foreach h (`seq -12 $time_window_hours 9`)
        # hours UTC after applying offset from 00z
        # used in output filename
        # could be any date yyyymmdd
        set hh=""
        foreach n (`seq 0 3 $extra_time_window_hours`)
            set hh=$hh`date -u --date "20080601 ${h}hours +${n}hours" +%H`
        end
        set batch=$TMPDIR/$fill.$line$barb$hh.pbs

        cat<<END>$batch
#!/bin/csh
#PBS -A NMMM0021
#PBS -N $fill$hh
#PBS -e $TMPDIR/$fill.$hh.err
#PBS -o $TMPDIR/$fill.$hh.out
#PBS -q casper
#PBS -l walltime=06:00:00
#PBS -l select=1:ncpus=1:mem=4GB

setenv TMPDIR /glade/scratch/$USER/temp
mkdir -p $TMPDIR

module reset
module load python
ncar_pylib 20201220_casper
module load ncl
END


        # These are the "days-to-subtract" to accomodate different starting hours. Subtracting a whole number of days does not change the hour UTC.
        set start_at_00z=`echo "$h/24"|bc` 
        set start_at_03z=`echo "($h+21)/24"|bc`
        set start_at_06z=`echo "($h+18)/24"|bc`
        set start_at_09z=`echo "($h+15)/24"|bc`
        set start_at_12z=0
        set start_at_15z=`echo "($h+9)/24"|bc`
        set start_at_18z=`echo "($h+6)/24"|bc`
        set start_at_21z=`echo "($h+3)/24"|bc`

        # ALL - different start-hours. 
        set desc=all_LTC
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Ivan      2004 `date -u --date "20040917 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Katrina   2005 `date -u --date "20050830 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Rita      2005 `date -u --date "20050925 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Florence  2018 `date -u --date "20180917 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gustav    2008 `date -u --date "20080902 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Irma      2017 `date -u --date "20170910 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Frances   2004 `date -u --date "20040907 ${h}hours -${start_at_03z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Jeanne    2004 `date -u --date "20040927 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Ike       2008 `date -u --date "20080914 ${h}hours -${start_at_00z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Lili      2002 `date -u --date "20021004 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Lee       2011 `date -u --date "20110905 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Hermine   2016 `date -u --date "20160902 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Alberto   2006 `date -u --date "20060614 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gaston    2004 `date -u --date "20040831 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Cindy     2005 `date -u --date "20050707 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gabrielle 2001 `date -u --date "20010915 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Bill      2003 `date -u --date "20030702 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Debby     2012 `date -u --date "20120625 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Claudette 2003 `date -u --date "20030716 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Humberto  2007 `date -u --date "20070914 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Isabel    2003 `date -u --date "20030919 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gordon    2018 `date -u --date "20180906 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Edouard   2008 `date -u --date "20080806 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Dennis    1999 `date -u --date "19990905 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Hanna     2002 `date -u --date "20020915 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Arlene    2005 `date -u --date "20050612 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Grace     2003 `date -u --date "20030901 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch


        # Pre-landfall 
        #        strong_LTC_many_tornadoes_prelandfall - different start-hours. 
        set desc=strong_LTC_many_tornadoes_prelandfall
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Ivan      2004 `date -u --date "20040915 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Katrina   2005 `date -u --date "20050829 ${h}hours -${start_at_03z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Rita      2005 `date -u --date "20050923 ${h}hours -${start_at_06z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Florence  2018 `date -u --date "20180913 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gustav    2008 `date -u --date "20080901 ${h}hours -${start_at_00z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Irma      2017 `date -u --date "20170909 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Frances   2004 `date -u --date "20040904 ${h}hours -${start_at_03z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Jeanne    2004 `date -u --date "20040925 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Ike       2008 `date -u --date "20080912 ${h}hours -${start_at_00z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch

        # Pre-landfall 
        #        weak-to-intermediate_LTC_many_tornadoes_prelandfall - all sequences start with 12z
        set desc=weak-to-intermediate_LTC_many_tornadoes_prelandfall
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Lili      2002 `date -u --date "20021003 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Lee       2011 `date -u --date "20110903 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Hermine   2016 `date -u --date "20160901 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Alberto   2006 `date -u --date "20060612 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gaston    2004 `date -u --date "20040829 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Cindy     2005 `date -u --date "20050706 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gabrielle 2001 `date -u --date "20010913 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Bill      2003 `date -u --date "20030630 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Debby     2012 `date -u --date "20120624 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch


        # Pre-landfall 
        #        no_or_few_tornadoes_prelandfall
        set desc=no_or_few_tornadoes_prelandfall
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Claudette 2003 `date -u --date "20030715 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Humberto  2007 `date -u --date "20070913 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Isabel    2003 `date -u --date "20030918 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gordon    2018 `date -u --date "20180904 ${h}hours -${start_at_06z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Edouard   2008 `date -u --date "20080805 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Dennis    1999 `date -u --date "19990904 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Hanna     2002 `date -u --date "20020914 ${h}hours -${start_at_06z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Arlene    2005 `date -u --date "20050611 ${h}hours -${start_at_00z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Grace     2003 `date -u --date "20030831 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch


        # Post landfall - all sequences start with 12z
        #        strong_LTC_many_tornadoes
        set desc=strong_LTC_many_tornadoes
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Ivan      2004 `date -u --date "20040917 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Katrina   2005 `date -u --date "20050830 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Rita      2005 `date -u --date "20050925 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Florence  2018 `date -u --date "20180917 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gustav    2008 `date -u --date "20080902 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Irma      2017 `date -u --date "20170910 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Frances   2004 `date -u --date "20040907 ${h}hours -${start_at_03z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Jeanne    2004 `date -u --date "20040927 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Ike       2008 `date -u --date "20080914 ${h}hours -${start_at_00z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch

        # Post landfall - all sequences start with 12z
        #        weak-to-intermediate_LTC_many_tornadoes
        set desc=weak-to-intermediate_LTC_many_tornadoes
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Lili      2002 `date -u --date "20021004 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Lee       2011 `date -u --date "20110905 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Hermine   2016 `date -u --date "20160902 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Alberto   2006 `date -u --date "20060614 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gaston    2004 `date -u --date "20040831 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Cindy     2005 `date -u --date "20050707 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gabrielle 2001 `date -u --date "20010915 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Bill      2003 `date -u --date "20030702 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Debby     2012 `date -u --date "20120625 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch

        # Post landfall - all sequences start with 12z
        #        no_or_few_tornadoes
        set desc=no_or_few_tornadoes
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Claudette 2003 `date -u --date "20030716 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Humberto  2007 `date -u --date "20070914 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Isabel    2003 `date -u --date "20030919 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Gordon    2018 `date -u --date "20180906 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Edouard   2008 `date -u --date "20080806 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Dennis    1999 `date -u --date "19990905 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Hanna     2002 `date -u --date "20020915 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Arlene    2005 `date -u --date "20050612 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Grace     2003 `date -u --date "20030901 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch

        # Tornadoes well inland
        set desc=tornadoes_well_inland
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Jeanne  2004 `date -u --date "20040928 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Cindy   2005 `date -u --date "20050708 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Erin    2007 `date -u --date "20070819 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Fay     2008 `date -u --date "20080827 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Hermine 2010 `date -u --date "20100909 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Isaac   2012 `date -u --date "20120902 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Bill    2015 `date -u --date "20150620 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Harvey  2017 `date -u --date "20170901 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Nate    2017 `date -u --date "20171009 ${h}hours +${extra_hours}hours" +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch

        # Tornadoes near the Coast - all, but Irene start at 12z
        set desc=tornadoes_near_coast
        set a=$stormlistdir/$desc.${hh}z.txt
        if (-e $a) rm $a
        foreach extra_hours (`seq 0 3 $extra_time_window_hours`)
            echo Andrea    2013 `date -u --date "20130607 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
            echo Jeanne    2004 `date -u --date "20040927 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
            echo Gabrielle 2001 `date -u --date "20010915 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
            echo Isidore   2002 `date -u --date "20020926 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
            echo Floyd     1999 `date -u --date "19990916 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
            echo Bret      1999 `date -u --date "19990823 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
            echo Alberto   2006 `date -u --date "20060614 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
            echo Irene     2011 `date -u --date "20110828 ${h}hours -${start_at_12z}days +${extra_hours}hours" +%Y%m%d%H` >> $a
            echo Irma      2017 `date -u --date "20170911 ${h}hours +${extra_hours}hours"                      +%Y%m%d%H` >> $a
        end
        echo python ~ahijevyc/bin/NARR_composite.py --fill $fill $lineargs $barbargs $a $normalize_str --ofile $odir/$desc.$fill.$line$barb${hh}z.png --netcdf $odir/nc/$desc.$fill.${hh}z.nc >> $batch



        echo qsub $batch \; sleep 5
    end

end
