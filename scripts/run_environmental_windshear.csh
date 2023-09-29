#!/bin/csh

foreach f (strong_LTC_many_tornadoes weak-to-intermediate_LTC_many_tornadoes some_tornadoes)
    foreach p ( _prelandfall "" )
         python environmental_windshear.py $TMPDIR/$f$p.0003060912151821z.txt
     end
end
 
# remove 3 out of 20 cases from no_or_few_tornadoes that are not in corresponding prelandfall list.
python environmental_windshear.py $TMPDIR/no_or_few_tornadoes_prelandfall.0003060912151821z.txt
cat $TMPDIR/no_or_few_tornadoes.0003060912151821z.txt|\
    grep -v "Bertha 2002"|\
    grep -v "Bertha 2020"|\
    grep -v "Tammy 2005" > $TMPDIR/same17
python environmental_windshear.py $TMPDIR/same17

python environmental_windshear.py $TMPDIR/all_LTC_prelandfall.0003060912151821z.txt
cat $TMPDIR/all_LTC.0003060912151821z.txt|\
    grep -v "Bertha 2002"|\
    grep -v "Bertha 2020"|\
    grep -v "Tammy 2005" > $TMPDIR/same69
python environmental_windshear.py $TMPDIR/same69
