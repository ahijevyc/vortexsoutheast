#!/bin/csh

foreach f (strong_LTC_many_tornadoes weak-to-intermediate_LTC_many_tornadoes some_tornadoes no_or_few_tornadoes all_LTC)
    foreach p ( _prelandfall "" )
         python environmental_windshear.py $TMPDIR/$f$p.0003060912151821z.txt
     end
 end
