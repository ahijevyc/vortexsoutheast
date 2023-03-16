#!/bin/csh
conda activate 
set repo=/glade/scratch/ahijevyc/vortexsoutheast
cd $repo/output
foreach f (../categories/*.txt)
    python ../scripts/track_torn_composite.py $f 
    python ../scripts/track_torn_composite.py $f --onetrackcolor
end
