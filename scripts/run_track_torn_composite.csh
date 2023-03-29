#!/bin/csh

conda activate 

set repo=/glade/scratch/ahijevyc/vortexsoutheast
cd $repo/output

# regular and rejected tracks
foreach f (../categories/*.txt ../categories/rejects/*.txt)
    python ../scripts/track_torn_composite.py $f 
    python ../scripts/track_torn_composite.py $f --onetrackcolor
end
