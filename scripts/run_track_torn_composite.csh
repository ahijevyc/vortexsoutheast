#!/bin/csh

conda activate 

set repo=/glade/scratch/ahijevyc/vortexsoutheast
cd $repo/output

# regular and rejected tracks
foreach f (../categories/*.txt ../categories/rejects/*.txt)
    python ../scripts/track_torn_composite.py $f --nolegend
    python ../scripts/track_torn_composite.py $f --onetrackcolor --ofile `basename $f txt`colorbytrack.legend.ps
    python ../scripts/track_torn_composite.py $f --nolegend --onetrackcolor --ofile `basename $f txt`colorbytrack.ps
end
