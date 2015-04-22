=======================================================================================
          About:
=======================================================================================
This is a tool for inserting unused contigs for reference-assisted assembly.
It works with Ragout: https://github.com/fenderglass/Ragout


=======================================================================================
          Input data:
=======================================================================================
- Ragout workdir
- Ragout recipe


=======================================================================================
          Output:
=======================================================================================
- Ragout_workdir/unused_contigs_workdir/contigs_between_blocks.txt
    For every used block size contains list of block pairs
    with contigs names and distances to this blocks.
- Ragout_workdir/unused_contigs_workdir/contigs_coords.txt
    Contains information about scaffold name, contig name and
    approximate coordinates for insertion.


=======================================================================================
          Usage:
=======================================================================================
usage: insert_unused_contigs.py -b <path to bwa> -w <path to ragout workdir> -r <ragout recipe>

Required arguments:
  -b BWA, --bwa BWA     
                        path to bwa
  -w WORKDIR, --workdir WORKDIR
                        ragout workdir
  -r RECIPE, --recipe RECIPE
                        ragout recipe

Optional arguments:
  -h, --help            
                        show help message and exit

