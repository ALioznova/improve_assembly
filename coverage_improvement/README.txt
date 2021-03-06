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
- Ragout_workdir/unused_contigs_workdir/scaffolds.fasta
    New scaffolds.
- Ragout_workdir/unused_contigs_workdir/scaffolds_as_contigs.txt
    Contains result of the tool: for every scaffold contains 
    contigs ids and gaps between them, also for each contig
    contains a flag indicating whether it was inserted by the
    tool.
- Ragout_workdir/unused_contigs_workdir/unused_contigs.fasta
     Contains contigs which weren't used by Ragout in fasta
     format.
- Ragout_workdir/unused_contigs_workdir/<ref_name>_aligned.sam
    Contains alignment of unused contigs to reference genome
    ref_name.
The rest files are available in debug mode.
- Ragout_workdir/unused_contigs_workdir/scaffolds_new.links
    File with new links for evaluation.
- Ragout_workdir/unused_contigs_workdir/new_contigs_coords.txt
    Coords of inserted contigs.
- Ragout_workdir/unused_contigs_workdir/contigs_between_blocks.txt
    For every used block size contains list of block pairs
    with contigs names and distances to this blocks.
- Ragout_workdir/unused_contigs_workdir/contigs_coords_info.txt
    Contains information about scaffold name, contig name and
    approximate coordinates for insertion.
- Ragout_workdir/unused_contigs_workdir/scaffolds_as_blocks.txt
    For every block size contains Ragout scaffolds (input) as
    list of block ids and coords.
- Ragout_workdir/unused_contigs_workdir/new_contigs_coords.txt
    Contains information about new coordinates of every contig.
    It contains name of contig, coordinat of begin and length
    of contig.


=======================================================================================
          Usage:
=======================================================================================
usage: insert_unused_contigs.py -b <path to bwa> -w <path to ragout workdir> -r <ragout recipe>
                                [-h] [-d]

Required arguments:
  -b BWA, --bwa BWA     
                        path to bwa
  -w WORKDIR, --workdir WORKDIR
                        ragout workdir
  -r RECIPE, --recipe RECIPE
                        ragout recipe

Optional arguments:
  -d, --debug
                        add debugging files to output
  -h, --help            
                        show help message and exit

