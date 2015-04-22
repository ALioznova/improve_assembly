=======================================================================================
          About:
=======================================================================================
This is a tool for misassembly detection without reference.


=======================================================================================
          System requirements:
=======================================================================================
This tool can be run on Unix or Mac platform.
You need to install:
- Python 2.6 or 2.7
(visit https://www.python.org/download/)
- Biopython
(visit http://biopython.org/wiki/Download)
and one of two aligners: 
- bwa (preferable) 
(visit http://sourceforge.net/projects/bio-bwa/files/)
or
- bowtie2
(visit http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2)


=======================================================================================
          Input data:
=======================================================================================
- file with contigs (format: *.fasta or *.fa)
- files with right and left reads (format: *.fastq)


=======================================================================================
          Output:
=======================================================================================
- out_dir/result.txt
    The main file with results. Contains information about 
    positions where large number of reads ends and begins, 
    fragments with large insert size, 
    fragments with low or high coverage.
- out_dir/beg_and_end_reads.txt
    Detailed information about positions 
    where large number of reads ends and begins.
- out_dir/large_insert_size.txt
    Detailed information about fragments with large insert size.
- out_dir/atypical_coverage.txt
    Detailed information about intervals with low or high coverage.


=======================================================================================
          Usage:
=======================================================================================
usage: find_misassemblies.py --out_dir OUT_DIR --contigs CONTIGS --reads_1 READS_1 --reads_2 READS_2 
                [--bowtie BOWTIE] [--bwa BWA] [--aligner {bowtie,bwa}]
                [-h] [--version]

Required arguments:
  --out_dir OUT_DIR, -o OUT_DIR
                        name of output directory
  --contigs CONTIGS, -c CONTIGS
                        file with contigs
  --reads_1 READS_1, -r1 READS_1
                        file with left reads
  --reads_2 READS_2, -r2 READS_2
                        file with right reads

Optional arguments:
  -h, --help            
                        show help message and exit
  --version, -v         
                        show program's version number and exit
  --bowtie BOWTIE, -bow BOWTIE
                        path to bowtie-2 dir
  --bwa BWA, -bwa BWA   
                        path to bwa dir
  --aligner {bowtie,bwa}, -a {bowtie,bwa}
                        method to align


