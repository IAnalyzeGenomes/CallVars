This folder should contains below human reference genome files.

genome.dict
genome.fa
genome.fa.amb
genome.fa.ann
genome.fa.bwt
genome.fa.fai
genome.fa.pac
genome.fa.sa
GenomeSize.xml

These reference genome files can be downloaded in 2bit format using below link.
http://hgdownload.cse.ucsc.edu/gbdb/hg19/
The utility program, twoBitToFa (available from the kent src tree), can be used to extract .fa file(s) from this file.  
A pre-compiled version of this command line tool can be found at: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
Converting 2bit to fa --> ./twoBitToFa hg19.2bit hg19.fa
renaming --> mv hg19.fa genome.fa
indexing .fa file --> bwa index genome.fa
Creating .dict file --> gatk CreateSequenceDictionary –R genome.fa
