############################################################
BED tools
############################################################

intersectBed:
	Only for Bed format. Other format should be added.

bedsAnnotation: (need update for tab)
	Only for Bed format. Tab format should be added. Input might be bowtie format.

############################################################
Hi-C tools
############################################################

endsMapAbility:
	Ends map ability of the Resitriction Enzyme fragments. The mapability is 1 over # of best hits.

genomeToREFrags:
	Split genome file into Resitriction Enzyme fragments.

seqToContact:
	Bowtie mapped files to RE fragment contacts.	

############################################################
Fastq tools
############################################################

fastqFilter:
	Filter fastq format files according to average Q scores.

fastqTrimmer:
	Trim fastq format files according to average Q scores in given window size.




