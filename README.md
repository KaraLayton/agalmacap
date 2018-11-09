# agalmacap


## Usage:
```
python3 agalmacap.py --param parameterfile.txt
```

## Description:
This pipeline automates the steps required to generate exon alignments from Agalma amino acid (AA) data. These exon alignments are DNA sequences of homologous genes identified from an Agalma transcriptome search that are cut at the intron sites of the closest reference genome. 



## Input Files:

+ Agalma supermatrix and partition file. Output of the Agalma phylogeny pipeline.
+ Transcriptome assemblies of the taxa included in the Agalma analysis. Note all of the transcriptome assembly files must have the name of the sequences in the supermatrix. The file endings must be '.fas'
+ Parameter file. Place all the values before the '#'. See Param-AgalmaCapBlank.txt 
+ Refseq genome files in a folder. The genomic.fna, cds_from_genomic.fna files and translated_CDS.faa are required. Download via FTP from here: ftp.ncbi.nih.gov/genomes/refseq/ Example:
```
curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna.gz
curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_genomic.fna.gz
curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_translated_cds.faa.gz
```

## Required:
Python3 Libraries:
+ Pool
+ BioPython
+ regex
+ subprocess
+ Pathlib
+ textwrap
+ multiprocessing

Programs in Path:
+ tBlastx
+ Exonerate
+ Mafft
