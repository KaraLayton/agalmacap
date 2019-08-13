# agalmacap


## Usage:
```
python3 agalmacap.py --param parameterfile.txt
```

## Description:
This pipeline automates the steps required to generate exon alignments from DNA or AA gene alignment data. These exon alignments are DNA sequences of homologous genes cut at the intron sites of the provided reference genome. These exon alignments can aid in more efficient exon-capture bait design where baits do not span intron boundaries. This workflow is non-model organism friendly and relatedness of the reference genome is relaxed (ie same family-order). Agalmacap was innitially designed to take the Agalma output files to generate intron alignments, but has since been upgraded to allow more input file formats. 


## Input Files:

+ DNA or AA alignments of homologous genes from a group of transcriptomes. Alternatively an Agalma supermatrix and partition file can be used as inputs. If the reference genome is distant from the ingroups, working with AA data may provide better results.
+ Transcriptome assemblies of the taxa included in the gene alignments or Agalma analysis. All of the transcriptome assembly files must have the name of the sequences in the gene alignments or supermatrix. The file endings must be '.fas'
+ Parameter file. Type all the values to the left of the '#'. See Param-AgalmaCapBlank.txt 
+ Refseq genome files in a folder. The genomic.fna, cds_from_genomic.fna files and translated_CDS.faa are required. Download via FTP from here: ftp.ncbi.nih.gov/genomes/refseq/ Example:
```
curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_cds_from_genomic.fna.gz
curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_genomic.fna.gz
curl -O ftp.ncbi.nih.gov/genomes/refseq/invertebrate/Aplysia_californica/latest_assembly_versions/GCF_000002075.1_AplCal3.0/GCF_000002075.1_AplCal3.0_translated_cds.faa.gz
```

## Required:
Python3 Libraries:

+ BioPython
+ regex
+ subprocess
+ Pathlib
+ textwrap
+ multiprocessing

Programs in Path:
+ tBlastx (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
+ Exonerate (https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
+ Mafft (https://mafft.cbrc.jp/alignment/software/)
