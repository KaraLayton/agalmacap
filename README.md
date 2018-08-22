# agalmacap

#########TODO#########
<!-- 
Usage:
```
python3 TxtmFishing.py --param parameters.txt 
```

Pulls orthologous genes from transcriptomes and trims them to target exons. Output files are formated to be used as a reference 'target file' for HybPiper.
Picks the best e-value by iterating through blast e-values until a minimum number of transcriptomes with hits is reached. These are specified in the parameter file.



## Input Files:

1) Transcriptome files formated as a genbank query databases in a folder.
2) Full CDS of target genes as AA.fa
3) Concatenated target exons of genes as AA.fa
4) Text file of transcriptome names.
5) Text file of a list of the gene names.
6) Parameter file. See Param-default.txt

## Required:
Python3 Libraries:
+ Pool
+ BioPython
+ regex

Programs in Path:
+ tBlastx
+ Exonerate
+ Mafft
 -->