#!/usr/bin/env python3

import re
from pathlib import Path
import sys

usage = """
Usage: VulgarityFilter.py --in exonerate_outfile.txt

Parses the output of exonerate with: --ryo '%qi\t%pi\t%qas\t%V\tEND\n'
output: inputname.fa of a fasta with each exon listed as a seperate sequence
example exonerate command:
exonerate --model est2genome Test44.fasta /Users/josec/Desktop/NudiSilicoTest/Exonerate/acl_ref_AplCal3.0_chrUn.fa -Q DNA -T DNA --showvulgar F --showalignment F --percent 90 --verbose 0 --ryo '%qi\t%pi\t%qas\t%V\tEND\n' --fsmmemory 20G --bestn 1 > exonerate_outfile.txt


This chops up the sequences in the exonerate_outfile into exons and saves them
as a fasta with headers:
>Gene_exonNumber
EXONSEQUENCE

"""


def openfile(infile):
    # Read exonerate output and make list of output for each target (gene)
    # The results are stored in target_list
    with open(infile, 'Ur') as inhandle:
        target_list = inhandle.read().split('END')
    return target_list


def parser(target):
    # This takes the results and parses them into the target_dict
    target_dict = {}
    tlist = [x.replace('\n', '') for x in target.split('\t')]
    target_dict['Target'] = tlist[0]
    target_dict['Percent'] = tlist[1]
    target_dict['cds'] = tlist[2]
    # Extract exon lengths from vulgar output, then add to target_dict
    # Each exon length is preceded by 'M '
    # Exon lengths are all in sequential list.
    vulgar_raw = tlist[3]
    vlist = re.findall(r'M\s\d*', vulgar_raw)
    vulgar = [int(x.replace('M ', '')) for x in vlist]
    target_dict['Vulgar'] = vulgar
    return target_dict


def splitter(cds, vlist):
    # Cuts up the cds into exon sequences and stores them in a list
    # This is done from the vulgar list of exon lengths
    # cds is the sequence of concatenated exons
    counter = 0
    exonseq_list = []
    for v in vlist:
        exonseq_list.append(cds[counter:counter + v])
        counter += v
    return exonseq_list


def writer(target_dict, gene_out_folder):
    # Writes the exon sequences into the outfile.fa
    # Each exon sequence is named by the order in the target (gene)
    outstring = ""
    vlist = target_dict['Vulgar']
    cds = target_dict['cds']
    try:
        gene_raw = target_dict['Target'].split("cds_")[1].split("_")
        gene = f"{gene_raw[0]}_{gene_raw[1]}"
    except KeyError:
        gene = target_dict['Target']
    exonseq_list = splitter(cds, vlist)
    out_path = str(gene_out_folder/f"{gene}.fas")
    with open(out_path, 'w') as out_handle:
        for exon_number, exonseq in enumerate(exonseq_list):
            outstring += f">{gene}_{exon_number}\n{exonseq}\n"
        out_handle.write(outstring)
    return outstring


def vulgarity_filter(infile,gene_out_folder):
    target_list = openfile(infile)
    for target in target_list:
        # check if target file is blank
        if len(target) > 2:
            target_dict = parser(target)
            # Filter to remove low ID hits
            if float(target_dict['Percent']) < 98:
                continue
            outstring = writer(target_dict, Path(gene_out_folder))
    return outstring


def main():
    # Engine for the program.
    args = sys.argv[1:]
    # args=['--in','/Users/josec/Desktop/NudiPreBait/NudiSilicoTest/test444out.txt'] #For testing
    # Print usage if no input is put in by user
    if not args:
        print(usage)
        sys.exit(1)
    if args[0] == '--in':
        infile = args[1]
    vulgarity_filter(infile)

    return


if __name__ == '__main__':
    main()
