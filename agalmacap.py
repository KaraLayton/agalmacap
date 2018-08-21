#!/usr/bin/env python3

import re
import subprocess
import sys

from AgalmaAA2dna import agalmaaa2txtmdna, codex_file_reader, blaster
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from multiprocessing import Pool
from pathlib import Path
from VulgarityFilter import vulgarity_filter


def raxpartition_reader(part_file):
    """Parse Raxml style partition file.

    Args:
        part_file: partition file output by agalma. In Raxml format. Looks like this:
        "WAG, 23621 = 1-440"

    Returns:
        A python list of the start and stop coordinates of every gene
        in the supermatirx (concatenated gene alignments)
    """
    partition_strings = []
    partitions = []
    with open(part_file) as f:
        for line in f:
            line = line.strip()
            line = line.split(" = ")
            partition_strings.append(line[1])
    for part in partition_strings:
        m = re.search("(.*)-(.*)", part)
        if m:
            start = int(m.groups()[0])
            end = int(m.groups()[1])
        partitions.append((start, end))
    return partitions


def alncutter(partitions, aln_file, aln_type, genealn_outdir):
    """Divides the aln_file into sub-alignments based on positions in partition list
    Writes gene alignments to genealn_outdirectory as fasta files.
    Ignores sequences made entirely of gaps

    Args:
        partitions: A python list of the start and stop coordinates of every gene
        aln_file: Path string the alignment file to subsample partitions from
        aln_type: Biopython alignments type (eg. 'fasta' or 'phyllip relaxed')
        genealn_outdir: The path to write sub-alignments to.
    """
    records = list(SeqIO.parse(aln_file, aln_type))
    for i, (start, end) in enumerate(partitions):
        out_path = str(Path(genealn_outdir)/f"{str(Path(aln_file).stem)}_{i}.fas")
        with open(out_path, 'w') as out:
            for record in records:
                sequence = record.seq[start:end]
                # Skip seqs made only of gaps
                if len(set(sequence)) > 1:
                    out.write(f">{record.id}\n")
                    out.write(f"{record.seq[start:end]}\n")
    return


def consensus_generator(aln_folder, aln_type, consensus_file):
    """Creates a simple consensus sequence (70% majority) of each alignment
     in aln_folder and writes them to the consensus_file

    Args:
        aln_folder: The folder with alignments to make consensus seqs of
        aln_type: Biopython alignments type (eg. 'fasta' or 'phyllip relaxed')
        consensus_file: path to outfile (eg. 'consensus_seqs.fas')
    """
    alns = [x for x in Path(aln_folder).iterdir() if x.is_file() and not x.stem.startswith('.')]
    with open(Path(consensus_file), 'w') as out_file:
        for file in alns:
            aln = AlignIO.read(str(file), aln_type)
            summary_align = AlignInfo.SummaryInfo(aln)
            consensus = summary_align.dumb_consensus()
            out_file.write(f'>{file.stem}\n{consensus}\n')
    return


def run_command(command, command_out_path=None):
    """Runs command with subprocess library.
    Prints error if commmand did not run.
    if outpath given it will print output there, otherwise this function returns it. """
    command_output = subprocess.getstatusoutput(command)
    if command_output[0]:
        print(f"{command} Failed")
        return None
    if command_out_path:
        with open(command_out_path, "w") as out_handle:
            out_handle.write(command_output[1])
    return command_output[1]


def generate_partitions(substring_fasta, aln_file):
    """

    Args:
        substring_fasta:
        aln_file:

    Out:
        partitions:


    substrings must be in order
    '_cds_' must be in fasta header of substring concat sequence(cds_string) """

    for record in SeqIO.parse(aln_file, 'fasta'):
        if '_cds_' in record.name:
            cds_string = str(record.seq).upper()
    exon_starts = []  # list of the start index of the exons.
    for sub_record in SeqIO.parse(substring_fasta, 'fasta'):
        exon_seq = str(sub_record.seq).upper()
        last_exon_string = exon_seq  # save this for later finding the end of the gene in the aln

        # find longest substring of exonseq in cds_string by cutting one letter off
        # the end of the exon_seq until there is a match.
        # Do this step because of gaps in the alignment result in non-exact matches of the substring
        while exon_seq not in cds_string:
            exon_seq = exon_seq[:-1]
        exon_start = cds_string.index(exon_seq)
        exon_starts.append(exon_start)

    # find the end of the last index
    while last_exon_string not in cds_string:
        last_exon_string = last_exon_string[1:]
    end_of_gene = len(last_exon_string) + cds_string.index(last_exon_string) + 1
    exon_starts.append(end_of_gene)

    partitions = []
    # skip last item in exon_starts list because this is the index of the end of the gene
    for i, start in enumerate(exon_starts[:-1]):
        partitions.append((start, exon_starts[i+1]-1))
    return partitions


def intronerator(threads, genome_fasta, target_cds):
    """Runs exonerate to find exon boundry sites on target genes given a reference genome.
      The output is a "roll your own format" string from exonerate that will be parsed by XXXX"""
    exonerate_command_prefix = f'exonerate --model est2genome -q {target_cds} -t {genome_fasta} -Q DNA -T DNA --showvulgar F --showalignment F --softmasktarget T --verbose 0 --ryo \"%qi\\t%pi\\t%qas\\t%V\\tEND\\n\" --fsmmemory 20G --bestn 1 --querychunktotal {threads} --querychunkid '
    exonerate_commands = [exonerate_command_prefix+str(i+1) for i in range((int(threads)))]
    p = Pool(threads)
    introterate_out = p.map(run_command, exonerate_commands)
    return introterate_out


def fasta_subseter(subset_list, in_fasta_path, out_fasta_path):
    """Given a list of fasta headers, it will extract the seqs in that list from in_fasta_path
    and write them to out_fasta_path
    """
    input_dict = SeqIO.to_dict(SeqIO.parse(in_fasta_path, "fasta"))
    out_records = []
    keys = input_dict.keys()
    for name in subset_list:
        for key in keys:
            if name in key:
                out_records.append(input_dict[key])
    with open(out_fasta_path, "w") as out_handle:
        SeqIO.write(out_records, out_handle, 'fasta')
    return


def fetch_refseqcds(refseqcds, codex_file, outfile):
    to_fetch = codex_file_reader(codex_file).values()
    # print(to_fetch)
    fasta_subseter(subset_list=to_fetch, in_fasta_path=refseqcds, out_fasta_path=outfile)
    return


def lessgappy_maffter(fasta_file_path, out_file_path, num_threads):
    """ Runs mafft alignment on fasta file. Options result in
    an alignment with very few gaps (simmilar to muscle).
    Sequences are checked for reverse complementing
    and the output string removes the R_ prefix that mafft adds to
    sequences that were reverse complemented
    """
    command = f"mafft --reorder --adjustdirection --leavegappyregion --op 3.0 --ep 1.0 --maxiterate 1000 --retree 1 --genafpair --thread {num_threads} --quiet --auto {fasta_file_path}"
    mafft_output = run_command(command)
    if mafft_output:
        with open(out_file_path, "w") as out_handle:
                out_handle.write(mafft_output.replace("_R_", ""))
        return
    else:
        return


def cds_loci_merger(codex_file, cds_file, dnabyloci_folder, dnacdsbyloci_folder):
    """ For every line of the codex file this will pull the CDS sequence from
    cds_file and combine it with the gene file from dnabyloci_folder.
    Results are writen to dnacdsbyloci_folder

    """
    Path(dnacdsbyloci_folder).mkdir(exist_ok=True)
    to_fetch = codex_file_reader(codex_file).values()
    cds_records = SeqIO.parse(cds_file, "fasta")
    for refseqID in to_fetch:
        loci_file_path = Path(dnabyloci_folder)/f"{refseqID}.fas"
        loci_cds_file_path = Path(dnacdsbyloci_folder)/f"{refseqID}.fas"
        records_to_write = SeqIO.parse(str(loci_file_path), "fasta")
        if not records_to_write:
            # Skip entry if there are no sequences to write in the DNAbyLoci folder
            break
        else:
            for cds in cds_records:
                if refseqID in cds.name:
                    # print(refseqID,cds.name)
                    with open(loci_cds_file_path, 'w') as out_handle:
                        SeqIO.write(cds, out_handle, 'fasta')
                        SeqIO.write(records_to_write, out_handle, 'fasta')
                        break

    return


def cut_genealns_to_exonalns(codex_file, cds_exon_folder, cds_loci_folder, exonaln_foler):
    """ For each loci in the codex, this will cut the loci alignment by exons
    and save it to exonaln_folder as Loci_Exonnumber.fas
    """
    good_loci = codex_file_reader(codex_file=codex_file).values()
    for loci in good_loci:
        cds_exon_file = Path(cds_exon_folder)/f"{loci}.fas"
        loci_aln_file = Path(cds_loci_folder)/f"{loci}.fas"
        if cds_exon_file.exists() and loci_aln_file.exists():
            loci_partitions = generate_partitions(substring_fasta=str(cds_exon_file), aln_file=str(loci_aln_file))
            alncutter(loci_partitions, aln_file=str(loci_aln_file), aln_type='fasta', genealn_outdir=exonaln_foler)
    return


def exon_aln_filter(exonaln_folder, filtered_exonaln_folder, min_exon_length=0, min_taxoncov=0):
    """
        
    """
    Path(filtered_exonaln_folder).mkdir(exist_ok=True)
    exon_files = [str(exon_file) for exon_file in Path(exonaln_folder).iterdir() if '.fa' in exon_file.suffix]
    exons2write = []
    for exon_file in exon_files:
        records = list(SeqIO.parse(exon_file, 'fasta'))
        if len(records) > 0:
            if len(records) > min_taxoncov \
                        and len(records[0].seq) > min_exon_length:
                exons2write.append(exon_file)
                outfile_stem = Path(exon_file).name
                out_path = Path(filtered_exonaln_folder)/outfile_stem
                with open(out_path, 'w') as out_handle:
                    SeqIO.write(records, out_handle, 'fasta')
    print(len(exons2write))

    return


aglama_partition_file = ''
agalma_supermatrix_file = ''


# partitions = raxpartition_reader(part_file='/Users/josec/Desktop/exoncap/corals/Txtms/100.supermatrix.partition.txt')
# alncutter(partitions=partitions, aln_file='/Users/josec/Desktop/exoncap/corals/Txtms/100.supermatrix.fa', aln_type='fasta', genealn_outdir='/Users/josec/Desktop/exoncap/corals/AgalmaAAaln2')

# consensus_generator('/Users/josec/Desktop/exoncap/Mygal/alns','fasta','/Users/josec/Desktop/exoncap/Mygal/reports/test.fa')

# Genrate Codex
# blaster(threads=8, query='/Users/josec/Desktop/exoncap/Mygal/reports/test.fa',blastdb_prefix='/Users/josec/Desktop/exoncap/Mygal/RefSeq/GCF_000365465.2_Ptep_2.0_protein.faa',blast_out_path='/Users/josec/Desktop/exoncap/Mygal/blasttest.txt')

# # DONT run agalmaaa2txtmdna multiple times because it appends to sequence files in DNAbyLoci, not overwrites
# agalmaaa2txtmdna(codex_file='/Users/josec/Desktop/exoncap/Mygal/blasttest.txt',alnaa_folder='/Users/josec/Desktop/exoncap/Mygal/alns/',txtm_folder='/Users/josec/Desktop/exoncap/Mygal/Assemblies/')



# # Generate Gene alignments with refseq and minimal gaps
# fetch_refseqcds(refseqcds='/Users/josec/Desktop/exoncap/Mygal/RefSeq/GCF_000365465.2_Ptep_2.0_cds_from_genomic.fna', codex_file='/Users/josec/Desktop/exoncap/Mygal/blasttest.txt', outfile='/Users/josec/Desktop/exoncap/Mygal/refseqcds_loci.fas')
# cds_loci_merger(codex_file='/Users/josec/Desktop/exoncap/Mygal/blasttest.txt', cds_file='/Users/josec/Desktop/exoncap/Mygal/refseqcds_loci.fas', dnabyloci_folder='/Users/josec/Desktop/exoncap/Mygal/alns/DNAbyLoci',dnacdsbyloci_folder='/Users/josec/Desktop/exoncap/Mygal/alns/DNAcdsbyLoci')
# to_aln_files = [x for x in Path('/Users/josec/Desktop/exoncap/Mygal/alns/DNAcdsbyLoci').iterdir() if '.fa' in x.suffix]
# aln = Path("/Users/josec/Desktop/exoncap/Mygal/LociAlns")
# for loci in to_aln_files:
#     lessgappy_maffter(fasta_file_path=loci, out_file_path=aln/f"{loci.name}",num_threads=7)

# def ugly():
    # good_loci = codex_file_reader(codex_file='/Users/josec/Desktop/exoncap/Mygal/blasttest.txt').values()
    # with open('/Users/josec/Desktop/exoncap/Mygal/refseqcds_intronerator_out.txt','r') as outhandle1:
    #     intronerator_out_string = outhandle1.read()
    #     intronerator_out = intronerator_out_string.split("END")
    # intronerator_out_path = '/Users/josec/Desktop/exoncap/Mygal/refseqcds_intronerator_out_filtered.txt'

    # with open(intronerator_out_path,'w') as out_handle:
    #     for item in intronerator_out:
    #         tname = item.split('\t')[0]
    #         gene_raw = tname.split("cds_")[1].split("_")
    #         gene = f"{gene_raw[0]}_{gene_raw[1]}"
    #         if gene in good_loci:
    #             out_handle.write(item+"END")
    # return
# ugly()


# This step takes a whiiiiiiiiiiiiile!! 3 days
# intronerator_out = intronerator(threads=8, genome_fasta='/Users/josec/Desktop/exoncap/Mygal/RefSeq/GCF_000365465.2_Ptep_2.0_genomic.fna', target_cds='/Users/josec/Desktop/exoncap/Mygal/refseqcds_loci.fas')
# intronerator_out_path = '/Users/josec/Desktop/exoncap/Mygal/refseqcds_intronerator_out_filtered.txt'
# with open(intronerator_out_path,'w') as out_handle:
#     for item in intronerator_out:
#         out_handle.write(item)
# VulgarityFilter.vulgarity_filter(intronerator_out_path)

# cut_genealns_to_exonalns(codex_file='/Users/josec/Desktop/exoncap/Mygal/blasttest.txt', cds_exon_folder='/Users/josec/Desktop/exoncap/Mygal/Exons', cds_loci_folder='/Users/josec/Desktop/exoncap/Mygal/LociAlns', exonaln_foler='/Users/josec/Desktop/exoncap/Mygal/ExonAlns')

# exon_aln_filter(exonaln_folder='/Users/josec/Desktop/exoncap/Mygal/ExonAlns', filtered_exonaln_folder='/Users/josec/Desktop/exoncap/Mygal/ExonAlns_Filtered', min_exon_length=80, min_taxoncoverage=4, min_PID=0)


