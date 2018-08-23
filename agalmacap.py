#!/usr/bin/env python3

import re
import subprocess

from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from multiprocessing import Pool
from pathlib import Path
from textwrap import dedent

from AgalmaAA2dna import agalmaaa2txtmdna, codex_file_reader, blaster
from VulgarityFilter import vulgarity_filter

usage = """
"""


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
    with open(consensus_file, 'w') as out_file:
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
    """From the VulgarityFilter output run on the RefSeq genome we now have
    a RefSeqCDS cut into exons. From the lessgappy_maffter output we have the
    RefSeqCDS aligned to the transcriptome data. This fuction generates partitions
    based on the exons in the substring_fasta so that we can later cut the
    gene alignment into exon alignments.
    Args:
        substring_fasta: path to VulgarityFilter output. The exon seqs for each RefSeqCDS in order
        aln_file: path to lessgappy_mafter output. The RefSeqCDS sequence must be in the alignment
                    and it must have '_cds_' in the fasta header. Should have it from
                    Genbank.

    Out:
        partitions: a python list of the start and stop positions of each exon in
                    the aln_file.
"""

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


def intronerator(num_threads, genome_fasta, target_cds):
    """Runs exonerate to find exon boundry sites on target genes given a reference genome.
      The output is a "roll your own format" string from exonerate that will be parsed by XXXX"""
    exonerate_command_prefix = dedent(f'''
                                        exonerate
                                        --model est2genome
                                        -q {target_cds}
                                        -t {genome_fasta}
                                        -Q DNA
                                        -T DNA
                                        --showvulgar F
                                        --showalignment F
                                        --softmasktarget T
                                        --verbose 0
                                        --ryo \"%qi\\t%pi\\t%qas\\t%V\\tEND\\n\"
                                        --fsmmemory 20G
                                        --bestn 1
                                        --querychunktotal {num_threads}
                                        --querychunkid
                                        ''').replace('\n', ' ')
    exonerate_commands = [exonerate_command_prefix+str(i+1) for i in range((int(num_threads)))]
    p = Pool(num_threads)
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
    """Find the CDS sequence in refseqcds for each gene in the codex_file.
    sequences are saved to outfile.
    """
    to_fetch = codex_file_reader(codex_file).values()
    fasta_subseter(subset_list=to_fetch, in_fasta_path=refseqcds, out_fasta_path=outfile)
    return


def lessgappy_maffter(fasta_file_path, out_file_path, num_threads):
    """ Runs mafft alignment on fasta file. Options result in
    an alignment with very few gaps (simmilar to muscle).
    Sequences are checked for reverse complementing
    and the output string removes the R_ prefix that mafft adds to
    sequences that were reverse complemented
    """
    command = dedent(f'''
                    mafft --reorder
                    --adjustdirection
                    --leavegappyregion
                    --op 3.0
                    --ep 1.0
                    --maxiterate 1000
                    --retree 1
                    --genafpair
                    --thread {num_threads}
                    --quiet
                    --auto {fasta_file_path}
                    ''').replace('\n', ' ')
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
                    with open(loci_cds_file_path, 'w') as out_handle:
                        SeqIO.write(cds, out_handle, 'fasta')
                        SeqIO.write(records_to_write, out_handle, 'fasta')
                        break
    return


def cut_genealns_to_exon_alns(codex_file, cds_exon_folder, cds_loci_folder, exonaln_foler):
    """ For each loci in the codex, this will cut the loci alignment by exons
    and save it to exonaln_folder as Loci_ExonNumber.fas
    """
    good_loci = codex_file_reader(codex_file=codex_file).values()
    for loci in good_loci:
        cds_exon_file = Path(cds_exon_folder)/f"{loci}.fas"
        loci_aln_file = Path(cds_loci_folder)/f"{loci}.fas"
        if cds_exon_file.exists() and loci_aln_file.exists():
            loci_partitions = generate_partitions(
                                                substring_fasta=str(cds_exon_file),
                                                aln_file=str(loci_aln_file))
            alncutter(
                    loci_partitions, aln_file=str(loci_aln_file),
                    aln_type='fasta', genealn_outdir=exonaln_foler)
    return


def aln_filter(aln_folder, filtered_aln_folder, min_exon_length=0, min_taxoncov=0):
    """Reads all of the fasta sequences in the aln_folder and saves them to the filtered_aln_folder
    if they are longer than the min_exon_length AND have more sequences in the alignment than
    min_taxoncov. The filtered_aln_folder is created if it does not exist.
    """
    Path(filtered_aln_folder).mkdir(exist_ok=True)
    efls = [str(exon_file) for exon_file in Path(aln_folder).iterdir() if '.fa' in exon_file.suffix]
    exons2write = []
    for exon_file in efls:
        records = list(SeqIO.parse(exon_file, 'fasta'))
        if len(records) > 0:
            if len(records) >= min_taxoncov \
                        and len(records[0].seq) >= min_exon_length:
                exons2write.append(exon_file)
                outfile_stem = Path(exon_file).name
                out_path = Path(filtered_aln_folder)/outfile_stem
                with open(out_path, 'w') as out_handle:
                    SeqIO.write(records, out_handle, 'fasta')
    return


def param_reader(paramfile_path):
    """Reads a parameter file and returns parameters in a dictionary.
    numbers are stored as python integers. Everything else stored as strings.
    """
    param = {}
    keys = [
            'projectwd',
            'aglama_partition_file',
            'agalma_supermatrix_file',
            'transcriptome_folder',
            'refseqcds_path',
            'refseqgenome_path',
            'refseqprotein_path',
            'supermatrix_aln_type',
            'threads',
            'min_taxoncov']
    with open(paramfile_path, 'r') as param_handle:
        values = [p.split('#')[0].strip() for p in param_handle.readlines()]
    param = dict(zip(keys, values))
    # Turn numbers in parameters into python integer types
    for key in param:
        try:
            param[key] = int(param[key])
        except ValueError:
            pass
    return param


def outdir_creator(param):

    basepath = Path(param['projectwd'])

    # make scratch folder for intermediate files
    scratch = basepath/'Agalmacap_scratch_files'
    scratch.mkdir(exist_ok=True)
    aa_loci = scratch/'1_aa_aln'            # Agalma supermatrix cut into gene alignments. Each agalma gene is now called a loci
    fil_aln_loci = scratch/'2_aa_faln'      # AA loci alingments with min_taxoncov filtered applied
    cons_aa_loci = scratch/'3_aa_cons.fas'  # Consensus seqs of each AA loci alignment
    codex_path = '4_codex.txt'              # Consensus of agalma loci blasted to RefSeq AA. Used to annotate each loci as a RefseqProtein ID
    aa_loci_txtm = scratch/''               # AA sequences from gene alignments grouped by transcriptome
    dna_loci = scratch/'5_dna_loci'         # DNA sequences of each loci
    refseqcds_loci = scratch/'6_refseqcds_loci.fas'             # All of the RefSeq CDS that were in the codex (ie all the CDS that agalma identified as homologous loci)
    dnacds_loci = scratch/''                # DNA sequences of each loci with the RefSeqCDS added
    intronerator_out_path = scratch/''      # Exonerate output file listing locations of introns to be parsed by vulgarityfilter
    cds_exon_folder = scratch/''            # Each Refseq CDS loci cut into exons

    # make directory for output of pipeline
    alignments = basepath/'Agalmacap_output_alignments'
    alignments.mkdir(exist_ok=True)
    loci_alns = alignments/''                # Alignments of each loci with the RefSeqCDS as the first sequence. 
    exon_alns = alignments/''                # Loci alignments cut by the intron boundries of the RefSeqCDS genome.
    exon_alns_filtered = alignments/''       # Exon alignments after the min length and min taxon coverage filters are applied

    keys = [
            'aa_loci',
            'fil_aln_loci',
            'cons_aa_loci',
            'codex_path',
            'aa_loci_txtm',
            'dna_loci',
            'refseqcds_loci',
            'dnacds_loci',
            'intronerator_out_path',
            'cds_exon_folder',
            'loci_alns',
            'exon_alns',
            'exon_alns_filtered']
    path_values = [
                    aa_loci,
                    fil_aln_loci,
                    cons_aa_loci,
                    codex_path,
                    aa_loci_txtm,
                    dna_loci,
                    refseqcds_loci,
                    dnacds_loci,
                    intronerator_out_path,
                    cds_exon_folder,
                    loci_alns,
                    exon_alns,
                    exon_alns_filtered]
    values = []
    for path in path_values:
        if not path.suffix:
            path.mkdir(exist_ok=True)
        values.append(str(path))
    outdir = dict(zip(keys, values))

    return outdir


def pipeline():
    param = param_reader('/Users/josec/Desktop/exoncap/Mygal/AgalmacapTesting/Param-test1.txt')
    outdir = outdir_creator(param)

    partitions = raxpartition_reader(part_file=param['aglama_partition_file'])

    alncutter(
            partitions=partitions, aln_file=param['agalma_supermatrix_file'],
            aln_type=param['supermatrix_aln_type'], genealn_outdir=outdir['aa_loci'])

    aln_filter(
            aln_folder=outdir['aa_loci'], filtered_aln_folder=outdir['fil_aln_loci'],
            min_exon_length=0, min_taxoncov=param['min_taxoncov'])

    consensus_generator(
                        aln_folder=outdir['fil_aln_loci'], aln_type='fasta',
                        consensus_file=outdir['cons_aa_loci'])
    # Genrate Codex
    blaster(
            num_threads=param['threads'], query=outdir['cons_aa_loci'],
            blastdb_prefix=param['refseqprotein_path'], blast_out_path=outdir['codex_path'])

    # DONT run agalmaaa2txtmdna multiple times because it appends to sequence files in DNAbyLoci
    # ie it does not overwrite.
    agalmaaa2txtmdna(
                    codex_file=outdir['codex_path'],
                    alnaa_folder='/Users/josec/Desktop/exoncap/Mygal/alns/',
                    txtm_folder=param['transcriptome_folder'],
                    loci_dna_out_folder=outdir['dna_loci'])
    # Generate Gene alignments with refseq and minimal gaps
    fetch_refseqcds(
                    refseqcds=param['refseqcds_path'],
                    codex_file=outdir['codex_path'],
                    outfile=outdir['refseqcds_loci'])
    # cds_loci_merger(codex_file=codex_path, cds_file='/Users/josec/Desktop/exoncap/Mygal/refseqcds_loci.fas', dnabyloci_folder='/Users/josec/Desktop/exoncap/Mygal/alns/DNAbyLoci',dnacdsbyloci_folder='/Users/josec/Desktop/exoncap/Mygal/alns/DNAcdsbyLoci')
    # to_aln_files = [x for x in Path('/Users/josec/Desktop/exoncap/Mygal/alns/DNAcdsbyLoci').iterdir() if '.fa' in x.suffix]
    # aln = Path("/Users/josec/Desktop/exoncap/Mygal/loci_alns")
    # for loci in to_aln_files:
    #     lessgappy_maffter(fasta_file_path=loci, out_file_path=aln/f"{loci.name}",num_threads=threads)
    return

   
pipeline()

# # This step takes a whiiiiiiiiiiiiile!! 3 days
# intronerator_out = intronerator(num_threads=threads, genome_fasta=refseqgenome_path, target_cds='/Users/josec/Desktop/exoncap/Mygal/refseqcds_loci.fas')
# intronerator_out_path = '/Users/josec/Desktop/exoncap/Mygal/refseqcds_intronerator_out_filtered.txt'
# with open(intronerator_out_path, 'w') as out_handle:
#     for item in intronerator_out:
#         out_handle.write(item)
# VulgarityFilter.vulgarity_filter(intronerator_out_path)

# cut_genealns_to_exon_alns(codex_file=codex_path, cds_exon_folder='/Users/josec/Desktop/exoncap/Mygal/Exons', cds_loci_folder='/Users/josec/Desktop/exoncap/Mygal/loci_alns', exonaln_foler='/Users/josec/Desktop/exoncap/Mygal/exon_alns')

# aln_filter(aln_folder='/Users/josec/Desktop/exoncap/Mygal/exon_alns', filtered_aln_folder='/Users/josec/Desktop/exoncap/Mygal/exon_alns_filtered', min_exon_length=80, min_taxoncov=4)
