#!/usr/bin/env python3

import re
import subprocess
import sys

from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from multiprocessing import Pool
from pathlib import Path
from textwrap import dedent
# from timeit import default_timer as timer

from AgalmaAA2dna import agalmaaa2txtmdna, codex_file_reader, blaster
from VulgarityFilter import vulgarity_filter

usage = """
python3 agalmacap.py --param parameterfile.txt
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


def alncutter(partitions, aln_file, aln_type, genealn_outdir,write_ref):
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


def codex_generator(num_threads, query, blastdb_prefix, codex_path):
    """ Blasts the AA consensus of each aglama loci to the Refseqprotein
    This function removes any agalma genes that blast to the same
    RefSeq gene. This ensures that we do not design baits that target multiple
    portions of the genome."""
    blast_output = blaster(num_threads, query, blastdb_prefix).split("\n")
    raw_list = []
    for line in blast_output:
        try:
            raw_list.append(line.split()[1])
        except IndexError:
            pass
    with open(codex_path, 'w') as codex_handle:
        for line in blast_output:
            try:
                if raw_list.count(line.split()[1]) == 1:
                    codex_handle.write(line+"\n")
            except IndexError:
                pass
    return


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


def gene_aligner(dnacdsbyloci_folder, locialn_folder, num_threads):
    """Creates a gene alignment for each fasta file in the
    dnacdsbyloci_folder. The alignment is simmilar to muscle
    output but run with mafft. See lessgappy_mafter for more info.
    """
    to_aln_files = [x for x in Path(dnacdsbyloci_folder).iterdir() if '.fa' in x.suffix]
    for loci in to_aln_files:
        lessgappy_maffter(
                            fasta_file_path=loci,
                            out_file_path=Path(locialn_folder)/f"{loci.name}",
                            num_threads=num_threads)
    return


def clean_exonerate(genome_fasta, target_cds, intronerator_out_path, cds_exon_folder, num_threads):
    """
    """
    intronerator_out = intronerator(
                                    num_threads=num_threads,
                                    genome_fasta=genome_fasta,
                                    target_cds=target_cds)
    with open(intronerator_out_path, 'w') as out_handle:
        for item in intronerator_out:
            out_handle.write(item)
    vulgarity_filter(intronerator_out_path, cds_exon_folder)
    return


def cut_genealns_to_exon_alns(codex_file, cds_exon_folder, loci_alns, exonaln_foler):
    """ For each loci in the codex, this will cut the loci alignment by exons
    and save it to exonaln_folder as Loci_ExonNumber.fas
    """
    good_loci = codex_file_reader(codex_file=codex_file).values()
    for loci in good_loci:
        cds_exon_file = Path(cds_exon_folder)/f"{loci}.fas"
        loci_aln_file = Path(loci_alns)/f"{loci}.fas"
        if cds_exon_file.exists() and loci_aln_file.exists():
            loci_partitions = generate_partitions(
                                                substring_fasta=str(cds_exon_file),
                                                aln_file=str(loci_aln_file))
            alncutter(
                    loci_partitions, aln_file=str(loci_aln_file),
                    aln_type='fasta', genealn_outdir=exonaln_foler)
    return


def aln_filter(aln_folder, filtered_aln_folder, write_ref, min_exon_length=0, min_taxoncov=0):
    """Reads all of the fasta sequences in the aln_folder and saves them to the filtered_aln_folder
    if they are longer than the min_exon_length AND have more sequences in the alignment than
    min_taxoncov. The filtered_aln_folder is created if it does not exist.
    If write_ref='True' the reference CDS sequence will be retained in final alingment.
    """
    Path(filtered_aln_folder).mkdir(exist_ok=True)
    efls = [str(exon_file) for exon_file in Path(aln_folder).iterdir() if '.fa' in exon_file.suffix]
    for exon_file in efls:
        records = list(SeqIO.parse(exon_file, 'fasta'))
        if len(records) > 0:
            if len(records) >= min_taxoncov \
                        and len(records[0].seq) >= min_exon_length:
                outfile_stem = Path(exon_file).name
                out_path = Path(filtered_aln_folder)/outfile_stem
                if write_ref == 'True':
                    with open(out_path, 'w') as out_handle:
                        SeqIO.write(records, out_handle, 'fasta')
                elif write_ref == 'False':
                    recs_to_write = [rec for rec in records if '_cds_' not in rec.name]
                    with open(out_path, 'w') as out_handle:
                        SeqIO.write(recs_to_write, out_handle, 'fasta')
                else:
                    print('Error reading write_ref variable')
                    sys.exit()

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
            'min_taxoncov',
            'min_exon_length',
            'write_ref']
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
    try:
        scratch.mkdir(exist_ok=False)
    except FileExistsError:
        print(f"delete folder: {scratch} ")
        # NOT SURE IF THIS MATTERS. NEED TO TEST
        # scratch2 = basepath/'OLD_Agalmacap_scratch_files'
        # scratch.replace(scratch2)
        # scratch.mkdir(exist_ok=False)
    aa_loci = scratch/'1_aa_aln'            # Agalma supermatrix cut into gene alignments. Each agalma gene is now called a loci
    fil_aln_loci = scratch/'2_aa_faln'      # AA loci alingments with min_taxoncov filtered applied
    cons_aa_loci = scratch/'3_aa_cons.fas'  # Consensus seqs of each AA loci alignment
    codex_path = scratch/'4_codex.txt'              # Consensus of agalma loci blasted to RefSeq AA. Used to annotate each loci as a RefseqProtein ID
    dna_loci = scratch/'5_dna_loci'         # DNA sequences of each loci
    refseqcds_loci = scratch/'6_refseqcds_loci.fas'             # All of the RefSeq CDS that were in the codex (ie all the CDS that agalma identified as homologous loci)
    dnacds_loci = scratch/'7_dnacds_loci'                # DNA sequences of each loci with the RefSeqCDS added
    intronerator_out_path = scratch/'8_intronerator_out.txt'      # Exonerate output file listing locations of introns to be parsed by vulgarityfilter
    cds_exon_folder = scratch/'9_cds_exon_folder'            # Each Refseq CDS loci cut into exons

    # make directory for output of pipeline
    alignments = basepath/'Agalmacap_output_alignments'
    alignments.mkdir(exist_ok=True)
    loci_alns = alignments/'loci_alns'                # Alignments of each loci with the RefSeqCDS as the first sequence. 
    exon_alns = alignments/'exon_alns'                # Loci alignments cut by the intron boundries of the RefSeqCDS genome.
    exon_alns_filtered = alignments/'exon_alns_filtered'       # Exon alignments after the min length and min taxon coverage filters are applied

    keys = [
            'aa_loci',
            'fil_aln_loci',
            'cons_aa_loci',
            'codex_path',
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
        if not Path(path).suffix:
            path.mkdir(exist_ok=True)
        values.append(str(path))
    outdir = dict(zip(keys, values))

    return outdir


def pipeline(param):
    param = param_reader(param)
    outdir = outdir_creator(param)

    partitions = raxpartition_reader(part_file=param['aglama_partition_file'])

    alncutter(
                partitions=partitions,
                aln_file=param['agalma_supermatrix_file'],
                aln_type=param['supermatrix_aln_type'],
                genealn_outdir=outdir['aa_loci'])

    aln_filter(
                aln_folder=outdir['aa_loci'],
                filtered_aln_folder=outdir['fil_aln_loci'],
                min_exon_length=0,
                min_taxoncov=param['min_taxoncov'],
                write_ref=True)

    consensus_generator(
                        aln_folder=outdir['fil_aln_loci'],
                        aln_type='fasta',
                        consensus_file=outdir['cons_aa_loci'])

    codex_generator(
                    num_threads=param['threads'],
                    query=outdir['cons_aa_loci'],
                    blastdb_prefix=param['refseqprotein_path'],
                    codex_path=outdir['codex_path'])

    agalmaaa2txtmdna(
                    codex_file=outdir['codex_path'],
                    alnaa_folder=outdir['fil_aln_loci'],
                    txtm_folder=param['transcriptome_folder'],
                    loci_dna_out_folder=outdir['dna_loci'],
                    num_threads=param['threads'])

    fetch_refseqcds(
                    refseqcds=param['refseqcds_path'],
                    codex_file=outdir['codex_path'],
                    outfile=outdir['refseqcds_loci'])

    cds_loci_merger(
                    codex_file=outdir['codex_path'],
                    cds_file=outdir['refseqcds_loci'],
                    dnabyloci_folder=outdir['dna_loci'],
                    dnacdsbyloci_folder=outdir['dnacds_loci'])

    gene_aligner(
                dnacdsbyloci_folder=outdir['dnacds_loci'],
                locialn_folder=outdir['loci_alns'],
                num_threads=param['threads'])

    clean_exonerate(
                    genome_fasta=param['refseqgenome_path'],
                    target_cds=outdir['refseqcds_loci'],
                    intronerator_out_path=outdir['intronerator_out_path'],
                    cds_exon_folder=outdir['cds_exon_folder'],
                    num_threads=param['threads'])

    cut_genealns_to_exon_alns(
                            codex_file=outdir['codex_path'],
                            cds_exon_folder=outdir['cds_exon_folder'],
                            loci_alns=outdir['loci_alns'],
                            exonaln_foler=outdir['exon_alns'])

    aln_filter(
                aln_folder=outdir['exon_alns'],
                filtered_aln_folder=outdir['exon_alns_filtered'],
                min_exon_length=param['min_exon_length'],
                min_taxoncov=param['min_taxoncov'],
                write_ref=param['write_ref'])
    return


# pipeline(param='/Users/josec/Desktop/exoncap/Mygal/AgalmacapTesting/Param-test1.txt')


def main():
    # Engine for the program.
    args = sys.argv[1:]
    # Print usage if no input is put in by user
    if not args:
        print(usage)
        sys.exit(1)
    if args[0] == '--param':
        param = args[1]

    # start = timer()
    pipeline(param=param)
    # end = timer()
    # print(end - start)

    return


if __name__ == '__main__':
    main()
