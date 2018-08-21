#!/usr/bin/env python3

import re
import subprocess
import sys

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


def dna_or_aa(seq_file):
    """Reads 5 lines of fasta seq_file and guesses if DNA or RNA.
    Exits if invalid residues are present"""
    records = SeqIO.parse(seq_file, "fasta")
    first5 = list(next(records) for _ in range(4))
    seqset = set(''.join([str(record.seq) for record in first5]))
    dna = set("NATGC")
    nondna = seqset-dna
    if len(nondna) == 0:
        alphabet = 'dna'
    elif len(nondna) > 4 and '-' not in nondna:
        alphabet = 'aa'
    elif 'X' in nondna or '-' in nondna:
        print(f"ERROR:{seq_file} has '-' or 'X' in it. Please fix and run again")
        sys.exit(1)
    else:
        print(f"Error reading {seq_file} ")
        sys.exit(1)
    return alphabet


def blaster(threads, query, blastdb_prefix, blast_out_path=None):
    """Runs a blast search and returns the names of the top hit sequences from the blast database
    (min eval = 1e-50)

    Args:
        threads: the number of threads
        query: path to query file for blast search
        blastdb_prefix: path to blast database file. Does not include blast extension.
        Creates a blast database if does not exist
    Returns:
        blast results as a string. Saves string to text file in blast_out_path if specified
    """

    qalphabet = dna_or_aa(query)
    dbalphabet = dna_or_aa(blastdb_prefix)
    if dbalphabet == 'aa':
        blast_suffix, dbtype = '.pin', 'prot'
        if qalphabet == 'dna':
            blast_prog = 'tblastx'
        elif qalphabet == 'aa':
            blast_prog = 'blastp'
    if dbalphabet == 'dna':
        blast_suffix, dbtype = '.nin', 'nucl'
        if qalphabet == 'dna':
            blast_prog = 'blastn'
        elif qalphabet == 'aa':
            blast_prog = 'tblastn'
    # Create blast database if there is none
    if not Path(f'{blastdb_prefix}{blast_suffix}').exists():
        dbcommand = f"makeblastdb -in {blastdb_prefix} -dbtype {dbtype}"
        run_command(dbcommand)
    # Run Blast
    bcommand = f"{blast_prog} -query {query} -db {blastdb_prefix} -evalue 1e-50 -outfmt '6 qseqid sseqid' -max_target_seqs 1 -max_hsps 1 -num_threads {threads}"
    blast_output = run_command(bcommand, command_out_path=blast_out_path)
    return blast_output


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


def fetchseq(names_tofetch, seqfile):
    """Searches for a list of names in fasta seqfile and
     returns them as a biopython SeqIO records"""
    seq_dict = SeqIO.to_dict(SeqIO.parse(seqfile, "fasta"))
    fetched_records = [seq_dict[name] for name in names_tofetch]
    return fetched_records


def groupagalmabytxtm(codex_dict, alnaa_folder, txtm):
    """Saves the sequence that has the same header as txtm from each
    alignment in alnaa_folder to a file named "AAbytxtm/<txtm>.fas"
    Strips the gaps from the AA sequences and removes seqs shorter than 50AA
    See agalmaaa2txtmdna for more info"""
    txtm_aaseq_dict = {}
    for agalID, RefseqProtID in codex_dict.items():
        agal_aln_path = str(Path(alnaa_folder)/f"{agalID}.fas")
        agal_aln = SeqIO.index(agal_aln_path, 'fasta')
        try:
            aaseq = str(agal_aln[txtm].seq).replace('-', '')
            if len(aaseq) > 50:  # Causes errors with blast if too short
                txtm_aaseq_dict[RefseqProtID] = aaseq
        except KeyError:
            txtm_aaseq_dict[RefseqProtID] = None
    txtm_aa_out_folder = Path(alnaa_folder)/"AAbytxtm"
    txtm_aa_out_folder.mkdir(exist_ok=True)
    with open(txtm_aa_out_folder/f"{txtm}.fas", 'w') as txtm_aa_out_handle:
        for RefseqProtID, aaseq in txtm_aaseq_dict.items():
            if aaseq:
                txtm_aa_out_handle.write(f">{RefseqProtID}\n{aaseq}\n")
    return


def blastbytxtm(alnaa_folder, txtm_folder, txtm):
    """Each amino acid sequence is blasted to the DNA transcriptome assembly.
    The results are written to the  txtm_dna_out_folder in fasta format.
    The fasta headers are '>txtm-loci' to be legible for cat_by_gene()
    See agalmaaa2txtmdna for more info"""
    query_path = Path(alnaa_folder)/"AAbytxtm"
    blastdb_path = str(Path(txtm_folder)/f"{txtm}.fas")
    blast_out = str(query_path/f"{txtm}_recipblast.txt")
    # This step takes 10 min.
    blaster(threads=7, query=str(query_path/f"{txtm}.fas"), blastdb_prefix=blastdb_path, blast_out_path=blast_out)
    # store blast results as dictionary
    rblast_dict = {}
    with open(blast_out, 'r') as rblastout_handle:
        rblast_lines = [line.split() for line in rblastout_handle.readlines()]
        for line in rblast_lines:
            rblast_dict[line[1]] = line[0]
    # fetch each sequences from the transcriptome and write them to a new file send to cat_by_gene()
    blast_hit_records = fetchseq(names_tofetch=rblast_dict.keys(), seqfile=blastdb_path)
    txtm_dna_out_folder = Path(alnaa_folder)/"DNAbytxtm"
    txtm_dna_out_folder.mkdir(exist_ok=True)
    dna_out_path = str(txtm_dna_out_folder/f"{txtm}_DNA.fas")
    with open(dna_out_path, 'w') as dna_out_handle:
        for blast_hit_record in blast_hit_records:
            string = f">{txtm}-{rblast_dict[blast_hit_record.name]}\n{blast_hit_record.seq}\n"
            dna_out_handle.write(string)
    return txtm_dna_out_folder


def cat_by_gene(txtm_list, loci_list, in_path, file_ending, out_path):
    """This concatenates sequence files for all txtmanisms by loci.
    Sequence files must have fasta heading of '>txtm-loci'
    """
    txtm_dict = {}
    for txtm in txtm_list:
        seq_dict = SeqIO.to_dict(SeqIO.parse(str(in_path/f"{txtm}{file_ending}"), "fasta"))
        txtm_dict[txtm] = seq_dict
    for loci in loci_list:
        out_file_path = f"{out_path}/{loci}.fas"
        with open(out_file_path, 'a') as out_handle:
            for txtm in txtm_list:
                record_id = f"{txtm}-{loci}"
                try:
                    out_handle.write(f">{txtm}\n{txtm_dict[txtm][record_id].seq}\n")
                except KeyError:
                    pass
    return


def agalmaaa2txtmdna(codex_file, alnaa_folder, txtm_folder):
    """Starting from Agalma AA loci alignments, this will blast each AA sequence to
    the transcriptome it came from and will save the results by loci as DNA fasta files.

    Args:
        codex_file: textfile of each loci as a line: agalma_gene_id \t RefSeqCDS_topBlastHit
        alnaa_folder: path to output of  alncutter() on agalma supermatrix.
                        The loci AA alignments with agalma fasta headers
        txtm_folder: folder with transcriptomes from agalma.
                        The names of the transcriptome files must match the fasta headers in the
                        AA alignments (whatever agalma sp was designated in the catalog)
    Out:
        written to alnaa_folder/DNAbyLoci/RefSeqProtIDxx.fas

    """
    txtms = []
    for file in Path(txtm_folder).iterdir():
        if not file.stem.startswith('.'):
            if 'fa' in file.suffix:
                txtms.append(file.stem)
    # Create a dictionary to translate agalma ID's to the top hit refseq Protein ID
    codex_dict = codex_file_reader(codex_file)
    # Make new fasta of agalma output grouped by transcriptome, strips gaps
    for txtm in txtms:
        groupagalmabytxtm(codex_dict, alnaa_folder, txtm)
        txtm_dna_out_folder = blastbytxtm(alnaa_folder, txtm_folder, txtm)
    # writes sequences with fasta header of txtm and the fasta file name is the loci name
    loci_dna_out_folder = Path(alnaa_folder)/"DNAbyLoci"
    loci_dna_out_folder.mkdir(exist_ok=True)
    cat_by_gene(txtm_list=txtms, loci_list=codex_dict.values(), in_path=txtm_dna_out_folder, file_ending="_DNA.fas", out_path=loci_dna_out_folder)
    return


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


def codex_file_reader(codex_file):
    """The codex dictionary that this function returns is the blast result of the
    agalma alignment reference to the reference genome CDS.
    The dicitonary key is the agalma gene alignment ID and the
    dictionary value is the RefSeq Protein ID from the reference genome.

    This function removes any agalma genes that blast to the same
    RefSeq gene. This ensures that we do not design baits that target multiple
    portions of the genome.
    """
    raw_list = []
    bad_list = []
    codex_dict = {}
    with open(codex_file) as f:
        for line in f:
            raw_list.append(line.split()[1])

    for item in raw_list:
        if raw_list.count(item) > 1:
            if item not in bad_list:
                bad_list.append(item)

    with open(codex_file) as f:
        for line in f:
            (agalID, RefseqProtID) = line.split()
            if RefseqProtID not in bad_list:
                codex_dict[agalID] = RefseqProtID
    return codex_dict


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


