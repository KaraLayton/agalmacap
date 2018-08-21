#!/usr/bin/env python3

import subprocess
import sys

from Bio import SeqIO
from pathlib import Path


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


def fetchseq(names_tofetch, seqfile):
    """Searches for a list of names in fasta seqfile and
     returns them as a biopython SeqIO records"""
    seq_dict = SeqIO.to_dict(SeqIO.parse(seqfile, "fasta"))
    fetched_records = [seq_dict[name] for name in names_tofetch]
    return fetched_records


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
