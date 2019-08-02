from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from HLAQuant.mhc_G_domain import mhc_G_domain
import tempfile
import json
import pandas as pd
import subprocess
import os
import sys


class HLAQuant:
    """HLAQuant provides a pipeline for allele specific expression of HLA genes.
    User must provide the number of threads to run Salmon with, as well as two files,
    one containing a sample id mapped to alleles, the other containing a sample id mapped
    to the relevant FASTQ files.
    """

    def __init__(self, threads, out_dir):
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
        self.IMGT_URL = ("ftp://ftp.ebi.ac.uk/pub/databases/ipd/",
                         "imgt/hla/hla_nuc.fasta")
        self.salmon_dir = os.path.join(
            self.package_directory, 'data/salmon-0.11.2-linux_x86_64/bin/salmon')
        self.imgt_loc = os.path.join(
            self.package_directory, 'data/hla_nuc.fasta')
        self.threads = threads
        if out_dir[-1] != "/":
            self.out_dir = out_dir + "/"
        else:
            self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir, mode=0o777)

    def run_cmd(self, cmd, input_string=''):
        """Run the cmd with input_string as stdin and return output.

        cmd -- The shell command to run
        input_string -- String to pass in via stdin
        """
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             stdin=subprocess.PIPE, stderr=subprocess.PIPE,
                             universal_newlines=True, close_fds=True,
                             env=dict(os.environ, my_env_prop='value'),
                             shell=True)

        out, stderr = p.communicate(input=input_string)
        if p.returncode:
            raise Exception('Cmd {} failed: {}'.format(cmd, stderr))

        return out

    def load_IMGT_db(self):
        """Loads the IMGT db from FASTA file into list of records"""
        with open(self.imgt_loc, 'r') as infile:
            records = list(SeqIO.parse(infile, "fasta"))
        return records

    def get_IMGT_db(self):
        """Pulls the IMGT database from their FTP server in FASTA format"""
        imgt_cmd = (" ").join(["wget", self.IMGT_URL, ">",
                               self.imgt_loc])
        self.run_cmd(imgt_cmd)

    def build_index(self, sample):
        """Runs the Salmon indexing command

        sample -- the sample id
        """
        salmon_index = (" ").join([self.salmon_dir, "index", "-p", self.threads, "-t",
                                   self.out_dir + sample + ".fasta", "-i",
                                   self.out_dir + sample + "_index"])
        self.run_cmd(salmon_index)

    def quantify(self, sample, paired, fq_files):
        """Runs the Salmon quant command.

        sample -- the sample id
        paired -- whether or not input is paired-end or single-end
        fq_files -- the fastq files used for quantification
        """
        if paired:
            fq1, fq2 = fq_files
            salmon_quant = (" ").join([self.salmon_dir, "quant", "-p", self.threads,
                                       "-i", self.out_dir + sample + "_index", "-l", "A",
                                       "-1", fq1, "-2", fq2, "--validateMappings",
                                       "--rangeFactorizationBins", "4"
                                       "--minScoreFraction", "1", "-o", self.out_dir + sample])
        else:
            fq = fq_files
            salmon_quant = (" ").join([self.salmon_dir, "quant", "-p", self.threads,
                                       "-i", self.out_dir + sample + "_index", "-l", "A",
                                       "-r", fq, "--validateMappings",
                                       "--rangeFactorizationBins", "4",
                                       "--minScoreFraction", "1", "-o", self.out_dir + sample])
        self.run_cmd(salmon_quant)

    def get_g_dom(self, seq):
        """Returns the groove domain of an input HLA sequence"""
        gd = mhc_G_domain(seq)
        g_dom = gd.get_g_domain()
        if g_dom:
            return g_dom
        else:
            return None

    def prep_index(self, sample_id, allele_seqs):
        """Takes haplotype FASTA file and retrieves the
        associated g domains to be used for salmon indexing

        hapl_f -- FASTA format file containing alleles from IMGT
        """
        g_doms = []
        for allele in allele_seqs:
            g_dom_seq = self.get_g_dom(str(allele.seq))
            if g_dom_seq:
                record = SeqRecord(g_dom_seq, id=allele.id,
                                   description=allele.description)
                g_doms.append(record)
            else:
                g_doms.append(allele)
        SeqIO.write(g_doms, self.out_dir + sample_id + ".fasta", "fasta")
        return self.out_dir + sample_id + ".fasta"

    def build_haplotype(self, alleles, imgt_seqs):
        allele_seqs = []
        for allele in alleles:
            for seq in imgt_seqs:
                if allele in seq.description:
                    seq.id = seq.description.split(" ")[1]
                    allele_seqs.append(seq)
                    break
        return allele_seqs

    def parse_alleles(self, tsv_f):
        """Parses tab separated input containing donors and their
        HLA alleles

        tsv_f -- name of input file
        """
        allele_dict = {}
        donors = []
        with open(tsv_f, 'r') as infile:
            for line in infile:
                line_list = line.rstrip().split("\t")
                donor = line_list[0]
                donors.append(donor)
                alleles = line_list[1:]
                allele_dict[donor] = alleles

        # This should be {'sample_id': ["HLA-A(1)",...,"HLA-DRB1(1)"]}
        return allele_dict, donors

    def parse_fastq(self, fastq_f, paired):
        """Parses tab separated input containing donors and their
        FASTQ files. Interprets if single or paired end data

        fastq_f -- name of input file
        """
        try:
            fastq_dict = {}
            with open(fastq_f, 'r') as infile:
                for line in infile:
                    line_list = line.rstrip().split("\t")
                    donor = line_list[0]
                    if paired:
                        fastq_dict[donor] = (line_list[1], line_list[2])
                    else:
                        fastq_dict[donor] = line_list[1]
            return fastq_dict
        except:
            print(
                "Invalid FASTQ file. Please check the README for information on valid formats.", file=sys.stderr)
            sys.exit()

    def run_pipeline(self, sample_file, fastq_file, paired=False):
        """Main driver for the pipeline

        paired -- whether input data is paired end or not
        infile -- input file containing sample ids and their alleles or FASTQ
                  locations
        """
        allele_dict, donors = self.parse_alleles(sample_file)
        fastq_dict = self.parse_fastq(fastq_file, paired)
        print("Loading the IMGT database")
        print("-------------------------")
        imgt_seqs = self.load_IMGT_db()

        for donor in donors:
            print("Processing sample: " + donor)
            hapl_f = self.build_haplotype(allele_dict[donor], imgt_seqs)
            self.prep_index(donor, hapl_f)
            self.build_index(donor)
            fq_files = fastq_dict[donor]
            self.quantify(donor, paired, fq_files)
            print("Done!")
