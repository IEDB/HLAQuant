from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from HLAQuant.mhc_G_domain import mhc_G_domain
import tempfile
import json
import pandas as pd
import subprocess
import os


class HLAQuant:

    def __init__(self, sample_file=None, fastq_file=None, threads=1):
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
        self.IMGT_URL = ("ftp://ftp.ebi.ac.uk/pub/databases/ipd/",
                         "imgt/hla/hla_nuc.fasta")
        self.imgt_loc = os.path.join(
            self.package_directory, 'data/hla_nuc.fasta')
        self.threads = threads
        if sample_file:
            self.samples = sample_file
        if fastq_file:
            self.fastq_file = fastq_file

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
        with open(self.imgt_loc, 'r') as infile:
            records = list(SeqIO.parse(infile, "fasta"))
        return records

    def get_IMGT_db(self):
        """Pulls the IMGT database from their FTP server in FASTA format"""
        imgt_cmd = (" ").join(["wget", self.IMGT_URL, ">",
                               self.imgt_loc])
        self.run_cmd(imgt_cmd)

    def build_index(self, sample_id):
        # Run Salmon index
        salmon_index = (" ").join(["salmon", "index", "-p", self.threads, "-t",
                                   sample_id + ".fasta", "-i",
                                   sample + "_index"])
        self.run_cmd(salmon_index)

    def quantify(self, sample):
        # Run Salmon quant
        salmon_quant = (" ").join(["salmon", "quant", "-p", self.threads,
                                   ])
        self.run_cmd(salmon_quant)

    def get_g_dom(self, seq):
        # Get the groove domain of a sequence
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
                record = SeqRecord(Seq(g_dom_seq), id=allele.id,
                                   description=allele.description)
                g_doms.append(record)
            else:
                g_doms.append(allele)
        SeqIO.write(g_doms, sample_id + ".fasta", "fasta")
        return sample_id + ".fasta"

    def build_haplotype(self, alleles, imgt_seqs):
        allele_seqs = []
        for allele in alleles:
            for seq in imgt_seqs:
                if allele in seq.description:
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
        fastq_dict = {}
        with open(fastq_f, 'r') as infile:
            for line in infile:
                line_list = line.rstrip().split("\t")
                donor = line_list[0]
                if paired:
                    fastq_dict[donor] = (line_list[1], line_list[2])
                else:
                    fastq_dict[donor] = line_list[1]

    def run_pipeline(self, sample_file, fastq_file, paired=False):
        """Main driver for the pipeline

        paired -- whether input data is paired end or not
        infile -- input file containing sample ids and their alleles or FASTQ
                  locations
        """
        allele_dict, donors = self.parse_alleles(sample_file)
        fastq_dict = self.parse_fastq(fastq_file, paired)
        imgt_seqs = self.load_IMGT_db()

        for donor in donors:
            hapl_f = self.build_haplotype(allele_dict[donor], imgt_seqs)
            index_f = self.prep_index(donor, hapl_f)
            self.build_index(donor, index_f)
            self.quantify(donor)
