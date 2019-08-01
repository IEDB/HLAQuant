from Bio import SeqIO
from mhc_G_domain import mhc_G_domain
import tempfile
import json
import pandas as pd
import subprocess
import os


class HLAQuant:
    def __init__(self, sample_file=None, threads=1):
        # Initialization of some basic stuff
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
        self.IMGT_URL = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta"
        self.imgt_loc = os.path.join(
            self.package_directory, 'data/hla_nuc.fasta')
        self.threads = threads
        if sample_file:
            self.samples = sample_file

    def run_cmd(self, cmd, input_string=''):
        """Run the cmd with input_string as stdin and return output.

        cmd -- The shell command to run
        input_string -- String to pass in via stdin
        """
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True,
                             env=dict(os.environ, my_env_prop='value'), shell=True)

        out, stderr = p.communicate(input=input_string)
        if p.returncode:
            raise Exception('Cmd {} failed: {}'.format(cmd, stderr))

        return out

    def load_IMGT_db(self):
        with open(self.imgt_loc, 'r') as infile:
            records = list(SeqIO.parse(infile))
        return records

    def get_IMGT_db(self):
        """Pulls the IMGT database from their FTP server in FASTA format"""
        imgt_cmd = (" ").join(["wget", self.IMGT_URL, ">",
                               self.imgt_loc])
        self.run_cmd(imgt_cmd)

    def build_index(self, sample):
        # Run Salmon index
        salmon_index = (" ").join(["salmon", "index", "-p", self.threads, "-t",
                                   sample + ".fasta", "-i", sample + "_index"])
        self.run_cmd(salmon_index)

    def quantify(self, sample):
        # Run Salmon quant
        salmon_quant = (" ").join(["salmon", "index", "-p", self.threads, "-t",
                                   sample + ".fasta", "-i", sample + "_index"])
        self.run_cmd(salmon_quant)

    def get_g_dom(self, seq):
        # Get the groove domain of a sequence
        gd = mhc_G_domain(seq)
        g_dom = gd.get_g_domain()
        if g_dom:
            return g_dom
        else:
            return None

    def prep_index(self, hapl_f):
        """Takes haplotype FASTA file and retrieves the
        associated g domains to be used for salmon indexing
        
        hapl_f -- FASTA format file containing alleles from IMGT
        """
        # Get G-domains
        # Add G-domain to template (unique name)
        print("placeholder")

    def build_haplotype(self, alleles, imgt_seqs):
        allele_seqs = []
        for allele in alleles:

    def parse_json(self, json_f):
        """Parses json input containing donors and their
        HLA alleles

        json_f -- name of input file
        """
        with open(json_f, 'r') as jsf:
            try:
                donor_dict = json.load(jsf)
            except ValueError:
                print("Error loading JSON file, see README for proper format")

        # This should be {'sample_id': ["HLA-A(1)",...,"HLA-DRB1(1)"]}
        return donor_dict

    def parse_tsv(self, tsv_f):
        """Parses tab separated input containing donors and their
        HLA alleles

        tsv_f -- name of input file
        """
        donor_dict = {}
        with open(tsv_f, 'r') as tsvf:
            for line in tsvf:
                line_list = line.rstrip().split("\t")
                donor = line_list[0]
                alleles = line_list[1:]
                donor_dict[donor] = alleles

        # This should be {'sample_id': ["HLA-A(1)",...,"HLA-DRB1(1)"]}
        return donor_dict

    def run_pipeline(self, json=False, typed=True, infile):
        """Main driver for the pipeline
        
        json -- whether input is in JSON format, default is False
        typed -- whether or not input file has alleles or FASTQ files
        infile -- input file containing sample ids and their alleles or FASTQ locations
        """
        if typed:
            if json:
                donor_dict = parse_json(infile)
            else:
                donor_dict = parse_tsv(infile)
        else:
            print("HLA Typing is unsupported in the prototype")
        
        imgt_seqs = self.get_IMGT_db()

        for donor in donor_dict.keys():
            hapl_f = build_haplotype(donor_dict[donor], imgt_seqs)
            index_f = prep_index(hapl_f)
            build_index(donor)
            quantify(donor)
    """
    def hla_type(self):
        # Run HLA-HD
        print("placeholder")
    """