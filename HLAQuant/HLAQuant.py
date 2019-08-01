# Need a Bio import at some point
from mhc_G_domain import mhc_G_domain
import json
import pandas as pd
import subprocess
import os


class HLAQuant:
    def __init__(self, sample_file=None, threads=1):
        # Initialization of some basic stuff
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
        self.IMGT_URL = "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta"
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

    def get_IMGT_db(self):
        """Pulls the IMGT database from their FTP server in FASTA format"""
        imgt_cmd = (" ").join(["wget", self.IMGT_URL], ">",
                              self.package_directory, "/data/")
        self.run_cmd(imgt_cmd)

    def get_g_dom(self, seq):
        # Get the groove domain of a sequence
        gd = mhc_G_domain(seq)
        g_dom = get_g_domain()
        if g_dom:
            return g_dom
        else:
            return None

    def hla_type(self):
        # Run HLA-HD
        print("placeholder")

    def prep_index(self):
        # Get G-domains
        # Add G-domain to template (unique name)
        print("placeholder")

    def build_index(self, sample):
        # Run Salmon index
        salmon_index = (" ").join(["salmon", "index", "-p", self.threads, "-t",
                                   sample + ".fasta", "-i", sample + "_index"])
        run_cmd(salmon_index)

    def quantify(self):
        # Run Salmon quant
        salmon_quant = (" ").join(["salmon", "index", "-p", self.threads, "-t",
                                   sample + ".fasta", "-i", sample + "_index"])        

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
