import pandas as pd
import subprocess
import os

class HLAQuant:
    def __init__(self, sample_file=None):
        #Initialize stuff for the HLAQuant class
        self.package_directory = os.path.dirname(os.path.abspath(__file__))
        if sample_file:
            self.samples = sample_file

    def run_cmd(self, cmd, input_string=''):
        """
        Run the cmd with input_string as stdin and return output.

        @param cmd: The shell command to run
        @param input_string: String to pass in via stdin
        """
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                             stderr=subprocess.PIPE, universal_newlines=True, close_fds=True,
                             env=dict(os.environ, my_env_prop='value'), shell=True)

        out, stderr = p.communicate(input=input_string)
        if p.returncode:
            raise Exception('Cmd {} failed: {}'.format(cmd, stderr))

        return out