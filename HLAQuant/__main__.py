import argparse
import sys
from HLAQuant import HLAQuant

logo = r"""
 _    _ _               ____                    _   
| |  | | |        /\   / __ \                  | |  
| |__| | |       /  \ | |  | |_   _  __ _ _ __ | |_ 
|  __  | |      / /\ \| |  | | | | |/ _` | '_ \| __|
| |  | | |____ / ____ \ |__| | |_| | (_| | | | | |_ 
|_|  |_|______/_/    \_\___\_\\__,_|\__,_|_| |_|\__|
"""
print(logo)

prsr = argparse.ArgumentParser(
    prog='HLAQuant', description='Quantify HLA expression')
prsr.add_argument('-t', help="Number of threads to be used by Salmon", type=str, 
                  metavar='# threads', required=True)
prsr.add_argument('-hla', help="Tab separated file containing sample_id and then HLA alleles",
                  type=str, metavar='hla_file', required=True)
prsr.add_argument('-fastq', help="Tab separated file containing sample_id and then FASTQ files",
                  type=str, metavar='fastq_file', required=True)
prsr.add_argument('-o', help="Specify output directory to place results in", type=str, 
                  metavar='out_dir', required=True)
args = prsr.parse_args()

hq = HLAQuant.HLAQuant(args.t, args.o)
hq.run_pipeline(args.hla, args.fastq)
