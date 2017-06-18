import sys
import argparse
import os
from datetime import datetime
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype

delim = "\t"

parser = argparse.ArgumentParser("Generate a final report from  GTC files")
parser.add_argument("manifest", help="BPM manifest file")
parser.add_argument("gtc_file", help="Directory containing GTC files")
parser.add_argument("output_file", help="Location to write report")

args = parser.parse_args()

try:
    manifest = BeadPoolManifest(args.manifest)
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

if os.path.isfile(args.output_file):
    sys.stderr.write("Output file already exists, please delete and re-run\n")
    sys.exit(-1)

with open(args.output_file, "w") as output_handle:
    
    samples = []
    gtc_file = args.gtc_file
    if gtc_file.lower().endswith(".gtc"):
    	samples.append(gtc_file)	
    
    
    output_handle.write(delim.join(["Name", "Chr", "Position", "G.type", "LogR_Ratio", "B_AlleleFreq"]) + "\n")
    for gtc_file in samples:
        sys.stderr.write("Processing " + gtc_file + "\n")
        gtc_file = os.path.join(gtc_file)
        gtc = GenotypeCalls(gtc_file)
        genotypes = gtc.get_genotypes()
        #plus_strand_genotypes = gtc.get_base_calls_plus_strand(manifest.snps, manifest.ref_strands)
        #forward_strand_genotypes = gtc.get_base_calls_forward_strand(manifest.snps, manifest.source_strands)
        logr_ratios = gtc.get_logr_ratios()
        ballele_freqs = gtc.get_ballele_freqs()
        

        assert len(genotypes) == len(manifest.names)
        for (name, chrom, map_info,  genotype, logr_ratio, ballele_freq) in zip(manifest.names, manifest.chroms, manifest.map_infos, genotypes, logr_ratios, ballele_freqs):
            output_handle.write(delim.join([name, chrom, str(map_info), code2genotype[genotype], str(logr_ratio), str(ballele_freq)]) + "\n")
