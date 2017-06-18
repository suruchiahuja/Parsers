import sys
import argparse
import os
import subprocess
#import vcf
from datetime import date
from IlluminaBeadArrayFiles import GenotypeCalls, BeadPoolManifest, code2genotype

delim = " "

class ManifestWithReference(object):
    def __init__(self, manifest_path):
        self.chroms = []
        self.names = []
        self.positions = []
        self.ref_alleles = []
        self.snps = []
        self.source_strands = []
	self.ballele_freqs = []
	self.logr_ratios = []

        with open(manifest_path) as f:
            for line in f:
                line = line.rstrip()
                tokens = line.split(' ')
                self.chroms.append(tokens[0])
                self.positions.append(tokens[1])
                self.names.append(tokens[2])
                self.ref_alleles.append(tokens[3])
                self.snps.append(tokens[4])
                self.source_strands.append(int(tokens[5]))

def get_alt_allele(genotypes, ref):
    for g in genotypes:
        alleles = list(g)
        if alleles[0] != ref:
            return alleles[0]
        elif alleles[1] != ref:
            return alleles[1]
    return '.'
    
def get_sample_calls(sample_ids,genotypes, ref, alt):
    samples = []
    for (sample, g) in zip(sample_ids, genotypes):
        alleles = list(g)
        if len(alleles) == 1: #no call
            samples.append('./.')
        elif alleles[0] == alleles[1]:
            if ref in alleles:
                samples.append('0/0')
            elif alt in alleles:
                samples.append('1/1')
            else:
                samples.append('./.')

        else:
            samples.append('0/1')
    
    return samples


def write_vcf_record(writer, chromosome, position, name, ref, genotypes,
        sample_indexes, sample_ids):
    alt = get_alt_allele(genotypes, ref)
    if alt == 'I' or alt == 'D':
        return
    
    samples = get_sample_calls(sample_ids, genotypes, ref, alt)
    record = "\t".join([chromosome, str(position), name, ref, alt, '.','.','.','GT']) + "\t" + "\t".join(samples) + "\n"
    writer.write( record)

def write_vcf_header(writer, sample_ids):
    writer.write("##fileformat=VCFv4.2\n")
    d = date.today().strftime("%Y%m%d")
    writer.write("##fileDate=" + d + "\n")
    writer.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n")
    writer.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + "\t".join(sample_ids) + "\n") 

def create_vcf_file(vcf_path, sample_ids, gtc_files, manifest):
    sample_indexes = {}
    for idx, sample in enumerate(sample_ids):
        sample_indexes[sample] = idx

    with open(vcf_path, 'w') as writer:
        def get_genotypes(gtc_file):
            gtc = GenotypeCalls(gtc_file)
            return gtc.get_base_calls_forward_strand(manifest.snps,
                    manifest.source_strands)
        gtc_genotypes = list(map(lambda g: get_genotypes(g), gtc_files))

        write_vcf_header(writer, sample_ids)

        for data in zip(manifest.chroms,
                manifest.positions, manifest.names, manifest.ref_alleles,
                *gtc_genotypes):
            #out_list = [chromosome, position, name, ref]
            (chromosome, position, name, ref) = data[0:4]
            genotypes = data[4:]
            if chromosome != '0':
                write_vcf_record(writer, chromosome, int(position), name, ref, genotypes, sample_indexes, sample_ids)
            
            

#        for (sample, gtc_file) in zip(sample_ids, gtc_files):
#            gtc = GenotypeCalls(gtc_file)
#            genotypes = gtc.get_base_calls_forward_strand(manifest.snps,
#                    manifest.source_strands)
#            assert len(genotypes) == len(manifest.names)
#
#            for (chromosome, position, name, genotype) in zip(manifest.chroms,
#                    manifest.map_infos, manifest.names, genotypes):
#                allele_1 = 0
#                allele_2 = 0
#                if (genotype != "-"):
#                    genotypes = list(genotype)
#                    allele_1 = genotypes[0]
#                    allele_2 = genotypes[1]
#                lgen.write(delim.join([str(sample), str(sample),
#                    name, str(allele_1), str(allele_2)]) + "\n")
#

parser = argparse.ArgumentParser("Generate a vcf file from a directory of GTC files")
parser.add_argument("--group", default='test', help="base name for output files")
parser.add_argument("manifest", help="manifest file with reference allele")
parser.add_argument("gtcs", nargs='*',help="list of gtc files to process")

args = parser.parse_args()
samples = [os.path.splitext(os.path.basename(gtc))[0] for gtc in args.gtcs]
base = args.group

try:
    manifest = ManifestWithReference(args.manifest)
except:
    sys.stderr.write("Failed to read data from manifest\n")
    sys.exit(-1)

vcf_path = base + '.vcf'
create_vcf_file(vcf_path, samples, args.gtcs, manifest)

#lgen_path = base + '.lgen'
#create_lgen_file(lgen_path, samples, args.gtcs, manifest)

