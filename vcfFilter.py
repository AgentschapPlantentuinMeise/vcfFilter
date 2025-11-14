#!usr/bin/python3

###############################
# Written by Sander de Backer #
# Crop Wild Relatives Group   #
# Meise Botanic Garden        #
# 2023                        #
###############################

import os, sys, argparse, subprocess, glob, time, statistics
from datetime import datetime
from argparse import RawTextHelpFormatter

#############
# ARGUMENTS #
#############
# Create an ArgumentParser object
parser = argparse.ArgumentParser(description = 'Filter VCF file based on different criteria.\nAlleles are automatically restricted to being biallelic.\nAnd non-variants are removed after applying all other filters.', formatter_class=RawTextHelpFormatter)

# VCF
parser.add_argument('--prefix', type = str, default = 'VCF',
                        help = 'Prefix of the VCF file.')

# VCF filter options
parser.add_argument('--minmeanDP', default = 5, type = int,
                        help = 'The minimum mean read depth supporting a SNP variant.\nDefault = 5')
parser.add_argument('--minAD', default = 0, type = int,
                        help = 'The minimum allele depth for a REF or ALT allele at a given position.\nDefault = 0')
parser.add_argument('--maf', default = 0.0, type = float,
                        help = 'Minimum minor allele frequency. SNPs with a lower minimal minor allele frequency are discarded.\nDefault = 0.0')
parser.add_argument('--minQ', default = 20, type = int,
                        help = 'The minimum quality of the variant site.\nDefault = 20')
parser.add_argument('--minGQ', default = 0, type = int,
                        help = 'The minimum genotype quality of the call.\nDefault = 0')
parser.add_argument('-c', '--completeness', type = float, default = 0.0,
                        help = 'SNPs with a lower completeness level are discarded.\nDefault = 0')

# Additional options
parser.add_argument('-v', '--verbose', action = 'store_true', 
                        help = 'Print additional information to the console.\n*Potentially slows down certain parts of the analysis down*\nWARNING: This will literally flood the console...')

# Parse arguments to a dictionary
args = vars(parser.parse_args())

#############
# FUNCTIONS #
#############

def restrict_biallelic(i = args['prefix']):
    '''Filter VCF based in allelicity of calls.'''
    vcfin = open(i + '.vcf')
    vcfout = open('temp.vcf', 'w+')
    # Start counter for total and kept SNPs
    total, kept = 0, 0
    # Select SNPs by looking at ALT column
    for line in vcfin:
        if line.startswith('##') or line.startswith('#CHROM'):
            print(line, end = '', file = vcfout)
        else:
            total += 1
            # Specify line format as tab-seperated and divide
            line = line.rstrip().split('\t')
            # ALT allele in column 5
            variant = line[4]
            vlen = len(variant)
            if variant != '.' and vlen == 1:
                kept += 1
                print('\t'.join(e for e in line), file = vcfout)
            else:
                if v:
                    print('{} removed, not a SNP'.format(variant))
    sys.stdout.write('\n{}\tAfter removing non-SNPs, kept {} out of a possible {} sites.\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), kept, total))

def vcf_filter():
    '''Filter the VCF file based on genotype calls.'''
    vcfin = open('temp.vcf')
    vcfout = open('temp_filtered.vcf', 'w+')
    sys.stdout.write('\n{}\tFiltering genotype calls with criteria minmeanDP={}, minAD={}, maf={}, minQ={}, and minGQ={} ....\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), args['minmeanDP'], args['minAD'], args['maf'], args['minQ'], args['minGQ']))
    # Specify FORMAT field tags for genotype calls
    tags = ['GT', 'AD', 'DP', 'GQ']
    # Start count for total variants and number kept after filter
    total, kept = 0, 0
    # Copy header lines and iterate over variant lines
    for line in vcfin:
        if line.startswith('##') or line.startswith('#CHROM'):
            print(line, end = '', file = vcfout)
        else:
            # Total variants increases with 1 with each line
            total += 1
            # Specify line format as tab-seperated and divide
            line = line.rstrip().split('\t')
            # For site quality
            qual = float(line[5])
            if qual > args['minQ']:
                # Start count for reference allele and alternate allele
                ref_num, alt_num = 0, 0
                #Find index for GT, AD, DP, GQ and PL tags in FORMAT field (=ninth column, or line[8])
                criteria = {tag:line[8].split(':').index(tag) for tag in tags}
                # Observe each genotype (gt) for an increasing j with each loop
                for j, gt in enumerate(line[9:]):
                    # Split the genotype call on the colon character
                    gt = gt.split(':')
                    # If statement based on genotype length (./. = 3, all else is longer)
                    if len(gt) > 3:
                        #Check if the total read depth and the genotype quality are not lower than the predefined minima.
                        if int(gt[criteria['DP']]) >= args['minmeanDP'] and float(gt[criteria['GQ']]) >= args['minGQ']:
                            # Count reference alleles and alternate alleles
                            AD_ref = (gt[criteria['AD']].split(','))[0]
                            AD_alt = (gt[criteria['AD']].split(','))[1]
                            # Check if reference or alternate allele is present and the allele depths are not lower than the predefined minima
                            if ('0' not in gt[criteria['GT']] or int(AD_ref) >= args['minAD']) and ('1' not in gt[criteria['GT']] or int(AD_alt) >= args['minAD']):
                                # Add to reference or alternate allele count for minor allele frequency later on
                                ref_num += gt[criteria['GT']].split('/').count('0')
                                alt_num += gt[criteria['GT']].split('/').count('1')
                                # Define the normal format of a vcf genotype call and specify as i-th column corresponding to sample position (column10=[9], then add i)
                                call = '{}:{}:{}:{}'.format(gt[criteria['GT']], gt[criteria['AD']], gt[criteria['DP']], gt[criteria['GQ']])
                                line[j + 9] = str(call)
                            else:
                                line[j + 9] = './.'
                        else:
                            line[j + 9] = './.'
                # Only retain sites that are still polymorphic after filtering and with a minor allele frequency not lower than the predefined minimum.
                if ref_num > 0 and alt_num > 0 and min(ref_num, alt_num)/(ref_num + alt_num) >= args['maf']:
                    print('\t'.join(e for e in line), file = vcfout)
                    kept += 1
    sys.stdout.write('\n{}\tAfter filtering genotype calls, kept {} out of a possible {} sites.\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), kept, total))

def vcf_env(i = args['prefix'], c = args['completeness']):
    '''Remove non-variant sites.'''
    vcfin = open('temp_filtered.vcf')
    vcfout = open(i + '_filtered.vcf', 'w+')
    sys.stdout.write('\n{}\tExcluding non-variant calls....\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    total, kept = 0, 0
    for line in vcfin:
        if line.startswith('##') or line.startswith('#CHROM'):
            print(line, end = '', file = vcfout)
        elif not line.startswith('##'):
            total += 1
            line = line.rstrip().split('\t')
            # Start count for each genotype
            ref, het, alt, missing = 0, 0, 0, 0
            for i, genotype in enumerate(line[9:]):
                genotype = genotype.split(':')[0].split('/')
                if all(x == '0' for x in genotype):
                    ref += 1
                elif all(x == '.' for x in genotype):
                    missing += 1
                elif all(x == '1' for x in genotype):
                    alt += 1
                else:
                    het += 1
                line[9 + i] = '/'.join(e for e in genotype)
            if ref + missing < len(line[9:]) and alt + missing < len(line[9:]) and (missing / (missing + ref + alt + het)) <= (1 - c):
                print('\t'.join(e for e in line), file = vcfout)
                kept += 1
    sys.stdout.write('\n{}\tAfter removing non-variants, kept {} out of a possible {} sites.\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), kept, total))

############
# ANALYSIS #
############

restrict_biallelic()
vcf_filter()
vcf_env()

for temp in glob.glob('temp*vcf'):
    subprocess.run(['rm', '-r', temp])

sys.stdout.write('\n**DONE!\n')
