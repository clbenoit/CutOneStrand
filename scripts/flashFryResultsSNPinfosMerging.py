#!/usr/bin/python3

import re
import pandas as pd
import vcf
import argparse

parser = argparse.ArgumentParser(description='select SNPs with intersting sequence modifications on PAM sites ')
parser.add_argument('--outfile', nargs='?', help='final file to store results in')
parser.add_argument('--mutlist', nargs='?', help='final file to store results in')
parser.add_argument('--flashfry', nargs='?', help='final file to store results in')
args = parser.parse_args()

# open these files 

snp_infos = pd.read_csv(args.mutlist, sep='\t')

snp_infos["contig"] =  snp_infos["mutation_ID"] + "_" + snp_infos['REF'] + "/" + snp_infos['ALT']
#print(snp_infos["contig"])

flashry_scores = pd.read_csv(args.flashfry, sep='\t')

final_results = snp_infos.merge(flashry_scores,how= 'inner',left_on = 'contig',right_on = 'contig')

final_results = final_results[['PAM_sgRNA','cas9_cutting_sequence','mutation_ID','AF','mutation_position','ref_sequence','alternate_sequence','mutation_effect','start','stop','target','context','overflow','orientation','Moreno-Mateos2015OnTarget',	'Doench2014OnTarget','DoenchCFD_maxOT',	'DoenchCFD_specificityscore','dangerous_GC','dangerous_polyT','dangerous_in_genome','Hsu2013','basesDiffToClosestHit','closestHitCount','0-1-2-3-4_mismatch','otCount']]

### 'cas9_cutting_sequence' et cas particulier quand PAM répété plusieurs fois dans l'un des deux alleles

#final_results = final_results.sort_values(by=['DoenchCFD_maxOT'], ascending=True)
final_results = final_results.sort_values(by=['mutation_ID'], ascending=True)
final_results = final_results.loc[final_results['start'] == 10] 

final_results.to_csv(args.outfile, sep ='\t', index = False)





	
