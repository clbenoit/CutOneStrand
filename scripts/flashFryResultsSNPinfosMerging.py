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



#-----BEGIN PGP MESSAGE-----
#Version: Keybase OpenPGP v2.0.76
#Comment: https://keybase.io/crypto

#wcBMA4K4B8OVkl7jAQf/VbZgjDSBvofjTpDZGewMp6lIek+amHLT47zexq1RSR3L
#z9ljmk4O7ZwO2KYah6a9rNh0dwqrQQ+8uuz+32MWEuHPuX51BKdoGPkkXjTbsou7
#raess8oTz9Vu44RWNtq2UAsOqcwc49JErxMcUBgMKwuMmmrvPUzqWOosSTm0yqlf
#lMSCoQl/8OaAmgTpZhuN/mMxph+hP8NMNkc7WUA32vWbyXBOMjxB3BZ6cAS6906e
#oT0POJU1WtWHehQZr487kj3thmibWs23EAt7rxdM2EZQ+FE2Y2Djv0qcHwLTfUsp
#h/JFUagttg6we/EWE8ZvMLBca5T08+mA5eLGsC0U3NLCPAFhDRbDtMkE64VriBbB
#g179Q3XCz32RJbPY+2/TX+6JmUGND9DZ5HReRO8MdH8jlsBz3YmQjjcssb7Jo7As
#6Ps1rq+eZbkALolsEvuFFnERmnRtzJ3VBLWbmirnClXFRuyydcj8jz4Qdck1IaLU
#RftT5dSLzAoPSvy6BPvN3EQqNJfrUrackIiP1QNZvcIanjeIzUanYWcENvD0148P
#AIHon+9UxXI36clJQ3pTeRAs33IOx2M9X7yDXeWXs9a2rH/+LPX2zX2GEKrA1yj2
#NdPKXacOP1g4nCBxwgF3tpZaI7Lcyfnsep63agZBN40zmNGThZlI0kyxqRw3RLQv
#0rdvAYANX/1yJEM8rjsld+4rPWXd84hJF0eOlWvLTG0WryQzm1dc6f7mXvDeELQg
#vHOdYV3Xp+SpYVqdju3zTZLG31Ff1vWUIuFG1O3FBBzY/Gp4u6/TznG4OjHkvsbd
#sUzw9qipAzGeoVuqIXjHRRQGwjR4QT6YkHmg2mpBdTAXa/0ZX6LWDlAJruBZOiqg
#tC2Pwsg+FckaKE9uLa/IV7R9K3yUxqGvABLyV8zEoN2Bo/Iz9qP8LvUZ+ji2IdkE
#nam5nJYZxJlt+SKAoXCz/UX9iNzmRVM52v2ojFpoDAj7NZyU5JqtYLakZAK4SkVx
#1JmefoUpmaRl+q8P0mN1C7lDcWiZPwLxEzvp/a2Dl/wfT5TllFehxWTCWrusoWef
#kcdV2hLuxyUNX8ZUmSyxwEhh4u+gYylyG1b8FwdMSq/h2iSOvx1MyAQeSC6xt/lN
#LrE1pPdMGwpj5m2BSd0JK2jOJyzo9LIGAifLPKgB8VzCQkAE9Wf0XcGJKCa0SITC
#va07k/uPbrEjMN9cUXTxo98EO/YMkPULYAuYOVEjiEEA2hrpkUrvO4fhjzZaA4uA
#4Ke0NUbaD96OPc/La5YSfq6vWE/0489BGeNAh83j0wMXefzywbrjmi0nqqIp+Xym
#UxTGAQZueBrTibQT2ZVOcEi+EJH9k5jVPyC2ntxl
#=V0Lk
#-----END PGP MESSAGE-----

	
