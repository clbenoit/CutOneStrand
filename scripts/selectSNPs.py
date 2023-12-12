#!/usr/bin/python3

import re
import pandas as pd
import vcf
import argparse
import os

parser = argparse.ArgumentParser(description='select SNPs with intersting sequence modifications on PAM sites ')
parser.add_argument('--fasta', nargs='?', help='fasta file containing SNP alleles and their flanking sequences')
parser.add_argument('--vcf', nargs='?', help='candidates in vcf format')
parser.add_argument('--mutlist', nargs='?', help='List of interesting mutations according to selected PAM')
parser.add_argument('--candidatesfasta', nargs='?', help='fasta files containing genimic cintext for intesrting sites')
parser.add_argument('--cas', nargs='?', help='cas9 to scan sites for')
parser.add_argument('--gene', nargs='?', help='gene')
args = parser.parse_args()

pam_dict = {}
pam_dict["spcas9ngg"] = ["GG", "CC"]
#pam_dict["sa_21"] = ["G(A/G)(A/G)T","A(T/C)(T/C)C"]
#pam_dict["cpf1"] = ["TTT(A/C/G)","(T/C/G)AAA"]
#print(pam_dict)
#print(args.cas)
#print(type(args.cas))

pam_transformed = re.sub("\/","|",re.sub("\(","[",re.sub("\)","]",pam_dict[args.cas][0])))
complementary_pam_transformed = re.sub("\/","|",re.sub("\(","[",re.sub("\)","]",pam_dict[args.cas][1])))

#print(pam_transformed)
#print(complementary_pam_transformed)

## Load fasta file into dictionnary
motifs = dict()
with open(args.fasta) as f:
	lines = f.readlines()
	for i in range(0, len(lines)):
		s = lines[i].strip()
		if s[0] == '>':
			key = s[1:]
		else:
			motifs[key] = s

#print(motifs)

results_table = pd.DataFrame(columns = ["mutation_position","ref_sequence","alternate_sequence","mutation_effect","cas9_cutting_sequence","ALT","REF"])

for sequence_id, sequence in motifs.items():

	centered_sequence = re.findall(pattern = ".\[.*\].", string = sequence, flags=re.I)[0]
	ref_sequence = re.sub("\/.*\]","",re.sub("\[","",centered_sequence))
	alternate_sequence = re.sub("\[.*\/","",re.sub("\]","",centered_sequence))
	cinq_prime_flanking_sequence = re.sub("\[.*","",sequence)
	trois_prime_flanking_sequence = re.sub(".*\]","",sequence)

	pos = re.sub("chr19:","",sequence_id)
	if ( (bool(re.search(pam_transformed, ref_sequence)) is False) & (bool(re.search(pam_transformed, alternate_sequence)) is True) ):
		# initialize data of lists.
		data = {'mutation_position': [pos],
		        'ref_sequence': [ref_sequence],
			'alternate_sequence' : [alternate_sequence],
			'mutation_effect' : ['creates_cutting_site'],
			'cas9_cutting_sequence' : [cinq_prime_flanking_sequence + alternate_sequence + trois_prime_flanking_sequence[0:10]],
			"REF" : [ref_sequence[1:-1]],
			'ALT' : [alternate_sequence[1:-1]],
			'PAM_sgRNA' : [cinq_prime_flanking_sequence + alternate_sequence ]
		}
		results_table = pd.concat((results_table,pd.DataFrame(data)),axis = 0)

	elif ( (bool(re.search(complementary_pam_transformed, ref_sequence)) is False) & (bool(re.search(complementary_pam_transformed, alternate_sequence)) is True) ):
		data = {'mutation_position': [pos],
		        'ref_sequence': [ref_sequence],
			'alternate_sequence' : [alternate_sequence],
			'mutation_effect' : ['creates_cutting_site'],
			'cas9_cutting_sequence' : [cinq_prime_flanking_sequence[-10:] + alternate_sequence + trois_prime_flanking_sequence ],
			"REF" : [ref_sequence[1:-1]],
			'ALT' : [alternate_sequence[1:-1]],
			'PAM_sgRNA' : [ alternate_sequence + trois_prime_flanking_sequence ]
		}
		results_table = pd.concat((results_table,pd.DataFrame(data)),axis = 0)

	elif ( (bool(re.search(pam_transformed, ref_sequence)) is True) & (bool(re.search(pam_transformed, alternate_sequence)) is False) ):
		data = {'mutation_position': [pos],
		        'ref_sequence': [ref_sequence],
			'alternate_sequence' : [alternate_sequence],
			'mutation_effect' : ['deleted_cutting_site'],
			'cas9_cutting_sequence' : [cinq_prime_flanking_sequence + ref_sequence + trois_prime_flanking_sequence[0:10]],
			"REF" : [ref_sequence[1:-1]],
			'ALT' : [alternate_sequence[1:-1]],
			'PAM_sgRNA' : [cinq_prime_flanking_sequence + ref_sequence ]
		}
		results_table = pd.concat((results_table,pd.DataFrame(data)),axis = 0)	

	elif ( (bool(re.search(complementary_pam_transformed, ref_sequence)) is True) & (bool(re.search(complementary_pam_transformed, alternate_sequence)) is False) ):
		data = {'mutation_position': [pos],
		        'ref_sequence': [ref_sequence],
			'alternate_sequence' : [alternate_sequence],
			'mutation_effect' : ['deleted_cutting_site'],
			'cas9_cutting_sequence' : [cinq_prime_flanking_sequence[-10:] + ref_sequence + trois_prime_flanking_sequence ],
			"REF" : [ref_sequence[1:-1]],
			'ALT' : [alternate_sequence[1:-1]],
			'PAM_sgRNA' : [ref_sequence + trois_prime_flanking_sequence ]
		}
		results_table = pd.concat((results_table,pd.DataFrame(data)),axis = 0)

results_table['mutation_position'] = results_table['mutation_position'].astype('string')
results_table['ALT'] = results_table['ALT'].astype('string')


results_table.to_csv(os.path.dirname(args.candidatesfasta) + "/" +"resultstableraw.tsv", sep ='\t', index = False)

name = os.path.basename(args.fasta).replace("_SNPS_CANDIDATES_hg38_WTH_FLANKING_SEQUENCES.fasta", "")
vcf = vcf.Reader(open(os.path.dirname(args.candidatesfasta) + '/' + name + '_SNPS_CANDIDATES_hg38.vcf', 'r'))
candidates_vcf = pd.DataFrame(columns = ["mutation_position","mutation_ID","AF"])
for record in vcf:
	#record = {"mutation_position" : record.POS, "mutation_ID" : record.ID + "_" + record.REF + "/" + record.ALT, "AF" : record.INFO["AF"]}
	record = {"mutation_position" : record.POS, "mutation_ID" : record.ID, "AF" : record.INFO["AF"], "ALT" : record.ALT}
	candidates_vcf = pd.concat((candidates_vcf,pd.DataFrame(record)), axis= 0)	

candidates_vcf["mutation_position"] = candidates_vcf["mutation_position"].astype('string')
candidates_vcf["ALT"] = candidates_vcf["ALT"].astype('string')

results_table = results_table.merge(candidates_vcf,how= 'inner',left_on = ['mutation_position',"ALT"], right_on = ['mutation_position',"ALT"])

results_table.to_csv(args.mutlist, sep ='\t', index = False)


with open(args.candidatesfasta, "w") as fasta_file:
	for index, row in results_table.iterrows():
#		print(row["mutation_position"],row["cas9_cutting_sequence"])
		#fasta_file.write(">" + row["mutation_ID"] + "_" + row["mutation_position"] + '\n' + row["cas9_cutting_sequence"] + "\n")
		if ((row["mutation_ID"] is not None) & (row["cas9_cutting_sequence"]  is not None)):
			fasta_file.write(">" + row["mutation_ID"] + "_" + row['ref_sequence'][1:-1] + "/" + row['alternate_sequence'][1:-1]  + '\n' + row["cas9_cutting_sequence"] + "\n")
			#fasta_file.write(">" + row["mutation_ID"] + '\n' + row["cas9_cutting_sequence"] + "\n")
#



	
