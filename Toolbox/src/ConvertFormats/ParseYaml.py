#!/usr/bin/env python

import yaml, csv, sys

input_file = sys.argv[1]
output_dir = sys.argv[2]

def usage():
    print (
        'usage: ParseYaml.py <file_name.yaml> <output_dir/prefix> \n'
        'Insert file to be processed, and the desired directory and prefix \n')
    sys.exit(2)


def main(input_file, output_dir):

	with open(input_file, 'r') as f:
        	doc = yaml.load(f)
	print(input_file)
	bin = input_file.strip('.yaml').split("/")[-1]

	output_file = output_dir + bin +".out"

	with open( output_file, 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t')
		#write header
		writer.writerow(["bin","id","target","prot_len","contig","contig_len","strand","score","status","reason","seq_upstream","seq_downstream",'dna_end', 'dna_start', 'inframe_stopcodons', 'mismatchlist', 'nucl_end', 'nucl_start', 'overlap', 'prot_end', 'prot_start', 'seq', 'seqshifts', 'translation', 'type', 'undeterminedlist'])
		for i in sorted(doc.keys()):
			for elem in range(len(doc[i])):
				prot_l = doc[i][elem]['prot_len']
				target = doc[i][elem]['target'].replace(" ", "_")
				target_l = doc[i][elem]['target_len']
				strand = doc[i][elem]['strand']
				score = doc[i][elem]['score']
				status = doc[i][elem]['status'].replace(" ", "_")
				reason = doc[i][elem]['reason']
				if reason is None:
					reason = "None"
				else:
					reason = reason.replace(" ","_")
				up = doc[i][elem]['upstream']
				down = doc[i][elem]['downstream']
				id = doc[i][elem]['ID']
				general = [bin, id, i, prot_l, target, target_l, strand, score, status, reason, up, down ]
				for type_gen in doc[i][elem]['matchings']:
					if type_gen["type"] == "intron" or type_gen["type"] == "intron?" or type_gen["type"] == "gap":
						type_gen["prot_start"] = "None"
						type_gen['prot_end'] = "None"
						type_gen['nucl_end'] = "None"
						type_gen['seqshifts'] = "None"
						type_gen['mismatchlist'] = "None"
						type_gen['undeterminedlist'] = "None"
						type_gen['inframe_stopcodons'] = "None"
						type_gen['translation'] = "None"
						type_gen['overlap'] = "None"	
					vals = sorted(type_gen.keys())
					new_entries = list()
					curr_entry = str()       
					for entry in vals:
						curr_entry = type_gen[entry] 
						if type(curr_entry) is list:	
							len_list = len(curr_entry)
							if len_list == 0:
								new_entries.append("None")
							elif len_list == 1:
								new_entries.append(curr_entry[0])
							else: 
								new_entries.append(",".join(map(str,curr_entry)))
						elif type(curr_entry) is dict:
							new_entries.append(",".join('{} : {}'.format(key, value) for key, value in curr_entry[0].items()).replace(" ", ""))
						elif curr_entry is None:
							new_entries.append('None')
						else:
							new_entries.append(curr_entry)
					writer.writerow(general + new_entries)

if __name__ == "__main__":
	main(input_file, output_dir)
