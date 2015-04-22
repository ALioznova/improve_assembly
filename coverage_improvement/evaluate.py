from sets import Set
import numpy
import os
import sys
import argparse
import subprocess
import numpy
from Bio import SeqIO
from Bio.Seq import Seq

#TODO:
# find correct chromosome (scaffolds and target)

def build_alignment_bwa(bwa_path, data_name, ref_file, contigs_file):
	subprocess.call([os.path.join(bwa_path, "bwa"), "index", ref_file])
	with open(data_name + "_aligned.sam", "w") as sam_file:
		subprocess.call([os.path.join(bwa_path, "bwa"), "mem", ref_file, contigs_file], stdout = sam_file)
	os.remove(ref_file + ".amb")
	os.remove(ref_file + ".ann")
	os.remove(ref_file + ".bwt")
	os.remove(ref_file + ".pac")
	os.remove(ref_file + ".sa")
	return data_name + "_aligned.sam"

def get_alignment_length(cigar): 
	if cigar == "*":
		return 0
	cigar_parsed = []
	i = 0
	count = 0
	while i < len(cigar):
		while (cigar[i + count]).isdigit():
			count += 1
		cigar_parsed.append((int(cigar[i:i+count]), cigar[i+count]))
		i += count + 1
		count = 0
	end = 0
	for (num, act) in cigar_parsed:
		if (act == "M") or (act == "D") or (act == "N")  or (act == "X")  or (act == "="):
			end += num
	return end

def parse_sam_record(line):
	if not line.startswith('$') and not line.startswith('@'):
		line = line.strip().split()
		record = {
		'QNAME' : line[0],
		'FLAG'  : int(line[1]),
		'RNAME' : line[2],
		'POS'   : int(line[3]),
#		'MAPQ'  : int(line[4]),
		'CIGAR' : line[5],
#		'RNEXT' : line[6],
#		'PNEXT' : int(line[7]),
#		'TLEN'  : int(line[8]),
#		'SEQ'   : line[9],
#		'QUAL'  : line[10],
#		'optional' : []
		}
#		for optional in line[11:]:
#			record['optional'].append(optional.split(':'))
		return record

def sequence_unmapped(record):
	# 0x4 segment unmapped, 0x8 next segment in the template unmapped
	return (record ['FLAG'] & int(0x4)) or (record ['FLAG'] & int(0x8))

def supplementary_alignment(record):
	# 0x800 supplementary alignment
	return ((record ['FLAG']) & int(0x800))

def seq_reverse_complemented(record):
	# 0x10 reverse complemented
	return ((record ['FLAG']) & int(0x10))

def process_sam_file(sam_file, output_file_name, target_len):
	number = 0
	mapped_number = 0
	print
	f_out = open(output_file_name, 'w')
	for line in open(sam_file):
		record = parse_sam_record(line)
		if record:
			number += 1
			if sequence_unmapped(record) or supplementary_alignment(record):
				continue
			mapped_number += 1
			begin = record['POS']
			end = record['POS'] + get_alignment_length(record['CIGAR'].strip())
			ref_name = record['RNAME']
			seq_name = record['QNAME']
			strand = not seq_reverse_complemented(record)
			f_out.write(ref_name + ':\t' + (["-", "+"][strand]) + seq_name + '\t' + str(begin) + '\n')
			rc_begin = target_len[ref_name] - end
			f_out.write(ref_name + ':\t' + (["-", "+"][not strand]) + seq_name + '\t' + str(rc_begin) + '\n')
	f_out.close()
	print "Total number of sequences is: ", number
	print "Number of mapped sequences is: ", mapped_number

def process_ref_file(ref_file):
	handle = open(ref_file, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	target_len = {}
	for elem in records:
		target_len[elem.name] = len(elem.seq)
	return target_len

if __name__ == "__main__":
	if len(sys.argv) == 1:
		print "Usage:", sys.argv[0], "-b <path to bwa> -u <unused contigs> -t <target genome> -c <contigs coords>"
		print "Please use the --help option to get more usage information."
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Evaluate insertion unused contigs.')
	parser.add_argument("-b", "--bwa", help="path to bwa", required=True)
	parser.add_argument("-u", "--unused", help="unused contigs", required=True)
	parser.add_argument("-t", "--target", help="target genome", required=True)
	parser.add_argument("-c", "--coords", help="contigs coords", required=True)

	args = parser.parse_args()
	bwa_path = args.bwa
	unused_contigs_path = args.unused
	target_path = args.target
	contigs_coords_path = args.coords

	data_name = os.path.join(os.path.dirname(unused_contigs_path), 'unused_to_target')
	sam_file = build_alignment_bwa(bwa_path, data_name, target_path, unused_contigs_path)
	target_len = process_ref_file(target_path)
	output_file_name = os.path.join(os.path.dirname(unused_contigs_path), 'true_coords.txt')
	process_sam_file(sam_file, output_file_name, target_len)
	
	print
	print '=============================='
	print 'Result can be found in', output_file_name

