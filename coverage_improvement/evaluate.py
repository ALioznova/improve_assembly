from sets import Set
import numpy
import os
import sys
import argparse
import subprocess
import numpy
from Bio import SeqIO
from Bio.Seq import Seq

def build_alignment_bwa(bwa_path, data_name, ref_file, contigs_file):
	subprocess.call([os.path.join(bwa_path, "bwa"), "index", ref_file])
	with open(data_name + "_aligned.sam", "w") as sam_file:
		subprocess.call([os.path.join(bwa_path, "bwa"), "mem", ref_file, contigs_file, '-a'], stdout = sam_file)
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
	return bool(record ['FLAG'] & int(0x4)) or (record ['FLAG'] & int(0x8))

def supplementary_alignment(record):
	# 0x800 supplementary alignment
	return bool((record ['FLAG']) & int(0x800))

def seq_reverse_complemented(record):
	# 0x10 reverse complemented
	return bool((record ['FLAG']) & int(0x10))

def process_sam_file_contigs(sam_file, target_len):
	number = 0
	mapped_number = 0
	alignment = {}
	print
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
			rc_begin = target_len[ref_name] - end
			rc_end = target_len[ref_name] - begin
			if not alignment.has_key((["-", "+"][strand]) + seq_name):
				alignment[(["-", "+"][strand]) + seq_name] = []
			alignment[(["-", "+"][strand]) + seq_name].append((ref_name, begin, end, False))
			if not alignment.has_key((["-", "+"][not strand]) + seq_name):
				alignment[(["-", "+"][not strand]) + seq_name] = []
			alignment[(["-", "+"][not strand]) + seq_name].append((ref_name, rc_begin, rc_end, True))
	print "Total number of sequences is: ", number
	print "Number of mapped sequences is: ", mapped_number
	return alignment

def process_ref_file(ref_file):
	handle = open(ref_file, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	target_len = {}
	for elem in records:
		target_len[elem.name] = len(elem.seq)
	return target_len

def process_contigs_coords_names(contigs_coords_path):
	insertion = {}
	scaff_name = None
	for line in open(contigs_coords_path):
		if line.startswith('>'):
			scaff_name = (line[1:]).strip()
			continue
		(contig_name, coords, cont_len, is_new) = line.split()
		if is_new == '-':
			continue
		coords = int(coords)
		if not insertion.has_key(contig_name):
			insertion[contig_name] = []
		insertion[contig_name].append((scaff_name, coords))
	return insertion

def compare_alignment_and_insertion(alignment, insertion, scaff_len, output_file_name):
	f_out = open(output_file_name, 'w')
	for (contig_name, coords_list) in insertion.iteritems():
		f_out.write(contig_name + '\t')
		f_out.write('tool:')
		for (scaff_name, coords) in coords_list:
			f_out.write(scaff_name + '\t' + str(coords) + '\t')
		f_out.write('\n' + contig_name + '\t')
		f_out.write('real:')
		if not alignment.has_key(contig_name):
			f_out.write('-\n')
		else:
			for ((ref_name, begin, end, rc)) in alignment[contig_name]:
				f_out.write(ref_name + ' ' + (["[  ]", "[RC]"][rc]) + '\t' + str(begin) + '\t')
			f_out.write('\n')
		f_out.write('\tdiff:\t')
		if not alignment.has_key(contig_name):
			f_out.write('-\n')
		else:
			for ((ref_name, begin, end, rc)) in alignment[contig_name]:
				for (scaff_name, coords) in coords_list:
					f_out.write([str(begin - coords), str(scaff_len[scaff_name] + begin - coords)][begin - coords < 0] + '\t')
			f_out.write('\n')
	f_out.close()

def process_contigs_coords_coords(contigs_coords_path, scaff_len):
	insertion = {}
	scaff_name = None
	for line in open(contigs_coords_path):
		if line.startswith('>'):
			scaff_name = (line[1:]).strip()
			insertion[scaff_name] = [None] * scaff_len[scaff_name]
			continue
		(contig_name, coords, cont_len, is_new) = line.split()
		coords = int(coords)
		cont_len = int(cont_len)
		assert insertion[scaff_name][coords] == None
		if is_new == '-':
			is_new = False
		elif is_new == '+':
			is_new = True
		insertion[scaff_name][coords] = (contig_name, cont_len, is_new)
	return insertion

def parse_cigar(cigar, beg_pos):
	if cigar == "*":
		return None
	cigar_parsed = []
	i = 0
	count = 0
	while i < len(cigar):
		while (cigar[i + count]).isdigit():
			count += 1
		cigar_parsed.append((int(cigar[i:i+count]), cigar[i+count]))
		i += count + 1
		count = 0
	matching_coords_list = []
	target_pos = beg_pos
	scaff_pos = 0
	for (num, act) in cigar_parsed:
		if act == 'M' or act == '=':
			matching_coords_list.append((target_pos, scaff_pos, num))
			target_pos += num
			scaff_pos += num
		elif act == 'I' or act == 'S' or act == 'H':
			scaff_pos += num
		elif act == 'D' or act == 'N':
			target_pos += num
		elif act == 'P':
			pass
		elif act == 'X':
			target_pos += num
			scaff_pos += num
		else:
			return None
	return matching_coords_list

def process_sam_file_scaffolds(sam_file_scaffolds):
	alignment = {}
	for line in open(sam_file_scaffolds):
		record = parse_sam_record(line)
		if record:
			if sequence_unmapped(record) or supplementary_alignment(record):
				continue
			begin = record['POS']
			ref_name = record['RNAME']
			seq_name = record['QNAME']
			strand = not seq_reverse_complemented(record)
			matching_list = parse_cigar(record['CIGAR'].strip(), begin)
			if not alignment.has_key(seq_name):
				alignment[seq_name] = []
			alignment[seq_name].append((ref_name, strand, matching_list))
	return alignment

def compare_alignment_and_insertion_sam_scaffolds(alignment_scaffolds, insertion_coords, evaluation_scaffolds_filename, scaff_len):
	contigs_insertion = {}
	f_out = open(evaluation_scaffolds_filename, 'w')
	for (scaff_name, scaff_alignment_list) in alignment_scaffolds.iteritems():
		cur_scaff_len = scaff_len[scaff_name]
		contigs_insertion[scaff_name] = [None] * len(insertion_coords[scaff_name])
		for (ref_name, scaff_strand, matching_list) in scaff_alignment_list:
			for (ref_pos, scaff_pos, matching_len) in matching_list:
				scaff_beg = scaff_pos
				scaff_end = scaff_pos + matching_len
				if not scaff_strand:
					scaff_beg = cur_scaff_len - (scaff_pos + matching_len)
					scaff_end = cur_scaff_len - scaff_pos
				contigs_insertion[scaff_name][scaff_beg] = (ref_pos, ref_name, scaff_end, scaff_strand)
		f_out.write('>' + scaff_name + '\n')
		for i in xrange(cur_scaff_len):
			if contigs_insertion[scaff_name][i] != None:
				(ref_pos, ref_name, scaff_end, scaff_strand) = contigs_insertion[scaff_name][i]
				scaff_beg = i
				f_out.write(ref_name + '\t' + str(ref_pos) + '\tscaff: ' + str(scaff_beg) + '-' + str(scaff_end) + '\t' + (["[RC]", "[  ]"][scaff_strand]) + '\n')
				contig_beg = i
				while contig_beg >= 0 and insertion_coords[scaff_name][contig_beg] == None:
					contig_beg -= 1
				(contig_name, cont_len, is_new) = insertion_coords[scaff_name][contig_beg]
				if contig_beg + cont_len >= scaff_beg:
					f_out.write('\t' + contig_name + '\t' + str(contig_beg) + '-' + str(contig_beg + cont_len) + '\t' + (["old", "new"][is_new]) + '\n')
				contig_beg += 1
				while contig_beg < cur_scaff_len and contig_beg <= scaff_end:
					if insertion_coords[scaff_name][contig_beg] != None:
						(contig_name, cont_len, is_new) = insertion_coords[scaff_name][contig_beg]
						f_out.write('\t' + contig_name + '\t' + str(contig_beg) + '-' + str(contig_beg + cont_len) + '\t' + (["old", "new"][is_new]) + '\n')
					contig_beg += 1
				
	f_out.close()

def build_alignment_mummer(mummer_path, data_name, ref_file, qry_file):
	subprocess.call([os.path.join(mummer_path, "nucmer"), '-p', data_name + '_nucmer', ref_file, qry_file])
	with open(data_name + '_nucmer.coords', "w") as coords_file:
		subprocess.call([os.path.join(mummer_path, "show-coords"), '-qldc' , data_name + '_nucmer.delta'], stdout = coords_file)
	return data_name + '_nucmer.coords'

def process_coords_file_scaffolds(coords_file_scaffolds):
	alignment = {}
	with open(coords_file_scaffolds) as fp:
		fp.readline()
		fp.readline()
		fp.readline()
		fp.readline()
		fp.readline()
		while True:
			line = fp.readline()
			if line == '':
				break
			info = (line.strip()).split()
			s1 = int(info[0])
			e1 = int(info[1])
			s2 = int(info[3])
			e2 = int(info[4])
			len1 = int(info[6])
			len2 = int(info[7])
			idy = float(info[9])
			lenr = int(info[11])
			lenq = int(info[12])
			covr = float(info[14])
			covq = float(info[15])
			strand_r = [True, False][int(info[17]) == -1]
			strand_q = [True, False][int(info[18]) == -1]
			ref_name = info[19]
			scaff_name = info[20]
			if not alignment.has_key(scaff_name):
				empty_list = []
				for i in xrange(lenq):
					empty_list.append([])
				alignment[scaff_name] = empty_list
			alignment[scaff_name][min(s2, e2)].append((min(s1, e1), ref_name, max(s2, e2), [False, True][strand_r == strand_q]))
	return alignment

def compare_alignment_and_insertion_coords_scaffolds(alignment_scaffolds, insertion_coords, evaluation_scaffolds_filename):
	f_out = open(evaluation_scaffolds_filename, 'w')
	for (scaff_name, contigs_insertion) in alignment_scaffolds.iteritems():
		f_out.write('>' + scaff_name + '\n')
		for i in xrange(len(contigs_insertion)):
			for (ref_pos, ref_name, scaff_end, scaff_strand) in contigs_insertion[i]:
				scaff_beg = i
				f_out.write(ref_name + '\t' + str(ref_pos) + '\tscaff: ' + str(scaff_beg) + '-' + str(scaff_end) + '\t' + (["[RC]", "[  ]"][scaff_strand]) + '\n')
				contig_beg = i
				while contig_beg >= 0 and insertion_coords[scaff_name][contig_beg] == None:
					contig_beg -= 1
				(contig_name, cont_len, is_new) = insertion_coords[scaff_name][contig_beg]
				if contig_beg + cont_len >= scaff_beg:
					f_out.write('\t' + contig_name + '\t' + str(contig_beg) + '-' + str(contig_beg + cont_len) + '\t' + (["old", "new"][is_new]) + '\n')
				contig_beg += 1
				while contig_beg < len(contigs_insertion) and contig_beg <= scaff_end:
					if insertion_coords[scaff_name][contig_beg] != None:
						(contig_name, cont_len, is_new) = insertion_coords[scaff_name][contig_beg]
						f_out.write('\t' + contig_name + '\t' + str(contig_beg) + '-' + str(contig_beg + cont_len) + '\t' + (["old", "new"][is_new]) + '\n')
					contig_beg += 1
	f_out.close()

if __name__ == "__main__":
	if len(sys.argv) == 1:
		print "Usage:", sys.argv[0], "-b <path to bwa> -u <unused contigs> -t <target genome> -c <new contigs coords> -s <scaffolds>"
		print "Please use the --help option to get more usage information."
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Evaluate insertion unused contigs.')
	parser.add_argument("-b", "--bwa", help="path to bwa", required=True)
	parser.add_argument("-m", "--mummer", help="path to mummer", required=True)
	parser.add_argument("-u", "--unused", help="unused contigs", required=True)
	parser.add_argument("-t", "--target", help="target genome", required=True)
	parser.add_argument("-c", "--coords", help="new contigs coords", required=True)
	parser.add_argument("-s", "--scaff", help="scaffolds", required=True)

	args = parser.parse_args()
	bwa_path = args.bwa
	mummer_path = args.mummer
	unused_contigs_path = args.unused
	target_path = args.target
	contigs_coords_path = args.coords
	scaffolds_path = args.scaff

	work_dir = os.path.join(os.path.dirname(unused_contigs_path), 'evaluation')
	if not os.path.exists(work_dir):
		os.mkdir(work_dir)

	data_name_contigs = os.path.join(work_dir, 'unused_to_target')
	sam_file_contigs = build_alignment_bwa(bwa_path, data_name_contigs, target_path, unused_contigs_path)
	target_len = process_ref_file(target_path)
	alignment_contigs = process_sam_file_contigs(sam_file_contigs, target_len)
	insertion_names = process_contigs_coords_names(contigs_coords_path)
	evaluation_contigs_filename = os.path.join(work_dir, 'evaluation_contigs_result.txt')
	scaff_len = process_ref_file(scaffolds_path)
	compare_alignment_and_insertion(alignment_contigs, insertion_names, scaff_len, evaluation_contigs_filename)

	data_name_scaffolds = os.path.join(work_dir, 'scaffolds_to_target')
	coords_file_scaffolds = build_alignment_mummer(mummer_path, data_name_scaffolds, target_path, scaffolds_path)
	insertion_coords = process_contigs_coords_coords(contigs_coords_path, scaff_len)
	alignment_scaffolds = process_coords_file_scaffolds(coords_file_scaffolds)
	evaluation_scaffolds_filename = os.path.join(work_dir, 'evaluation_scaffolds_result.txt')
	compare_alignment_and_insertion_coords_scaffolds(alignment_scaffolds, insertion_coords, evaluation_scaffolds_filename)
	
	print
	print '=============================='
	print 'Result can be found in', evaluation_contigs_filename

