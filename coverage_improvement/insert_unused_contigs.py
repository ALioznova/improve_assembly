from sets import Set
import numpy
import os
import sys
import argparse
import subprocess
import numpy
from Bio import SeqIO
from Bio.Seq import Seq

def parse_recipe(ragout_recipe):
	fasta_pathes = {}
	target_name = None
	blocks = []
	base_name = os.path.dirname(ragout_recipe)
	with open(ragout_recipe) as rcp:
		for line in rcp:
			if (not line.strip()) or line.startswith('#'):
				continue
			elif line.startswith('.target'):
				target_name = line.split('=')[1].strip()
			elif line.startswith('.blocks'):
				blocks = [int(b) for b in (line.split('=')[1]).split(',')]
			elif not line.startswith('.'):
				if len(line.split('=')) == 2:
					[name_fasta, path] = [e.strip() for e in line.split('=')]
					if name_fasta.endswith('.fasta'):
						name = name_fasta.split('.')[0]
						if not os.path.isabs(path):
							path = os.path.normpath(os.path.join(base_name, path))
						fasta_pathes[name] = path
	return (target_name, blocks, fasta_pathes)

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

def process_ref_file(ref_file, ref_name):
	print
	print "Parsing reference...", ref_name
	handle = open(ref_file, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	coverage = {}
	for elem in records:
		name = ref_name + '.' + elem.name
		seq = elem.seq
		empty_list = []
		for i in xrange(len(seq)):
			empty_list.append([])
		coverage[name] = empty_list
		print name, "\tlen=", len(empty_list)
	print "Found", len(coverage), "ref seq in", ref_name
	return coverage

def sequence_unmapped(record):
	# 0x4 segment unmapped, 0x8 next segment in the template unmapped
	return (record ['FLAG'] & int(0x4)) or (record ['FLAG'] & int(0x8))

def supplementary_alignment(record):
	# 0x800 supplementary alignment
	return ((record ['FLAG']) & int(0x800))

def seq_reverse_complemented(record):
	# 0x10 reverse complemented
	return ((record ['FLAG']) & int(0x10))

def process_sam_file(ref_file, sam_file, reference_name):
	alignment = process_ref_file(ref_file, reference_name)
	number = 0
	mapped_number = 0
	print
	print "Reading alignment to", reference_name, "..."
	for line in open(sam_file):
		record = parse_sam_record(line)
		if record:
			number += 1
			if sequence_unmapped(record) or supplementary_alignment(record):
				continue
			mapped_number += 1
			begin = record['POS']
			end = record['POS'] + get_alignment_length(record['CIGAR'].strip())
			ref_name = reference_name + '.' + record['RNAME']
			seq_name = record['QNAME']
			strand = not seq_reverse_complemented(record)
			alignment[ref_name][begin].append((end, seq_name, strand))
	print "Total number of sequences is: ", number
	print "Number of mapped sequences is: ", mapped_number
	return alignment

def process_contigs_file(file_name):
	all_contigs_length = {}
	handle = open(file_name, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	for elem in records:
		all_contigs_length[elem.id] = len(elem.seq)
	return all_contigs_length

def process_links_file(file_name, all_contigs_length):
	used_contigs_set = Set()
	scaffolds_as_contig_pairs = {}
	with open(file_name) as fp:
		cur_scaff_name = None
		for line in fp:
			if line.startswith("ragout-scaffold-"):
				cur_scaff_name = line.strip()
				scaffolds_as_contig_pairs[cur_scaff_name] = []
			elif (not (line.startswith("--"))) and (line.startswith('+') or line.startswith('-')):
				[contig_1, contig_2, gap] = line.split()[:3]
				gap = int(gap)
				if contig_1.startswith('+'):
					cont1_strand = True
				elif contig_1.startswith('-'):
					cont1_strand = False
				cont1_name = contig_1[1:]
				if contig_2.startswith('+'):
					cont2_strand = True
				elif contig_2.startswith('-'):
					cont2_strand = False
				cont2_name = contig_2[1:]
				used_contigs_set.add(cont1_name)
				used_contigs_set.add(cont2_name)
				scaffolds_as_contig_pairs[cur_scaff_name].append((cont1_name, cont1_strand, cont2_name, cont2_strand, gap))
	scaffolds_as_contigs = {}
	for (scaff_name, cont_pairs) in scaffolds_as_contig_pairs.iteritems():
		scaffold_len = 0
		for (cont1_name, cont1_strand, cont2_name, cont2_strand, gap) in cont_pairs:
			scaffold_len += all_contigs_length[cont1_name]
			scaffold_len += gap
		scaffold_len += all_contigs_length[cont2_name]
		scaffolds_as_contigs[scaff_name] = [None] * scaffold_len
		current_pos_in_scaff = 0
		for (cont1_name, cont1_strand, cont2_name, cont2_strand, gap) in cont_pairs:
			scaffolds_as_contigs[scaff_name][current_pos_in_scaff] = (cont1_name, cont1_strand, all_contigs_length[cont1_name])
			current_pos_in_scaff += all_contigs_length[cont1_name]
			current_pos_in_scaff += gap
		scaffolds_as_contigs[scaff_name][current_pos_in_scaff] = (cont2_name, cont2_strand, all_contigs_length[cont2_name])
	return (used_contigs_set, scaffolds_as_contigs)

def write_unused_contigs(workdir_name, all_contigs_path, unused_contigs_set):
	work_dir = os.path.join(workdir_name, "unused_contigs_workdir")
	if not os.path.exists(work_dir):
		os.mkdir(work_dir)
	unused_contigs_path = os.path.join(work_dir, "unused_contigs.fasta")
	handle = open(all_contigs_path, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	handle.close()
	sequences = []
	for elem in records:
		if elem.id in unused_contigs_set:
			sequences.append(elem)
	output_handle = open(unused_contigs_path, "w")
	SeqIO.write(sequences, output_handle, "fasta")
	output_handle.close()
	return unused_contigs_path

def process_blocks_coords_file(file_name):
	seq_id_to_name = {}
	name_to_seq_id = {}
	seq_as_blocks = {}
	with open(file_name) as fp:
		fp.readline()
		while True:
			line = fp.readline()
			if line.startswith('--'):
				break
			[seq_id, size, description] = line.split()[:3]
			seq_id = int(seq_id)
			size = int(size)
			seq_id_to_name[seq_id] = (description, size)
			name_to_seq_id[description] = seq_id
		for (seq_id, (description, size)) in seq_id_to_name.iteritems():
			seq_as_blocks[seq_id] = [None] * size
		while True:
			line = fp.readline().strip()
			if line == '':
				break
			if line.startswith('--'):
				continue
			elif line.startswith('Block #'):
				cur_block = int(line.split('#')[1])
				fp.readline()
			else:
				[seq_id, strand, start, end, length] = line.split()[:5]
				seq_id = int(seq_id)
				start = int(start)
				end = int(end)
				length = int(length)
				if strand == '+':
					strand = True
				elif strand == '-':
					strand = False
					(start, end) = (end, start)
				assert seq_as_blocks[seq_id][start] == None
				seq_as_blocks[seq_id][start] = (cur_block, strand, length)
	return (seq_id_to_name, name_to_seq_id, seq_as_blocks)

def reverse_seq_as_blocks(contig_as_blocks):
	contig_len = len(contig_as_blocks)
	reverse_contig_as_blocks = [None] * contig_len
	for i in xrange(len(contig_as_blocks)):
		if contig_as_blocks[i] != None:
			(cur_block, strand, block_len) = contig_as_blocks[i]
			reverse_contig_as_blocks[contig_len - i - block_len] = (cur_block, not strand, block_len)
	return reverse_contig_as_blocks

def get_scaffolds_as_blocks(scaffolds_as_contigs, name_to_seq_id, seq_as_blocks):
	scaffolds_as_blocks = {}
	for (name, seq) in scaffolds_as_contigs.iteritems():
		scaffolds_as_blocks[name] = [None] * len(seq)
		for i in xrange(len(seq)):
			if seq[i] != None:
				(contig_name, cont_strand, cont_len) = seq[i]
				contig_seq_id = name_to_seq_id[target_name + '.' + contig_name]
				contig_as_blocks = seq_as_blocks[contig_seq_id]
				if cont_strand == False:
					contig_as_blocks = reverse_seq_as_blocks(contig_as_blocks)
				for j in xrange(len(contig_as_blocks)):
					scaffolds_as_blocks[name][i+j] = contig_as_blocks[j]
	return scaffolds_as_blocks

def get_unused_contigs_aligmnent(work_dir, unused_contigs_path, fasta_pathes, target_name, bwa_path):
	sam_files = {}
	unused_contigs_alignment = {}
	for ref_name in fasta_pathes.iterkeys():
		if ref_name == target_name:
			continue
		ref_path = fasta_pathes[ref_name]
		data_name = os.path.join(work_dir, ref_name)
		sam_file = build_alignment_bwa(bwa_path, data_name, ref_path, unused_contigs_path)
		sam_files[ref_name] = sam_file
	for ref_name in fasta_pathes.iterkeys():
		if ref_name == target_name:
			continue
		ref_path = fasta_pathes[ref_name]
		sam_file = sam_files[ref_name]
		alignment = process_sam_file(ref_path, sam_file, ref_name)
		for (name, seq) in alignment.iteritems():
			unused_contigs_alignment[name] = seq
	return unused_contigs_alignment

def get_neighbour_blocks_to_unused_contigs(unused_contigs_alignment, name_to_seq_id, seq_as_blocks):
	neighbour_blocks = {}
	for (name, seq) in unused_contigs_alignment.iteritems():
		for i in xrange(len(seq)):
			if seq[i] == []:
				continue
			for (cont_end, cont_name, cont_strand) in seq[i]:
				(left_block, right_block) = (None, None)
				(left_dist, right_dist) = (-1, -1)
				seq_id = name_to_seq_id[name]
				pos_in_ref = i
				while (pos_in_ref >= 0) and (seq_as_blocks[seq_id][pos_in_ref] == None):
					pos_in_ref -= 1
				if pos_in_ref >= 0:
					(cur_block, block_strand, block_length) = seq_as_blocks[seq_id][pos_in_ref]
					(left_block, left_strand) = (cur_block, block_strand)
					left_dist = i - (pos_in_ref + block_length)
				else:
					continue
				contig_aligned_to_block = False
				for j in xrange(i, cont_end):
					if seq_as_blocks[seq_id][j] != None:
						contig_aligned_to_block = True
						break
				if pos_in_ref + block_length >= i:
					contig_aligned_to_block = True
				if contig_aligned_to_block:
					continue
				pos_in_ref = cont_end
				while (pos_in_ref < len(seq_as_blocks[seq_id])) and (seq_as_blocks[seq_id][pos_in_ref] == None):
					pos_in_ref += 1
				if pos_in_ref < len(seq_as_blocks[seq_id]):
					(cur_block, block_strand, block_length) = seq_as_blocks[seq_id][pos_in_ref]
					(right_block, right_strand) = (cur_block, block_strand)
					right_dist = pos_in_ref - cont_end
				else:
					continue
				contig_aligned_to_block = False
				if cont_end >= pos_in_ref:
					contig_aligned_to_block = True
				if contig_aligned_to_block:
					continue
				if not neighbour_blocks.has_key((left_block, left_strand, right_block, right_strand)):
					neighbour_blocks[(left_block, left_strand, right_block, right_strand)] = []
				neighbour_blocks[(left_block, left_strand, right_block, right_strand)].append((cont_name, cont_strand, left_dist, right_dist))
	return neighbour_blocks

def find_blocks_to_insert(scaffolds_as_blocks, neighbour_blocks):
	blocks_to_insert = {}
	for (name, seq) in scaffolds_as_blocks.iteritems():
		prev_block = None
		prev_strand = None
		prev_block_coord = -1
		cur_pos_in_seq = -1
		for elem in seq:
			cur_pos_in_seq += 1
			if elem != None:
				(cur_block, cur_strand, length) = elem
				if prev_block == None:
					prev_block = cur_block
					prev_strand = cur_strand
					prev_block_coord = cur_pos_in_seq + length
					continue
				(left_coord, right_coord) = (prev_block_coord, cur_pos_in_seq)
				contig_list = []
				if neighbour_blocks.has_key((prev_block, prev_strand, cur_block, cur_strand)):
					contig_list = neighbour_blocks[(prev_block, prev_strand, cur_block, cur_strand)]
				elif neighbour_blocks.has_key((cur_block, not cur_strand, prev_block, not prev_strand)):
					contig_list = neighbour_blocks[(cur_block, not cur_strand, prev_block, not prev_strand)]
					for i in xrange(len(contig_list)):
						(contig_name, contig_strand, left_dist, right_dist) = contig_list[i]
						contig_list[i] = (contig_name, not contig_strand, right_dist, left_dist)
				if len(contig_list) != 0:
					blocks_to_insert[(prev_block, prev_strand, cur_block, cur_strand)] = (contig_list, (left_coord, right_coord, name))
				prev_block = cur_block
				prev_strand = cur_strand
				prev_block_coord = cur_pos_in_seq + length
	return blocks_to_insert

def output_contigs_between_blocks(contigs_between_blocks, contigs_between_blocks_file_name, all_contigs_length):
	f_contigs_between_blocks = open(contigs_between_blocks_file_name, 'w')
	for (block_size, res) in contigs_between_blocks.iteritems():
		f_contigs_between_blocks.write('--------------------\n')
		f_contigs_between_blocks.write('Block size = ' + str(block_size) + '\n')
		if len(res) == 0:
			f_contigs_between_blocks.write('Nothing\n')
		else:
			for (block_pair, (contig_list, (left_coord, right_coord, scaff_name))) in res.iteritems():
				(left, left_strand, right, right_strand) = block_pair
				f_contigs_between_blocks.write('\nBetween blocks ' + (["-", "+"][left_strand]) + str(left) + ' and block ' + (["-", "+"][right_strand]) + str(right) + ':\n')
				f_contigs_between_blocks.write('Block ' + (["-", "+"][left_strand]) + str(left) + ' (block_end: ' + str(left_coord) + ')' + ' in ' + scaff_name + '\n')
				f_contigs_between_blocks.write('Block ' + (["-", "+"][right_strand]) + str(right) + ' (block_beg: ' + str(right_coord) + ')' + ' in ' + scaff_name + '\n')
				for (contig_name, contig_strand, left_dist, right_dist) in contig_list:
					f_contigs_between_blocks.write('\t' + (["-", "+"][contig_strand]) + contig_name + '\tlen:' + str(all_contigs_length[contig_name]) + '\t\tleft_dist: ' + str(left_dist) + '\tright_dist: ' + str(right_dist) + '\n')
	f_contigs_between_blocks.close()

def output_contigs_coords(contigs_between_blocks, contigs_coords_file_name, all_contigs_length):
	f_contigs_coords = open(contigs_coords_file_name, 'w')
	for (block_size, res) in contigs_between_blocks.iteritems():
		if len(res) != 0:
			for (block_pair, (contig_list, (left_coord, right_coord, scaff_name))) in res.iteritems():
				(left, left_strand, right, right_strand) = block_pair
				for (contig_name, contig_strand, left_dist, right_dist) in contig_list:
					coord_to_insert_left = left_coord + left_dist
					coord_to_insert_right = right_coord - right_dist - all_contigs_length[contig_name]
					if coord_to_insert_left > coord_to_insert_right:
						(coord_to_insert_left, coord_to_insert_right) = (coord_to_insert_right, coord_to_insert_left)
					f_contigs_coords.write(scaff_name + '\t' + (["-", "+"][contig_strand]) + contig_name + '\t\t' + str(coord_to_insert_left) + '-' + str(coord_to_insert_right) + '\n')
	f_contigs_coords.close()

def output_scaffolds_as_blocks(scaffolds_as_blocks, scaffolds_as_blocks_filename, block_size):
	f_scaffolds_as_blocks = open(scaffolds_as_blocks_filename, 'a')
	f_scaffolds_as_blocks.write('==============================\n')
	f_scaffolds_as_blocks.write('Block size: ' + block_size + '\n')
	f_scaffolds_as_blocks.write('==============================\n')
	for (name, seq) in scaffolds_as_blocks.iteritems():
		f_scaffolds_as_blocks.write('>' + name + '\n')
		f_scaffolds_as_blocks.write('beg\tname\tend\n')
		for i in xrange(len(seq)):
			if seq[i] != None:
				(cur_block, cur_strand, length) = seq[i]
				f_scaffolds_as_blocks.write(str(i) + '\t' + (["-", "+"][cur_strand]) + str(cur_block) + '\t' + str(i + length) + '\n')
	f_scaffolds_as_blocks.write('\n')
	f_scaffolds_as_blocks.close()


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print "Usage:", sys.argv[0], "-b <path to bwa> -w <path to ragout workdir> -r <ragout recipe> "
		print "Please use the --help option to get more usage information."
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Insert unused contigs.')
	parser.add_argument("-b", "--bwa", help="path to bwa", required=True)
	parser.add_argument("-w", "--workdir", help="ragout workdir", required=True)
	parser.add_argument("-r", "--recipe", help="ragout recipe", required=True)

	args = parser.parse_args()
	bwa_path = args.bwa
	ragout_workdir_path = args.workdir
	ragout_recipe = args.recipe

	if os.path.exists(os.path.join(ragout_workdir_path, "scaffolds_refined.links")):
		scaffolds_links_path = os.path.join(ragout_workdir_path, "scaffolds_refined.links")
	else:
		scaffolds_links_path = os.path.join(ragout_workdir_path, "scaffolds.links")

	(target_name, blocks, fasta_pathes) = parse_recipe(ragout_recipe)
	all_contigs_length = process_contigs_file(fasta_pathes[target_name])
	(used_contigs_set, scaffolds_as_contigs) = process_links_file(scaffolds_links_path, all_contigs_length)
	unused_contigs_set = (Set(all_contigs_length.keys())).difference(used_contigs_set)
	unused_contigs_path = write_unused_contigs(ragout_workdir_path, fasta_pathes[target_name], unused_contigs_set)
	work_dir = os.path.dirname(unused_contigs_path)
	unused_contigs_alignment = get_unused_contigs_aligmnent(work_dir, unused_contigs_path, fasta_pathes, target_name, bwa_path)

	print
	print "used contigs number\t", len(used_contigs_set)
	print "unused contigs number\t", len(unused_contigs_set)

	scaffolds_as_blocks_filename = os.path.join(work_dir, 'scaffolds_as_blocks.txt')
	if os.path.exists(scaffolds_as_blocks_filename):
		os.remove(scaffolds_as_blocks_filename)

	contigs_between_blocks = {}
	for b in blocks:
		block_size = str(b)
		if os.path.exists(os.path.join(ragout_workdir_path, "sibelia-workdir")):
			blocks_coords_path = os.path.join(ragout_workdir_path, "sibelia-workdir", block_size, "blocks_coords.txt")
		elif os.path.exists(os.path.join(ragout_workdir_path, "maf-workdir")):
			blocks_coords_path = os.path.join(ragout_workdir_path, "maf-workdir", block_size, "blocks_coords.txt")
		(seq_id_to_name, name_to_seq_id, seq_as_blocks) = process_blocks_coords_file(blocks_coords_path)
		scaffolds_as_blocks = get_scaffolds_as_blocks(scaffolds_as_contigs, name_to_seq_id, seq_as_blocks)
		neighbour_blocks = get_neighbour_blocks_to_unused_contigs(unused_contigs_alignment, name_to_seq_id, seq_as_blocks)

		blocks_to_insert = find_blocks_to_insert(scaffolds_as_blocks, neighbour_blocks)
		contigs_between_blocks[b] = blocks_to_insert
		
		output_scaffolds_as_blocks(scaffolds_as_blocks, scaffolds_as_blocks_filename, block_size)
		
	contigs_between_blocks_file_name = os.path.join(work_dir, 'contigs_between_blocks.txt')
	contigs_coords_file_name = os.path.join(work_dir, 'contigs_coords.txt')
	output_contigs_between_blocks(contigs_between_blocks, contigs_between_blocks_file_name, all_contigs_length)
	output_contigs_coords(contigs_between_blocks, contigs_coords_file_name, all_contigs_length)
	print
	print '=============================='
	print 'Result can be found in', contigs_between_blocks_file_name

