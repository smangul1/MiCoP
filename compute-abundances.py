import argparse
import operator
import random
import sys
import time


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Compute abundance estimations for species in a sample.')
    parser.add_argument('refdb', help='Reference database file. Required.')
    parser.add_argument('sam', help='.sam file or file with list of SAM files to process. Required.')
    parser.add_argument('--abundance_cutoff', type=float, default=0.0, help='Organism abundance to count it as present.')
    parser.add_argument('--min_map', type=int, default=100, help='Minimum bases mapped to count a hit.')
    parser.add_argument('--max_ed', type=int, default=1, help='Maximum edit distance from a reference to count a hit.')
    parser.add_argument('--normalize', type=bool, choices=[True,False], default=True,
    					help='Normalize species abundance by genome length or not. Default: True')
    parser.add_argument('--output', default='abundances.txt', help='Output abundances file. Default: abundances.txt')
    #parser.add_argument('--paired', action='store_true', default=False, help='Use if reads are paired end.')
    parser.add_argument('--pct_id', type=float, default=0.99, help='Minimum percent identity from reference to count a hit.')
    parser.add_argument('--read_cutoff', type=int, default=100, help='Number of reads to count an organism as present.')
    parser.add_argument('--taxon', choices=['genus','species','strain'], default='species',
    					help='genus/species/strain level of taxonomic classification')
    args = parser.parse_args()
    return args


def find_taxid(tag):
        if not '|' in tag:
                return tag
        else:
                splits = tag.split('|')
                for sp in splits:
                        if '_' in sp:
                                return sp
        return tag


def ids2len(refdb):
	genlens, curtag = {}, ''
	with(open(refdb, 'r')) as ref:
		for line in ref:
			if line.startswith('>'):
				curtag = find_taxid(line.strip().split(' ')[0][1:])
				if curtag in genlens:
					print 'Warning: TaxID ' + curtag + ' occurrs twice in reference'
				else:
					genlens[curtag] = 0
			else:
				genlens[curtag] += len(line.strip())
	return genlens


def spe2ids(args):
	level = ['','genus','species','strain'].index(args.taxon)
	if level == 3:
		level = 999  # if strain, keep entire organism name
	spe2ids = {}
	with(open(args.refdb, 'r')) as ref:
		for line in ref:
			if not line.startswith('>'):
				continue
			splits = line.strip().split(' ')
			taxid = find_taxid(splits[0][1:])
			spe = ' '.join(splits[1:1+level])
			if spe.endswith(','):
				spe = spe[:-1]
			if spe in spe2ids:
				spe2ids[spe].append(taxid)
			else:
				spe2ids[spe] = [taxid]
	return spe2ids


def filter_line(args, splits):  # determine whether to filter this line out
    cigar = splits[5]
    if cigar == '*':
    	return True
    '''if 'M' not in cigar:
            return True
    loc = cigar.find('M')
    try:
            mapped = int(cigar[loc-2:loc])
    except:
            return True
    if mapped < args.min_map:
            return True'''
    matched_len, total_len, cur = 0, 0, 0
    for ch in cigar:
    	if not ch.isalpha():  # ch is a number
    		cur = (cur * 10) + int(ch)
    	else:  # ch is a letter
    		if ch == 'M' or ch == '=':
    			matched_len += cur
    		total_len += cur
    		cur = 0
    if matched_len < args.min_map:
    	return True
    edit_distance = int(splits[11][5:])
    if edit_distance > args.max_ed:
    	return True
    #if float(mapped) / float(mapped + edit_distance) < args.pct_id:
    if float(matched_len) / float(total_len) < args.pct_id:
    	return True
    return False  # if read passes quality checks, don't filter


def compute_abundances(args, samfile, genlens, spe2id):
	infile = open(samfile, 'r')
	ids, ids2abs, spe2abs = [], {}, {}
	prev_read_num, prev_tag, prev_count, ignore = '', '', 0.0, False
	multimapped, ids2reads, read_ordering = {}, {}, []
	lc = 0

	print 'Reading sam file ' + samfile
	for line in infile:
		lc += 1
		if lc % 1000000 == 0:
			print 'Done reading ' + str(lc) + ' lines of sam file'
		if line.startswith('@'):
			continue
		splits = line.strip().split('\t')
		if filter_line(args, splits):
			continue

		tag = find_taxid(splits[2])
		if tag == '*':
			continue
		#read_num = int(splits[0])
		read_num = splits[0]
		read_ordering.append(read_num)
		if read_num == prev_read_num and tag == prev_tag:
			pass
			#if prev_count < 2.0 and args.paired == True:
			#	prev_count += 1.0
		elif read_num == prev_read_num and tag != prev_tag:
			ignore = True
			strnum = str(prev_read_num)
			if strnum not in multimapped:
				multimapped[strnum] = [prev_tag]
			else:
				multimapped[strnum].append(prev_tag)
			prev_tag = tag
		else:
			if not(prev_read_num == '' or ignore == True):
				if prev_tag not in ids2reads:
					ids2reads[prev_tag] = [prev_read_num]
				else:
					ids2reads[prev_tag].append(prev_read_num)
				if prev_tag not in ids:
					ids.append(prev_tag)
				if prev_tag in ids2abs:
					ids2abs[prev_tag] += prev_count
				else:
					ids2abs[prev_tag] = prev_count
			elif ignore == True:
				multimapped[prev_read_num].append(prev_tag)
	                prev_read_num = read_num
	                prev_tag = tag
	                prev_count = 1.0
	                ignore = False
	infile.close()
	print 'Done reading sam file.'

	if not(prev_read_num == '' or ignore == True):
		if prev_tag not in ids2reads:
	                ids2reads[prev_tag] = [prev_read_num]
	        else:
	                ids2reads[prev_tag].append(prev_read_num)
	        if prev_tag not in ids:
	                ids.append(prev_tag)
	       	if prev_tag in ids2abs:
	                ids2abs[prev_tag] += prev_count
	       	else:
	                ids2abs[prev_tag] = prev_count
	elif ignore == True:
		multimapped[read_num].append(prev_tag)

	print 'Deleting species/reads with insufficient evidence...'
	for spe in spe2id.keys():
	        if not (spe in spe2abs):
	                spe2abs[spe] = 0.0
	        for taxid in spe2id[spe]:
	                if taxid in ids2abs:
	                        spe2abs[spe] += ids2abs[taxid]
	for spe in spe2abs:
		if spe2abs[spe] < args.read_cutoff:
			spe2abs[spe] = 0.0
			for taxid in spe2id[spe]:
				if taxid in ids2abs:
					ids2abs[taxid] = 0
	for taxid in ids:
		if taxid in ids2abs and ids2abs[taxid] == 0:
			del ids2abs[taxid]
	print 'Done deleting species/reads.'

	print 'Assigning multimapped reads...'
	added = {}
	for read in multimapped.keys():
		randnum = random.random()
		options, total = {}, 0.0
		for taxid in multimapped[read]:
			if taxid not in ids2abs or taxid in options:
				continue
			ab = ids2abs[taxid]
			total += ab
			options[taxid] = ab
		if len(options.keys()) == 0:
			continue
		for key in options.keys():
			ab = options[key] / total
			if ab >= randnum:
				if key in added:
					added[key] += 1.0 #2.0
				else:
					added[key] = 1.0 #2.0
				break
			else:
				randnum -= ab
			
	for key in added.keys():
		val = added[key]
		ids2abs[key] += val
	print 'Multimapped reads assigned.'

	if args.normalize:
		for taxid in ids2abs.keys():
			ids2abs[taxid] /= genlens[taxid]  # normalize by genome length

	total_ab = 0.0
	for taxid in ids2abs.keys():
		total_ab += float(ids2abs[taxid])

	for taxid in ids2abs.keys():
		ids2abs[taxid] = float(ids2abs[taxid]) * 100.0 / total_ab  # normalize abundances

	spe2abs = {}
	for spe in spe2id.keys():
		if not (spe in spe2abs):
			spe2abs[spe] = 0.0
		for taxid in spe2id[spe]:
			if taxid in ids2abs:
				spe2abs[spe] += ids2abs[taxid]

	return spe2abs, ids2abs


def main():
	args = parseargs()
	if args.pct_id > 1.0 or args.pct_id < 0.0:
		print 'Error: --pct_id must be between 0.0 and 1.0, inclusive'
		sys.exit()
	samfiles = []
	if args.sam.endswith('.sam'):
		samfiles.append(args.sam)
	else:
		with(open(args.sam, 'r')) as filenames:
			for line in filenames:
				samfiles.append(line.strip())
	genlens = ids2len(args.refdb)  # maps NCBI taxID to genome length for normalization
	spe2id = spe2ids(args)  # maps NCBI taxIDs to species

	spe2abs, ids2abs = {}, {}
	for sam in samfiles:
		s2a, i2a = compute_abundances(args, sam, genlens, spe2id)
		for spe in s2a:
			if spe not in spe2abs:
				spe2abs[spe] = s2a[spe]
			else:
				spe2abs[spe] += s2a[spe]
		for taxid in i2a:
			if taxid not in ids2abs:
				ids2abs[taxid] = i2a[taxid]
			else:
				ids2abs[taxid] += i2a[taxid]

	for spe in spe2abs:
		spe2abs[spe] /= len(samfiles)
	for taxid in ids2abs:
		ids2abs[taxid] /= len(samfiles)
	sorted_ids2abs = sorted(ids2abs.items(), key=operator.itemgetter(1), reverse=True)
	sorted_spe2abs = sorted(spe2abs.items(), key=operator.itemgetter(1), reverse=True)

	outfile = open(args.output, 'w')
	print 'Writing genome and species abundances...'
	outfile.write('Abundances by Species:\n')
	for line in sorted_spe2abs:
		desc = line[0]
		if line[1] > args.abundance_cutoff:
			outfile.write(str(desc) + '\t' + str(line[1]) + '\n')

	outfile.write('\n\nAbundances by NCBI Taxonomic ID:\n')
	for line in sorted_ids2abs:
		desc = line[0]
		if '|' in desc:
			desc = desc.split('|')[0]
		if line[1] > args.abundance_cutoff:
			outfile.write(str(desc) + '\t' + str(line[1]) + '\n')
	outfile.close()
	print 'Done.'


if __name__ == '__main__':
	main()
#
