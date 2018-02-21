import argparse, shlex, subprocess, sys


bwa_exec = 'bwa.kit/bwa'
virus_ind = 'data/viral-refseq.fna'
fungi_ind = 'data/fungi-refseq.fna'


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Run BWA mem for viruses or fungi.')
    parser.add_argument('reads', help='Reads to run BWA on. Required.')
    parser.add_argument('--fungi', action='store_true', help='Run BWA on fungi. This or --virus required.')
    parser.add_argument('--virus', action='store_true', help='Run BWA on viruses. This or --fungi required.')
    parser.add_argument('--output', default='alignments.sam', help='Where to output BWA results to.')
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    if (args.virus and args.fungi) or not(args.virus or args.fungi):
        print 'Error: must specify either --virus or --fungi.'
        sys.exit()

    if args.virus:
        cmd = ' '.join([bwa_exec, 'mem', '-a', virus_ind, args.reads])
    else:  # args.fungi
        cmd = ' '.join([bwa_exec, 'mem', '-a', fungi_ind, args.reads])
    with(open(args.output, 'w')) as outfile:
        subprocess.call(shlex.split(cmd), stdout=outfile)


if __name__ == '__main__':
	main()
#
