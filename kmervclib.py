import os, shutil, logging, sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np
from scipy import special
from scipy import stats

CWD = os.getcwd()
MERGED_TABLE_COLUMN_COUNT = 9	# extra columns for distinct extract starts/ends per variant merged variants
MINIMUM_REQUIRED_FLANKING_REGION = 2	# minumum number of bases required on either side of indel extraction region
SCRIPT_PATH = os.path.realpath(sys.argv[0]).rsplit('/', 1)[0]	# path to folder holding real path of script
REFERENCES_FOLDER = os.path.join(SCRIPT_PATH, 'references')	# path to folder holding auxiliary reference files
GENOME_FASTA = os.path.join(REFERENCES_FOLDER, 'grch38_canonical_chrs_chrM.fa')	# path to genome reference fasta file
GENOME_JELLYFISH = None	# Initiliaze in setup_output_directory
VARIANT_COL_NAME = 'Variant'
VARIANT_TYPE_COL_NAME = 'Variant_Type'
WILDTYPE_SEQUENCE_COL_NAME = 'Wildtype_Sequence'
WILDTYPE_CONTROL_COUNT_COL_NAME = 'Wildtype_Control_Count'
WILDTYPE_TEST_COUNT_COL_NAME = 'Wildtype_Test_Count'
WILDTYPE_KMER_COUNT_COL_NAME = 'Wildtype_Kmer_Count'
MUTATION_SEQUENCE_COL_NAME = 'Mutation_Sequence'
MUTATION_CONTROL_COUNT_COL_NAME = 'Mutation_Control_Count'
MUTATION_TEST_COUNT_COL_NAME = 'Mutation_Test_Count'
MUTATION_KMER_COUNT_COL_NAME = 'Mutation_Kmer_Count'
REJECT_NULL_CONTROL_COL_NAME = 'Reject_Null_Control'
REJECT_NULL_TEST_COL_NAME = 'Reject_Null_Test'
CONTROL_P_VALUE_COL_NAME = 'Control_P_Val'
TEST_P_VALUE_COL_NAME = 'Test_P_Val'
VARIANT_CALL_COL_NAME = 'Variant_Call'
WILDTYPE_UNIQUE_COUNT_COL_NAME = 'Wildtype_Unique_Count'
MUTATION_ZERO_COUNT_COL_NAME = 'Mutation_Unique_Count'
SUFFICIENT_DATA_MINIMUM_KMERS = 5	# minimum count of unique kmers required for to analyze variant call
SEQUENCE_ERROR_PROBABILITY = 0.01
ALPHA = 0.01



def setup_logger():
	# Setup basic logging configuration
	logging.basicConfig(
			filename='application.log',
			level=logging.DEBUG,
			format='%(asctime)s [%(levelname)s] %(message)s',
			datefmt='%Y-%m-%d %H:%M:%S')


def get_arg_parser():
	# Get arguments passed in by user and return args when using compare function
	parser = ArgumentParser(description='kmerVC Command Line Argument Parser')
	parser.add_argument('subcommand', choices=['compare', 'assess'])

	is_compare = sys.argv[1] == 'compare'

	vc_group = parser.add_mutually_exclusive_group(required=is_compare)
	vc_group.add_argument('-v', '--vcf', dest='vcf_input', help='Input vcf file')
	vc_group.add_argument('-b', '--bed', dest='bed_input', help='Input ved file')

	fastq_group = parser.add_argument_group('fastq_group', 'Fastq input files')
	fastq_group.add_argument('-t1', '--test1', dest='test_fastq1', help='Fastq file 1 from test sample')
	fastq_group.add_argument('-t2', '--test2', dest='test_fastq2',  help='Fastq file 2 from test sample')
	fastq_group.add_argument('-c1', '--control1', dest='control_fastq1', help='Fastq file 1 from control sample')
	fastq_group.add_argument('-c2', '--control2', dest='control_fastq2', help='Fastq file 2 from control sample')

	jellyfish_group = parser.add_argument_group('jellyfish_group', 'Jellyfish input files')
	jellyfish_group.add_argument('-j1', '--jellyfish_test', dest='jellyfish_test', help='Jellyfish file of test input')
	jellyfish_group.add_argument('-j2', '--jellyfish_control', dest='jellyfish_control', help='Jellyfish file of control input')

	parser.add_argument('-k', '--kmer_size', dest='kmer_size', required=is_compare, type=int, help='Size of kmer to use for analysis')
	parser.add_argument('-o', '--output_name', dest='output_name', required=is_compare, help='Output file directory name')

	parser.add_argument('-fi', '--reference_genome_fasta', dest='reference_genome_fasta', help='Reference genome fasta file to use, if different than default')
	parser.add_argument('-d', '--delimiter', dest='delimiter', help='Delimiter to use in final analysis table. Pick from TAB, COMMA, SPACE')
	parser.add_argument('-m', '--multiple_mutations', action='store_true', help='Flag indicating whether consecutive mutations in overlapping regions should be considered in combination.')
	parser.add_argument('-r', '--rna', action='store_true', help='Flag indicating if doing RNA analysis')
	parser.add_argument('-poi', '--poisson', action='store_true', help='Flag indicating if using doing poisson distribution for variant analysis')
	parser.add_argument('-a', '--alpha', dest='alpha', default=0.01, type=float, required=not is_compare, help='Alpha value used in hypothesis testing')

	parser.add_argument('-p', '-penultimate', dest='penultimate', help='Penultimate final table to make assessments on.', required=not is_compare)
	return parser


def print_help_message():
	# Message to the user to indicate how to use kmerVC
	print('For help, call: python kmerVC.py -h')

def too_few_arguments():
	# Message to the user when too few arguments are provided
	print('kmerVC:\ttools for variant calling utiling k-mers.')
	print('usage:\tpython kmerVC.py <subcommand> [options]\n')
	print_help_message()
	logging.info('Too few arguments.')
	exit_program()

def exit_program():
	logging.info('Exiting program...')
	raise SystemExit()


#############################################################
#				General Helper Methods						#
#############################################################

def reverse_complement(kmer):
	#Takes an input kmer and returns reverse complement.
	complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	return ''.join([complement[base] for base in reversed(kmer)])

def dict_from_kmer_counts(filename):
	# Takes as input a jellyfish count file with lines of the format:
	# kmer count
	# and creates and returns a python dictionary
	dataframe = pd.read_csv(filename, sep=' ', header=None, names=['Kmer', 'Count'])
	return dataframe.set_index('Kmer')['Count'].to_dict()

#############################################################
#	 Helper Methods Used Specifically for Compare Function	#
#############################################################

def create_jellyfish(args):
	def create_command(fastq_1, fastq_2, is_control_fastq):
		if fastq_1 is None: return None, None
		output_filename = '{}_control_{}mer.jf'.format(args.output_name, args.kmer_size) if is_control_fastq \
				else '{}_test_{}mer.jf'.format(args.output_name, args.kmer_size)
		command = 'jellyfish count -m {} -s 100M -t 24 -C -o {}'.format(
				args.kmer_size,
				output_filename)
		if fastq_1.endswith('.gz'):
			command = '{} <(zcat {})'.format(command, fastq_1)
			if fastq_2 is not None: command = '{} <(zcat {})'.format(command, fastq_2)
		else:
			if fastq_2 is None: command = '{} {}'.format(command, fastq_1)
			else: command = '{} -F 2 {} {}'.format(command, fastq_1, fastq_2)
		return command, output_filename

	logging.info('No jellyfish count files provided. Now creating jellyfish files from fastq input.')
	if args.test_fastq1 is None:
		logging.error('Fastq files not provided. Cannot proceed in creating jellyfish files.')
		print('error:\tfastq files not provided. At least test_fastq1 is required.')
		exit_program()
	test_command, test_output_filename = create_command(args.test_fastq1, args.test_fastq2, False)
	control_command, control_output_filename = create_command(args.control_fastq1, args.control_fastq2, True)
	os.system(test_command)
	if control_command is not None: os.system(control_command)
	return test_output_filename, control_output_filename


def setup_output_directory(output_folder_name, reference_genome_fasta_file, kmer_size, alpha):
	global GENOME_JELLYFISH, GENOME_FASTA, ALPHA
	logging.info('Setting up output directory.')
	if os.path.isdir(output_folder_name): shutil.rmtree(output_folder_name)	# remove existing output folder if it exists
	output_path = os.path.join(CWD, output_folder_name)
	os.mkdir(output_path)
	ALPHA = alpha
	if reference_genome_fasta_file:
		GENOME_FASTA = reference_genome_fasta_file	# update genome fasta given user provided input
		filename, file_extension = os.path.splitext(GENOME_FASTA)
		basename = os.path.basename(filename)
		expected_genome_jellyfish = '{}_{}mer.jf'.format(basename, kmer_size)
		if not os.path.isfile(expected_genome_jellyfish):
			logging.info('Genome jellyfish does not exist.')
			print('error:\tGenome jellyfish does not exist. Please create and place in the same directory as your input reference genome file. Call this command in the same directory as your reference genome file to create:\n\tjellyfish count -m {} -s 100M -t 24 -C -o {} {}'.format(kmer_size, '{}_{}mer.jf'.format(basename, kmer_size), '{}.fa'.format(basename)))
			exit_program()
		else:
			GENOME_JELLYFISH = expected_genome_jellyfish
	else:
		GENOME_JELLYFISH = os.path.join(REFERENCES_FOLDER, 'grch38_canonical_chrs_chrM_{}mer.jf'.format(kmer_size))
	if not os.path.isfile(GENOME_JELLYFISH):
		logging.info('Genome jellyfish file does not exist. Please create with the name grch38_canonical_chrs_chrM_{}mer.jf and place it in the references folder'.format(kmer_size))
		print('error:\tGenome jellyfish file does not exist. Please create with the name grch38_canonical_chrs_chrM_{}mer.jf and place it in the references folder'.format(kmer_size))
		exit_program()



def convert_to_bed(vcf_input):
	# Converts a vcf file to a bed file
	def fix_deletions(deletions):
		# modify respective columns of deletion entries to match bed format
		deletions_lengths = deletions[3].str.len() - deletions[4].str.len()
		deletions[1] = deletions[1] + 1	# offset due to indexing of indels
		deletions[2] = deletions[2] + deletions_lengths
		deletions[3] = deletions[3].str[1:]
		deletions[4] = '-'
		return deletions

	def fix_insertions(insertions):
		# modify respective columns of inserton entries to match bed format
		insertions[3] = '-'
		insertions[4] = insertions[4].str[1:]
		return insertions

	variants_dataframe = pd.read_csv(vcf_input, sep='\t', comment='#', header=None)
	variants_dataframe.insert(1, 'start', variants_dataframe[1].copy() - 1)
	variants_dataframe.drop([2, 5, 6], axis=1, inplace=True)
	variants_dataframe.rename(columns={'start':1, 1:2, 7:5}, inplace=True)
	deletions = fix_deletions(variants_dataframe.loc[variants_dataframe[3].str.len() > variants_dataframe[4].str.len(), :].copy())
	insertions = fix_insertions(variants_dataframe.loc[variants_dataframe[4].str.len() > variants_dataframe[3].str.len(), :].copy())
	variants_dataframe.update(deletions)
	variants_dataframe.update(insertions)
	variants_dataframe[1] = variants_dataframe[1].astype(int)
	variants_dataframe[2] = variants_dataframe[2].astype(int)
	return variants_dataframe


def make_variant_bed_file(variants_dataframe, args):
	variants_dataframe['original_start'] = variants_dataframe['start'].copy()
	variants_dataframe['original_end'] = variants_dataframe['end'].copy()
	extract_start = variants_dataframe['start'].copy() - args.kmer_size + 1
	extract_end = variants_dataframe['start'].copy() + args.kmer_size
	variants_dataframe['start'] = extract_start
	variants_dataframe['end'] = extract_end

	# handle difference in extraction range for deletions
	deletions = variants_dataframe.loc[variants_dataframe['mutation'] == '-', :].copy()
	deletion_lengths = deletions['wildtype'].str.len()
	deletions['end'] += deletion_lengths
	deletions['start'] -= 1
	deletions['end'] -= 1
	variants_dataframe.update(deletions)

	# handle difference in extraction range for insertions
	insertions = variants_dataframe.loc[variants_dataframe['wildtype'] == '-', :].copy()
	insertions['start'] -= 1
	insertions['end'] -= 1
	variants_dataframe.update(insertions)

	variants_dataframe['start'], variants_dataframe['end'] = variants_dataframe['start'].astype(int), variants_dataframe['end'].astype(int)
	variants_dataframe['original_start'], variants_dataframe['original_end'] = variants_dataframe['original_start'].astype(int), variants_dataframe['original_end'].astype(int)
	output_path = os.path.join(CWD, args.output_name)
	variants_dataframe.to_csv(os.path.join(output_path, 'variants.bed'), sep='\t', index=False, header=False)
	os.system('sort -V {} >{}'.format(os.path.join(output_path, 'variants.bed'), os.path.join(output_path, 'sorted_variants.bed')))
	os.system('bedtools merge -c 2,3,4,5,6,7 -o collapse -i {} >{}'.format(os.path.join(output_path, 'sorted_variants.bed'), os.path.join(output_path, 'merged_variants.bed')))
	variants_file = os.path.join(output_path, 'merged_variants.bed') if args.multiple_mutations else os.path.join(output_path, 'sorted_variants.bed')
	new_variants_dataframe = pd.read_csv(variants_file, sep='\t', header=None)
	if new_variants_dataframe.shape[1] == MERGED_TABLE_COLUMN_COUNT:
		new_variants_dataframe.rename(columns={
			0: 'chromosome',
			1: 'start',
			2: 'end',
			3: 'all_starts',
			4: 'all_ends',
			5: 'wildtype',
			6: 'mutation',
			7: 'all_original_starts',
			8: 'all_original_ends'}, inplace=True)
		new_variants_dataframe['all_starts'] = new_variants_dataframe['all_starts'].apply(str)
		new_variants_dataframe['all_ends'] = new_variants_dataframe['all_ends'].apply(str)
		new_variants_dataframe['all_original_starts'] = new_variants_dataframe['all_original_starts'].apply(str)
		new_variants_dataframe['all_original_ends'] = new_variants_dataframe['all_original_ends'].apply(str)
	else:
		new_variants_dataframe.rename(columns={
			0: 'chromosome',
			1: 'start',
			2: 'end',
			3: 'wildtype',
			4: 'mutation',
			5: 'original_start',
			6: 'original_end'}, inplace=True)
	return new_variants_dataframe


def filter_variants(variants_dataframe, args):
	single_mutations, two_mutations, multiple_mutations = variants_dataframe, pd.DataFrame(), pd.DataFrame()
	if variants_dataframe.shape[1] == MERGED_TABLE_COLUMN_COUNT:
		single_mutations = variants_dataframe.loc[
				variants_dataframe['all_starts'].str.count(',') == 0,
				['chromosome', 'start', 'end', 'wildtype', 'mutation', 'all_original_starts', 'all_original_ends']]
		single_mutations.rename(columns={'all_original_starts': 'original_start', 'all_original_ends': 'original_end'}, inplace=True)
		multiple_mutations = variants_dataframe.loc[variants_dataframe['all_starts'].str.count(',') >= 1, :]
	single_mutations = single_mutations.loc[
			(single_mutations['wildtype'].str.len() <= args.kmer_size - MINIMUM_REQUIRED_FLANKING_REGION) &
			(single_mutations['mutation'].str.len() <= args.kmer_size - MINIMUM_REQUIRED_FLANKING_REGION), :]
	return single_mutations, multiple_mutations


# Split column given by column_name in dataframe on ',' delimiter. Return pandas dataframe with split entities per each column
def get_split_columns(dataframe, column_name):
	return dataframe[column_name].str.split(',', expand=True)

def filter_multiple_mutations(single_mutations, multiple_mutations, kmer_size, consider_multiple_mutations):
	two_mutations = pd.DataFrame(columns=multiple_mutations.columns.values)
	if not consider_multiple_mutations: return single_mutations, two_mutations	# no filtering necessary since merged entries are irrelevant
	new_single_mutations, new_two_mutations = [], []
	for index, row in multiple_mutations.iterrows():	# row iteration required due to variable mutation merge lengths
		wildtypes, mutations = get_split_columns(multiple_mutations, 'wildtype'), get_split_columns(multiple_mutations, 'mutation')
		starts, ends = get_split_columns(multiple_mutations, 'all_starts'), get_split_columns(multiple_mutations, 'all_ends')
		original_starts, original_ends = get_split_columns(multiple_mutations, 'all_original_starts'), get_split_columns(multiple_mutations, 'all_original_ends')
		for col in range(len(row['all_starts'].split(','))):
			new_single_mutation = {
					'chromosome': row['chromosome'],
					'start': int(row['all_starts'].split(',')[col]),
					'end': int(row['all_ends'].split(',')[col]),
					'wildtype': row['wildtype'].split(',')[col],
					'mutation': row['mutation'].split(',')[col],
					'original_start': row['all_original_starts'].split(',')[col],
					'original_end': row['all_original_ends'].split(',')[col]
				}
			new_single_mutations.append(new_single_mutation)
		for col in range(len(row['all_starts'].split(',')) - 1):
			starts_split, ends_split = [int(x) for x in row['all_starts'].split(',')], [int(x) for x in row['all_ends'].split(',')]
			original_starts_split, original_ends_split = [int(x) for x in row['all_original_starts'].split(',')], [int(x) for x in row['all_original_ends'].split(',')]
			wildtypes_split, mutations_split = row['wildtype'].split(','), row['mutation'].split(',')
			if (starts_split[col + 1] - starts_split[col]) + 1 > kmer_size or \
				starts_split[col + 1] == starts_split[col] or \
				len(wildtypes_split[col]) > len(mutations_split[col]) and \
				starts_split[col] + len(wildtypes_split[col]) > starts_split[col + 1]:
					continue # if kmer does not encompass both mutations or if deletion and first region crosses into second, skip
			new_two_mutation = {
					'chromosome': row['chromosome'],
					'start': int(row['all_starts'].split(',')[col]),
					'end': int(row['all_ends'].split(',')[col + 1]),
					'all_starts': '{},{}'.format(starts_split[col], starts_split[col + 1]),
					'all_ends': '{},{}'.format(ends_split[col], ends_split[col + 1]),
					'wildtype': '{},{}'.format(wildtypes_split[col], wildtypes_split[col + 1]),
					'mutation': '{},{}'.format(mutations_split[col], mutations_split[col + 1]),
					'all_original_starts': '{},{}'.format(original_starts_split[col], original_starts_split[col + 1]),
					'all_original_ends': '{},{}'.format(original_ends_split[col], original_ends_split[col + 1])
				}
			new_two_mutations.append(new_two_mutation)
	if not len(new_single_mutations) == 0:
		single_mutations = single_mutations.append(new_single_mutations, ignore_index=True)
	two_mutations = pd.DataFrame(new_two_mutations, columns=multiple_mutations.columns.values)
	return single_mutations, two_mutations


def extract_variant_sequences(single_mutations, two_mutations, args):
	def get_output_name(filename):
		return os.path.join(output_path, filename)
	output_path = os.path.join(CWD, args.output_name)
	if not single_mutations.empty:
		single_mutations['fasta_id'] = single_mutations.apply(
			lambda x: '_'.join(x[['chromosome', 'original_start', 'original_end', 'wildtype', 'mutation', 'start', 'end']].values.astype(str).tolist()),
			axis=1)
	if not two_mutations.empty:
		two_mutations['fasta_id'] = two_mutations.apply(
			lambda x: '_'.join(x[['chromosome', 'all_original_starts', 'all_original_ends', 'wildtype', 'mutation', 'all_starts', 'all_ends']].values.astype(str).tolist()),
			axis=1)
	single_mutations.drop(['original_start', 'original_end', 'wildtype', 'mutation'], axis=1, inplace=True)
	single_variant_bed = get_output_name('single_variants.bed')
	single_variant_fasta = get_output_name('single_variants.fa')
	single_mutations.to_csv(single_variant_bed, sep='\t', index=False, header=False)
	os.system('bedtools getfasta -fi {} -fo {} -bed {} -nameOnly'.format(GENOME_FASTA, single_variant_fasta, single_variant_bed))
	two_variant_bed = get_output_name('two_variants.bed')
	two_variant_fasta = get_output_name('two_variants.fa')
	if args.multiple_mutations:
		two_mutations.drop(['all_starts', 'all_ends', 'all_original_starts', 'all_original_ends', 'wildtype', 'mutation'], axis=1, inplace=True)
		two_mutations.to_csv(two_variant_bed, sep='\t', index=False, header=False)
		os.system('bedtools getfasta -fi {} -fo {} -bed {} -nameOnly'.format(GENOME_FASTA, two_variant_fasta, two_variant_bed))
	return single_variant_fasta, two_variant_fasta


def make_mutation_sequences(single_variant_fasta, two_variant_fasta, kmer_size, output_name, multiple_mutations):
	def modified_single_variant_sequences(fasta_entry_name, sequence):
		chromosome, original_start, original_end, wildtype, mutation, start, end = fasta_entry_name.split('_')
		is_insertion, is_deletion = wildtype == '-', mutation == '-'
		if is_insertion:
			wildtype_sequence = sequence[1:]
			mutation_sequence = sequence[1 + len(mutation): kmer_size] + mutation + sequence[kmer_size: -len(mutation)]
		elif is_deletion:
			wildtype_sequence = sequence[len(wildtype): -len(wildtype) - 1]
			mutation_sequence = sequence[: kmer_size] + sequence[kmer_size + len(wildtype): - 1]
		else:	# is_substitution
			wildtype_sequence = sequence
			mutation_sequence = sequence[: kmer_size - 1] + mutation + sequence[kmer_size: ]
		return (wildtype_sequence, mutation_sequence)

	def modified_two_variant_sequences(fasta_entry_name, sequence):
		chromosome, original_start, original_end, wildtype, mutation, start, end = fasta_entry_name.split('_')
		def split_two(entry): return entry.split(',')
		starts, ends = [int(x) for x in split_two(start)], [int(x) for x in split_two(end)]
		wildtypes, mutations = split_two(wildtype), split_two(mutation)
		first_is_insertion, first_is_deletion = wildtypes[0] == '-', mutations[0] == '-'
		second_is_insertion, second_is_deletion = wildtypes[1] == '-', mutations[1] == '-'
		variants_distance = starts[1] - starts[0]
		first_variant_mutation_start, second_variant_mutation_end = kmer_size - 1, variants_distance + kmer_size - 1
		first_variant_wildtype_start, second_variant_wildtype_end = kmer_size - 1, variants_distance + kmer_size - 1
		if second_is_insertion:
			semi_mutation_sequence = sequence[: second_variant_mutation_end + 1] + mutations[1] + sequence[second_variant_mutation_end + 1: ]
			second_variant_mutation_end += 1 + len(mutations[1])
			second_variant_wildtype_end += 1
		elif second_is_deletion:
			semi_mutation_sequence = sequence[: second_variant_mutation_end] + sequence[second_variant_mutation_end + len(wildtype[1]): ]
			second_variant_mutation_end -= (len(wildtype[1]) - 1)
			second_variant_wildtype_end += 1
		else:	# second_is_substitution
			semi_mutation_sequence = sequence[: second_variant_mutation_end] + mutations[1] + sequence[second_variant_mutation_end + 1: ]

		if first_is_insertion:
			mutation_sequence = semi_mutation_sequence[: kmer_size] + mutations[0] + semi_mutation_sequence[kmer_size: ]
			second_variant_mutation_end += len(mutations[0])
		elif first_is_deletion:
			mutation_sequence = semi_mutation_sequence[: kmer_size] + semi_mutation_sequence[kmer_size + len(wildtype[0]): ]
			second_variant_mutation_end -= len(wildtype[0])
			first_variant_mutation_start -= 1
			first_variant_wildtype_start -= 1
		else:	# first_is_substitution
			mutation_sequence = semi_mutation_sequence[: kmer_size - 1] + mutations[0] + semi_mutation_sequence[kmer_size: ]

		wildtype_extract_start, wildtype_extract_end = second_variant_wildtype_end - kmer_size + 1, first_variant_wildtype_start + kmer_size
		mutation_extract_start, mutation_extract_end = second_variant_mutation_end - kmer_size + 1, first_variant_mutation_start + kmer_size
		return(sequence[wildtype_extract_start: wildtype_extract_end], mutation_sequence[mutation_extract_start: mutation_extract_end])

	def open_and_read_file(variant_fasta, modification_function):
		with open(variant_fasta, 'r') as f:
			variant_wildtype_sequences = {
					entry.rstrip().split('\n')[0]: entry.rstrip().split('\n')[1]
					for entry in f.read().split('>')[1:]
				}
		variant_wildtype_and_mutation_sequences = {
				name: modification_function(name, sequence)
				for name, sequence in variant_wildtype_sequences.items()
			}
		return variant_wildtype_and_mutation_sequences

	def write_output(variant_sequences, output_file, writing_mutation):
		idx = 1 if writing_mutation else 0
		for variant, sequences in variant_sequences.items():
			output_file.write('>{}\n{}\n'.format(variant, sequences[idx]))

	single_variant_sequences = open_and_read_file(single_variant_fasta, modified_single_variant_sequences)
	output_path = os.path.join(CWD, output_name)
	wildtype_path, mutation_path = os.path.join(output_path, 'wildtype_extracted_sequences.fa'), os.path.join(output_path, 'mutation_extracted_sequences.fa')
	wildtype_output, mutation_output = open(wildtype_path, 'w+'), open(mutation_path, 'w+')
	write_output(single_variant_sequences, wildtype_output, False)
	write_output(single_variant_sequences, mutation_output, True)
	if multiple_mutations:
		two_variant_sequences = open_and_read_file(two_variant_fasta, modified_two_variant_sequences)
		write_output(two_variant_sequences, wildtype_output, False)
		write_output(two_variant_sequences, mutation_output, True)
	wildtype_output.close(); mutation_output.close()
	return wildtype_path, mutation_path


def query_sequences(wildtype_sequence_file, mutation_sequence_file, output_name, jellyfish_test, jellyfish_control):
	def query_kmer_frequencies(sequence_file, output_filename, jellyfish_file):
		output_filepath = os.path.join(CWD, output_name, output_filename)
		os.system('jellyfish query -s {} -o {} {}'.format(sequence_file, output_filepath, jellyfish_file))
		return output_filepath

	test_wildtype_file = query_kmer_frequencies(wildtype_sequence_file, 'wildtype_counts_in_test_sample.txt', jellyfish_test)
	test_mutation_file = query_kmer_frequencies(mutation_sequence_file, 'mutation_counts_in_test_sample.txt', jellyfish_test)
	genome_wildtype_file = query_kmer_frequencies(wildtype_sequence_file, 'wildtype_counts_in_reference_genome.txt', GENOME_JELLYFISH)	# reference genome counts
	genome_mutation_file = query_kmer_frequencies(mutation_sequence_file, 'mutation_counts_in_reference_genome.txt', GENOME_JELLYFISH)	# reference genome counts
	control_wildtype_file, control_mutation_file = None, None
	if jellyfish_control:
		control_wildtype_file = query_kmer_frequencies(wildtype_sequence_file, 'wildtype_counts_in_control_sample.txt', jellyfish_control)
		control_mutation_file = query_kmer_frequencies(mutation_sequence_file, 'mutation_counts_in_control_sample.txt', jellyfish_control)

	return (test_wildtype_file, test_mutation_file, control_wildtype_file, control_mutation_file, genome_wildtype_file, genome_mutation_file)

def construct_variant_information_table(kmer_frequency_files, kmer_size, using_control_sample, wildtype_path, mutation_path):
	#TODO refactor this function
	test_wildtype_file, test_mutation_file, control_wildtype_file, control_mutation_file, genome_wildtype_file, genome_mutation_file = kmer_frequency_files
	wildtype_in_genome, mutation_in_genome = dict_from_kmer_counts(genome_wildtype_file), dict_from_kmer_counts(genome_mutation_file)
	wildtype_in_test, mutation_in_test = dict_from_kmer_counts(test_wildtype_file), dict_from_kmer_counts(test_mutation_file)
	if using_control_sample:
		wildtype_in_control, mutation_in_control = dict_from_kmer_counts(control_wildtype_file), dict_from_kmer_counts(control_mutation_file)

	wildtype_sequences_file, mutation_sequences_file = open(wildtype_path, 'r'), open(mutation_path, 'r')
	sequences = { entry.rstrip().split('\n')[0]: [entry.rstrip().split('\n')[1]] for entry in wildtype_sequences_file.read().split('>')[1:] }
	mutation_sequences = [entry.rstrip().split('\n') for entry in mutation_sequences_file.read().split('>')[1:]]

	for entry in mutation_sequences:	# add to list in sequences dictionary
		sequences[entry[0]].append(entry[1])
	wildtype_sequences_file.close(); mutation_sequences_file.close()

	kmer_counts_per_variant = {}
	for variant in sequences:
		wildtype_sequence, mutation_sequence = sequences[variant]
		wildtype_kmers = set([wildtype_sequence[i : i + kmer_size] for i in range(len(wildtype_sequence) - kmer_size + 1)])
		mutation_kmers = [mutation_sequence[i : i + kmer_size] for i in range(len(mutation_sequence) - kmer_size + 1)]
		mutation_kmers = [kmer for kmer in mutation_kmers if kmer not in wildtype_kmers]	# remove kmer in mutation_kmers that exists in wildtpye_kmers

		if using_control_sample:
			wildtype_out = [ (kmer, wildtype_in_control[kmer], wildtype_in_test[kmer]) if kmer in wildtype_in_test
							else (kmer, wildtype_in_control[reverse_complement(kmer)], wildtype_in_test[reverse_complement(kmer)])
							for kmer in wildtype_kmers ]
			mutation_out = [ (kmer, mutation_in_control[kmer], mutation_in_test[kmer]) if kmer in mutation_in_test
							else (kmer, mutation_in_control[reverse_complement(kmer)], mutation_in_test[reverse_complement(kmer)])
							for kmer in mutation_kmers ]
		else:
			wildtype_out = [ (kmer, None, wildtype_in_test[kmer]) if kmer in wildtype_in_test
							else (kmer, None, wildtype_in_test[reverse_complement(kmer)])
							for kmer in wildtype_kmers ]
			mutation_out = [ (kmer, None, mutation_in_test[kmer]) if kmer in mutation_in_test
							else (kmer, None, mutation_in_test[reverse_complement(kmer)])
							for kmer in mutation_kmers ]
		kmer_counts_per_variant[variant] = (wildtype_out, mutation_out)

	wildtype_sequences_file, mutation_sequences_file = open(wildtype_path, 'r'), open(mutation_path, 'r')
	wildtype_sequences = { entry.rstrip().split('\n')[0]: entry.rstrip().split('\n')[1] for entry in wildtype_sequences_file.read().split('>')[1:] }
	mutation_sequences = { entry.rstrip().split('\n')[0]: entry.rstrip().split('\n')[1] for entry in mutation_sequences_file.read().split('>')[1:] }
	wildtype_sequences_file.close(); mutation_sequences_file.close()
	all_entries = []

	def get_median(counts): return 0 if len(counts) == 0 else round(np.median(counts), 3)

	for variant in kmer_counts_per_variant:
		wildtype_unique, mutation_zero = [], []
		wildtype_kmers, mutation_kmers = kmer_counts_per_variant[variant]
		for kmer_info in wildtype_kmers:
			kmer = kmer_info[0] if kmer_info[0] in wildtype_in_genome else reverse_complement(kmer_info[0])
			if wildtype_in_genome[kmer] == 1: wildtype_unique.append((kmer_info[1], kmer_info[2]))
		for kmer_info in mutation_kmers:
			kmer = kmer_info[0] if kmer_info[0] in mutation_in_genome else reverse_complement(kmer_info[0])
			if mutation_in_genome[kmer] == 0: mutation_zero.append((kmer_info[1], kmer_info[2]))
		# in_control lists composed of Nones when not using control_sample
		wildtype_in_control, wildtype_in_test = [count[0] for count in wildtype_unique], [count[1] for count in wildtype_unique]
		mutation_in_control, mutation_in_test = [count[0] for count in mutation_zero], [count[1] for count in mutation_zero]

		wildtype_control_median, wildtype_test_median, mutation_control_median, mutation_test_median = 0, 0, 0, 0
		is_sufficient_data = True if len(wildtype_unique) >= SUFFICIENT_DATA_MINIMUM_KMERS and len(mutation_zero) >= SUFFICIENT_DATA_MINIMUM_KMERS else False
		if is_sufficient_data:
			wildtype_in_test, mutation_in_test = [float(x) for x in wildtype_in_test], [float(x) for x in mutation_in_test]
			wildtype_test_median, mutation_test_median = get_median(wildtype_in_test), get_median(mutation_in_test)
			if using_control_sample:
				wildtype_in_control, mutation_in_control = [float(x) for x in wildtype_in_control], [float(x) for x in mutation_in_control]
				wildtype_control_median, mutation_control_median = get_median(wildtype_in_control), get_median(mutation_in_control)

		#all_wildtype_in_control_counts = None if len(wildtype_in_control) == 0 else '|'.join([str(count) for count in wildtype_in_control])
		#all_mutation_in_test_counts = None if len(mutation_in_test) == 0 else '|'.join([str(count) for count in mutation_in_test])
		new_entry = {
				VARIANT_COL_NAME: variant,
				WILDTYPE_SEQUENCE_COL_NAME: wildtype_sequences[variant],
				WILDTYPE_TEST_COUNT_COL_NAME: wildtype_test_median,
				WILDTYPE_KMER_COUNT_COL_NAME: len(wildtype_sequences[variant]) - (kmer_size - 1),
				MUTATION_SEQUENCE_COL_NAME: mutation_sequences[variant],
				MUTATION_TEST_COUNT_COL_NAME: mutation_test_median,
				MUTATION_KMER_COUNT_COL_NAME: len(mutation_sequences[variant]) - (kmer_size - 1),
				WILDTYPE_UNIQUE_COUNT_COL_NAME: len(wildtype_unique),
				MUTATION_ZERO_COUNT_COL_NAME: len(mutation_zero),
		}
		if using_control_sample:
			new_entry[WILDTYPE_CONTROL_COUNT_COL_NAME] = wildtype_control_median
			new_entry[MUTATION_CONTROL_COUNT_COL_NAME] = mutation_control_median
		all_entries.append(new_entry)

	all_columns = [
		VARIANT_COL_NAME,
		WILDTYPE_SEQUENCE_COL_NAME,
		WILDTYPE_CONTROL_COUNT_COL_NAME,
		WILDTYPE_TEST_COUNT_COL_NAME,
		WILDTYPE_KMER_COUNT_COL_NAME,
		MUTATION_SEQUENCE_COL_NAME,
		MUTATION_CONTROL_COUNT_COL_NAME,
		MUTATION_TEST_COUNT_COL_NAME,
		MUTATION_KMER_COUNT_COL_NAME,
		WILDTYPE_UNIQUE_COUNT_COL_NAME,
		MUTATION_ZERO_COUNT_COL_NAME]
	if not using_control_sample:
		all_columns = [column for column in all_columns if column not in [WILDTYPE_CONTROL_COUNT_COL_NAME, MUTATION_CONTROL_COUNT_COL_NAME]]

	variant_info_dataframe = pd.DataFrame(all_entries, columns=all_columns)
	return variant_info_dataframe


def hypothesis_test(variant_info_dataframe, kmer_size, using_poisson, using_control_sample, output_name):
	# Poisson distribution pmf with mean mu for k in { 0, 1, 2, ... }
	def poisson_pmf(k, mu): return np.exp(k * np.log(mu) - mu - special.gammaln(k + 1))

	# P(X >= k) for poisson distributed random variable X with mean mu
	def poisson_cdf(k, mu): return 1.0 - sum(poisson_pmf(k_, mu) for k_ in xrange(0, k))

	# P(X >= k) for binomial distributed random variable X with probability prob
	#def binomial_cdf(k, n, p): return 1.0 - stats.binom.cdf(k - 1, n, p)

	def binomialpmf(k, n, p):
		# Binomial distribution pmf with mean n*p on set {0,1,...,n} at k in {0,1,...,n}
		p = np.exp(k*np.log(p) + (n-k)*np.log(1-p) + special.gammaln(n+1) - special.gammaln(k+1) - special.gammaln(n-k+1))
		return p

	def binomial_cdf(k, n, p):
	# P(X >= k) for X Binomial-distributed with mean n*p on set {0,1,...,n}
		return 1 - sum(binomialpmf(k_, n, p) for k_ in range(k))

	sample_count = variant_info_dataframe.shape[0]
	alpha = ALPHA / sample_count # bonferroni correction
	if using_control_sample:
		control_variant_indices, test_variant_indices = [], []
		control_normal_indices = []
		control_p_values, test_p_values = {}, {}
		# Select control variants. selects the entries for which we accept the hypothesis M_N = 0
		for index, row in variant_info_dataframe.iterrows():
			wildtype, mutation = row[VARIANT_COL_NAME].split('_')[3:5]
			is_indel = wildtype == '-' or mutation == '-'
			indel_size = len(wildtype.replace('-', '')) - len(mutation.replace('-', ''))
			power = 1 # + 0.25 * (indel_size - 1) if is_indel else 1
			sequence_error_probability = np.power(SEQUENCE_ERROR_PROBABILITY, power)
			total_control_kmer_count = row[WILDTYPE_CONTROL_COUNT_COL_NAME] + row[MUTATION_CONTROL_COUNT_COL_NAME]
			mutation_count_mean_estimate = total_control_kmer_count * np.power(1 - sequence_error_probability, kmer_size - 1) * (sequence_error_probability / 3)

			if np.int(row[MUTATION_CONTROL_COUNT_COL_NAME]) > 0: # Need > 0 trials to perform binomial test
				p_value = poisson_cdf(np.int(row[MUTATION_CONTROL_COUNT_COL_NAME]), np.int(mutation_count_mean_estimate)) if using_poisson else \
					binomial_cdf(np.int16(row[MUTATION_CONTROL_COUNT_COL_NAME]), np.int16(total_control_kmer_count), sequence_error_probability)
				if p_value < alpha: control_variant_indices.append(index)
				else: control_normal_indices.append(index)
				control_p_values[index] = p_value
			else:
				control_p_values[index] = '-'

		# Among the entries for which we accept the first hypothesis, we select the entries for which we reject the hypothesis M_T = 0
		for index, row in variant_info_dataframe.iterrows():
			wildtype, mutation = row[VARIANT_COL_NAME].split('_')[3:5]
			is_indel = wildtype == '-' or mutation == '-'
			indel_size = len(wildtype.replace('-', '')) - len(mutation.replace('-', ''))
			power = 1 # + 0.25 * (indel_size - 1) if is_indel else 1
			sequence_error_probability = np.power(SEQUENCE_ERROR_PROBABILITY, power)
			total_test_kmer_count = row[WILDTYPE_TEST_COUNT_COL_NAME] + row[MUTATION_TEST_COUNT_COL_NAME]
			mutation_count_mean_estimate = total_test_kmer_count * np.power(1 - sequence_error_probability, kmer_size - 1) * (sequence_error_probability / 3)
			if np.int(row[MUTATION_TEST_COUNT_COL_NAME]) > 0: # Need > 0 trials to perform binomial test
				p_value = poisson_cdf(np.int(row[MUTATION_TEST_COUNT_COL_NAME]), np.int(mutation_count_mean_estimate)) if using_poisson else \
					binomial_cdf(np.int16(row[MUTATION_TEST_COUNT_COL_NAME]), np.int16(total_test_kmer_count), sequence_error_probability)
				if p_value < alpha: test_variant_indices.append(index)
				test_p_values[index] = p_value
			else:
				test_p_values[index] = '-'

		reject_null_control_variants, reject_null_test_variants, variant_calls= [], [], []
		control_p_value_list, test_p_value_list = [], []
		for index in range(variant_info_dataframe.shape[0]):
			reject_null_control_value = True if index in control_variant_indices else False
			reject_null_test_value = True if index in test_variant_indices else False
			variant_call_value = not reject_null_control_value and reject_null_test_value # Need reject_null_control = False and reject_null_test = True
			variant_call_value = True if variant_call_value else False
			reject_null_control_variants.append(reject_null_control_value)
			reject_null_test_variants.append(reject_null_test_value)
			variant_calls.append(variant_call_value)
			control_p_value_list.append(control_p_values[index])
			test_p_value_list.append(test_p_values[index])

		variant_info_dataframe.insert(1, VARIANT_CALL_COL_NAME, variant_calls)
		variant_info_dataframe.insert(2, REJECT_NULL_CONTROL_COL_NAME, reject_null_control_variants)
		variant_info_dataframe.insert(3, CONTROL_P_VALUE_COL_NAME, control_p_value_list)
		variant_info_dataframe.insert(4, REJECT_NULL_TEST_COL_NAME, reject_null_test_variants)
		variant_info_dataframe.insert(5, TEST_P_VALUE_COL_NAME, test_p_value_list)
		output_path = os.path.join(CWD, output_name)
		variant_info_dataframe.to_csv(os.path.join(output_path, 'variant_info_summary.txt'), sep='\t', index=False)

	else:
		test_variant_indices, test_normal_indices, test_p_values = [], [], {}
		# Select control variants. selects the entries for which we accept the hypothesis M_N = 0
		for index, row in variant_info_dataframe.iterrows():
			wildtype, mutation = row[VARIANT_COL_NAME].split('_')[3:5]
			is_indel = wildtype == '-' or mutation == '-'
			indel_size = len(wildtype.replace('-', '')) - len(mutation.replace('-', ''))
			power = 1 # + 0.25 * (indel_size - 1) if is_indel else 1
			sequence_error_probability = np.power(SEQUENCE_ERROR_PROBABILITY, power)
			total_test_kmer_count = row[WILDTYPE_TEST_COUNT_COL_NAME] + row[MUTATION_TEST_COUNT_COL_NAME]
			mutation_count_mean_estimate = total_test_kmer_count * np.power(1 - sequence_error_probability, kmer_size - 1) * (sequence_error_probability / 3)

			if np.int(row[MUTATION_TEST_COUNT_COL_NAME]) > 0: # Need > 0 trials to perform binomial test
				p_value = poisson_cdf(np.int(row[MUTATION_CONTROL_COUNT_COL_NAME]), np.int(mutation_count_mean_estimate)) if using_poisson else \
					binomial_cdf(np.int(row[MUTATION_TEST_COUNT_COL_NAME]), np.int(total_test_kmer_count), sequence_error_probability)
				if p_value < alpha: test_variant_indices.append(index)
				else: test_variant_indices.append(index)
				test_p_values[index] = p_value
			else:
				test_p_values[index] = '-'


		reject_null_test_variants, variant_calls = [], []
		test_p_value_list = []
		for index in range(variant_info_dataframe.shape[0]):
			reject_null_test_value = True if index in test_variant_indices else False
			variant_call_value = reject_null_test_value
			reject_null_test_variants.append(reject_null_test_value)
			variant_calls.append(variant_call_value)
			test_p_value_list.append(test_p_values[index])

		variant_info_dataframe.insert(1, VARIANT_CALL_COL_NAME, variant_calls)
		variant_info_dataframe.insert(2, REJECT_NULL_TEST_COL_NAME, reject_null_test_variants)
		variant_info_dataframe.insert(3, TEST_P_VALUE_COL_NAME, test_p_value_list)
		output_path = os.path.join(CWD, output_name)
		variant_info_dataframe.to_csv(os.path.join(output_path, 'variant_info_summary.txt'), sep='\t', index=False)

	return variant_info_dataframe

def create_variant_call_summary_table(variant_call_info_dataframe, command_args, kmer_frequency_files, output_name, original_dataframe, using_control_sample):
	def estimate_sequencing_coverage():	# TODO: Use only unique kmers
		test_wildtype_file, test_mutation_file, control_wildtype_file, control_mutation_file, genome_wildtype_file, genome_mutation_file = kmer_frequency_files
		wildtype_in_genome = dict_from_kmer_counts(genome_wildtype_file)
		test_counts_dataframe = pd.read_csv(test_wildtype_file, sep=' ', header=None, names=['Kmer', 'Count'])
		test_counts_dataframe['Genome_Count'] = test_counts_dataframe['Kmer'].apply(lambda kmer: wildtype_in_genome[kmer] if kmer in wildtype_in_genome else wildtype_in_genome[reverse_complement(kmer)])
		test_counts_dataframe = test_counts_dataframe[test_counts_dataframe['Genome_Count'] == 1]
		test_sequencing_coverage = test_counts_dataframe['Count'].median()
		control_sequencing_coverage = None
		if using_control_sample:
			control_counts_dataframe = pd.read_csv(control_wildtype_file, sep=' ', header=None, names=['Kmer', 'Count']) if control_wildtype_file != None else None
			control_counts_dataframe['Genome_Count'] = control_counts_dataframe['Kmer'].apply(lambda kmer: wildtype_in_genome[kmer] if kmer in wildtype_in_genome else wildtype_in_genome[reverse_complement(kmer)]) if control_wildtype_file != None else None
			control_counts_dataframe = control_counts_dataframe[control_counts_dataframe['Genome_Count'] == 1] if control_wildtype_file != None else None
			control_sequencing_coverage = control_counts_dataframe['Count'].median() if control_wildtype_file != None else None
		return control_sequencing_coverage, test_sequencing_coverage
	control_sequencing_coverage, test_sequencing_coverage = estimate_sequencing_coverage()
	validated_call_minimum_required_count = stats.binom.ppf(1 - ALPHA, test_sequencing_coverage, SEQUENCE_ERROR_PROBABILITY)	# TODO control/test coverage usage distinction
	def is_sufficient_data(row): return row[WILDTYPE_UNIQUE_COUNT_COL_NAME] >= SUFFICIENT_DATA_MINIMUM_KMERS and row[MUTATION_ZERO_COUNT_COL_NAME] >= SUFFICIENT_DATA_MINIMUM_KMERS
	def is_snp_affected(row): return is_sufficient_data(row) and row[WILDTYPE_CONTROL_COUNT_COL_NAME] < validated_call_minimum_required_count
	def is_validated(row): return is_sufficient_data(row) and not is_snp_affected(row) and row[VARIANT_CALL_COL_NAME]
	def is_germline(row): return row[REJECT_NULL_CONTROL_COL_NAME]
	def is_multiple(row): return len(row[VARIANT_COL_NAME].split(',')) > 1
	variant_types = []
	for index, row in variant_call_info_dataframe.iterrows():
		if using_control_sample:
			if is_snp_affected(row): variant_type = 'SNP'
			elif not is_sufficient_data(row): variant_type = 'Insufficient'
			elif is_validated(row): variant_type = 'Validated'
			elif is_germline(row): variant_type = 'Germline'
			else: variant_type = 'Invalidated'
			if is_multiple(row): variant_type = '{}_Multiple'.format(variant_type)
		else:
			if not is_sufficient_data(row): variant_type = 'Insufficient'
			elif is_sufficient_data(row) and row[VARIANT_CALL_COL_NAME]: variant_type = 'Validated'
			else: variant_type = 'Invalidated'
			if is_multiple(row): variant_type = '{}_Multiple'.format(variant_type)
		variant_types.append(variant_type)
	variant_call_info_dataframe.insert(1, VARIANT_TYPE_COL_NAME, variant_types)
	variant_type_counts = variant_call_info_dataframe[VARIANT_TYPE_COL_NAME].value_counts()
	variant_call_info_dataframe['Variant'] = variant_call_info_dataframe['Variant'].apply(lambda x: '_'.join(x.split('_')[:-2]))
	variant_call_info_dataframe.drop_duplicates(inplace=True, subset='Variant')
	variant_call_info_dataframe = pd.merge(variant_call_info_dataframe, original_dataframe, on='Variant', how='outer')
	variant_call_info_dataframe = variant_call_info_dataframe.set_index(['Variant'])
	multiple_variants_dataframe = variant_call_info_dataframe.drop(original_dataframe['Variant'])	# Since dropped in indexing
	variant_call_info_dataframe = variant_call_info_dataframe.drop(multiple_variants_dataframe.index)
	variant_call_info_dataframe = variant_call_info_dataframe.reindex(original_dataframe['Variant'])
	variant_call_info_dataframe = variant_call_info_dataframe.append(multiple_variants_dataframe)
	sample_count = variant_call_info_dataframe.shape[0]
	with open(os.path.join(CWD, '{}_penultimate_variant_summary_table.txt'.format(output_name)), 'w+') as output:
		penultimate_variant_call_info_dataframe = variant_call_info_dataframe.drop(
			[VARIANT_TYPE_COL_NAME,
			REJECT_NULL_CONTROL_COL_NAME,
			REJECT_NULL_TEST_COL_NAME,
			VARIANT_CALL_COL_NAME], axis=1)
		penultimate_variant_call_info_dataframe.to_csv(output, sep='\t', index=True)
	with open(os.path.join(CWD, '{}_variant_summary_table.txt'.format(output_name)), 'w+') as output:
		output.write('Command:{}\n'.format(' '.join(['python'] + command_args)))
		if using_control_sample: variant_type_counts.to_csv(output, sep=':', index=True, header=False)
			#output.write('Control_Estimated_Sequencing_Coverage:{}\n'.format(control_sequencing_coverage))
		#output.write('Test_Estimated_Sequencing_Coverage:{}\n'.format(test_sequencing_coverage))
		#variant_type_counts.to_csv(output, sep=':', index=True)
		output.write('Alpha:{}\n'.format(ALPHA))
		output.write('Adjusted alpha value with bonferroni correction by n={}\n'.format(sample_count))
		variant_call_info_dataframe.to_csv(output, sep='\t', index=True)

def validate_mutations(penultimate_dataframe, alpha, output_name, command_args):
	# Given df, get Reject null control/test and give validations
	# TODO - account for single sample case
	reject_null_control_variants, reject_null_test_variants, variant_calls= [], [], []
	control_p_value_list, test_p_value_list = [], []
	for index, row in penultimate_dataframe.iterrows():
		reject_null_control_value = False
		reject_null_test_value = False
		control_p_val, test_p_val = row[[CONTROL_P_VALUE_COL_NAME, TEST_P_VALUE_COL_NAME]]
		if control_p_val != '-' and float(control_p_val) < alpha:
			reject_null_control_value = True
		if test_p_val != '-' and float(test_p_val) < alpha:
			reject_null_test_value = True
		variant_call_value = not reject_null_control_value and reject_null_test_value # Need reject_null_control = False and reject_null_test = True
		variant_call_value = True if variant_call_value else False
		reject_null_control_variants.append(reject_null_control_value)
		reject_null_test_variants.append(reject_null_test_value)
		variant_calls.append(variant_call_value)
	penultimate_dataframe.insert(1, VARIANT_CALL_COL_NAME, variant_calls)
	penultimate_dataframe.insert(2, REJECT_NULL_CONTROL_COL_NAME, reject_null_control_variants)
	penultimate_dataframe.insert(4, REJECT_NULL_TEST_COL_NAME, reject_null_test_variants)

	def is_sufficient_data(row): return row[WILDTYPE_UNIQUE_COUNT_COL_NAME] >= SUFFICIENT_DATA_MINIMUM_KMERS and row[MUTATION_ZERO_COUNT_COL_NAME] >= SUFFICIENT_DATA_MINIMUM_KMERS
	def is_snp_affected(row): return is_sufficient_data(row)  and row[WILDTYPE_CONTROL_COUNT_COL_NAME] < 2 # TODO- CHANGE hardcodoing to validated_call_minimum_required_count
	def is_validated(row): return is_sufficient_data(row) and not is_snp_affected(row) and row[VARIANT_CALL_COL_NAME]
	def is_germline(row): return row[REJECT_NULL_CONTROL_COL_NAME]
	def is_multiple(row): return len(row[VARIANT_COL_NAME].split(',')) > 1
	variant_types = []
	for index, row in penultimate_dataframe.iterrows():
		if is_snp_affected(row): variant_type = 'SNP'
		elif not is_sufficient_data(row): variant_type = 'Insufficient'
		elif is_validated(row): variant_type = 'Validated'
		elif is_germline(row): variant_type = 'Germline'
		else: variant_type = 'Invalidated'
		if is_multiple(row): variant_type = '{}_Multiple'.format(variant_type)
		variant_types.append(variant_type)

	penultimate_dataframe.insert(1, VARIANT_TYPE_COL_NAME, variant_types)
	variant_type_counts = penultimate_dataframe[VARIANT_TYPE_COL_NAME].value_counts()
	sample_count = penultimate_dataframe.shape[0]
	with open(os.path.join(CWD, '{}_variant_summary_table.txt'.format(output_name)), 'w+') as output:
		output.write('Command:{}\n'.format(' '.join(['python'] + command_args)))
		variant_type_counts.to_csv(output, sep=':', index=True)
		output.write('Alpha:{}\n'.format(ALPHA))
		output.write('Adjusted alpha value with bonferroni correction by n={}\n'.format(sample_count))
		penultimate_dataframe.to_csv(output, sep='\t', index=False)


#############################################################
#	 					Set Up util File					#
#############################################################
setup_logger()
