#!/usr/bin/env python
from kmervclib import *

def run_kmervc():
	parser = get_arg_parser()
	args = parser.parse_args()
	if args.subcommand == 'compare':
		logging.info('Performing compare subcommand. Calling perform_compare')
		perform_compare(args)
	elif args.subcommand == 'assess':
		logging.info('Performing assess subcommand. Calling perform_assess')
		perform_assess(args)

def perform_compare(args):
	jellyfish_test, jellyfish_control = args.jellyfish_test, args.jellyfish_control
	if not args.jellyfish_test and not args.jellyfish_control:  # jellyfish not provided
		jellyfish_test, jellyfish_control = create_jellyfish(args)
	setup_output_directory(args.output_name, args.reference_genome_fasta, args.kmer_size, args.alpha)
	variants_dataframe, original_dataframe = construct_variants_dataframe(args)	# consolidates vcf/bed into one standard dataframe
	variants_dataframe = make_variant_bed_file(variants_dataframe, args)	# handles sorting/merging/writing new variants dataframe
	single_mutations, multiple_mutations = filter_variants(variants_dataframe, args)
	single_mutations, two_mutations = filter_multiple_mutations(single_mutations, multiple_mutations, args.kmer_size, args.multiple_mutations)
	single_variant_fasta, two_variant_fasta = extract_variant_sequences(single_mutations, two_mutations, args)
	wildtype_sequence_file, mutation_sequence_file = make_mutation_sequences(single_variant_fasta, two_variant_fasta, args.kmer_size, args.output_name, args.multiple_mutations)
	kmer_frequency_files = query_sequences(wildtype_sequence_file, mutation_sequence_file, args.output_name, jellyfish_test, jellyfish_control)
	variant_info_dataframe = construct_variant_information_table(kmer_frequency_files, args.kmer_size, jellyfish_control, wildtype_sequence_file, mutation_sequence_file)
	variant_call_info_dataframe = hypothesis_test(variant_info_dataframe, args.kmer_size, args.poisson, jellyfish_control, args.output_name)
	create_variant_call_summary_table(variant_call_info_dataframe, sys.argv, kmer_frequency_files, args.output_name, original_dataframe, jellyfish_control)

def perform_assess(args):
	penultimate_dataframe = pd.read_csv(args.penultimate, sep='\t')
	validate_mutations(penultimate_dataframe, args.alpha, args.penultimate.split('_')[0], sys.argv)

def construct_variants_dataframe(args):
	# Construct a dataframe of the input vcf/bed file, retaining columns of the chromosome name,
	#	bed start and end positions, wildtype sequence bases, and mutation sequence bases
	if args.vcf_input: variants_dataframe = convert_to_bed(args.vcf_input)
	else: variants_dataframe = pd.read_csv(args.bed_input, sep='\t', comment='#', header=None)
	original_dataframe = variants_dataframe.copy()
	original_dataframe['Variant'] = original_dataframe.apply(lambda x: '_'.join(x[[0, 1, 2, 3, 4]].values.astype(str).tolist()), axis=1)
	original_dataframe.drop([0, 1, 2, 3, 4], inplace=True, axis=1)
	variants_dataframe = variants_dataframe.filter([0, 1, 2, 3, 4])
	variants_dataframe.rename(columns={0: 'chromosome', 1: 'start', 2: 'end', 3: 'wildtype', 4: 'mutation'}, inplace=True)
	original_dataframe.drop_duplicates(inplace=True, subset='Variant')
	return variants_dataframe, original_dataframe

if __name__ == '__main__':
	setup_logger()
	run_kmervc()
