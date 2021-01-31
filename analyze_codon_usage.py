### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.25 ###

__usage__ = """
	python analyze_codon_usage.py\n
	--in <FULL_PATH_TO_INPUT_FILE>
	--out<FULL_PATH_TO_OUTPUT_FILE>
	
	optional:
	--exp <FULL_PATH_TO_EXPRESSION_FILE>
	--name <NAME>
	--top100 <ANALYZE_ONLY_TOP100>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

from operator import itemgetter
import numpy as np
import sys

# --- end of imports --- #


def load_genetic_code():
	"""! @brief return standard genetic code """
	
	genetic_code = {	'CTT': 'L',
									'ATG': 'M',
									'AAG': 'K',
									'AAA': 'K',
									'ATC': 'I',
									'AAC': 'N',
									'ATA': 'I',
									'AGG': 'R',
									'CCT': 'P',
									'ACT': 'T',
									'AGC': 'S',
									'ACA': 'T',
									'AGA': 'R',
									'CAT': 'H',
									'AAT': 'N',
									'ATT': 'I',
									'CTG': 'L',
									'CTA': 'L',
									'CTC': 'L',
									'CAC': 'H',
									'ACG': 'T',
									'CCG': 'P',
									'AGT': 'S',
									'CAG': 'Q',
									'CAA': 'Q',
									'CCC': 'P',
									'TAG': '*',
									'TAT': 'Y',
									'GGT': 'G',
									'TGT': 'C',
									'CGA': 'R',
									'CCA': 'P',
									'TCT': 'S',
									'GAT': 'D',
									'CGG': 'R',
									'TTT': 'F',
									'TGC': 'C',
									'GGG': 'G',
									'TGA': '*',
									'GGA': 'G',
									'TGG': 'W',
									'GGC': 'G',
									'TAC': 'Y',
									'GAG': 'E',
									'TCG': 'S',
									'TTA': 'L',
									'GAC': 'D',
									'TCC': 'S',
									'GAA': 'E',
									'TCA': 'S',
									'GCA': 'A',
									'GTA': 'V',
									'GCC': 'A',
									'GTC': 'V',
									'GCG': 'A',
									'GTG': 'V',
									'TTC': 'F',
									'GTT': 'V',
									'GCT': 'A',
									'ACC': 'T',
									'TTG': 'L',
									'CGT': 'R',
									'TAA': '*',
									'CGC': 'R'
								}
	return genetic_code


def load_all_seqs_from_multiple_fasta_file( filename ):
	"""! @brief load all sequences from multiple fasta file """
	
	data = {}
	
	with open( filename, "r" ) as f:
	 	header = f.readline().strip()[1:].split(' ')[0]
		line = f.readline()
		seq = []
		while line:
			if line[0] == '>':
				data.update( { header: "".join( seq ) } )
				header = line.strip()[1:].split(' ')[0]
				seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		data.update( { header: "".join( seq ) } )
	return data


def analyze_codon_usage( CDS, codons, expression ):
	"""! @brief analyze codon usage in representative CDS """
	
	codon_usage = {}
	exp_codon_usage = {}
	for codon in codons:
		codon_usage.update( { codon: 0 } )
		exp_codon_usage.update( { codon: 0 } )
	
	# --- run this analysis if expression data is not available --- #
	if expression == {}:
		for key in CDS.keys():
			seq = CDS[ key ].upper()
			if len( seq ) % 3 == 0:
				if not 'N' in seq:
					blocks = [ seq[ i:i+3 ] for i in range( 0, len( seq ), 3 ) ]
					for block in blocks:
						codon_usage[ block ] += 1
				else:
					print "ERROR: N in coding sequence!"
			else:
				print "ERROR: CDS length is not multiple of 3!"
	
	# --- run this analysis if expression data is available --- #
	else:
		for key in CDS.keys():
			seq = CDS[ key ].upper()
			if len( seq ) % 3 == 0:
				if not 'N' in seq:
					blocks = [ seq[ i:i+3 ] for i in range( 0, len( seq ), 3 ) ]
					for block in blocks:
						exp_codon_usage[ block ] += expression[ key ]
						codon_usage[ block ] += 1
				else:
					print "ERROR: N in coding sequence!"
			else:
				print "ERROR: CDS length is not multiple of 3!"
	return codon_usage, exp_codon_usage


def save_codon_usage( codon_usage, codon_usage_exp, genetic_code, output_file ):
	"""! @brief write codon usage results into output file """
	
	# --- invert genetic code dictionary --- #
	codons_per_aa = {}
	for key in genetic_code.keys():
		try:
			codons_per_aa[ genetic_code[ key ] ].append( key )
		except KeyError:
			codons_per_aa.update( { genetic_code[ key ]: [ key ] } )
	
	# --- calculate relative codon usage without expression data --- #
	if codon_usage_exp == {}:
		all_data = []
		for aa in codons_per_aa:
			counts_per_codon = []
			for codon in codons_per_aa[ aa ]:
				counts_per_codon.append( codon_usage[ codon ] )
			for codon in codons_per_aa[ aa ]:
				try:
					all_data.append( { 'aa': aa, 'codon': codon, 'counts': codon_usage[ codon ],
													'frequency': codon_usage[ codon ] / float( sum( counts_per_codon ) ),
													'global_freq': codon_usage[ codon ] / float( sum( codon_usage.values() ) )
													} )
				except ZeroDivisionError:
					all_data.append( { 'aa': aa, 'codon': codon, 'counts': codon_usage[ codon ], 'frequency': "n/a", 'global_freq': "n/a" } )
	
	# --- calculate relative codon usage with expression data --- #
	else:
		all_data = []
		for aa in codons_per_aa:
			counts_per_codon = []
			exp_counts_per_codon = []
			for codon in codons_per_aa[ aa ]:
				counts_per_codon.append( codon_usage[ codon ] )
				exp_counts_per_codon.append( codon_usage_exp[ codon ] )
			for codon in codons_per_aa[ aa ]:
				try:
					all_data.append( { 'aa': aa, 'codon': codon, 'counts': codon_usage[ codon ],
													'frequency': codon_usage[ codon ] / float( sum( counts_per_codon ) ),
													'global_freq': codon_usage[ codon ] / float( sum( codon_usage.values() ) ),
													'exp': codon_usage_exp[codon]  / float( sum( exp_counts_per_codon ) ),
													'exp_global': codon_usage_exp[ codon ] / float( sum( codon_usage_exp.values() ) )
												} )
				except ZeroDivisionError:
					all_data.append( { 'aa': aa, 'codon': codon, 'counts': codon_usage[ codon ],
													'frequency': "n/a",
													'global_freq': "n/a",
													'exp': "n/a",
													'exp_global': "n/a"
												} )
	
	# --- write output file --- #
	with open( output_file, "w" ) as out:
		if codon_usage_exp == {}:
			out.write( 'AminoAcid\tCodon\tCounts\tFrequency\tGlobalFrequency\n' )
			for entry in sorted( all_data, key=itemgetter( 'aa', 'codon' ) ):
				out.write( "\t".join( map( str, [ entry['aa'], entry['codon'], entry['counts'], entry['frequency'], entry['global_freq'] ] ) ) + '\n' )
		else:
			out.write( 'AminoAcid\tCodon\tCounts\tFrequency\tGlobalFrequency\tExpressionBased\tExpressionBasedGlobalFrequency\n' )
			for entry in sorted( all_data, key=itemgetter( 'aa', 'codon' ) ):
				out.write( "\t".join( map( str, [ entry['aa'], entry['codon'], entry['counts'], entry['frequency'], entry['global_freq'], entry['exp'], entry['exp_global'] ] ) ) + '\n' )
	return all_data


def load_gene_expression( expression_file ):
	"""! @brief load gene expression values per gene """
	
	expression = {}
	with open( expression_file, "r" ) as f:
		f.readline()
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			expression.update( { parts[0]: np.median( map( float, parts[1:] ) ) } )
			line = f.readline()
	return expression


def construct_codon_string( entry ):
	"""! @brief construct string of one codon """
	
	try:
		cblock = [ 	entry['codon'],
							entry['aa'],
							str( round( entry['exp'], 2 ) ) + " " * (4-len( str( round( entry['exp'], 2 ) ) ) ),	#complete string to three characters
							" " * ( 5 - len( str( round( entry['exp_global']*1000, 1 ) ) ) ) + str( round( entry['exp_global']*1000, 1 ) ),
							"(" + " "*( 6-len( str( entry['counts'] )  ) ) +  str( entry['counts'] ) + ")"
						]
	except KeyError:
		cblock = [ 	entry['codon'],
							entry['aa'],
							str( round( entry['frequency'], 2 ) ) + " " * (4-len( str( round( entry['exp'], 2 ) ) )),	#complete string to three characters
							str( round( entry['global_freq']*1000, 1 ) ) + " " * ( 5 - len( str( round( entry['global_freq']*1000, 1 ) ) ) ),
							"(" + " "*( 6-len( str( entry['counts'] )  ) ) + str( entry['counts'] ) + ")"
						]
	return " ".join( cblock )


def GC_calculations( CDS ):
	"""! GC content calculations """
	
	nt1 = []
	nt2 = []
	nt3 = []
	for key in CDS.keys():
		seq = CDS[ key ]
		if len( seq ) % 3 == 0:
			if not 'N' in seq:
				blocks = [ seq[ i:i+3 ] for i in range( 0, len( seq ), 3 ) ]
				for block in blocks:
					nt1.append( block[0] )
					nt2.append( block[1] )
					nt3.append( block[2] )
	
	gc1 = 100.0 * ( nt1.count( "G" ) + nt1.count( "C" ) ) / len( nt1 )
	gc2 = 100.0 * ( nt2.count( "G" ) + nt2.count( "C" ) ) / len( nt2 )
	gc3 = 100.0 * ( nt3.count( "G" ) + nt3.count( "C" ) ) / len( nt3 )
	gc = ( gc1+gc2+gc3 ) / 3.0
	return gc, gc1, gc2, gc3


def construct_10D_table( table_file, codon_usage, number_of_CDS, name, gc, gc1, gc2, gc3 ):
	"""! @brief construct 10D output format """
	
	total_codons = 0
	codon_dict = {}
	for entry in codon_usage:
		total_codons += entry['counts']
		codon_dict.update( { entry['codon']: entry } )
	
	codon_orders = [ ]
	for nt1 in [ "T", "C", "A", "G" ]:
		for nt2 in [ "T", "C", "A", "G" ]:
			for nt3 in [ "T", "C", "A", "G" ]:
				codon_orders.append( nt1+nt3+nt2 )
	codon_blocks = [ codon_orders[i:i + 4] for i in xrange(0, len( codon_orders ), 4) ]
	#blocks of four codons = one line
	
	output_lines = []
	for idx, block in enumerate( codon_blocks ):
		if idx > 0 and idx % 4 == 0:
			output_lines.append( "" )	#inserts an empty line after four codon lines
		new_line = []
		for codon in block:
			new_line.append( construct_codon_string( codon_dict[ codon ] ) )
		output_lines.append( "  ".join( new_line ) )
	
	prefix1 = '<html><head>\n<meta http-equiv="Content-Type" content="text/html; charset=windows-1252">\n<title>Codon usage table</title>\n</head>\n<body bgcolor="#F0F0F0">\n'
	prefix2 = '<strong><i>Cyanidioschyzon merolae strain 10D </i>[gbpln]: ' + str( number_of_CDS) + ' CDS\'s (' + str( total_codons ) + ' codons)</strong>\n\n'
	prefix3 = '<hr size="1" align="LEFT">\nfields: [triplet] [amino acid] [fraction] [frequency: <strong>per thousand</strong>] ([number])\n<hr size="1" align="LEFT">\n<pre>'
	suffix1 = '</pre>\n<hr size="1" align="LEFT">\nCoding GC ' + str( gc )+ '%\n1st letter GC ' + str( gc1 )+ '%\n2nd letter GC ' + str( gc2 )+ '%\n3rd letter GC ' + str( gc3 )+ '%<br>\n\n<strong>'
	suffix2 = 'Genetic code 1: Standard</strong>\n<hr size="1" align="LEFT">\n\n<strong>Format:</strong><br>\n</body>\n</html>'
	
	with open( table_file, "w" ) as out:
		out.write( prefix1 + prefix2 + prefix3 + "\n".join( output_lines ) + "\n" + suffix1 + suffix2 )


def reduce_to_top100( CDS, expression ):
	"""! @brief select the top100 transcripts """
	
	exp_per_transcript = []
	for key in expression.keys():
		exp_per_transcript.append( { 'ID': key, 'val': expression[ key ] } )
	
	top100 = sorted( exp_per_transcript, key=itemgetter('val') )[::-1][:100]
	
	top_cds = {}
	for each in top100:
		try:
			top_cds.update( { each['ID']: CDS[ each['ID'] ] } )
		except KeyError:
			pass
	return top_cds


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	input_cds_file = arguments[ arguments.index( '--in' )+1 ]
	output_file = arguments[ arguments.index( '--out' )+1 ]
	
	top100 = False
	expression_status = False
	if '--exp' in arguments:
		expression_file = arguments[ arguments.index( '--exp' )+1 ]
		expression_status = True
		if '--top100' in arguments:
			top100 = True
	
	if '--name' in arguments:
		name = arguments[ arguments.index( '--exp' )+1 ]
	else:
		name = "xxx"
	
	# --- get genetic code --- #
	genetic_code = load_genetic_code( )	#add genetic code as dictionary to avoid additional file import
	
	# --- get sequences --- #
	CDS = load_all_seqs_from_multiple_fasta_file( input_cds_file )
	
	# --- get expression --- #
	if expression_status:
		expression = load_gene_expression( expression_file )
		if top100:
			CDS = reduce_to_top100( CDS, expression )
		error_counter = 0
		for ID in CDS.keys():
			try:
				expression[ ID ]
			except:
				error_counter += 1
				expression.update( { ID: 0 } )
		print "number of genes missing in expression data set: " + str( error_counter )
		codon_usage, codon_usage_exp = analyze_codon_usage( CDS, genetic_code.keys(), expression )
		codon_usage_info = save_codon_usage( codon_usage, codon_usage_exp, genetic_code, output_file )
		
	else:
		codon_usage, codon_usage_exp = analyze_codon_usage( CDS, genetic_code.keys(), {} )
		codon_usage_info = save_codon_usage( codon_usage, {}, genetic_code, output_file )
	
	# --- generate 10D output format file --- #
	output_10D_table_file = output_file + ".10D.html"
	gc, gc1, gc2, gc3 = GC_calculations( CDS )
	construct_10D_table( output_10D_table_file, codon_usage_info, len( CDS.keys() ), name, gc, gc1, gc2, gc3 )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
