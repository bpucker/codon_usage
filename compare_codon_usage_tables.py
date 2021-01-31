### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.25 ###

__usage__ = """
	python compare_codon_usage_tables.py\n
	--in1 <INPUT_FILE1>
	--in2 <INPUT_FILE2>
	--out<FULL_PATH_TO_OUTPUT_FOLDER>
	
	optional:
	--name1 <NAME1>
	--name2 <NAME2>
	
	bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

from operator import itemgetter
import numpy as np
import sys, os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# --- end of imports --- #

def load_content( html_file ):
	"""! @brief load HTML content """
	
	with open( html_file, "r" ) as f:
		content = f.read()
	lines = content.split('<pre>')[1].split('</pre>')[0].split('\n')
	data = {}
	for line in lines:
		if len( line ) >0:
			blocks = line.split(')  ')
			for block in blocks:
				block = block.replace( "  ", " " ).replace( "  ", " " ).replace( "  ", " " ).replace( "  ", " " )
				parts = block.split(' ')
				data.update( { parts[1] + "_" + parts[0]: { 'codon': parts[0], 'aa': parts[1], 'freq': float( parts[2] ) } } )
	return data


def comparison( data1, data2, figfile, name1, name2 ):
	"""! @brief compare both datasets  """
	
	
	fig, ax = plt.subplots( figsize=(8, 3) )
	
	x_values = []
	y_values = []
	labels = []
	for idx, key in enumerate( sorted( data1.keys() ) ):
		x_values.append( idx+1 )
		y_values.append( data1[ key ]['freq'] )
		labels.append( key )
	ax.bar( x_values, y_values, width=0.2, color="blue", tick_label=labels )
	
	x_values = []
	y_values = []
	labels = []
	for idx, key in enumerate( sorted( data2.keys() ) ):
		x_values.append( idx+1.3 )
		y_values.append( data2[ key ]['freq'] )
		labels.append( key )
	ax.bar( x_values, y_values, width=0.2, color="red", tick_label=labels )
	
	
	my_legend = [ 	mpatches.Patch( color='blue', label=name1 ),
								mpatches.Patch( color='red', label=name2 )
							]
	
	ax.legend( handles=my_legend, bbox_to_anchor=( 0.5, 1.05 ), loc="center", ncol=2 )
	
	ax.set_xlim( 0, len( x_values ) + 1 )
	ax.set_ylabel( "codon usage" )
	
	ax.tick_params(axis='both', which='major', labelsize=7, rotation=90)
	ax.tick_params(axis='both', which='minor', labelsize=7, rotation=90)
	
	plt.subplots_adjust( left=0.05, right=0.99, top=0.9, bottom=0.15 )
	
	fig.savefig( figfile, dpi=300 )


def main( arguments ):
	"""! @brief runs all parts of this script """
	
	input_file1 = arguments[ arguments.index( '--in1' )+1 ]
	input_file2 = arguments[ arguments.index( '--in2' )+1 ]
	output_folder = arguments[ arguments.index( '--out' )+1 ]
	
	if not os.path.exists( output_folder ):
		os.makedirs( output_folder )
	
	if '--name1' in arguments:
		name1 = arguments[ arguments.index( '--name1' )+1 ]
	else:
		name1 = "name1"
	
	if '--name2' in arguments:
		name2 = arguments[ arguments.index( '--name2' )+1 ]
	else:
		name2 = "name2"
	
	content1 = load_content( input_file1 )
	content2 = load_content( input_file2 )
	
	figfile = output_folder + name1 + "_vs_" + name2 + "comparison.pdf"
	comparison( content1, content2, figfile, name1, name2 )


if '--in1' in sys.argv and '--out' in sys.argv and '--in2' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
