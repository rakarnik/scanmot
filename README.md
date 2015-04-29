*scanmot* and *scantop* are programs that scan sequences for instances of motifs.

## scanmot
Finds the highest scoring instances of the motif in the input sequences.

Usage:
`scanmot <ace_file> <motif_number> <sequence_file> <output_file> <max_sites>`
- ace_file: AlignACE output file with motifs
- motif_number:	Motif number within ace_file to use
- sequence_file: File with sequences to scan in FASTA format
- output_file: Output file with motif instances
- max_sites: Maximum number of sites to output

## scantop
Finds the top-scoring instance(s) of the motif in each sequence.

Usage:
`scantop <ace_file> <motif_number> <sequence_file> <background_file> <site_count> <output_file>`
- ace_file: AlignACE output file with motifs
- motif_number:	Motif number within ace_file to use
- sequence_file: File with sequences to scan (FASTA format)
- background_file: File with sequences used to build background model (FASTA format)
- site_count: Number of top scoring instances to output per sequence
- output_file: Output file with motif instances


