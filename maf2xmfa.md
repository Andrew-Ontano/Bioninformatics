# maf2xmfa tool
This tool will convert an maf file of a single sequence/chromosome (produced from cactus-hal2maf from a minigraph-cactus
chromosome sub-problem) and output a mauve-viewer and Geneious compatible xmfa.

## Usage
`maf2xmfa.py [-h] -i INPUT_MAF -o OUTPUT_XMFA [-v] [-g]`<br>
*maf2xmfa: a tool for converting maf to xmfa.*<br>
*optional arguments:*<br>
  *-h, --help*<br>Output of these arguments.

  *-i INPUT_MAF, --input INPUT_MAF*<br>Path to input multiple alignment formatted file.

  *-o OUTPUT_XMFA, --output OUTPUT_XMFA*<br>Path to output xmfa file.

  *-v, --verbose*<br>Output progress to the terminal.

  *-g, --gappy*<br>Allows for gappy alignment blocks to be generated in the xmfa. Gappy alignments have headers with "$i:0-0 +
  $Sequence.$Chrom", with a sequence of "N" the length of the alignment where a sequence was absent from the alignment.
