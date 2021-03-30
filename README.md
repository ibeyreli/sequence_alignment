# Sequence Alignment
Aim: Given two DNA sequences in a single FASTA-formatted file, obtain
both global and local alignments using naive and affine gap penalties
where the match score is 2, mismatch score is -3, and gap opening and 
extentions penaties are user inputs.

# allalign.py

allalign is a program that performs sequence alignment in 4 modes:
1. global: Needleman-Wunsch with naive gap scoring
2. local: Smith-Waterman with naiveve gap scoring
3. aglobal: Needleman-Wunsch with affine gap scoring
4. alocal: Smith-Waterman with affine gap scoring 

## Installation

allalign uses Python 3 standard libraries.
If you have installed Python 3 before, no additional package is needed.
If not, you may install Python 3 from "https://www.python.org/downloads/".
Make sure you have selected the appropriate distribution for your operating system.

## Usage

You can run the program using terminal.

The syntax:

"python allalign.py --mode global --input sequences.fa --gapopen Wg"
"python allalign.py --mode aglobal --input sequences.fa --gapopen Wg --gapext(optional) We"
"python allalign.py --mode local --input sequences.fa --gapopen Wg"
"python allalign.py --mode alocal --input sequences.fa --gapopen Wg --gapext(optional) We"

--mode				: Mode selection
					Appropriate selections are:
					- global for Needleman-Wunsch with naive gap scoring
					- local for Smith-Waterman with naiveve gap scoring
					- aglobal for Needleman-Wunsch with affine gap scoring
                    - alocal for Smith-Waterman with affine gap scoring 
					
--input			: Indicator that the next input is the input file which contains the sequences to be aligned	  
'sequences.fa' 	: File in the FASTA format, containing the text T in which the pattern will be searched
--gapopen		: Indicator that the next input is the gap opening penalty
Wg				: An integer denoting the gap opening penalty
--gapext		: (Optional)Indicator that the next input is the gap extention penalty for aglobal and alocal modes
We				: An integer denoting the gap extention penalty

## Examples:

hw3example.fa :
>my_first_sequence
TCGACCCAAGTAGGGAAAGAATATCAACACAAAGGCTCGAGAAGAGCCACC
CCATGAGCCACCGCATCTACCCCGTGCCCCAGCAAATTAAGAATAG
>another_sequence
TCGACCCATGTAGGGAAAGCATATCAATTTCACAAAGGCTCGAGAAGAGCC
ACATGAGCCACCGCATCTACCCCAGCAAATTAAGAAAAG

Input >>> python allalign.py --mode global --input hw3example.fa --gapopen -3

Output >>>
Score:120

>my_first_sequence TCGACCCAAGTAGGGAAAGAATATCAA---CACAAAGGCTCGAGAAGAGCCACCCCATGA
>another_sequence  TCGACCCATGTAGGGAAAGCATATCAATTTCACAAAGGCTCGAGAAGAGCCA---CATGA

>my_first_sequence GCCACCGCATCTACCCCGTGCCCCAGCAAATTAAGAATAG
>another_sequence  GCCACCGCATCTA-------CCCCAGCAAATTAAGAAAAG
