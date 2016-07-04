# This is a short-read genome assembler for Python 3.X
The actual code that can be run is a few directories down in /BMES544FinalProject/BMES544FinalProject/src/ and requires you have Biopython, PYaml, and numpy installed. The whole directory structure is a little messed up since this was written in VS Community and on a Linux cluster on different occasions. Make sure both the BMES544FinalProject.py script and the ASM\_CONFIG.yaml file are present in the same folder when running.

## Usage
There are two flags that can be used for this assembler:

> -o will specify how much of a nucleotide overlap you want to have between reads

> -g will specify the mode of operation that will prioritize contig length over contig number

> -h will specify the mode of operation that will prioritize contig number over contig length

## Required libraries:

Python: [Python Software Foundation](https://www.python.org)

Biopython: [Biopython](http://biopython.org/wiki/Download)

PYaml - pip install pyaml
