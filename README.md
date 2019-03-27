# consensus

Collapses multiple aligned sequences using iupac ambiguity bases. No base calling or anything. '-' and 'n' are ignored
unless found in every sequence in the alignment. Processes all fasta files in the input directory with .fasta or .fna
extension and outputs all the consensus sequences in one output fasta file.

## usage
```
usage: consensus.py [-h] [-i IN_DIRECTORY] [-o OUT_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -i IN_DIRECTORY, --in_directory IN_DIRECTORY
                        path to input directory
  -o OUT_FILE, --out_file OUT_FILE
                        path to output file
```
