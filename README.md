# consensus

Collapses multiple aligned sequences using iupac ambiguity bases. No base calling or anything. `-` and `n` are ignored unless found in every sequence in the alignment. By default, takes one fasta file and outputs one fasta file. In batch mode, processes all fasta files in the input directory with `.fasta` or `.fna` extension and outputs all the consensus sequences in one output fasta file.

## usage
```
usage: consensus.py [-h] [-i INPUT] [-o OUT_FILE] [--batch_mode]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input file, or input directory if using batch
                        mode
  -o OUT_FILE, --out_file OUT_FILE
                        path to output file
  --batch_mode          take an input directory instead of input file

```
