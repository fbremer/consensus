# consensus

`consensus.py` collapses multiple aligned sequences using iupac ambiguity bases, with the exception of `n` discussed below. No base calling or cutoffs or anything.

`-` and `n` are treated specially and similarly. By default both use an *all* rule, where they are ignored unless all sequences contain the symbol at that location. They each have a flag to change to an *any* rule, where if any of the symbol is found at a location, the consensus will take that symbol. If `--any_gap` and `--any_n` are set, gaps take precedence.

By default, `consensus.py` takes one fasta file and outputs one fasta file. In batch mode, it processes all files in the input directory with a `.fasta` or `.fna` extension as fasta's and outputs all the consensus sequences in one output fasta file.

## usage
```
usage: consensus.py [-h] [-i INPUT] [-o OUT_FILE] [--batch_mode] [--any_n]
                    [--any_gap]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        path to input file, or input directory if using batch
                        mode
  -o OUT_FILE, --out_file OUT_FILE
                        path to output file
  --batch_mode          take an input directory instead of input file
  --any_n               By default, if *all* bases at a location are 'n', the
                        consensus will be a 'n'. This flag changes that rule
                        to return an 'n' if *any* 'n' are found at that loc.
  --any_gap             Same as --any_n above. If --any_gap and --any_n are
                        set, gaps take precedence
```
