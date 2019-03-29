# strict_iupac_consensus

`consensus.py` collapses multiple aligned sequences into strict consensus sequences using IUPAC ambiguity bases. By strict, we mean that it does not use any conventions for determining whether an IUPAC code is introduced (e.g. minimum required occurrences of minor bases, such as in [Cavener 1987](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/15.4.1353), which is used by [Bio.motifs](http://biopython.org/DIST/docs/tutorial/Tutorial.html) in BioPython). If >1 base is present at a given site, `consensus.py` will introduce an ambiguity. By default, it will also ignore `-` and `n` characters if other bases are present at a given position (see below).  

By default, `consensus.py` takes one fasta file and outputs one fasta file. In batch mode, it processes all files in the input directory with a `.fasta` or `.fna` extension as fasta's and outputs all the consensus sequences in one output fasta file.

### `-` and `n` characters
`-` and `n` are treated specially and similarly. By default both use an **all** rule, where they are ignored unless all sequences contain the symbol at that location. They each have a flag to change to an **any** rule, where if any of the symbol is found at a location, the consensus will take that symbol. If both `--any_gap` and `--any_n` are set, gaps take precedence.


### usage
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


### example
Input: 
```
>ind1
NCNCNT
>ind2
GCNCTC
>ind3
GCNTT-
>ind4
GCATTT
```
Outputs:
```
>consensus (default options)
GCAYTY
>consensus (with --any_n)
NCNYNY
>consensus (with --any_gap)
GCAYT-
```

### credits
 * [Forest Bremer](https://github.com/Woods26): author
