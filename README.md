# dnalc

Local complexity estimation from DNA sequences

## Build

This program depends on [libdivsufsort](https://github.com/y-256/libdivsufsort) and
[SDSL](https://github.com/simongog/sdsl-lite).

To build dnalc, first obtain the source:

```
git clone https://github.com/EvolBioInf/dnalc.git
cd dnalc
```

If libdivsufsort is not installed, build it first:

```
make libdivsufsort
```

Simiarly, if SDSL is not installed, enter:

```
make sdsl
```

Finally, build dnalc:

```
make
```

The `dnalc` binary is now located in the `build` directory.

## The match complexity

dnalc implements a measure of complexity we call the match complexity,
MC. MC is a measure of non-repetitiveness of a given sequence. To
compute MC, a sequence gets factorized into consecutive unique
substrings, the match factors.  A substring is unique if the sequence
does not contain it at any other position. From left to right
characters are added to the current factor until it becomes a unique
string, then the factor is complete and the next one begins at the
next character. The process is continued until the whole sequence is
processed. For example, `CCCCGCTCTCCA` factorizes to
`CCC.CG.CTC.TC.C.A`. The factors are separated by dots, each factor
being unique with respect to the full sequence. The MC value for some
interval of the sequence (or the whole sequence) is obtained by
counting the factors within it and dividing the count by the expected
number for a random sequence. This expected number of match factors is
computed using equation (6) in **[1]**.

MC lies between 0 for sequences that are repeated exactly elsewhere
and 1 for random sequences. The definition of MC means that no value
can be less than 0. However, it is possible to obtain MC values
slighly larger than 1, as this maximum is only an expectation that
holds in the limit of long sequences.

## Usage

dnalc can be applied directly to FASTA files. If a file contains
multiple sequences, they are treated as one long sequence in the
factorization step.  The MC values are always calculated with respect
to the complete input sequence.

For example, to investigate the repetitiveness of chromosomes within a
genome, first create a FASTA file that contains the sequences of all
chromosomes. Note that keeping the chromosomes in separate files may
lead to different results, because possible repeats between
chromosomes cannot be detected. This feature of dnalc allows
repetitiveness to be investigated on various scales,
e. g. chromosome-wide vs. genome-wide - the "frame of reference" is
defined by input file.

If no parameters are specified, dnalc returns a single MC value for
the complete concatenated sequence. The `-n` parameter restricts dnalc
to specific regions. If the file contains multiple sequences `-n NAME`
selects the sequence with the given name. The name is a prefix of the FASTA
header of the sequence until either the first whitespace or 32 characters are reached.
It is also possible to restrict dnalc to a region of that sequence by using `-n
NAME:FROM-TO`.

**Examples**

```
# get MC value for the whole concatenated sequence:
dnalc seq.fa
# get MC value for the complete second sequence in the file:
dnalc -n chrZ seq.fa
# get MC value for the interval 12-345 of the second sequence in the file:
dnalc -n chrZ:12-345 seq.fa
```

### Batch modes
dnalc has two batch modes: *sliding window* or *list of queries*.
When given the sliding window size parameter `-w`, dnalc returns a
series of MC values. Alternatively, when given a file with intervals,
dnalc prints an MC value for each interval.

The `-w` parameter (optionally combined with `-n`) yields the values
of windows within the specified region. Using `-w` without `-n`
returns the window results for all the sequences in the file. The `-k`
parameter defines the step between windows and is 1/10 of the window
by default. In our experience, window lengths between 1000 and 100000
yield interesting results. Larger windows tend to average out the
changes in repetitiveness along the sequence. We recommend
experimenting with different window sizes. In the sliding window mode
the output columns are: *sequence number, window midpoint, MC value*.

Given a list of specific intervals of interest - rather than the a
contiguous sequence - the `-f` option specifies the file listing
these intervals, one per line, in the same format as for
`-n`. Note that `-f` is mutually exclusive with `-n` and `-w`. In the
`-f` mode the output columns are: *query number, MC value*.

**Examples:**

```
#produces a series of values for sliding windows of size 10000 and step 1000
dnalc seq.fa -w 10000
#produces a series of values for non-overlapping windows (k=w) of size 10000
dnalc seq.fa -w 10000 -k 10000
#produces a value for each query
dnalc seq.fa -f queries.txt
```

### FASTA vs index file
Raw sequence data can either be piped into dnalc, or passed as a
file. However, instead of raw sequence data, dnalc also accepts an
index to that sequence as input. Working with such a pre-computed
index can speed up the analysis of long sequences by orders of
magnitude. An index file is obtained by calling dnalc with the `-s`
flag and piping the output into a file: `dnalc seq.fa -s > seq.idx`.

Load an index file by using the `-i` flag; so if a file `seq.fa` was
transformed into the index `seq.idx`, use `dnalc -i seq.idx` instead
of `dnalc seq.fa`, everything else stays the same. In fact, whenever
dnalc is applied to a FASTA file, internally the index structure is
calculated on the fly and thrown away after the analysis.

To inspect an index file, use `dnalc -i someindex.idx -l`.  This
returns a list of all sequences indexed, in the same order as in the
input file. This also lists the possible arguments for the `-n`
parameter.

### Renaming
If you want to rename the sequences in the index (e.g. if the name deduced from
the FASTA header is not human readable), you can create a list of new names in a
text file, one name per line. Each name shall have at most 32 characters and be
unique. Then use the `-r` parameter to rename the sequences in an index. The
new names stored in the text file should be in the order the sequences appear in the
index (check with `-l`).

```
# rename sequences in existing index
dnalc -r new_names_file.txt some_seq.idx
```

### Gnuplot integration
The results of a sliding window analysis are best visualized. This can
be done using the script `dnalc_plot.sh`, which is part of the dnalc
repository. The script accepts dnalc output and visualizes it with
gnuplot. To use it, add the `-g` flag to dnalc and pipe the result
into the script:

```
dnalc -i seq.idx -n chrZ -w 10000 -g | ./dnalc_plot.sh
```

## References
**[1]** Estimating mutation distances from unaligned genomes.
Haubold, Pfaffelhuber, et al., Journal of Computational Biology,
Vol. 16, Number 10, 2009
