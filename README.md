# dnalc

Experimental local neighborhood complexity estimation for DNA,
based on match-length factorization (derived from Lempel-Ziv-factorization).

## Build

This program depends on [libdivsufsort](https://github.com/y-256/libdivsufsort).
If you want to work with very long sequences, you also need
[SDSL](https://github.com/simongog/sdsl-lite).

The build is straightforward. First obtain the source:

```
git clone https://github.com/EvolBioInf/dnalc.git
cd dnalc
```

If you do not have libdivsufsort installed, build it first:

```
make libdivsufsort
```

If you want to use the space-efficient SDSL structures (you most certainly do!):

```
make sdsl
```

Now build the tool itself:

```
make
```

If you don't plan to work with very long sequences (> 2^32), you can disable
SDSL in the Makefile or by invoking make like this: `make USE_SDSL=0`. Then a
faster 32-bit version will be built.

After the build is complete you can find the `dnalc` binary in the `build` directory.

## What is match-length complexity (MC) ?

MC is a measure of non-repetitiveness of a given sequence. First a sequence gets
factorized into consecutive unique strings, the match factors.
A substring is unique iff the sequence does not contain it at any other
position. From left to right characters are added to the current factor until
it becomes a unique string, then the factor is complete and the next one begins
at the next character. The process is continued until the whole sequence is
processed. For example, `CCCCGCTCTCCA` will be factorized to
`CCC.CG.CTC.TC.C.A`. The factors are separated by dots and as you can easily
verify that each factor is unique in the sequence.

The MC value for some interval of the sequence (or the whole sequence) is
obtained by counting the number of such factors within it and dividing it
by the expected number for a random sequence (the calculation of the expected
value is quite sophisticated, see **[1]** for details).

Usually real sequences are more repetitive than completely random ones, so that
there are less, but longer factors, therefore the quotient is smaller,
usually the value will be between 0 and 1. It is possible (but quite rare) to
obtain a larger value, hinting at an especially complex (i.e. non-repetitive)
region.

## Usage
dnalc can be used directly on FASTA files. If a file contains multiple
sequences, they are treated as one long sequence for the factorization step.
The MC values are always depending on the underlying sequence.

For example, to investigate the repetitiveness of chromosomes within a genome
you would create a FASTA file that contains the sequences of all chromosomes.
If you have the chromosomes in different files, you will get different results,
because non-local repeats (e.g. same sequence on different chromosomes) will
not be identified. Keep this in mind, as you might wish to look from different
scales (e.g. genome-wide or chromosome-wide) - the file is the "frame of
reference" for repetitivity.

If you do not specify any parameters, you will get a single MC value for
the complete concatenated sequence. With the `-n` parameter you can restrict
dnalc to a specific region - if the file contains multiple sequences,
with `-n NUM` you select the sequence. If you want to calculate for a part
of that sequence, you can specify it as `-n NUM:FROM-TO`.

**Examples**

```
# get MC value for the whole concatenated sequence:
dnalc seq.fa
# get MC value for the complete second sequence in the file:
dnalc -n 2 seq.fa
# get MC value for the interval 12-345 of second sequence in the file:
dnalc -n 2:12-345 seq.fa
```

### Batch modes
dnalc has two batch modes: *sliding window* or *list of queries*.
When you provide the sliding window size parameter `-w`, dnalc will return a
series of MC values. If you provide a file with queries, if will output the
MC values for the requested intervals.

If you set the `-w` parameter (optionally combined with `-n`), you will get the
values of windows within the specified region. If you use `-w` without
specifying anything with `-n`, you will get the window results for all the
sequences inside the file. The `-k` parameter defines the step between windows
and is 1/10 of the window as default setting. For realistic sequences a window
between 1000 and 100000 yields good results. Larger windows smooth the change of
repetitiveness along the sequence. You should experiment with different window
sizes. In the sliding window mode the output columns are: *sequence number,
offset, MC value*.

If you are not interested in the complexity landscape of your sequence, but
instead have a list of specific regions of interest, you can use the `-f`
parameter to specify the name of the file which contains a list of queries, one
per line, in the same format as used for `-n`. Note that `-f` is mutually
exclusive with `-n` and `-w`. In this mode the output columns are: *query
number, MC value*.

**Examples:**

```
#produces a series of values for windows of size 10000 and step 1000
dnalc seq.fa -w 10000
#produces a series of values for non-overlapping windows (k=w) of size 10000
dnalc seq.fa -w 10000 -k 10000
#produces a value for each query
dnalc seq.fa -f queries.txt
```

### FASTA vs index file
You can either pipe a sequence into dnalc, or pass it a FASTA file as an
argument. For reasonably small sequences (up to the size of small chromosomes)
this may be fine, but for serious usage you are strongly encouraged to create an
index file. The MC value calculation is very fast while the preparation step is
slow and always the same for a given file, so it is a good idea to do that
preparation just once. An index file is obtained by calling dnalc with the `-s`
flag and pipe the output into a file, like this: `dnalc seq.fa -s > seq.idx`.

To tell dnalc that the input is an index file you have to use the `-i` flag, so
if you had a file `seq.fa` and created the index `seq.idx`, you will need to
call `dnalc -i seq.idx` instead of `dnalc seq.fa`, everything else stays the
same. In fact, if you use dnalc on a FASTA file, internally the index
structure is calculated every time and thrown away afterwards.

You can get information about an index file with `dnalc -i someindex.idx -l`.
You will get a list of all contained sequences, in the same order like in the
FASTA file. Here you can also see the number you need to pass with the `-n`
parameter.

### Gnuplot integration
If you are using the sliding window mode, you probably want to plot the result.
For convenience, this repository provides the script `dnalc_plot.sh`, which can
accept dnalc output and visualize it with gnuplot. If you wish to use it, add
the `-g` flag to dnalc and pipe the result into the script, like this:

```
dnalc -i seq.idx -n 1 -w 10000 -g | ./dnalc_plot.sh
```

## References
**[1]** Estimating mutation distances from unaligned genomes.
Haubold, Pfaffelhuber, et al., Journal of Computational Biology, Vol. 16, Number 10, 2009
