# Naive-RoH
Calculate homozygosity in overlapping windows, from a vcf/bcf file.  There
exist some very clever ways to identify runs of homozygosity (e.g.
[using a hidden Markov model](https://samtools.github.io/bcftools/howtos/roh-calling.html)),
however the code here takes a naive approach---it merely counts homozygous
and heterozygous sites in the window based on the genotype call (FORMAT/GT
field).

Two components are provided here, (1) a C program, `hom_windows`, which prints
counts of homozygous sites within windows, and (2) a script
`plot_hom_windows.py` that produces a pdf figure of the homozygosity along a chromosome.

# Prerequisites
`hom_windows` uses **htslib** to parse vcf/bcf files.  The plotting script
requires **python** (tested with version **2.7.14** and **3.6.1**) and
**matplotlib** (tested with version **2.1.0**).

# Installation
Clone the git repository, then build with `make`.

# Usage
The recommended way to use `hom_windows` is to first split your vcf by
chromosome and it run separately on each chromosome file.  E.g. to split the
file `infile.vcf.gz` by chromosome (for an organism with 26 autosomes),
then run `hom_windows` and plot the result:
```
for c in $(seq 26) X; do
	chr="chr$c"
	bcftools view -O z -o ${chr}.vcf.gz -r ${chr} infile.vcf.gz
	hom_windows ${chr}.vcf.gz > hom_windows.${chr}.txt
	plot_hom_windows.py --title $chr -o $chr.pdf hom_windows.${chr}.txt
done
```

Additional control over the window size and the step size for moving along
chromosomes can be obtained via command line parameters to `hom_windows`:
```
hom_windows v1
usage: ./hom_windows [...] file.vcf

  -s INT         Move window along chromosomes in steps of INT bp [200000].
  -w INT         Output windows of size INT bp [5000000].
  -S STR[,...]   For a multi-sample vcf, specify the sample to use [].
                 Multiple samples may be specified, separated with a comma,
                 in which case loci not segregating among the samples are
                 counted as homozygous.
  -h INT[,...]   Ignore sites with depth higher than INT [1000].
                 If multiple samples are specified, comma separated max depths
                 must be specified for each sample.
  -l INT[,...]   Ignore sites with depth lower than INT [0].
                 If multiple samples are specified, comma separated min depths
                 must be specified for each sample.
```

The homozygosity from multiple individuals may be included in the same figure
by specifying multiple input files and labels for the figure legend.
```
plot_hom_windows.py \
	--title "Homozygosity along chr1" \
	-o chr1.pdf \
	-l "Individual 1,Individual 2 (inbred line)" \
	hom_windows.ind1.chr1.txt \
	hom_windows.ind2-inbred.chr1.txt
```

Further control over the scale and aspect ratio of the figure can obtained via
command line parameters to the script.
```
usage: plot_hom_windows.py [-h] [--wide] [--scale SCALE] [--title TITLE]
                           [-o OPDF] [-l LABELS]
                           infiles [infiles ...]

plot `hom_windows' output

positional arguments:
  infiles               input file

optional arguments:
  -h, --help            show this help message and exit
  --wide                plot widescreen ratio (16x9) [False]
  --scale SCALE         scale the plot [1.0]
  --title TITLE         plot title
  -o OPDF, --opdf OPDF  output filename [out.pdf]
  -l LABELS, --labels LABELS
                        comma separated list of labels to correspond with
                        input files
```
