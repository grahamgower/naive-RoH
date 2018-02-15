#!/usr/bin/env python

from __future__ import print_function
import sys
import os.path
import itertools

import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

import numpy as np

def parse_hombed(fn, skip=1, missing_threshold=1.0/3.0):
    with open(fn) as f:
        while skip > 0:
            # skip header
            next(f)
            skip -= 1
        for line in f:
            fields = line.split()
            chrom = fields[0]
            start, end, n, hom = map(int, fields[1:])
            if n < missing_threshold*(end-start):
                continue
            yield chrom, (start+end)/2, float(hom)/n

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="plot `hom_windows' output")
    parser.add_argument("--wide", action="store_true", default=False, help="plot widescreen ratio (16x9) [%(default)s]")
    parser.add_argument("--scale", type=float, default=1.0, help="scale the plot [%(default)s]")
    parser.add_argument("--title", type=str, help="plot title")
    parser.add_argument("-o", "--opdf", type=str, default="out.pdf", help="output filename [%(default)s]")
    parser.add_argument("-l", "--labels", type=str, help="comma separated list of labels to correspond with input files")
    parser.add_argument("infiles", nargs="+", help="input files")
    args = parser.parse_args()

    if args.labels:
        args.labels = args.labels.split(",")
        if len(args.labels) != len(args.infiles):
            print("Number of labels does not match the number of input files.".format(), file=sys.stderr)
            exit(1)
    else:
        args.labels = [os.path.bastname(fn) for fn in args.infiles]

    return args


if __name__ == "__main__":

    args = parse_args()

    pdf = PdfPages(args.opdf)
    if args.wide:
        fig_w, fig_h = plt.figaspect(9.0/16.0)
    else:
        fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(args.scale*fig_w, args.scale*fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    y_min = 1.0
    x_min = float('inf')
    x_max = 0.0

    # attempt to distinguish each of the plots in the figure
    colours = itertools.cycle(plt.get_cmap("tab10").colors)
    lstyles = itertools.cycle([':', '--', '-'])
    lwidths = itertools.cycle([2,2,1])

    for infile, label, col, ls, lw in zip(args.infiles, args.labels, colours, lstyles, lwidths):
        _, x, y = zip(*list(parse_hombed(infile)))

        ax1.plot(x, y, color=col, linestyle=ls, lw=lw, label=label)

        y_min = min(min(y), y_min)
        x_min = min(min(x), x_min)
        x_max = max(max(x), x_max)

    xticks = np.arange(0, 3e8, 2e7)
    xlabels = [str(int(xt/1e6))+" mbp" for xt in xticks]
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xlabels, size="small")

    ax1.set_xlim(x_min*0.85, x_max*1.01)
    ax1.set_ylim(y_min*.999, 1)

    ax1.set_ylabel("Homozygosity")
    if args.title:
        ax1.set_title(args.title)

    ax1.legend(loc="lower right")

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
