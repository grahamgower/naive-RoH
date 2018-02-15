#!/usr/bin/env python

from __future__ import print_function
import sys
import numpy as np

class Datum:
    def __init__(self, gen, ninds):
        self.gen = gen
        self.ninds = ninds
        self.segsites = None
        self.positions = None
        self.inds = []

def parse_ms(fn):
    data = {}
    genlist = iter([20,40,60,80,100])
    gen = None
    datum = None
    with open(fn) as f:
        for line in f:
            if line.startswith("//"):
                continue
            if line.startswith("#OUT:"):
                fields = line.split()
                gen = fields[1]
                if gen == "20:100":
                    gen = next(genlist)
                else:
                    gen = int(gen)
                ninds = int(fields[3])
                if datum is not None:
                    data[datum.gen] = datum
                datum = Datum(gen, ninds)
                continue
            if line.startswith("segsites:"):
                datum.segsites = int(line.split()[1])
                continue
            if line.startswith("positions:"):
                datum.positions = np.array(map(float, (line.split()[1:])))
                continue
            if datum is not None and datum.positions is not None:
                datum.inds.append(line.strip())
                continue
            

    if datum is not None:
        data[datum.gen] = datum

    return data

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: {} chr chrlen ms.txt".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    chrom = sys.argv[1]
    chrlen = int(sys.argv[2])
    ms_fn = sys.argv[3]

    win = 5*1000*1000
    step = 200*1000
    n_chunks = win / step

    data = parse_ms(ms_fn)

    #for gen in [20,40,60,80,100]:
    for gen in [40]:
        datum = data[gen]
        assert datum.ninds == len(datum.inds)/2, "mismatch for number of haplotypes for gen " + str(gen)
        for i in [0]: #range(0,datum.ninds):

            print("CHROM","START","END","N","HOM",sep="\t")

            hetlist = [0]*n_chunks

            g0 = datum.inds[i*2]
            g1 = datum.inds[i*2+1]
            chunk_start = win_start = 0
            window_het = 0
            hi = 0
            last_pos = -1
            for j, p in enumerate(datum.positions):
                pos = int(chrlen * p)
                if pos == last_pos:
                    continue
                last_pos = pos
                if pos > chunk_start + step:
                    window_het += hetlist[hi]
                    hi = (hi+1) % n_chunks

                    if pos > win_start + win:
                        win_end = min(win_start+win, chrlen)
                        n = win_end - win_start
                        print(chrom, win_start, win_end, n, n-window_het, sep="\t")
                        window_het -= hetlist[hi]
                        hetlist[hi] = 0
                        win_start += step
                    chunk_start += step

                if g0[j]!=g1[j]:
                    hetlist[hi] += 1

            # mop up
            if pos > win_start:
                if pos > chunk_start:
                    window_het += hetlist[hi]
                win_end = min(win_start+win, chrlen)
                n = win_end - win_start
                print(chrom, win_start, win_end, n, n-window_het, sep="\t")
