/*
 * Copyright (c) 2017 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <ctype.h>

#include <htslib/vcf.h>

#define HOM_WINDOWS_VERSION "1"

// This should be in htslib, but I couldn't find it.
static inline uint32_t bcf_hdr_id2ctglen(const bcf_hdr_t *hdr, int rid) { return hdr->id[BCF_DT_CTG][rid].val->info[0]; }

#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

typedef struct {
	int *min_depth;
	int *max_depth;
	int window;
	int step;
	char *samples;
	int n_samples;
} opt_t;

typedef struct block {
	uint32_t hom;
	uint32_t n;
} block_t;

int
hom_windows(opt_t *opt, char *vcf_fn)
{
	htsFile *fp;
	bcf_hdr_t *hdr;
	bcf1_t *rec;
	int ret;

	int n_dp = 0;
	int *dp = NULL;
	int ngt;
	int ngt_arr = 0;
	int32_t *gt_arr = NULL;

	fp = hts_open(vcf_fn, "r");
	if (fp == NULL) {
		ret = -1;
		goto err0;
	}

	hdr = bcf_hdr_read(fp);
	if (hdr == NULL) {
		ret = -2;
		goto err1;
	}

	if (bcf_hdr_nsamples(hdr) != 1) {
		int sret;

		if (!opt->samples) {
			fprintf(stderr, "%s: multisample bcf/vcf, but no samples specified with `-S'.\n", vcf_fn);
			ret = -3;
			goto err1;
		}

		if ((sret = bcf_hdr_set_samples(hdr, opt->samples, 0)) != 0) {
			if (sret < 0) {
				fprintf(stderr, "%s: couldn't process samples `%s'\n", vcf_fn, opt->samples);
				ret = -4;
				goto err1;
			} else if (sret > 0) {
				// One of the samples isn't in the vcf. Get a pointer to it.
				int i, j;
				char *s = NULL;
				for (i=0, j=1; i<strlen(opt->samples); i++) {
					if (opt->samples[i] == ',') {
						opt->samples[i] = '\0';
						j++;
					}
					if (sret == j)
						s = opt->samples + i;
				}

				if (s != NULL)
					fprintf(stderr, "%s: couldn't find sample `%s' in bcf/vcf\n", vcf_fn, s);
				else
					fprintf(stderr, "%s: one of the samples was not found, bcf_hdr_set_samples() returned %d\n", vcf_fn, sret);
				ret = -5;
				goto err1;
			}
		}
	}


	rec = bcf_init1();
	if (rec == NULL) {
		perror("bcf_init1: calloc");
		ret = -6;
		goto err2;
	}

	int32_t last_chrom = -1, win_start = 0, chunk_start = 0;
	int n_chunks = opt->window / opt->step;
	block_t *blocklist;
	block_t *b = blocklist;
	int bi = 0, bcount = 0;
	int window_n = 0, window_hom = 0;

	b = blocklist = calloc(n_chunks, sizeof(*blocklist));
	if (blocklist == NULL) {
		ret = -7;
		goto err3;
	}

	printf("CHROM\tSTART\tEND\tN\tHOM\n");

	while (bcf_read(fp, hdr, rec) >= 0) {

		// ignore all but REF or SNP calls
		if (!bcf_is_snp(rec) || bcf_get_info(hdr, rec, "INDEL"))
				break;

		if (bcf_get_format_int32(hdr, rec, "DP", &dp, &n_dp) < 0) {
			fprintf(stderr, "%s: missing FORMAT/DP field at %s:%d\n",
					vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
			ret = -8;
			goto err4;
		}

		int i;
		int skip = 0;
		for (i=0; i<n_dp; i++) {
			if (dp[i] < opt->min_depth[i] || dp[i] > opt->max_depth[i]) {
				skip = 1;
				break;
			}
		}
		if (skip)
			continue;

		if ((ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr)) < 1) {
                        fprintf(stderr, "%s: missing FORMAT/GT field at %s:%d\n",
                                        vcf_fn, bcf_seqname(hdr, rec), rec->pos+1);
                        ret = -9;
                        goto err4;
                }

		int j;
		int gt0 = -1;
		int hom = 1;
		int max_ploidy = ngt / opt->n_samples;

		for (i=0; i<opt->n_samples; i++) {
			int32_t *ptr = gt_arr + i*max_ploidy;
			int gt;
			for (j=0; j<max_ploidy; j++) {
				// if true, the sample has smaller ploidy
				if (ptr[j] == bcf_int32_vector_end)
					break;
				// missing allele
				if (bcf_gt_is_missing(ptr[j])) {
					skip = 1;
					break;
				}
				gt = bcf_gt_allele(ptr[j]);
				if (gt0 == -1)
					gt0 = gt;
				else if (gt != gt0) {
					hom = 0;
					break;
				}
			}
			if (skip)
				break;
		}
		if (skip)
			continue;

		if (rec->rid != last_chrom || rec->pos > chunk_start+opt->step) {
			window_n += b->n;
			window_hom += b->hom;

			if (rec->rid != last_chrom) {
				if (last_chrom != -1) {
					printf("%s\t%d\t%d\t%d\t%d\n",
							bcf_hdr_id2name(hdr, last_chrom),
							win_start,
							min(win_start+opt->window, bcf_hdr_id2ctglen(hdr,last_chrom)),
							window_n,
							window_hom);
				}
				bi = 0;
				b = blocklist;
				bcount = 0;
				window_n = 0;
				window_hom = 0;
				chunk_start = 0;
				win_start = 0;
				memset(blocklist, 0, n_chunks*sizeof(*blocklist));
				last_chrom = rec->rid;
				continue;

			}

			do {
				bi = (bi+1) % n_chunks;
				b = blocklist+bi;
				bcount += (bcount<n_chunks?1:0);
				
				if (rec->pos > win_start+opt->window) {
					printf("%s\t%d\t%d\t%d\t%d\n",
							bcf_hdr_id2name(hdr, last_chrom),
							win_start,
							min(win_start+opt->window, bcf_hdr_id2ctglen(hdr,last_chrom)),
							window_n,
							window_hom);

					window_n -= b->n;
					window_hom -= b->hom;
					assert(window_n >= 0);
					assert(window_hom >= 0);
					b->n = 0;
					b->hom = 0;
					win_start += opt->step;
				}
			} while (rec->rid == last_chrom && rec->pos > win_start+opt->window);

			chunk_start += opt->step;
		}

		b->n++;
		b->hom += hom;
	}

	// mop up
	if (rec->pos > win_start) {
		if (rec->pos > chunk_start) {
			window_n += b->n;
			window_hom += b->hom;
		}
		printf("%s\t%d\t%d\t%d\t%d\n",
			bcf_hdr_id2name(hdr, last_chrom),
			win_start,
			min(win_start+opt->window, bcf_hdr_id2ctglen(hdr,last_chrom)),
			window_n,
			window_hom);
	}

	ret = 0;


err4:
	if (dp)
		free(dp);
	if (gt_arr)
		free(gt_arr);
	free(blocklist);
err3:
	bcf_destroy1(rec);
err2:
	bcf_hdr_destroy(hdr);
err1:
	hts_close(fp);
err0:
	return ret;
}

unsigned long
parse_bp(char *s)
{
	unsigned long x = 0;

again:
	if (strlen(s) <= 1)
		return -1;

	switch (tolower(s[strlen(s)-1])) {
		case 'k':
			x = 1000;
			break;
		case 'm':
			x = 1000000;
			break;
		case 'b':
			if (x == -1)
				return -1;
			x = -1;
			s[strlen(s)-1] = 0;
			goto again;
		default:
			x = 1;
			break;
	}

	return x*strtoul(s, NULL, 0);
}

void
usage(opt_t *opt, char *argv0)
{
	fprintf(stderr, "hom_windows v%s\n", HOM_WINDOWS_VERSION);
	fprintf(stderr, "usage: %s [...] file.vcf\n", argv0);
	fprintf(stderr, "\n");
	fprintf(stderr, "  -s INT         Move window along chromosomes in steps of INT bp [%d].\n", opt->step);
	fprintf(stderr, "  -w INT         Output windows of size INT bp [%d].\n", opt->window);
	fprintf(stderr, "  -S STR[,...]   For a multi-sample vcf, specify the sample to use [].\n");
	fprintf(stderr, "                 Multiple samples may be specified, separated with a comma,\n");
	fprintf(stderr, "                 in which case loci not segregating among the samples are\n");
	fprintf(stderr, "                 counted as homozygous.\n");
	fprintf(stderr, "  -h INT[,...]   Ignore sites with depth higher than INT [%d].\n", 1000);
	fprintf(stderr, "                 If multiple samples are specified, comma separated max depths\n");
	fprintf(stderr, "                 must be specified for each sample.\n");
	fprintf(stderr, "  -l INT[,...]   Ignore sites with depth lower than INT [%d].\n", 0);
	fprintf(stderr, "                 If multiple samples are specified, comma separated min depths\n");
	fprintf(stderr, "                 must be specified for each sample.\n");
	exit(1);
}

int
main(int argc, char **argv)
{
	int i, c;
	int len_min_depth = 0, len_max_depth = 0;
	char *vcf_fn;
	opt_t opt;

	memset(&opt, 0, sizeof(opt));

	opt.window = 5*1000*1000;
	opt.step = 200*1000;
	opt.samples = NULL;
	opt.n_samples = 1;

	while ((c = getopt(argc, argv, "h:l:s:w:S:")) != -1) {
		switch (c) {
			case 'h':
				{
					unsigned long x;
					int n = 1;
					int *mem = NULL;
					char *argnext = optarg;
					do {
						if (*argnext == ',') {
							argnext++;
							n++;
						}
						x = strtoul(argnext, &argnext, 0);
						if (x < 0 || x > 10000) {
							fprintf(stderr, "max_depth[%d]=`%s' out of range\n", n, optarg);
							return -1;
						}
						mem = realloc(opt.max_depth, sizeof(*opt.max_depth)*n);
						if (mem == NULL) {
							perror("realloc");
							exit(1);
						}
						opt.max_depth = mem;
						opt.max_depth[n-1] = x;

					} while (*argnext != '\0');
					len_max_depth = n;
				}
				break;
			case 'l':
				{
					unsigned long x;
					int n = 1;
					int *mem = NULL;
					char *argnext = optarg;
					do {
						if (*argnext == ',') {
							argnext++;
							n++;
						}
						x = strtoul(argnext, &argnext, 0);
						if (x < 0 || x > 10000) {
							fprintf(stderr, "min_depth[%d]=`%s' out of range\n", n, optarg);
							return -1;
						}
						mem = realloc(opt.min_depth, sizeof(*opt.min_depth)*n);
						if (mem == NULL) {
							perror("realloc");
							exit(1);
						}
						opt.min_depth = mem;
						opt.min_depth[n-1] = x;

					} while (*argnext != '\0');
					len_min_depth = n;
				}
				break;
			case 's':
				opt.step = parse_bp(optarg);
				if (opt.step < 0 || opt.step > 100*1000*1000) {
					fprintf(stderr, "step=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 'w':
				opt.window = parse_bp(optarg);
				if (opt.window < 0 || opt.window > 100*1000*1000) {
					fprintf(stderr, "window=`%s' out of range\n", optarg);
					return -1;
				}
				break;
			case 'S':
				{
					opt.samples = optarg;
					for (i=0; i<strlen(opt.samples); i++) {
						if (opt.samples[i] == ',')
							opt.n_samples++;
					}
				}
				break;
			default:
				usage(&opt, argv[0]);
		}
	}

	if (argc-optind != 1) {
		usage(&opt, argv[0]);
	}

	if (opt.step > opt.window) {
		fprintf(stderr, "Error: step (%d) > window (%d) doesn't make sense.",
				opt.step, opt.window);
		return -1;
	}

	if (opt.window % opt.step != 0) {
		fprintf(stderr, "Error: step (%d) must divide window (%d) with no remainder.",
				opt.step, opt.window);
		return -1;
	}

	if (opt.min_depth == NULL) {
		opt.min_depth = calloc(opt.n_samples, sizeof(*opt.min_depth));
		if (opt.min_depth == NULL) {
			perror("calloc");
			exit(1);
		}
		for (i=0; i<opt.n_samples; i++) {
			opt.min_depth[i] = 0;
		}
	}
	if (opt.max_depth == NULL) {
		opt.max_depth = calloc(opt.n_samples, sizeof(*opt.max_depth));
		if (opt.max_depth == NULL) {
			perror("calloc");
			exit(1);
		}
		for (i=0; i<opt.n_samples; i++) {
			opt.max_depth[i] = 1000;
		}
	}

	if (len_min_depth != len_max_depth || len_min_depth != opt.n_samples) {
		fprintf(stderr, "Error: must have: number of samples (-S) == length of min_depth (-l) == length of max_depth (-h)\n");
		return -1;
	}

	vcf_fn = argv[optind];

	if (hom_windows(&opt, vcf_fn) < 0) {
		return -1;
	}

	free(opt.min_depth);
	free(opt.max_depth);

	return 0;
}

