# HGG multi-omics robustness audit

How robust are the conclusions of Wang et al. (2019) to the choice of data mining method? We re-analyse the
proteome, phosphoproteome and RNA-seq data from PDGFRA- vs NTRK1-driven high-grade glioma (mouse), comparing
the authors' methods with alternatives. Each claim ends with a verdict: robust, method-dependent, not
supported, or reproduced / not reproduced.

**Start here: [TEAM_GUIDE.md](TEAM_GUIDE.md)** (repository map, the paper's methods and claims, claim-by-claim
comparison, how to run things, open tasks).

**Credit:** Wang H. et al., *Nat Commun* 10:3718 (2019), doi:10.1038/s41467-019-11661-4. The authors' code is in
`JUMPn/` and `JUMPp/` (unmodified; github.com/hongwang198745/HGG_Source_Code). Their README is kept as
`README_original.md`.

| Folder or file | What it is |
|---|---|
| `analysis/C1_...` to `C11_...` | one folder per claim, each with a five-line `README.md` (Kept / Compared / Added / Result / Verdict) |
| `analysis/00_sample_key/` | sample and TMT channel mapping |
| `LOG.md` | pre-registrations, results and decisions, in order |
| `results/metrics.csv` | all quantitative results (claim, variant, metric, value, n) |
| `CLAUDE.md` | working rules, metric defaults and the known-issues list |
| `env/` | R session information |
| `data/external/` | not in git; the paper's supplementary files (see the guide, section 8) |

**Status:** all eleven claims analysed. Open: verdict table, figure fixes, report and slides
(guide, section 10).