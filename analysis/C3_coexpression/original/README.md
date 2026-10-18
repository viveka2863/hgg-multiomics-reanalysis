# C3 co-expression: original (control arm)

- **Kept:** Authors' wgcna_v1.0.R, paper settings (signed network, softPower 16, minModuleSize 100, deepSplit 2, merge height 0.2, raw intensities, 4,703 pre-selected DE proteins), run in scratch/JUMPnc with only 3 edits (workingDir, sizeGrWindow commented out, shipped outputs renamed).
- **Compared:** TBD
- **Added:** TBD
- **Result:** 4,703 rows in input, shipped and rerun. Rerun vs shipped ARI = 0.9993 (2 proteins differ, a blue/yellow swap), cluster sizes identical (turquoise 2341, blue 1592, yellow 307, green 356, black 107; shipped blue 1592/yellow 307, rerun same). Paper labels (WP-C1..5) vs rerun ARI = 0.9507 (same vs shipped); WP-C1 and WP-C2 each split across blue/yellow (60 and 22 proteins crossing), WP-C3 -> green (337/376) with 32 in blue, WP-C4 -> black (96/104), WP-C5 -> turquoise (2339/2344). Reproduced (>= 0.95) vs shipped output; paper-label agreement is just above the threshold.
- **Verdict:** TBD

Files: run_control.sh (exact steps), compare.R / compare.log, crosstab_*.csv, ari.csv, run.log.
