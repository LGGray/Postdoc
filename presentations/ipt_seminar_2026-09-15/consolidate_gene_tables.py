# Consolidate per-nucleus Allelome.PRO2 chrX locus tables (Allelome.PRO2_all_genes)
# into one long table: sample, cell_barcode, gene, start, end, A1_reads, A2_reads, total_reads
import os, sys, time
base = "/Users/graylachlan/cluster/OCM/Allelome.PRO2_all_genes"
out = sys.argv[1]
t0 = time.time()
with open(out, "w") as fo:
    fo.write("sample\tcell_barcode\tchr\tstart\tend\tgene\tA1_reads\tA2_reads\ttotal_reads\n")
    for s in ["9w", "78w", "Sham", "TAC"]:
        dirs = sorted(os.listdir(os.path.join(base, s)))
        n = 0
        for d in dirs:
            p = os.path.join(base, s, d, "locus_table.txt")
            if not os.path.exists(p):
                continue
            # dir name: <date>_<sample>_<barcode>_annotation_us_mm39_chrX.bed_1
            parts = d.split("_")
            bc = parts[4]
            with open(p) as f:
                next(f)
                for line in f:
                    x = line.rstrip("\n").split("\t")
                    fo.write(f"{s}\t{s}_{bc}\t{x[0]}\t{x[1]}\t{x[2]}\t{x[3]}\t{x[4]}\t{x[5]}\t{x[6]}\n")
            n += 1
            if n % 500 == 0:
                print(s, n, round(time.time() - t0), "s", flush=True)
        print("done", s, n, "cells", round(time.time() - t0), "s", flush=True)
print("FINISHED", round(time.time() - t0), "s")
