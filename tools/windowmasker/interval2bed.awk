# NCBI dustmasker/windowmasker `-outfmt interval` -> BED3 (chrom, start, end).
#
# ⛔ THE INTERVAL FORMAT IS 0-BASED AND INCLUSIVE AT BOTH ENDS, and this file spent its whole life
# converting it as if it were 1-based. MEASURED against dustmasker's own `-outfmt fasta` on the
# same sequence, which is the only ground truth that settles it:
#
#     -outfmt interval      873 - 880          1500 - 1800
#     lowercase runs        [873, 881)         [1500, 1801)      <- what was actually masked
#     what this emitted     872   880          1499   1800       <- one base LEFT, right width
#
# `s = $1 - 1; print c, s, $2` preserves the WIDTH, which is why nothing caught it: every coverage
# percentage, every `bedtools genomecov` column of the masking table, and
# scripts/verify_softmask_outputs.py (which compares the union against the lowercase runs of a
# FASTA masked from that same union) all agree with themselves. Only the POSITION is wrong -- by
# one base, at every dustmasker and windowmasker interval in every strain. Downstream that
# lower-cases one non-repeat base and leaves one repeat base uppercase at each boundary, and it
# puts `lc_classify` out of phase: a `(AT)n` array was reported as `(TA)n`, a `(CAGGT)n` unit as
# `(GTCAG)n`, and a pure polyA run at purity 997 instead of 1000.
#
# At a sequence start the old clamp turned the shift into a TRUNCATION: dustmasker reporting
# `0 - 199` for a 200 bp leading repeat became `0 199`, losing one masked base outright. With the
# correct conversion the clamp has nothing to do -- `$1` is already >= 0 -- so it is gone.
#
# (content for cols 4-6 is added downstream by lc_classify.py)
BEGIN { OFS = "\t" }
/^>/  { c = substr($1, 2); next }
/[0-9]+ *- *[0-9]+/ { gsub(/-/, " "); print c, $1, $2 + 1 }
