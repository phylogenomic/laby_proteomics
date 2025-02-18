## Columns in spreadsheets

RAW.T0-1: raw values for the proteomics dataset.

BATCH_CORRECT.T0-1: values after batch correction using SVA in Combat.	

T0_1: log2 values after VSN normalization (DEP).

T2_vs_T0_ratio: Average across T2 replicates - Average across T0 replicates. Also, log2(T2)-log2(T0), or log2(T2/T0).

T0_centered: Average across T0 replicates - Average across all timepoints and replicates. log2(T0)/log2(total.average).

T2_vs_T0_p.val: p.value calculated from limma (DEP)

T2_vs_T0_p.adj: adjusted p.value

T2_vs_T0_significant

T2_vs_T0_diff: similar to ratio, approximation by limma

T2_vs_T0_CI.L: Left CI of diff.

T2_vs_T0_CI.R: Right CI of diff.

FC_T2_vs_T0: 2^(ratio)

T2.reg: UP/DOWN/NO. If significant and above/below the threshold.
