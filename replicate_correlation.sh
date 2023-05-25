multiBigwigSummary bins \
 -b bigwig/KO-CCR6-8-21_FKDL202595257-1a.bw  bigwig/S_0527_KO_CCR6_H73WNDSX2_L1.bw bigwig/WT-CCR6-8-21_FKDL202595254-1a.bw bigwig/S_0527_WT_CCR6_H73WNDSX2_L1.bw \
 --labels  H3K4me3-CCR6-KO-rep1 H3K4me3-CCR6-KO-rep2 H3K4me3-CCR6-WT-rep1 H3K4me3-CCR6-WT-rep2 \
 -out plot/H3K4me3-CCR6_scores_per_bin.npz --outRawCounts plot/H3K4me3-CCR6_scores_per_bin.tab

plotCorrelation \
-in H3K4me3-CCR6_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o scatterplot_PearsonCorr_bigwigScores.png   \
--removeOutliers \
--outFileCorMatrix PearsonCorr_bigwigScores.tab


plotCorrelation \
-in H3K4me3-CCR6_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "spearman Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o scatterplot_spearmanCorr_bigwigScores.png   \
--removeOutliers \
--outFileCorMatrix spearmanCorr_bigwigScores.tab

