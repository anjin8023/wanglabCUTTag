tssbed=/home/liukuai/Documents/reference/ucsc_mm10_refseq/ucsc.mm10.refseq.tss.bed

computeMatrix scale-regions -S peakcor/bigwig/KO-CCR6-8-21_FKDL202595257-1a.bw \
                               peakcor/bigwig/S_0527_KO_CCR6_H73WNDSX2_L1.bw \
                               peakcor/bigwig/WT-CCR6-8-21_FKDL202595254-1a.bw \
                               peakcor/bigwig/S_0527_WT_CCR6_H73WNDSX2_L1.bw \
                              -R $tssbed  \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o heatmap/ccr6_matrix_gene.mat.gz \
                              -p 4


plotHeatmap -m heatmap/ccr6_matrix_gene.mat.gz -out heatmap/ccr6.pdf --sortUsing sum --samplesLabel "KO_rep1" "KO_rep2" "WT_rep1" "WT_rep2" --plotTitle "CCR6"  &
plotHeatmap -m heatmap/NKp46_matrix_gene.mat.gz -out heatmap/NKp46.pdf --sortUsing sum --samplesLabel "KO_rep1" "KO_rep2" "WT_rep1" "WT_rep2" --plotTitle "NKp46"  &
plotHeatmap -m heatmap/DN_matrix_gene.mat.gz -out heatmap/dn.pdf --sortUsing sum --samplesLabel "KO_rep1" "KO_rep2" "WT_rep1" "WT_rep2" --plotTitle "DN" &




multiBigwigSummary bins \
 -b bigwig/KO-CCR6-8-21_FKDL202595257-1a.bw  bigwig/S_0527_KO_CCR6_H73WNDSX2_L1.bw bigwig/WT-CCR6-8-21_FKDL202595254-1a.bw bigwig/S_0527_WT_CCR6_H73WNDSX2_L1.bw \
 --labels  H3K4me3-CCR6-KO-rep1 H3K4me3-CCR6-KO-rep2 H3K4me3-CCR6-WT-rep1 H3K4me3-CCR6-WT-rep2 \
 -out H3K4me3-CCR6_scores_per_bin.npz --outRawCounts H3K4me3-CCR6_scores_per_bin.tab




plotCorrelation \
-in H3K4me3-CCR6_scores_per_bin.npz \
--corMethod pearson --skipZeros \
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







multiBigwigSummary bins \
 -b bigwig/KO-DN-8-21_FKDL202595259-1a.bw  bigwig/S_0527_KO_DN_H73WNDSX2_L1.bw  bigwig/WT-DN-8-21_FKDL202595256-1a.bw bigwig/S_0527_WT_DN_H73WNDSX2_L1.bw \
 --labels  H3K4me3-DN-KO-rep1 H3K4me3-DN-KO-rep2 H3K4me3-DN-WT-rep1 H3K4me3-DN-WT-rep2 \
 -out plot/H3K4me3-DN_scores_per_bin.npz --outRawCounts plot/H3K4me3-DN_scores_per_bin.tab




plotCorrelation \
-in H3K4me3-DN_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o dn/scatterplot_PearsonCorr_bigwigScores.png   \
--removeOutliers \
--outFileCorMatrix dn/PearsonCorr_bigwigScores.tab



plotCorrelation \
-in H3K4me3-DN_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "spearman Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o dn/scatterplot_spearmanCorr_bigwigScores.png   \
--removeOutliers \
--outFileCorMatrix dn/spearmanCorr_bigwigScores.tab






multiBigwigSummary bins \
 -b bigwig/KO-NCR-8-21_FKDL202595258-1a.bw  bigwig/S_0527_KO_NKp46_H73WNDSX2_L1.bw  bigwig/WT-NCR-8-21_FKDL202595255-1a.bw bigwig/S_0527_WT_NKp46_H73WNDSX2_L1.bw \
 --labels  H3K4me3-NKP46-KO-rep1 H3K4me3-NKP46-KO-rep2 H3K4me3-NKP46-WT-rep1 H3K4me3-NKP46-WT-rep2 \
 -out plot/H3K4me3-NKP46_scores_per_bin.npz --outRawCounts plot/H3K4me3-NKP46_scores_per_bin.tab


plotCorrelation \
-in H3K4me3-NKP46_scores_per_bin.npz \
--corMethod pearson --skipZeros \
--plotTitle "Pearson Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o nkp46/scatterplot_PearsonCorr_bigwigScores.png   \
--removeOutliers \
--outFileCorMatrix nkp46/PearsonCorr_bigwigScores.tab



plotCorrelation \
-in H3K4me3-NKP46_scores_per_bin.npz \
--corMethod spearman --skipZeros \
--plotTitle "spearman Correlation of Average Scores Per Transcript" \
--whatToPlot scatterplot \
-o nkp46/scatterplot_spearmanCorr_bigwigScores.png   \
--removeOutliers \
--outFileCorMatrix nkp46/spearmanCorr_bigwigScores.tab