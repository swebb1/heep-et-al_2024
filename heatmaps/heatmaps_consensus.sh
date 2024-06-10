
computeMatrix-3.5.0 scale-regions -S bw_files/*.bw -R bed_files2/L1Md_A*bed bed_files2/L1Md_F*bed bed_files2/L1Md_anno/L1Md_T*bed -b 10 -o matrix_con.npz  --outFileNameMatrix matrix_con.tsv -m 5000 --upstream 2000 --downstream 2000 --sortRegions descend --smartLabels -p 24 --skipZeros

plotHeatmap -m matrix_con.npz -o consensus.heatmap.pdf --colorMap YlGnBu --whatToShow "plot, heatmap and colorbar" --plotFileFormat "pdf" -z "L1Md_A.O" "L1Md_A.Y" "L1Md_F.O" "L1Md_F.Y" "L1Md_T.O" "L1Md_T.Y" --heatmapWidth 6 --legendLocation upper-right -min 0 -x "Repetitive element" --startLabel "Start" --endLabel "End"
