
computeMatrix-3.5.0 scale-regions -S bw_files/*.bw -R transposon_annotation/bed_files/*filt.bed -b 10 -o matrix.npz  --outFileNameMatrix matrix.tsv -m 5000 --upstream 2000 --downstream 2000 --sortRegions descend --smartLabels -p 24 --skipZeros

plotHeatmap -m matrix.npz -o heatmap.pdf --colorMap YlGnBu --whatToShow "plot, heatmap and colorbar" --plotFileFormat "pdf" -z "IAPEy" "IAPEz" "L1Md_A" "L1Md_F" "L1Md_Gf" "L1Md_T" "MMERVK10C" --heatmapWidth 6 --legendLocation upper-right -min 0 -x "Repetitive element" --startLabel "Start" --endLabel "End"

