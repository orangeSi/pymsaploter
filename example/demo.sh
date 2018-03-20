muscle_x64 -in genes.fa -out genes.muscle.clwstrict.aln -clwstrict
python ../multi_alignment.ploter.py genes.muscle.clwstrict.aln aln nuc
convert -density 100 out.multil.svg out.multil.png
