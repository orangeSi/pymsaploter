#muscle_x64 -in genes.fa -out genes.muscle.clwstrict.aln -clwstrict
pre=out
python ../multi_alignment.ploter.py plot --align genes.muscle.clwstrict.aln --formats aln --types nuc --prefix $pre
convert  $pre.svg $pre.png
