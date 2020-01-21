######################
#MaxL tree build
######################
cat ref.fa queryAll.fa > full.fa
mafft --reorder --auto full.fa > full.align.fa
trimal -in full.align.fa -out full.trimal.fa -gt 0.3 -st 0.001
~/raxmlHPC-AVX-v8/raxml -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n tre -s full.trimal.fa

# reference tree based off of silva living tree project SSU v132
sed "s/'/\n/g" chlamydia.tre | grep -v ":" | grep -v "(" > ref.ids
sed -i 's/_.*//' ref.ids
seqtk subseq ~/refDB/silva_livingTree/LTPs132_SSU_aligned.fasta ref.ids > ref.fa
grep "Chlamydiae" ../sequence_taxonomy_table.16s.noUnassigned.txt | awk '{print $1}' | while read line; do grep -w $line ../rep_set.filt.fa -A 1; done > query.fa
# add in reads from Mitura 2017 paper
cat query.fa mitura2017.fa > queryAll.fa

######################
#EPA tree build
######################
#align your (unaligned) sequences to your curated reference alignment (not the trimmed one)
echo "Aligning sequences to reference alignment"
sina -i ref.fa --prealigned -o ref.arb
sina -i queryAll.fa -r ref.arb -o query.aligned.fa --fs-msc 0.01 --fs-full-len=10
cat ref.fa query.aligned.fa > queryPlusRef.align.fa
#positions with all gaps and those that do not cover the amplicon deleted manually in aliview
#clean up
sed -i 's/ .*//' queryPlusRef.trim.fa

#full maximum likelihood tree
/Users/mann/raxmlHPC-AVX-v8/raxml -f a -# 100 -m GTRCAT -n chlam.tre -s queryPlusRef.trim.fa -T 4 -p 12345 -x 12345

#build epa tree
echo "EPA tree building"
rm *epa.tre
/Users/mann/raxmlHPC-AVX-v8/raxml -f v -G 0.2 -m GTRCAT -n epa.tre -s queryPlusRef.trim.fa -t chlamydia.tre -T 4

#clean up tree so you can read into figtree
sed 's/QUERY___//g' RAxML_labelledTree.epa.tre | sed 's/\[I[0-9]*\]//g' > RAxML_placementTree.epa.tre

#######################
#Constraint tree build
#######################
rm *cons.tre
/Users/mann/raxmlHPC-AVX-v8/raxml -f a -N 100 -G 0.2 -m GTRCAT -n cons.tre -s queryPlusRef.trim.fa -g chlamydia.tre -T 4 -x 25734 -p 25793
