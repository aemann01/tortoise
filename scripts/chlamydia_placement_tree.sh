######################
#Reference tree build
######################
mafft --reorder --auto ref.fa > ref.align.fa
trimal -in ref.align.fa -out ref.trimal.fa -gt 0.3 -st 0.001
rm *ref.tre
/Users/mann/raxmlHPC-AVX-v8/raxml -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n ref.tre -s ref.trimal.fa

######################
#Constraint tree build
######################
#align your (unaligned) sequences to your curated reference alignment (not the trimmed one)
echo "Aligning sequences to reference alignment"
rm ref.arb
sina -i ref.align.fa --prealigned -o ref.arb
sina -i query.fa -r ref.arb -o query.aligned.fa --fs-msc 0.01 --fs-full-len=10
cat ref.align.fa query.aligned.fa > queryPlusRef.align.fa

#fix names
sed 's/:.*//' queryPlusRef.align.fa | sed 's/ .*//' > queryPlusRef.align.fix.fa
#run constraint tree
rm *cons.tre
/Users/mann/raxmlHPC-AVX-v8/raxml -f a -N 100 -G 0.2 -m GTRCAT -n cons.tre -s queryPlusRef.align.fix.fa -g RAxML_bestTree.ref.tre -T 4 -x 25734 -p 25793
