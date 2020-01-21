###########
#Chlamydia
###########
#reference tree based off of Mitura 2017 paper
grep "Chlamydia" ../assigntax/taxonomy.tsv | awk '{print $1}' | while read line; do grep -w $line ../rep_set.filt.fa -A 1; done > query_chlamydia.fa
cat ref_chlamydia.fa query_chlamydia.fa > all_chlamydia.fa
mafft --auto all_chlamydia.fa > all_chlamydia.align.fa
trimal -in all_chlamydia.align.fa -out all_chlamydia.trimal.fa -gt 0.99 -st 0.001 #trim to just the v4 region
raxmlHPC-PTHREADS-SSE3 -T 4 -m GTRCAT -c 25 -e 0.001 -p 31415 -f a -N 100 -x 02938 -n chlam.tre -s all_chlamydia.trimal.fa







