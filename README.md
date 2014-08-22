# This is my README
# To run the test examepl run the following command in the terminal from with in the vartab folder
$vartav/$ perl VarTab.pl -f 100 -b test-data/test.funseq.bed -a test-data/test.rmsk.bed.gz,test-data/hg19.gwas.bed.gz,test-data/hg19.CpG.bed.gz < test-data/test.vcf

time perl bin/VarTab.pl -a -v test-data/targeted/cases.chr6.vcf -o test-data/targeted/cases.chr6 --frequency 5 --maf1kg 0.01
