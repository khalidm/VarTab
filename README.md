# This is my README
# To run the test examepl run the following command in the terminal from with in the vartab folder
$vartav/$ perl VarTab.pl -f 100 -b test-data/test.funseq.bed -a test-data/test.rmsk.bed.gz,test-data/hg19.gwas.bed.gz,test-data/hg19.CpG.bed.gz < test-data/test.vcf

perl bin/VarTab.pl -f 100 -b test-data/test.funseq.bed -a test-data/hg19.rmsk.bed.gz,test-data/hg19.CpG.bed.gz,test-data/clinvar_20140303.bed.gz,test-data/gwascatalog.bed.gz,test-data/wgEncodeRegTfbsClusteredV2.cell_count.20130213.bed.gz,test-data/stam.125cells.dnaseI.hg19.bed.gz,test-data/hg19.polyphen.bed.gz,test-data/1kg.phase1.snp.bed.gz -m 0.01 -v test-data/test.vcf -o test-data/test -n

time perl bin/VarTab.pl -f 100 -a test-data/hg19.rmsk.bed.gz,test-data/hg19.CpG.bed.gz,test-data/clinvar_20140303.bed.gz,test-data/gwascatalog.bed.gz,test-data/wgEncodeRegTfbsClusteredV2.cell_count.20130213.bed.gz,test-data/stam.125cells.dnaseI.hg19.bed.gz,test-data/hg19.polyphen.bed.gz,test-data/1kg.phase1.snp.bed.gz,data/hg19.gerp.elements.bed.gz,data/hg19.ncrna_pgene.bed.gz,data/hg19.funseq.nc_sensitive.bed.gz,data/hg19.cadd_ph15.bed.gz,data/gencode.v7.promoter.sort.bed.gz -m 0.01 -v test-data/test.vcf -o test-data/test -n -k 0.0000005
