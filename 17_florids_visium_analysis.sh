### VISIUM ANALYSIS USING FLORIDA Transcriptome: 

## From florida_trans_lcm_reads_mapping_2021.sh
### Transcriptome is : 
## 
cd /home/jjm247/project/florida_transcriptome_dir/Assembled-transcriptomes/
florida=/home/jjm247/project/florida_transcriptome_dir/Assembled-transcriptomes/Trinity.fasta.reduced.cdscdhit_v2
### get longest iso
perl ~/programs/trinityrnaseq-v2.10.0/util/misc/get_longest_isoform_seq_per_trinity_gene.pl $florida > Trinity.fasta.reduced.cdscdhit_v2.longest.fasta
#### clean headers
sed 's/ .*//g' Trinity.fasta.reduced.cdscdhit_v2.longest.fasta > Trinity.fasta.reduced.cdscdhit_v2.longest.clenaheaders.fasta

######### original script

## make the exon gft file from transdecoder output: 
grep exon Trinity.fasta.transdecoder.gff3 > Trinity.fasta.transdecoder_exon.gff3
cat Trinity.fasta.transdecoder_exon.gff3 | sed 's/ID=/gene_id "/g' |sed 's/;Parent=/"; transcript_id "/g' | sed 's/$/"/g' > Trinity.fasta.transdecoder_exon.gtf
sed 's/\.p.*;/";/g' Trinity.fasta.transdecoder_exon.gtf > Trinity.fasta.transdecoder_exon2.gtf

### remove duplicates
cut -f1 Trinity.fasta.transdecoder_exon2.gtf | sort | uniq | while read line; do grep -m1 "$line" Trinity.fasta.transdecoder_exon2.gtf; done > Trinity.fasta.transdecoder_exon2_uniq.gtf
#### reduce gtf to longest cleanheadrs. 
grep '>' Trinity.fasta.reduced.cdscdhit_v2.longest.clenaheaders.fasta | sed 's/>//g' | while read line; do \
grep -m1 "$line." Trinity.fasta.transdecoder_exon2.gtf; \
done > Trinity.fasta.transdecoder_exon2_reduced.gtf
#############END OF ORIGINAL SCRIPT
################ TEST ONE (DIDN'T WORK)
awk 'BEGIN {FS="\t";OFS="\t"} {$7="+" ;print ;}' Trinity.fasta.transdecoder_exon2_reduced.gtf > Trinity.fasta.transdecoder_exon2_reduced_both_senses.gtf
awk 'BEGIN {FS="\t";OFS="\t"} {$7="-" ;print ;}' Trinity.fasta.transdecoder_exon2_reduced.gtf >> Trinity.fasta.transdecoder_exon2_reduced_both_senses.gtf
cat  Trinity.fasta.transdecoder_exon2_reduced_both_senses.gtf | tr ' ' '\t' > Trinity.fasta.transdecoder_exon2_reduced_both_senses2.gtf
############### ################ TEST ONE (DIDN'T WORK)

##### 
##### TEST 2 MAP TO CCM ONE CONTIG ONLY ########
### map only to the most highly expressed DE genes: 
# from fig2_table_CCM_DE, copies contigs names and will reduce the Trinity.fasta

cat ccm_target_genes.list | while read line; do grep "$line " Trinity.fasta -A1; done > Trinity.fasta.CCM.genes.fas
sed -i 's/ .*//g' Trinity.fasta.CCM.genes.fas
### GENE INFO FILE 
cat ccm_target_genes.list | while read line; do \
grep "$line.p" Trinity.fasta.transdecoder_exon2.gtf; \
done > Trinity.fasta.CCM.genes.gtf
###
INITY_DN11597_c0_g1_i2.p1

##### check how many has sense and anti sense ORFs
cut -f1,7 Trinity.fasta.CCM.genes.gtf | sort | uniq | cut -f1 | uniq -c | sort -nr | grep '2 '
## check from the lis of both senses, the ones that match the enzyme under study
grep TRINITY_DN25070_c1_g2_i1 trinotate_annotation_report.xls
### From the contigs with more than one CDS
Remove: 
#  INITY_DN25070_c1_g2_i1.p1 has two pyrophosphatase 
INITY_DN773_c0_g1_i11.p2
TRINITY_DN11597_c0_g1_i2.p1
INITY_DN7121_c0_g1_i4.p2
INITY_DN6907_c0_g1_i11.p2
INITY_DN25070_c1_g2_i1.p1
INITY_DN2073_c0_g1_i1.p1
TRINITY_DN1997_c3_g1_i3
INITY_DN11597_c0_g1_i2.p1
INITY_DN1108_c0_g1_i13.p2

## delete each: 
cat ccm_isoforms_remove.list | while read line;do sed -i "/$line/d" Trinity.fasta.CCM.genes.gtf; done
### check in the LCM and in rennatas
#### NEW SCRIPT TO INTRODUCE ALL THE SENSES::
### remove those duplicates with the same sense: 
cut -f1 Trinity.fasta.CCM.genes.gtf | sort | uniq -c | sort -nr | grep '2 '

TRINITY_DN91948_c0_g1_i1
TRINITY_DN3383_c0_g1_i3
TRINITY_DN282_c3_g1_i2
TRINITY_DN2102_c0_g1_i10
TRINITY_DN1747_c2_g1_i1
TRINITY_DN1529_c0_g1_i11
TRINITY_DN105_c0_g1_i1

cat remove_duplicates.list | while read line; do \
isoform=$(grep -m1 "$line" Trinity.fasta.CCM.genes.gtf | sed 's/.* //g'); \
sed -i "/$isoform/d" Trinity.fasta.CCM.genes.gtf; done


#### make reference: 
########
cd /home/jjm247/project/visium_may21
module load spaceranger
#rm -r florida_visium_transcriptome
spaceranger mkref --genome=florida_visium_ccm_genes_reference \
--fasta=/home/jjm247/project/florida_transcriptome_dir/Assembled-transcriptomes/Trinity.fasta.CCM.genes.fas \
--genes=/home/jjm247/project/florida_transcriptome_dir/Assembled-transcriptomes/Trinity.fasta.CCM.genes.gtf

##################################
mv florida_visium_ccm_genes_reference florida_visium_transcriptome
## the transcriptome is in 
/home/jjm247/project/visium_may21/florida_visium_transcriptome
## reads 
~/scratch60/visium_reads

### change scripts: 
imag='\/home\/jjm247\/project\/visium\/manually_aligned_images'
trans='\/home\/jjm247\/project\/visium_may21\/florida_visium_transcriptome \\'
for i in *sh; do \
sed "s/^images.*/images=$imag/" $i | sed "s/20april_folder/12may_dir/g" | sed "s/^--transcriptome.*/--transcriptome=$trans/" | \
sed "s/--fastqs=\/home\/jjm247\/project\/visium/--fastqs=\/home\/jjm247\/scratch60/" > $i.12may2021; done

##### going to run the scripts in scratch 60
cd /home/jjm247/scratch60/visium_analysis_12May21

for i in scripts_12may2021/*12may2021; do sbatch $i; done
### cp loupe files in the workind directiory and give name 
ls *out | sed 's/\.out//g' | sed 's/.*_//g' | while read line; do cp $line*12may_dir/outs/cloupe.cloupe $line-cloupe.cloupe;done

#### scp to external drive. 
ls 2*out | sed 's/\.out//g' | sed 's/.*_//g' | while read line; do cp $line*12may_dir/outs/cloupe.cloupe $line-cloupe.cloupe;done

