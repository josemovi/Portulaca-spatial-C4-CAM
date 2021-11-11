### Bash commands to reduce the Portulaca oleracea transcriptome from Gilman et al. (2021) to non-redundant contigs containing Coding Sequences and to functionally annotate sequences using Trinotate.

### Author: J.J. Moreno-Villena. Supplement Code for 'Spatial resolution of an integrated C4+CAM photosynthetic metabolism'

### TransDecoder identifies candidate coding regions within transcript sequences. Adapted from (https://github.com/TransDecoder/TransDecoder/wiki)
transdec_dir=directory_with_TransDecoder-v5.5.0
fasta=directory_with_Trinity.fas

## Step 1: Transdecoder: extract the long open reading frames
$transdec_dir/TransDecoder.LongOrfs -t $fasta/Trinity.fasta

## Step 2: Including homology searches as ORF retention criteria: BlastP Search
#Needs BLAST+
#Download protein database from uniprot: ~/uniprot/uniprot-reviewed-viridiplantae-A33090.fas
blastp -query longest_orfs.pep \
    -db ~/uniprot/uniprot-reviewed-viridiplantae-A33090.fas -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 12 > blastp1.outfmt6

## Step 3: predict the likely coding regions
$transdec_dir/TransDecoder.Predict -t $fasta/Trinity.fasta --retain_blastp_hits blastp1.outfmt6

### Reduce redundancy of codign sequences with CD-HIT. Modified from https://github.com/weizhongli/cdhit
outdir=dir_output_CD-HIT
cdhit_dir=dir_with_CD-HIT_program
fasta=dir_to_Trinity.fasta.transdecoder.cds

$cdhit/./cd-hit-est -i $fasta -o $workdir/Trinity.cds.cdhit.fas -c 0.95 -d 0 -M 48000 -T 12

### Subset mRNA contigs in the Trinity.fas to those sequences in CD-HIT output and clean sequence headers.

grep '>' Trinity.cds.cdhit.fas | sed 's/\.p.*//g' | sort | uniq | while read seq; do grep "$seq " Trinity.fasta -A1 ; done | sed 's/ .*//g' > Trinity.reduced.cdscdhit.fas

### Functional annotation with Trinotate and best blast hit
trinodir=dir_with_Trinotate_program
trinifas=dir_to_Trinity.fasta # Inital Trinity file
transpep=dir_to_Trinity.fasta.transdecoder.pep # Peptide file from the transdecoder analysis

# Programs installed needed: Trinity, BLAST+
# Trinity.fasta.gene_trans_map can be generated from the Fasta.file 

## run Blast and select the best hit per sequence
blastp -query $transpep -db $trinodir/uniprot_sprot.pep -num_threads 12 -outfmt 6 -evalue 1e-3 | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > best_single_hits.blastp

## Generate trinity transcripts to gene map
~/programs/trinityrnaseq-v2.10.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map

## Prepare files to run Trinotate
$trinodir/Trinotate $trinodir/Trinotate.sqlite init \
--gene_trans_map Trinity.fasta.gene_trans_map \ 
--transcript_fasta $trinifas \
--transdecoder_pep $transpep

## Run Trinotate using the Best blast hits
$trinodir/Trinotate $trinodir/Trinotate.sqlite LOAD_swissprot_blastp best_single_hits.blastp

## Create Trinotate report
$trinodir/Trinotate $trinodir/Trinotate.sqlite report -E 1e-3 > trinotate_annotation_report.xls

## GENERATE GENE TO TRANSCRIPTOME MAP WITH TRINOATE GENE ANNOTATIONS

cut -f1,2,7 ../kallisto_to_mRNA/trinotate_annotation_report.xls | sed 's/\^.*//g' > Trinity_reduced_gene_map_annotated.delim2

