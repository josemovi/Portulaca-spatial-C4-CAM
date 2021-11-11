### Bash commands to filter and trimm mRNA reads and use Kallisto pseudoaligner to cuantify reads counts per contig.

### Author: J.J. Moreno-Villena. Supplement Code for 'Spatial resolution of an integrated C4+CAM photosynthetic metabolism'

### Instructions to submit a job array and execute the script in parallel in a HPC with the 'Dead Simple Queue' tool.
## Save this script as 'trim_and_kallisto.sh'
## Generate list with library identifiers. Command: `ls *fastq | sed 's/_.*//g' | sort | uniq > list_lcm_libraries.txt`
## Create one execution command per line as: 

#cat list_lcm_libraries.txt | while read line; do \
#var1=$(echo "sed \"s/library/$line/g\" trim_and_kallisto.sh") ; \
#echo "eval \"\$($var1)\""; done > trim_and_jobs_list.txt

## Run the array of jobs

#dsq --job-file trim_and_jobs_list.txt -p scavenge -c 16 --mem-per-cpu 4g -t 500:00 --submit --max-jobs 15

## To run the command once, replace the word 'library' by the reads library ID. command: "sed 's/library/your_lib_ID/g' trim_and_kallisto.sh > new_trim_and_kallisto.sh"

### SCRIPT commands trim_and_kallisto.sh

readsdir=lcm_reads-dir # directorty containing read libraries

### The Illumina clip TruSeq3-PE-2.fa was modified to include primer sequences specific for SMARTer® Ultra® Low Input RNA for Illumina® sequencing (Takarabio labs) (remove #)
#>PrefixPE/1
#TACACTCTTTCCCTACACGACGCTCTTCCGATCT
#>PrefixPE/2
#GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
#>PE1
#TACACTCTTTCCCTACACGACGCTCTTCCGATCT
#>PE1_rc
#AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
#>PE2
#GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
#>PE2_rc
#AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
#>3SMARTCDSPrimerIIA_1
#AAGCAGTGGTATCAACGCAGAGTAC
#>3SMARTCDSPrimerIIA_2
#GTACTCTGCGTTGATACCACTGCTT
#>p-A
#AAAAAAAAAAAAAAAAAAAAAAAAA
#>p-T
#TTTTTTTTTTTTTTTTTTTTTTTTT

#### Execute Trimomatic, to quiality filter reads, remove Illumina adaptors and poly-A tails

tripath=path_to_Trimmomatic_directory

java -Xmx4G -jar $tripath/trimmomatic-0.39.jar PE -phred33 -threads 16 \
-trimlog logfile_library $readsdir/library_R1.fastq $readsdir/library_R2.fastq \
library_R1.PairedTrimmed.fastq \
library_R1.UnpairedTrimmed.fastq.gz \
library_R2.PairedTrimmed.fastq \
library_R2.UnpairedTrimmed.fastq.gz \
ILLUMINACLIP:$tripath/adapters/TruSeq3-PE-2.fa:3:30:10 SLIDINGWINDOW:5:20 LEADING:5 TRAILING:5 MINLEN:18 HEADCROP:10

### Execute Kallisto (needs to be installed in your system)

readsdir=path_to_lcm_reads-dir;
kallistoout=path_to_output_filtered_reads
transcript=path_to_reference_transcriptome.fasta; # e.g. Trinity.reduced.cdscdhit.fas

## Index transcriptome

kallisto index -i $transcript $transcript.kal

## Execute Kallisto pseudoaligner

kallisto quant -i $transcript.kal -o $kallistoout/library-dir -b 100 --threads 16 $readsdir/library_R1.PairedTrimmed.fastq $readsdir/library_R2.PairedTrimmed.fastq;

### GENERATE READS MAPPING STATISTICS TABLE

## in the folder where the reads are
for i in *R1.fastq; do echo $i $(( $(cat $i | wc -l) / 4 )); done > number_read_per_library.csv

# COUNT NUMBER OF READS, NUMBER OF FILTERED, NUMBER OF MAPPED BY KALLISTO
echo 'library n_raw_reads n_trimmed n_processed n_pseudoaligned n_unique p_pseudoaligned p_unique'> library_reads_mapping_stats_june2021.csv
cut -f1 -d' ' number_read_per_library.csv | while read i; do \
lib=$(echo $i | sed 's/_R.*//g'); \
echo $lib $(grep "^$lib\_R1" number_read_per_library.csv | cut -f2 -d' ') \
$(( $(cat $lib*R1.PairedTrimmed.fastq | wc -l) / 4 )) \
$(grep n_processed ../trimmed_transdec_cdhit_kallisto_bam/$lib-dir/run_info.json | sed 's/.* //g' | sed 's/,.*//g') \
$(grep n_pseudoaligned ../trimmed_transdec_cdhit_kallisto_bam/$lib-dir/run_info.json | sed 's/.* //g' | sed 's/,.*//g') \
$(grep n_unique ../trimmed_transdec_cdhit_kallisto_bam/$lib-dir/run_info.json | sed 's/.* //g' | sed 's/,.*//g') \
$(grep p_pseudoaligned ../trimmed_transdec_cdhit_kallisto_bam/$lib-dir/run_info.json | sed 's/.* //g' | sed 's/,.*//g') \
$(grep p_unique ../trimmed_transdec_cdhit_kallisto_bam/$lib-dir/run_info.json | sed 's/.* //g' | sed 's/,.*//g') \
; done > library_reads_mapping_stats_june2021.csv
