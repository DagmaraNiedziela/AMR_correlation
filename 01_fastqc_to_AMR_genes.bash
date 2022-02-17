#### Dagmara Niedziela AMR correlation E coli file analysis #### 

# Task 1 - fastQC 
./software/FastQC/fastqc ./Ecoli_data/fastq/ERR3333462_1.fastq.gz --outdir ./Ecoli_analysis/fastqc_results # Directory must exist before you specify it 
./software/FastQC/fastqc ./Ecoli_data/fastq/ERR3333462_2.fastq.gz --outdir ./Ecoli_analysis/fastqc_results 

# Task 2 - fastp 
./software/fastp --in1 ./Ecoli_data/fastq/ERR3333462_1.fastq.gz --in2 ./Ecoli_data/fastq/ERR3333462_2.fastq.gz --out1=./Ecoli_analysis/fastp_results/ERR3333462_R1.clean.fastq.gz --out2=./Ecoli_analysis/fastp_results/ERR3333462_R2.clean.fastq.gz -j ./Ecoli_analysis/fastp_results/ERR3333462.json -h ./Ecoli_analysis/fastp_results/ERR3333462.html -p -c --detect_adapter_for_pe --length_required 20 --qualified_quality_phred 20 --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

# Task 3 - SPAdes 
./software/SPAdes-3.13.1-Linux/bin/spades.py -k 21,33,55,77,99,127 --careful --only-assembler --pe1-1 ./Ecoli_analysis/fastp_results/ERR3333462_R1.clean.fastq.gz  --pe1-2 ./Ecoli_analysis/fastp_results/ERR3333462_R2.clean.fastq.gz -o ./Ecoli_analysis/SPAdes_results/ERR3333462
 
# Submit a SLURM script! 
cp samplejob.sh SPAdes_script.sh 
nano SPAdes_script.sh 

sbatch SPAdes_script.sh 

# Task 4 - Quast 

# First will need to move all SPAdes results to one folder 
mkdir Ecoli_analysis/SPAdes_results/all_contigs
cp ./Ecoli_analysis/SPAdes_results/ERR3333462/contigs.fasta ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333462.fasta 

python ./software/quast-5.0.2/quast.py --labels ERR3333462.fasta,ERR3333463.fasta,ERR3333464.fasta,ERR3333532.fasta,ERR3333533.fasta,ERR3333602.fasta,ERR3333603.fasta -o ./Ecoli_analysis/Quast_results ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333462.fasta  ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333463.fasta ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333464.fasta ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333532.fasta ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333533.fasta ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333602.fasta ./Ecoli_analysis/SPAdes_results/all_contigs/ERR3333603.fasta  # labels with commas, no space, final files with spaces! 

# Task 7 - abricate 
conda activate ./conda 
abricate --db resfinder --csv ./Ecoli_data/assembly/*.fasta > ./Ecoli_analysis/abricate_results/abricate_resfinder_all.csv 
abricate --db resfinder ./Ecoli_data/assembly/*.fasta > ./Ecoli_analysis/abricate_results/resfinder_results_all.tab
abricate --summary ./Ecoli_analysis/abricate_results/resfinder_results_all.tab > ./Ecoli_analysis/abricate_results/Resfinder_summary_all.tab 

# NEW ###### 
# Abricate - EcOH 
conda activate ./conda 
abricate --list # list databases 
# Get o and h serotypes 
abricate --db ecoh ./Ecoli_data/assembly/*.fasta > ./Ecoli_analysis/abricate_results/abricate_ecoh.tab 
abricate --db ecoh --mincov 80.0 --minid 95.0 ./Ecoli_data/assembly/*.fasta > ./Ecoli_analysis/abricate_results/abricate_ecoh_filtered.tab 
abricate --summary ./Ecoli_analysis/abricate_results/abricate_ecoh_filtered.tab > ./Ecoli_analysis/abricate_results/Ecoh_summary_all.tab 

# Mlst 

# Install mlst 
conda install -c conda-forge -c bioconda -c defaults mlst 

# Run code 
mlst --csv ./Ecoli_data/assembly/*.fasta  > ./Ecoli_analysis/abricate_results/mlst_all.csv 

# Phylogroups - ClermonTyping --- https://github.com/A-BN/ClermonTyping 
# Do it online here: http://clermontyping.iame-research.center/ 
# Make sure you reference: Beghain, J., Bridier-Nahmias, A., Le Nagard, H., Denamur, E. & Clermont, O. ClermonTyping: an easy-to-use and accurate in silico method for Escherichia genus strain phylotyping. Microbial Genomics (2018). doi:10.1099/mgen.0.000192 

# Install AMRFinder ###

conda activate ./conda

conda install -y -c bioconda ncbi-amrfinderplus 

amrfinder -u 
# This makes it download the database 
	
# Test on one strain 
amrfinder --plus -n ./Ecoli_data/assembly/ESC_ZA6363AA_AS.result.fasta ./Ecoli_data/assembly/ESC_ZA6364AA_AS.result.fasta -O Escherichia > ./Ecoli_analysis/AMRFinder/a2_strains_AMRfinder_results

# For slurm script 
for i in ./Ecoli_data/assembly/*.fasta; do sn=`echo $i | cut -c 23-36`;echo -e "amrfinder --plus -n $i -O Escherichia > ./Ecoli_analysis/AMRFinder/${sn}_AMRfinder_results" >> AMRFinder.sh ; done 

# Example line 
amrfinder --plus -n ./Ecoli_data/assembly/ESC_ZA6363AA_AS.result.fasta -O Escherichia > ./Ecoli_analysis/AMRFinder/ESC_ZA6363AA_A_AMRfinder_results

echo -n ./Ecoli_data/assembly/ESC_ZA6363AA_AS.result.fasta | wc -c 
# 50 characters 

# Changed the slurm script in nano but then also manually on my comp, and used dos2unix after transfer 
sbatch slurm_AMRFinder.sh
Submitted batch job 430374 

# Append file names - first copy the directory 
cp -r AMRFinder AMRFinder2
for f in file1 file2 file3; do sed -i "s/$/\t$f/" $f; done 

cp -r ./Ecoli_analysis/AMRFinder/ ./Ecoli_analysis/AMRFinder2/
for f in ls ./Ecoli_analysis/AMRFinder2/; do sed -i "s/$/\t$f/" $f; done 
for f in *; do sed -i "s/$/\t$f/" $f; done # This one was done in the folder 
# If the above does not work, get a list of file names and include

#For each file in the list, this will use sed to append to the end of each line a tab and the filename.
#Explanation:
#  Using the -i flag with sed to perform a replacement in-place, overwriting the file
#  Perform a substitution with s/PATTERN/REPLACEMENT/. In this example PATTERN is $, the end of the line, and REPLACEMENT is \t (= a TAB), and $f is the filename, from the loop variable. The s/// command is within double-quotes so that the shell can expand variables.

cd AMRFinder2
cat *results >> all_AMRFinder2.txt
# Rows with hearders are removed by hand in Excel 

# I can see ^M signs where I have the weird line breaks 
# Get rid of these 
sed -e "s/^M//" all_BacMet_EXP.txt > all_BacMet_EXP2.txt
sed -e "s/^M//" all_BacMet_predicted.txt > all_BacMet_predicted2.txt
# To get ^M in bash, hold CTRL and press V then M - this is the only way it works, not by copying and pasting the above 
# And voilla! I have a pretty file. Headings separate each file, and I can just add the file names manually 