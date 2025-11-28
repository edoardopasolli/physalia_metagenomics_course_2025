## physalia_metagenomics_course_2025
### Lecture 1 - Pre-processing

#### Get into the right place
```
cd /home/user<YOUR USER NAME>
```

#### Set up anaconda and check whether everything's visible !
```
DON'T INSTALL IT...
##wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh
##bash Anaconda3-2024.10-1-Linux-x86_64.sh

...WE HAVE ALREADY SET UP A VERSION
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
source ${path}/activate
```

#### With the following command, you get a list of all the environments that have been created for this course
```
conda info --envs
```
#### Do you see a list of different dedicated environments? They cover as many foundamental steps 
#### Raw data pre-processing on fastq example files "seq1.fastq" and "seq2.fastq" from https://github.com/biobakery/biobakery/wiki/kneaddata
```

##conda create -n trimmomatic -c bioconda trimmomatic ## DON'T DO IT. WE DID ALREADY
##conda create -n bowtie2 -c bioconda bowtie2 ## DON'T DO IT. WE DID ALREADY
##conda create -n samtools -c bioconda samtools ## DON'T DO IT. WE DID ALREADY

mkdir 1_pre-processing
cd 1_pre-processing

wget https://github.com/biobakery/kneaddata/files/4703820/input.zip
unzip -j input.zip input/seq1.fastq input/seq2.fastq
```

#### Define variable "s" with the sampleID
```
s="seq"
```

#### Run trimmomatic
```
source ${path}/activate trimmomatic

trimmomatic PE -threads 8 -phred33 -trimlog ${s}_trimmomatic.log ${s}1.fastq ${s}2.fastq \
${s}_filtered_1.fastq ${s}_unpaired_1.fastq ${s}_filtered_2.fastq ${s}_unpaired_2.fastq \
ILLUMINACLIP:${path}/../envs/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75

for i in *.fastq; do echo -ne "${i}\t"; cat "$i" | wc -l; done
```

#### Getting human genome and generate bowtie2 index
#### Getting the file GCF_009914755.1_T2T-CHM13v2.0.fna from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/
```
human_gen_path="/home/ubuntu/shotgun_course/human_genome/"
conda deactivate
source ${path}/activate bowtie2

## CREATE THE BOWTIE2 INDEX (DON'T RUN IT! IT TAKES A FEW HOURS TO BE EXECUTED)
## mkdir -p human_genome/
##bowtie2-build ${human_gen_path}GCF_009914755.1_T2T-CHM13v2.0.fna human_genome/GCF_009914755.1_T2T-CHM13v2.0

## MAPPING OF THE READS VS THE HUMAN GENOME
bowtie2 -x ${human_gen_path}GCF_009914755.1_T2T-CHM13v2.0 -1 ${s}_filtered_1.fastq -2 ${s}_filtered_2.fastq \
    -S ${s}.sam --very-sensitive-local -p 8

conda deactivate
source ${path}/activate samtools

samtools view -bS ${s}.sam > ${s}.bam
samtools view -b -f 12 -F 256 ${s}.bam > ${s}.bothunmapped.bam
samtools sort -n -m 5G -@ 2 ${s}.bothunmapped.bam -o ${s}.bothunmapped.sorted.bam
samtools fastq ${s}.bothunmapped.sorted.bam -1 >(gzip > ${s}_filtered.final_1.fastq.gz) -2 >(gzip > ${s}_filtered.final_2.fastq.gz) -0 /dev/null -s /dev/null -n # OUTPUT IN GZ FORMAT
#rm ${s}.sam; rm ${s}.bam; rm ${s}.bothunmapped.bam; rm ${s}.bothunmapped.sorted.bam ### IF YOU WANT TO REMOVE THE INTERMEDIATE FILES

for i in *.gz; do echo -ne "${i}\t"; zcat "$i" | wc -l; done
```
#### Did the preprocessing produce the same exact number of reads in R1 and R2 ?

### End of Lecture 1 - Pre-processing
### Lecture 2 - MetaPhlAn profiling
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
```

#### MetaPhlAn for taxonomic profiling (folder "2_metaphlan")
```
conda deactivate
source ${path}/activate

##conda create -n mpa -c conda-forge -c bioconda python=3.7 metaphlan=4.2.4 ## DON'T DO IT. WE DID ALREADY
source ${path}/activate mpa

mkdir 2_metaphlan
cd 2_metaphlan

```

#### Getting example files (6 fasta files) from https://github.com/biobakery/biobakery/wiki/metaphlan4
```
mpa_db="/home/ubuntu/shotgun_course/metaphlan_databases/"
db_version="mpa_vJan25_CHOCOPhlAnSGB_202503"

wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014476-Supragingival_plaque.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014494-Posterior_fornix.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014459-Stool.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014464-Anterior_nares.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014470-Tongue_dorsum.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014472-Buccal_mucosa.fasta.gz

s="SRS014476-Supragingival_plaque"
```

#### Let's have a look a the MetaPhlAn parameters
```
metaphlan -h
```

#### Let's now run MetaPhlAn
```
metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version}

s="SRS014494-Posterior_fornix"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version}
s="SRS014459-Stool"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version}
s="SRS014464-Anterior_nares"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version}
s="SRS014470-Tongue_dorsum"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version}
s="SRS014472-Buccal_mucosa"; metaphlan ${s}.fasta.gz --input_type fasta --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version}

merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
```

#### Getting another example file (the fastq file SRS013951.fastq.bz2) from here: http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS013951.fastq.bz2

s="SRS013951";

metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version}
metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}_sub.bowtie2.bz2 --samout ${s}_sub.sam.bz2 -o ${s}_sub_profile.txt \
    --profile_vsc --vsc_out ${s}_profile.vsc.txt \
    --stat_q 0.1 --nproc 8 --subsampling 1000 --db_dir ${mpa_db} --index ${db_version}
```

#### Generate heatmap with hclust2
```
## conda create -n hclust2 -c bioconda hclust2 python=2.7 ## DON'T DO IT. WE DID ALREADY
conda deactivate

source ${path}/activate hclust2
grep -E "s__|SRS" merged_abundance_table.txt | grep -v "t__" | sed "s/^.*|//g" | sed "s/SRS[0-9]*-//g" > merged_abundance_table_species.txt

hclust2.py \
-i merged_abundance_table_species.txt \
-o metaphlan4_abundance_heatmap_species.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.5 \
--log_scale \
--flabel_size 10 --slabel_size 10 \
--max_flabel_len 100 --max_slabel_len 100 \
--minv 0.1 \
--dpi 300
```

#### Characterization of the virome
```
## LET'S CONSIDER THIS PAPER (https://www.nature.com/articles/s41467-022-33395-6) AND IDENTIFY THE BIOPROJECT (PRJNA746322).
## WE COULD DOWNLOAD THE 68 METAGENOMES USING THE SRA EXPLORER (https://sra-explorer.info/)
## AND THEN PROCESS THEM USING THE COMMAND PROVIDED EARLIER.

## FOR NOW, JUST COPY THE PROFILES:

mkdir Bernard-RaichonL_2022
cd Bernard-RaichonL_2022

cp /home/ubuntu/shotgun_course/2_metaphlan/Bernard-RaichonL_2022/*_profile.txt ./
cp /home/ubuntu/shotgun_course/2_metaphlan/Bernard-RaichonL_2022/*_profile.vsc.txt ./

## AND GENERATE THE MERGED TABLE WITH THE VIRAL PROFILES:

merge_vsc_tables.py -g depth_of_coverage_mean *_profile.vsc.txt > Bernard-RaichonL_2022_profile.vsc.tsv
```

#### Processing of Nanopore long reads sequences
```
mkdir -p zymo
cd zymo
```

## Way n. 1: Downloading the Nanopore sequence data and running MetaPhlAn (Don't do it, it is time consuming)

```
## DOWNLOAD THE NANOPORE SEQUENCE DATA, MADE AVAILABLE HERE: https://github.com/LomanLab/mockcommunity
#wget https://nanopore.s3.climb.ac.uk/Zymo-GridION-EVEN-BB-SN.fq.gz

s="Zymo-GridION-EVEN-BB-SN"

#metaphlan ${s}.fq.gz --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt \
#    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version} --long_reads

#metaphlan ${s}.fq.gz --input_type fastq --mapout ${s}.split.bowtie2.bz2 --samout ${s}.split.sam.bz2 -o ${s}.split_profile.txt \
#    --stat_q 0.1 --nproc 8 --db_dir ${mpa_db} --index ${db_version} --split_reads
```

## Way n. 2: Copying the profiling files (generated by MetaPhlAn)
```
cp /home/ubuntu/shotgun_course/2_metaphlan/zymo/*_profile.txt ./
```

### End of Lecture 2 - MetaPhlAn profiling
### Lecture 3 - GraPhlAn
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
```

#### GraPhlAn for circular representations of taxonomic and phylogenetic trees (folder "3_graphlan")
```
conda deactivate
source ${path}/activate

## conda create -n graphlan -c bioconda graphlan ## DON'T DO IT. WE DID ALREADY
source ${path}/activate graphlan

mkdir 3_graphlan
cd 3_graphlan
```

#### Getting example files from https://github.com/biobakery/graphlan/tree/master/examples/guide
```
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/guide.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_0.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_1.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_2.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_3.txt

graphlan.py guide.txt step_0.png --dpi 300 --size 3.5
graphlan.py guide.txt step_0.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_0.txt guide.txt guide_1.xml
graphlan.py guide_1.xml step_1.png --dpi 300 --size 3.5
graphlan.py guide_1.xml step_1.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_1.txt guide_1.xml guide_2.xml
graphlan.py guide_2.xml step_2.png --dpi 300 --size 3.5
graphlan.py guide_2.xml step_2.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_2.txt guide_2.xml guide_3.xml
graphlan.py guide_3.xml step_3.png --dpi 300 --size 3.5
graphlan.py guide_3.xml step_3.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_3.txt guide_3.xml guide_4.xml
graphlan.py guide_4.xml step_4.png --dpi 300 --size 3.5 --pad 0.0
graphlan.py guide_4.xml step_4.svg --dpi 300 --size 3.5 --pad 0.0
```

#### Getting another example (PhyloPhlAn) from https://github.com/biobakery/graphlan/tree/master/examples/PhyloPhlAn
```
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/PhyloPhlAn/ppa_tol.xml 
wget -O annot_PhyloPhlAn.txt https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/PhyloPhlAn/annot.txt

graphlan_annotate.py --annot annot_PhyloPhlAn.txt ppa_tol.xml ppa_tol.annot.xml 
graphlan.py ppa_tol.annot.xml ppa_tol.png --dpi 200 --size 15 --pad 0.6
```

#### Getting another example (HMP_tree) from https://github.com/biobakery/graphlan/tree/master/examples/HMP_tree
```
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/HMP_tree/hmptree.xml
wget -O annot_HMP_tree.txt https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/HMP_tree/annot.txt

graphlan_annotate.py --annot annot_HMP_tree.txt hmptree.xml hmptree.annot.xml 
graphlan.py hmptree.annot.xml hmptree.png --dpi 150 --size 14 
```

### End of Lecture 3 - GraPhlAn

### Lecture 4 - StrainPhlAn 
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
```

#### StrainPhlAn for strain-level profiling (folder "4_strainphlan")
```
conda deactivate
source ${path}/activate
source ${path}/activate mpa

mkdir 4_strainphlan
cd 4_strainphlan
```

#### Getting example files (6 fastq files) from https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS013951.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS014613.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS019161.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS022137.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS055982.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS064276.fastq.bz2
```

## Way n. 1: Running MetaPhlAn 
```
mpa_db="/home/ubuntu/shotgun_course/metaphlan_databases/"
db_version="mpa_vJan25_CHOCOPhlAnSGB_202503"

## s="SRS013951"; metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --db_dir ${mpa_db} --index ${db_version}
## s="SRS014613"; metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --db_dir ${mpa_db} --index ${db_version}
## s="SRS019161"; metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --db_dir ${mpa_db} --index ${db_version}
## s="SRS022137"; metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --db_dir ${mpa_db} --index ${db_version}
## s="SRS055982"; metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --db_dir ${mpa_db} --index ${db_version}
## s="SRS064276"; metaphlan ${s}.fastq.bz2 --input_type fastq --mapout ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8 \
##     --db_dir ${mpa_db} --index ${db_version}
```

## Way n. 2: Copying the alignment files (generated by MetaPhlAn)
```
cp /home/ubuntu/shotgun_course/4_strainphlan/SRS013951.sam.bz2 ./
cp /home/ubuntu/shotgun_course/4_strainphlan/SRS014613.sam.bz2 ./
cp /home/ubuntu/shotgun_course/4_strainphlan/SRS019161.sam.bz2 ./
cp /home/ubuntu/shotgun_course/4_strainphlan/SRS022137.sam.bz2 ./
cp /home/ubuntu/shotgun_course/4_strainphlan/SRS055982.sam.bz2 ./
cp /home/ubuntu/shotgun_course/4_strainphlan/SRS064276.sam.bz2 ./

mpa_database=${mpa_db}"/"${db_version}".pkl"
sample2markers.py -i *.sam.bz2 -o ./ -n 8 -d ${mpa_database}

mkdir -p db_markers
## extract_markers.py -c t__SGB1877 -o db_markers/ -d ${mpa_database} ## TOO LONG,
## DO THIS INSTEAD:

cp /home/ubuntu/shotgun_course/4_strainphlan/db_markers/t__SGB1877.fna db_markers/
```

#### Getting a reference genome ("GCF000273725")
```
mkdir -p reference_genomes
wget -P reference_genomes/ http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/reference_genomes/G000273725.fna.bz2
```

#### Let's have a look at the StrainPhlAn params
```
strainphlan -h
```

#### Build the multiple sequence alignment and the phylogenetic tree:
```
mkdir -p output
strainphlan -s *.json.bz2 -m db_markers/t__SGB1877.fna -r reference_genomes/G000273725.fna.bz2 -o output -c t__SGB1877 -n 8 -d ${mpa_database}
```

#### Getting the metadata file ("metadata.txt")
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/metadata.txt
add_metadata_tree.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre -f metadata.txt -m subjectID --string_to_remove .fastq.bz2

conda deactivate
source ${path}/activate graphlan
${path}/../envs/mpa/bin/plot_tree_graphlan.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre.metadata -m subjectID
```

### Lecture 4 - StrainPhlAn 
### Lecture 5 - PanPhlAn 
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
```

#### PanPhlAn for Pangenome-based Phylogenomic Analysis (folder "5_panphlan")
```
conda deactivate
source ${path}/activate

## conda create -n panphlan -c bioconda panphlan ## DON'T DO IT. WE DID ALREADY
source ${path}/activate panphlan
## conda install -c conda-forge matplotlib ## DON'T DO IT. WE DID ALREADY

mkdir 5_panphlan
cd 5_panphlan
```

#### Getting fastq example files from https://github.com/SegataLab/panphlan/wiki/Tutorial-3_0
```
## wget https://www.dropbox.com/s/oi26jg0v7ktlavc/panphlan_tutorial_samples.tar.bz2
## tar -xvjf panphlan_tutorial_samples.tar.bz2 ## TOO LONG,
## DO THIS INSTEAD:

mkdir -p samples_fastq
cp /home/ubuntu/shotgun_course/5_panphlan/samples_fastq/CCMD34381688ST-21-0.fastq samples_fastq/

panphlan_download_pangenome.py -i Eubacterium_rectale -o ./

mkdir -p map_results

s="CCMD34381688ST-21-0"
#panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8 ## TOO LONG,
## DO THIS INSTEAD:

cp /home/ubuntu/shotgun_course/5_panphlan/map_results/${s}_erectale.tsv map_results/
```

#### The same script must be run for the other samples (fastq files)
```
##s="G78505"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
##s="G88884"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
##s="G88970"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
##s="G89027"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
##s="H2M514903"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
##s="HD-1"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
##s="T2D-063"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
##s="T2D-105"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8

#### WE CAN JUST COPY THEM...
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/G78505_erectale.tsv map_results/
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/G88884_erectale.tsv map_results/
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/G88970_erectale.tsv map_results/
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/G89027_erectale.tsv map_results/
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/H2M514903_erectale.tsv map_results/
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/HD-1_erectale.tsv map_results/
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/T2D-063_erectale.tsv map_results/
cp /home/ubuntu/shotgun_course/5_panphlan/map_results/T2D-105_erectale.tsv map_results/

panphlan_profiling.py -i map_results/ --o_matrix ./result_profile_erectale.tsv -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv --o_covplot ./erectale_covplot
panphlan_profiling.py -i map_results/ --o_matrix ./result_profile_erectale_annotation.tsv -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv --func_annot Eubacterium_rectale/panphlan_Eubacterium_rectale_annot.tsv --field 2
panphlan_profiling.py -i map_results/ --o_matrix ./result_profile_erectale_annotation_withref.tsv -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv --func_annot Eubacterium_rectale/panphlan_Eubacterium_rectale_annot.tsv --field 2 --add_ref
```

#### Post-processing to create heatmap and MDS from the PanPhlAn output results
We also use the file "metadata_erectale.txt" available in https://github.com/SegataLab/panphlan/wiki/Tutorial-3_0

See the R script and associated files in the subfolder [https://github.com/edoardopasolli/physalia_metagenomics_2025/tree/main/5_panphlan-Rcode](5_panphlan-Rcode)

### End of Lecture 5 - PanPhlAn 
### Lecture 6 - HUMAnN 4
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
```

#### HUMAnN for profiling the abundance of microbial metabolic pathways and other molecular functions (folder "6_humann")
```
conda deactivate
source ${path}/activate

## conda create -n humann -c bioconda python=3.7 ## DON'T DO IT. WE DID ALREADY
source ${path}/activate humann
## pip install humann==4.0.0a1 ## DON'T DO IT. WE DID ALREADY

mkdir 6_humann
```

#### Test the local HUMAnN environment
```
humann_test
humann_config
```

#### Let's look at the HUMAnN parameters !
```
humann -h
```

#### The HUMAnN databases
```
## humann_databases --download chocophlan full humann_dbs --update-config yes ## DON'T DO IT. WE DID ALREADY
## humann_databases --download uniref uniref90_ec_filtered_diamond humann_dbs --update-config yes ## DON'T DO IT. WE DID ALREADY
## humann_databases --download utility_mapping full humann_dbs --update-config yes ## DON'T DO IT. WE DID ALREADY
```

#### Getting example of fastq file from EBI
```
cd 6_humann
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/096/SRR15408396/SRR15408396.fastq.gz
```

#### Run HUMAnN to perform functional profiling
```
s="SRR15408396"

## NOW YOU CAN RUN:
## humann --input ${s}.fastq.gz --output ${s} --threads 8 --metaphlan ${path}/../envs/mpa-4.1/bin/ --metaphlan-options "--index mpa_vOct22_CHOCOPhlAnSGB_202403 -t rel_ab_w_read_stats"

## BUT IT TAKES A FEW HOURS... OR YOU CAN RUN:
mkdir -p ${s}
cp /home/ubuntu/shotgun_course/6_humann/${s}/* ${s}/
```

#### Manipulating HUMAnN output files
```
grep "|" ${s}/${s}_2_genefamilies.tsv | cut -f2 | paste -s -d+ - | bc

humann_regroup_table -i ${s}/${s}_2_genefamilies.tsv -o ${s}/${s}_level4ec.tsv --groups uniclust90_level4ec

humann_rename_table -i ${s}/${s}_level4ec.tsv -o ${s}/${s}_level4ec_named.tsv --names ec
```

#### Getting another fastq file from EBI and run HUMAnN
```
s="SRR15408398"

## SAME:
## wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/098/SRR15408398/SRR15408398.fastq.gz
## humann --input ${s}.fastq.gz --output ${s} --threads 8 --metaphlan ${path}/../envs/mpa-4.1/bin/ --metaphlan-options "--index mpa_vOct22_CHOCOPhlAnSGB_202403 -t rel_ab_w_read_stats"

## FOR NOW, RUN:
mkdir -p ${s}
cp /home/ubuntu/shotgun_course/6_humann/${s}/* ${s}/
```

#### Merge multiple outputs
```
mkdir -p merged

for s in SRR15408396 SRR15408398; do
   cp ${s}/${s}_2_genefamilies.tsv merged/
   cp ${s}/${s}_3_reactions.tsv merged/
   cp ${s}/${s}_4_pathabundance.tsv merged/
done

humann_join_tables -i merged -o merged/merged_2_genefamilies.tsv --file_name _2_genefamilies.tsv
humann_join_tables -i merged -o merged/merged_3_reactions.tsv --file_name _3_reactions.tsv
humann_join_tables -i merged -o merged/merged_4_pathabundance.tsv --file_name _4_pathabundance.tsv
```

#### Plotting stratified functions
We use the file "hmp_pathabund.pcl" (generated by HUMAnN v2) available in https://raw.githubusercontent.com/biobakery/biobakery/master/test_suite/biobakery_tests/data/humann/input/hmp_pathabund.pcl
```
mkdir -p hmp; cd hmp
wget https://raw.githubusercontent.com/biobakery/biobakery/master/test_suite/biobakery_tests/data/humann/input/hmp_pathabund.pcl

# conda install numpy scipy matplotlib ## DON'T DO IT. WE DID ALREADY

humann_barplot -i hmp_pathabund.pcl -o hmp_plot1.png --focal-metadata STSite --last-metadata STSite --focal-feature METSYN-PWY
humann_barplot -i hmp_pathabund.pcl -o hmp_plot2.png --focal-metadata STSite --last-metadata STSite --focal-feature METSYN-PWY --sort sum
humann_barplot -i hmp_pathabund.pcl -o hmp_plot3.png --focal-metadata STSite --last-metadata STSite --focal-feature METSYN-PWY --sort sum metadata --scaling logstack
```

### End of Lecture 6 - HUMAnN 3
### Lecture 7 - Metagenomic assembly
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
```

#### Megahit for de novo metagenomic assembly (folder "7_assembly")
```
# We can install Megahit in this way:
wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz

# Otherwise through the usual conda:
## conda create -n megahit -c bioconda megahit ## DON'T DO IT. NOT NECESSARY
```

#### The megahit paramaters:
```
${path}/../../MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -h
```

#### Getting example of fastq file from https://github.com/voutcn/megahit/wiki/An-example-of-real-assembly
```
mkdir 7_assembly
cd 7_assembly

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/SRR341725/SRR341725_[12].fastq.gz
```

#### Run megahit to generate metagenomic assembly
```
s="SRR341725"
## MEGAHIT WILL TAKE A FEW HOURS:
## ${path}/../../MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -o ${s}.megahit_asm -t 8

## FOR NOW WE CAN COPY THE RESULTS FROM MEGAHIT
mkdir -p ${s}.megahit_asm/

cp /home/ubuntu/shotgun_course/7_assembly/${s}.megahit_asm/final.contigs.fa  ${s}.megahit_asm/

## WE ALSO NEED TWO CUSTOM SCRIPTS:
cp /home/ubuntu/shotgun_course/7_assembly/filter_contigs.py ./
cp /home/ubuntu/shotgun_course/7_assembly/megahit2spades.py ./
```

#### Do some post-processing on the contigs file
```
source ${path}/activate megahit
## conda install -c conda-forge biopython ## DON'T DO IT. WE DID ALREADY
## conda install -c anaconda pandas ## DON'T DO IT. WE DID ALREADY

python megahit2spades.py ${s}.megahit_asm/final.contigs.fa ${s}.megahit_asm/contigs.fasta
python filter_contigs.py ${s}.megahit_asm/contigs.fasta ${s}.megahit_asm/contigs_filtered.fasta
python filter_contigs.py ${s}.megahit_asm/contigs.fasta ${s}.megahit_asm/contigs_filtered_50000.fasta -l 50000
```

#### Run Flye to perform assembly of long-read sequences
```
conda deactivate
source ${path}/activate

## conda create -n flye flye ## DON'T DO IT. WE DID ALREADY
source ${path}/activate flye

## flye --nano-raw ../2_metaphlan/zymo/Zymo-GridION-EVEN-BB-SN.fq.gz --out-dir zymo.flye_asm --meta --threads 8 ## TOO LONG,
## DO THIS INSTEAD:

cp -r /home/ubuntu/shotgun_course/7_assembly/zymo.flye_asm/ ./
```

### End of Lecture 7 - Metagenomic assembly
### Lecture 8 - MAG reconstruction
```
cd /home/user<YOUR USER NAME>
path="/home/ubuntu/shotgun_course/anaconda3course/bin/"
```

## Metabat2 for metagenomic binning (folder "8_MAG-reconstruction")
```
source ${path}/activate
 
## conda create -n metabat2 -c bioconda metabat2 ## DON'T DO IT. WE DID ALREADY
source ${path}/activate metabat2

## conda install -c bioconda bowtie2 ## DON'T DO IT. WE DID ALREADY
## conda install -c bioconda samtools ## DON'T DO IT. WE DID ALREADY

mkdir 8_MAG-reconstruction
cd 8_MAG-reconstruction
```

## Copy the raw reads and contigs generated in the previous tutorial (folder "7_assembly")
```
s="SRR341725"

cp ../7_assembly/SRR341725.megahit_asm/contigs_filtered.fasta ./
cp ../7_assembly/SRR341725_1.fastq.gz ./
cp ../7_assembly/SRR341725_2.fastq.gz ./
```

## Mapping of raw reads against contigs
```
bowtie2-build contigs_filtered.fasta contigs_filtered
## bowtie2 -x contigs_filtered -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -S ${s}.sam -p 8 2> ${s}.bowtie2.log ## TOO LONG
## samtools view -bS ${s}.sam > ${s}.bam ## TOO LONG
## samtools sort ${s}.bam -o sorted_${s}.bam ## TOO LONG

## COPY THE RESULT FOR NOW:
cp /home/ubuntu/shotgun_course/8_MAG-reconstruction/sorted_SRR341725.bam ./
```

## Run metabat2 to reconstruct metagenome-assembled genomes (MAGs)
```
jgi_summarize_bam_contig_depths --outputDepth ${s}_depth.txt sorted_${s}.bam 2> ${s}_depth.log
metabat2 -i contigs_filtered.fasta -a ${s}_depth.txt -o ${s}_bins/bin -m 1500 --unbinned -t 8 > ${s}_metabat2.log
```

## CheckM2 to estimate MAG quality
```
conda deactivate
source ${path}/activate

## conda create -n checkm2 -c bioconda checkm2 ## DON'T DO IT. WE DID ALREADY
source ${path}/activate checkm2

## LET'S DOWNLOAD THE DATABASE
## checkm2 database --download --path ../ ## DON'T DO IT. WE DID ALREADY

## WE CAN USE A COPY HERE:
checkm2_db="/home/ubuntu/shotgun_course/CheckM2_database/uniref100.KO.1.dmnd"
checkm2 testrun --database_path ${checkm2_db} --threads 8

## checkm2 predict -i SRR341725_bins -o SRR341725_checkm2 -x .fa --database_path ${checkm2_db} --threads 8 ## TOO LONG

## COPY THE RESULT FOR NOW:
cp -r /home/ubuntu/shotgun_course/8_MAG-reconstruction/SRR341725_checkm2 ./

awk -F'\t' '$2 > 50 && $3 < 5' SRR341725_checkm2/quality_report.tsv > SRR341725_checkm2/quality_report_filtered.tsv

mkdir -p ${s}_bins_filtered
cut -f1 SRR341725_checkm2/quality_report_filtered.tsv | while read -r value; do cp ${s}_bins/${value}.fa ${s}_bins_filtered/; done
```

## Run PhyloPhlAn to perform taxonomic assignment
```
conda deactivate

## conda create -n phylophlan -c anaconda biopython=1.73
source ${path}/activate phylophlan
## conda install -c bioconda phylophlan
```

## Let's have a look at PhyloPhlAn commands:
```
phylophlan_assign_sgbs -h
```

## Let's run the PhyloPhlAn taxonomic assignment tool
```
## WE'LL USE A COPY OF THIS DATABASE TO SPARE DOWNLOAD TIME
database_folder="/home/ubuntu/shotgun_course/phylophlan_databases/"

## WE RUN THE TAXONOMIC ASSIGNMENT
## phylophlan_assign_sgbs -i SRR341725_bins_filtered -o SRR341725_bins_filtered_phylophlan \ ## TOO LONG
##    -d SGB.Jun23 --database_folder ${database_folder} \
##    -n 1 --verbose --nproc 8 2>&1 | tee SRR341725_phylophlan.log

## COPY THE RESULT FOR NOW:
cp /home/ubuntu/shotgun_course/8_MAG-reconstruction/SRR341725_phylophlan.log ./
cp /home/ubuntu/shotgun_course/8_MAG-reconstruction/SRR341725_bins_filtered_phylophlan.tsv ./
```


#### Run Bakta for rapid genome annotation
```
conda deactivate

## conda create -n bakta -c conda-forge -c bioconda bakta  ## DON'T DO IT. WE DID ALREADY
source ${path}/activate bakta

## LET'S DOWNLOAD THE DATABASE
## bakta_db download --output ../bakta_db --type full ## DON'T DO IT. WE DID ALREADY
bakta_db="/home/ubuntu/shotgun_course/bakta_db/db"

cd SRR341725_bins_filtered
mkdir -p tmp
s="bin.8"
## bakta --db ${bakta_db} --min-contig-length 200 --prefix ${s} --output ${s}_bakta/ --compliant --meta --threads 8 --tmp-dir ./tmp ${s}.fa ## TOO LONG

## COPY THE RESULT FOR NOW:
cp -r /home/ubuntu/shotgun_course/8_MAG-reconstruction/SRR341725_bins_filtered/bin.8_bakta ./
```


### End of Lecture 8 - MAG reconstruction
### Lecture 9 - Statistical analysis
See the R script here in subfolder [https://github.com/edoardopasolli/physalia_metagenomics_course_2025/tree/main/9_statistical-analysis-Rcode](9_statistical-analysis-Rcode)
### End of Lecture 9 - Statistical analysis
