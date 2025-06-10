## Atlantic White‐cedar: Population Genetics Analysis Comparisons
### Prior to these analyses, SNPs were called by aligning reads to 3 separate references:

### De novo Reference
<details><summary> De novo - Assembled from GBS reads using the dDocent pipeline (ranbow + CD hit) </summary>

***



</details>

### Close Reference
<details><summary> Chamaecyparis obtusa - The closest reference available (HiFi-based) </summary>

***

<details><summary> 04. BWA: Alignment to the Chamaecyparis obtusa genome - 10/24/23 </summary>
<p>

1. I made a new directory for this analysis and created a directory to house my reference.
* filepath:
```
cd /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/
mkdir 03_C.obtusa_reference-aligned
cd 03_C.obtusa_reference-aligned
mkdir C.obtusa_reference
cd C.obtusa_reference
```

2. Then, I downloaded the Chamaecyparis obtusa genome from https://drive.google.com/drive/u/2/folders/1YKzyQJ3td4W3jAomRIeWizCf5D-1L-I5, and I secure-copied it from my desktop to Flora. The original paper is here: https://www-tandfonline-com.utk.idm.oclc.org/doi/full/10.1080/13416979.2023.2267304.

3. To get things ready for bwa and GATK, I unzipped and indexed the _C. obtusa_ genome with bwa and samtools. I also created a sequence dictionary with picard.
```
gunzip index COB_r1.0.fasta.gz
conda activate gatk # this environment includes my picard installation.
```
```
nano run_indexing.qsh
```
```
#!/bin/bash
#SBATCH --job-name=indexing
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=40G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=12:00:00
#SBATCH --mail-user=zsmith10@vols.utk.edu

INPUT_FILE=COB_r1.0.fasta
OUTPUT_FILE=$( basename $INPUT_FILE | sed 's/.fasta//')

module load bwa
module load samtools
module load gatk

bwa index $INPUT_FILE \
&& \
samtools faidx $INPUT_FILE \
&& \
picard CreateSequenceDictionary R=$INPUT_FILE O=$OUTPUT_FILE\.dict
```
```
sbatch run_indexing.qsh
```

4. Next, I created a main analysis directory with a bwa analysis directory in it and symbolically linked all of my trimmed read files.
```
cd /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/03_C.obtusa_reference-aligned
mkdir analysis
cd analysis
mkdir 04_bwa
cd 04_bwa
ln -s /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/01_T.plicata_reference-aligned/analysis/02_fastp/fastp_out/*fq.gz
```

5. In this case there was no need to remove the failed samples as they are already in an `excluded_samples` directory within my `fastp_out` directory.

6. I ran a task array on the ISAAC HPC to align all reads files to the _Chamaecyparis obtusa_ genome and output merged binary alignment map (BAM) files.
```
#!/bin/bash
#SBATCH --job-name=bwa_array
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=30G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --array=0-90
#SBATCH -o bwa_outs/%x_%A_%a.out
#SBATCH -e bwa_outs/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *r1.fq.gz))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# create an output file name
OUT=$(echo $ARRAY_FILE | sed 's/.r1.fq.gz//')

# define read 1 and read 2 from the array
r1=${ARRAY_FILE}
r2=`sed 's/.r1.fq.gz/.r2.fq.gz/' <(echo ${ARRAY_FILE})`

echo "Read 1: $r1"
echo "Read 2: $r2"


module load bwa
module load samtools

bwa mem -t 5 \
        /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/03_C.obtusa_reference-aligned/C.obtusa_reference/COB_r1.0.fasta \
        $r1 \
        $r2 \
        | samtools view -bSh \
        | samtools sort \
        -@ 5 -m 4G \
        -o $OUT\_sorted.bam

echo Files $r1 and $r2 were aligned by task number $SLURM_ARRAY_TASK_ID on $(date)
```

7. Then, I ran samtools flagstats to look at alignment statistics.
```
nano run_flagstats.qsh
```
```
#!/bin/bash
#SBATCH -J flagstats
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load samtools

for b in *.bam
do
        stats_out=$( basename $b | sed 's/.bam/.stats/g')

        samtools flagstat $b > $stats_out
done
```
```
sbatch run_flagstats.qsh
```

8. Then I collated the most useful statistics from the flagstats outputs into a single file.
```
nano collate_flagstats.sh
```
```
for file in *stats
do
        echo $file >> flagstat_summary.tsv &&\
        grep '+ 0 in total' $file >> flagstat_summary.tsv &&\
        grep '+ 0 mapped' $file >> flagstat_summary.tsv &&\
        grep '+ 0 primary mapped' $file >> flagstat_summary.tsv &&\
        grep '+ 0 supplementary' $file  >> flagstat_summary.tsv
done
```
```
bash collate_flagstats.sh
```

9. Finally, I curated the stats into a tsv file to copy into a spreadsheet in google drive [here](https://docs.google.com/spreadsheets/d/1KL36GsyMEJyrZzoKnBIFlJDbhraCc4ZqQv5sTebNm6w/edit#gid=250450826). 
```
nano curate_flagstats.sh
```
```
awk 'NR % 5 == 1' flagstat_summary.tsv > column1.tsv
awk 'NR % 5 == 2' flagstat_summary.tsv > column2.tsv
awk 'NR % 5 == 3' flagstat_summary.tsv > column3.tsv
awk 'NR % 5 == 4' flagstat_summary.tsv > column4.tsv
awk 'NR % 5 == 0' flagstat_summary.tsv > column5.tsv

paste column1.tsv column2.tsv column3.tsv column4.tsv column5.tsv > final_flagstat_summary.tsv
```
```
bash curate_flagstats.sh
```

</details>
</p>

***

<details><summary> 05. GATK: HaplotypeCaller </summary>
<p>

1. Create a GATK environment on ISAAC and install picard tools. **It is important to install this conda environment, but to call the program through spack by using the `module load gatk` command in order 1) avoid dependency issues, while 2) still ensuring that gatk is routed through the appropriate compiler on the ISAAC HPC system. I have had odd software conflicts/failures by attempting to run the version of GATK installed in the conda environment, as well as the natively installed GATK version outside of my conda environment. At some point, I will explore this odd interaction further, but this is a reliable fix for the moment.**

```
conda create gatk -c bioconda gatk4
conda install -c bioconda picard
conda activate gatk
```

2. Link the sorted BAMs from the 04_bwa directory in a new 05_gatk_haplotypecaller/ directory.
```
cd /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/03_C.obtusa_reference-aligned/analysis/
mkdir 04=5_gatk_haplotypecaller/
ln -s ../04_bwa/*sorted.bam .
```

6. Run the following array to add/replace read groups using picard.
```
nano run_picard_array.qsh
```
```
#!/bin/bash
#SBATCH --job-name=picard
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --array=0-90
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *_sorted.bam))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# create an output file name
OUT=`sed 's/_sorted.bam//' <(echo ${ARRAY_FILE})`

# define the bam file from the array
BAM=${ARRAY_FILE}

echo "BAM $BAM"
echo "OUT $OUT"

picard \
        AddOrReplaceReadGroups \
        I=${BAM} \
        O=${OUT}_sorted.RG.bam \
        RGSM=$OUT \
        RGLB=$OUT \
        RGPL=illumina \
        RGPU=$OUT

echo BAM file $BAM had reads groups added/replaced by picard-2.27.3 in Slurm Task ID $SLURM_ARRAY_TASK_ID on $(date).
```
```
sbatch run_picard_array.qsh
```
* Note 1: the array numbers were run in batches of 90 (Run 1: 0-90, Run 2: 91-180, Run 3: 181-269)

7. Index the sorted bam files using samtools.
```
nano run_samtools_index.qsh
```
```
#!/bin/bash
#SBATCH --job-name=samtools_index
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=6:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load samtools

for file in *_sorted.RG.bam; do samtools index $file; done
```
```
sbatch run_samtools_index.qsh
```

8. Run the following gatk HaplotypeCaller array to create `g.vcf` files.
```
nano run_haplotypecaller_array.qsh
```
```
#!/bin/bash
#SBATCH --job-name=GATK_HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=10G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --array=261-270
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *_sorted.RG.bam))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# define the BAM file from the array file
BAM=${ARRAY_FILE}

# create an output file name
OUT=$( basename $BAM | sed 's/_sorted.RG.bam//g' )


echo "BAM: $BAM"
echo "OUT: $OUT"

module load gatk

gatk HaplotypeCaller \
    -R /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/03_C.obtusa_reference-aligned/C.obtusa_reference/COB_r1.0.fasta \
    -I $BAM \
    -O ${OUT}.g.vcf \
    --native-pair-hmm-threads 10 \
    -ERC GVCF \

echo BAM file $BAM was processed by the GATK Haplotype Caller by Slurm Task ID $SLURM_ARRAY_TASK_ID on $(date)
```
```
sbatch run_haplotypecaller_array.qsh
```

9. To ensure the HaplotypeCaller jobs are completing, you can manually check the `.err` files with a text editor, which shows the entire task log, or use the following shortcut and ensure that your actual number of completed jobs matches with your expected number of completed jobs.
```
grep 'Traversal complete' *err | wc -l
```

10. Move array `.err` and `.out` files to a new directory to reduce clutter.
```
mkdir gatk_haplotypecaller_outs
mv *err *out gatk_haplotypecaller_outs
```


</details>
</p>

***

<details><summary> 06. GATK: GenomicsDBImport & GenotypeGVCFs </summary>
<p>

### 06. GATK:  GenomicsDBImport & GenotypeGVCFs

1. Make a new directory in the analysis directory to house this analysis.
```
cd ..
mkdir 06_gatk_genotype_gvcfs
cd 06_gatk_genotype_gvcfs
```

2. Link the `g.vcf` files into the new directory.
```
ln -s ../02_gatk/*g.vcf* .
```

3. Create a list of all `g.vcf` files to feed into the CombineGVCFs tool.
```
ls *g.vcf > gvcfs.list
```

4. Unfortunately, the CombineGVCFs tool did not perform well on this dataset, so I attempted to use GenomicsDBImport, another GATK module. First I prepared an intervals file to feed into the tool. The following script will print the name of each scaffold in the assembly, remove the `>` in front of the scaffold name, as well as remove any genomic coordinate data, such as ">29382484 **1871870 0 29346470+,...,29353346+**". Failure to remove these features will result in failure of the GenomicsDBImport tool as GATK's GenotypeGVCFs tool does not include any genomic coordinate data in the resulting `g.vcf` files--only the scaffold name, in the above case: **29382484**.
* Note: This file must end in a `.list` suffix to be properly read by GenomicsDBImport.

```
nano create_intervals_file.qsh
```
```
# Create the intervals.list file.
grep '>' /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/03_C.obtusa_reference-aligned/C.obtusa_reference/COB_r1.0.fasta | sed 's/>//' >> intervals.list

input_file="intervals.list"

while IFS= read -r line; do
  # Surround each line with backticks and print the result
  echo "$line" | cut -d' ' -f1 >> trimmed_intervals.list
done < "$input_file"

rm intervals.list
```
```
bash create_intervals_file.qsh
```

5. Now, run the GenomicsImportDB tool using the following script. The following tool requires quite a bit of nuance to run appropriately, and it does not produce a simple combined GVCF, instead it creates a database file that includes all of joint-call data.  
```
nano run_genomicsdb_import.qsh
```
```
#!/bin/bash
#SBATCH --job-name=GenomicsDBImport
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=300G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=condo-ut-genomics
#SBATCH --qos=genomics
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load gatk

gatk \
        --java-options "-Xms280G -Xmx280G -XX:ParallelGCThreads=2" \
        GenomicsDBImport \
        --genomicsdb-shared-posixfs-optimizations \
        --reference /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/03_C.obtusa_reference-aligned/C.obtusa_reference/COB_r1.0.fasta \
        --variant gvcfs.list \
        --genomicsdb-workspace-path Cthyoides_Cobtusa-ref_database \
        --intervals trimmed_intervals.list \
        --batch-size 50 \
        --merge-contigs-into-num-partitions 25 \
        --reader-threads 2
```
```
sbatch run_genomicsdb_import.qsh
```
* `--java-options "-Xmx280g -Xms280g"` Sets the heap memory for java to 460 Gb (~10% less than the total requested on ISAAC).
* `GenomicsDBImport` specifies the GATK tool.
* `--genomicsdb-shared-posixfs-optimizations` toggles optimizations for HPC systems, such as those run on lustre (like ISAAC).
* `--reference` designates the path to the reference sequence.
* `--variant gvcfs.list` designates the list of gvcf files.
* `--genomicsdb-workspace-path` designates the workspace output, aka the genomicsdb database directory. This will be your input file for GenotypeGVCFs.
* `--overwrite-existing-genomicsdb-workspace`. This tool fails if the database already exists. I enabled this due to several failed trial runs due to forgetting to remove the previous database file. Use at your own risk.
* `--intervals trimmed_intervals.list` designates the genomic intervals (scaffolds, in my case) in the reference to read over.
* `--batch-size 50` sets the batch size of files to process. 50 is the recommended value used in-house at the Broad Institute (GATK's developer). Using the default number (all files) caused me to exceed the SLURM walltime. More info here: https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines
* `--merge-contigs-into-num-partitions 25`. This flag merges smaller contigs (but supposedly not larger ones) to improve speed and efficiency of the tool. This setting was recommended by the Broad Institute.
* `--reader-threads 56`. This designates the number of parallel threads for GenomicsDBImport to run across.

7. Run the GenotypeGVCFs tool again on the combined GVCF file.
```
nano run_gatk_genotypegvcfs.qsh
```
```
#!/bin/bash
#SBATCH --job-name=GenotypeGVCFs_trimmed_intervals
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=30G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=condo-ut-genomics
#SBATCH --qos=genomics
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --time=144:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load gatk

gatk \
   --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=2" \
    GenotypeGVCFs \
   -R /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/03_C.obtusa_reference-aligned/C.obtusa_reference/COB_r1.0.fasta  \
   -V gendb://Cthyoides_Cobtusa-ref_database \
   --intervals trimmed_intervals.list \
   -O Cthyoides_Cobtusa-refaligned_intervals.vcf.gz
```
```
sbatch run_gatk_genotypegvcfs.qsh
```

</details>
</p>

***

</details>

### Distant Reference
<details><summary> Thuja plicata - The second closest reference (here, our "phylogenetically distant reference"; Illumina-based from a multigeneration-selfed individual) </summary>

***

<details><summary> 04. BWA: Alignment to the Thuja plicata genome - 7/14/23 </summary>
<p>

Testing the Thuja plicata genome:
1. Make a new directory for this analysis.
* filepath:
```
cd /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/01_Ct_reference-aligned/analysis
mkdir 01_bwa
cd 01_bwa
```

2. Download the Thuja plicata genome.
```
curl --cookie jgi_session=/api/sessions/06136b61cda6b143bb89e9fccea3f35a --output download.20230713.163729.zip -d "{\"ids\":{\"Phytozome-572\":[\"5f2da2419a211ae42a190397\",\"5f2da2419a211ae42a190396\",\"5f2da2419a211ae42a190398\",\"5f2da2409a211ae42a190389\",\"5f2da2409a211ae42a19038b\",\"5f2da2409a211ae42a19038a\",\"5f2da2409a211ae42a19038e\",\"5f2da2409a211ae42a19038f\",\"5f2da2409a211ae42a190390\",\"5f2da2409a211ae42a190393\",\"5f2da2419a211ae42a190395\",\"5f2da2409a211ae42a19038c\",\"5f2da2409a211ae42a190394\",\"5f2da2409a211ae42a19038d\",\"5f2da2409a211ae42a190391\",\"5f2da2409a211ae42a190392\",\"5f6bb3637a4cf8208a33e0d9\"]}}" -H "Content-Type: application/json" https://files.jgi.doe.gov/filedownload/
```

2. Symbolically link all trimmed read files.
```
ln -s ../../../CT_dDocentHPC/analysis/04_fastp/fastp_out/*.gz .
```

3. Unzip the T. plicata genome.
```
unzip download.20230713.163729.zip
```

4. Copy the T. plicata genome to the main working directory.
```
cp /lustre/isaac/proj/UTK0032/zsmith10/CT_thuja_genome/analysis/01_bwa/Phytozome/PhytozomeV13/Tplicata/v3.1/assembly/Tplicata_572_v3.fa.gz .
```

5. Index the T. plicata genome with bwa.
```
module load bwa
bwa index Tplicata_572_v3.fa.gz
```

6. Run a task array on the ISAAC HPC to align all reads files to the Thuja plicata and output merged binary alignment map (BAM) files.
```
#!/bin/bash
#SBATCH --job-name=CT_bwa_array_2
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --mem=30G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --array=181-270
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *r1.fq.gz))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# create an output file name
OUT=$(echo $ARRAY_FILE | sed 's/.r1.fq.gz//')

# define read 1 and read 2 from the array
r1=${ARRAY_FILE}
r2=`sed 's/.r1.fq.gz/.r2.fq.gz/' <(echo ${ARRAY_FILE})`

echo "Read 1: $r1"
echo "Read 2: $r2"


module load bwa
module load samtools

bwa mem -t 5 \
        Tplicata_572_v3.fa.gz \
        $r1 \
        $r2 \
        | samtools view -bSh \
        | samtools sort \
        -@ 10 -m 4G \
        -o $OUT\_sorted.bam

echo Files $r1 and $r2 were aligned by task number $SLURM_ARRAY_TASK_ID on $(date)
```

7. Then, I ran samtools flagstats to look at alignment statistics.
```
nano run_flagstats.qsh
```
```
#!/bin/bash
#SBATCH -J CT_denovo_flagstats
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load samtools

for b in *.bam
do
        stats_out=$( basename $b | sed 's/.bam/.stats/g')

        samtools flagstat $b > $stats_out
done
```
```
sbatch run_flagstats.qsh
```

8. Then I collated the most useful statistics from the flagstats outputs into a single file.
```
nano collate_flagstats.sh
```
```
for file in *stats
do
        echo $file >> flagstat_summary.tsv &&\
        grep '+ 0 in total' $file >> flagstat_summary.tsv &&\
        grep '+ 0 mapped' $file >> flagstat_summary.tsv &&\
        grep '+ 0 primary mapped' $file >> flagstat_summary.tsv &&\
        grep '+ 0 supplementary' $file  >> flagstat_summary.tsv
done
```
```
bash collate_flagstats.sh
```

9. Finally, I curated the stats into a tsv file to copy into a spreadsheet in google drive [here](https://docs.google.com/spreadsheets/d/1KL36GsyMEJyrZzoKnBIFlJDbhraCc4ZqQv5sTebNm6w/edit#gid=250450826). 
```
nano curate_flagstats.sh
```
```
awk 'NR % 5 == 1' flagstat_summary.tsv > column1.tsv
awk 'NR % 5 == 2' flagstat_summary.tsv > column2.tsv
awk 'NR % 5 == 3' flagstat_summary.tsv > column3.tsv
awk 'NR % 5 == 4' flagstat_summary.tsv > column4.tsv
awk 'NR % 5 == 0' flagstat_summary.tsv > column5.tsv

paste column1.tsv column2.tsv column3.tsv column4.tsv column5.tsv > final_flagstat_summary.tsv
```
```
bash curate_flagstats.sh
```



</details>
</p>

***

<details><summary> 05. GATK: HaplotypeCaller </summary>
<p>

1. Create a GATK environment on ISAAC and install picard tools. **It is important to install this conda environment, but to call the program through spack by using the `module load gatk` command in order 1) avoid dependency issues, while 2) still ensuring that gatk is routed through the appropriate compiler on the ISAAC HPC system. I have had odd software conflicts/failures by attempting to run the version of GATK installed in the conda environment, as well as the natively installed GATK version outside of my conda environment. At some point, I will explore this odd interaction further, but this is a reliable fix for the moment.**

```
conda create gatk -c bioconda gatk4
conda install -c bioconda picard
conda activate gatk
```

2. Link the sorted BAMs from the 01_bwa directory in a new 02_gatk directory.
```
cd ..
mkdir 02_gatk
ln -s ../01_bwa/*bam .
```

3. Copy the Thuja plicata reference to this directory.
```
cp ../01_bwa/Tplicata_572_v3.fa.gz .
```

4. Unzip the T. plicata reference genome.
```
gunzip Tplicata_572_v3.fa.gz
```

5. Reindex the genome with samtools.
```
module load samtools
samtools faidx Tplicata_572_v3.fa
picard CreateSequenceDictionary R=Tplicata_572_v3.fa O=Tplicata_572_v3.dict
```

6. Run the following array to add/replace read groups using picard.
```
nano run_picard_array.qsh
```
```
#!/bin/bash
#SBATCH --job-name=CT_GATK1
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --array=0-90
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *_sorted.bam))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# create an output file name
OUT=`sed 's/_sorted.bam//' <(echo ${ARRAY_FILE})`

# define the bam file from the array
BAM=${ARRAY_FILE}

echo "BAM $BAM"
echo "OUT $OUT"

picard \
        AddOrReplaceReadGroups \
        I=${BAM} \
        O=${OUT}_sorted.RG.bam \
        RGSM=$OUT \
        RGLB=$OUT \
        RGPL=illumina \
        RGPU=$OUT

echo BAM file $BAM had reads groups added/replaced by picard-2.27.3 in Slurm Task ID $SLURM_ARRAY_TASK_ID on $(date).
```
```
sbatch run_picard_array.qsh
```
* Note 1: the array numbers were run in batches of 90 (Run 1: 0-90, Run 2: 91-180, Run 3: 181-269)

7. Index the sorted bam files using samtools.
```
nano run_samtools_index.qsh
```
```
#!/bin/bash
#SBATCH --job-name=samtools_index
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=6:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load samtools

for file in *_sorted.RG.bam; do samtools index $file; done
```
```
sbatch run_samtools_index.qsh
```

8. Run the following gatk HaplotypeCaller array to create `g.vcf` files.
```
nano run_haplotypecaller_array.qsh
```
```
#!/bin/bash
#SBATCH --job-name=CT_HaplotypeCaller_array4
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=10G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --array=261-270
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *_sorted.RG.bam))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# define the BAM file from the array file
BAM=${ARRAY_FILE}

# create an output file name
OUT=$( basename $BAM | sed 's/_sorted.RG.bam//g' )


echo "BAM: $BAM"
echo "OUT: $OUT"

module load gatk

gatk HaplotypeCaller \
    -R /lustre/isaac/proj/UTK0032/zsmith10/CT_thuja_genome/analysis/02_gatk/Tplicata_572_v3.fa \
    -I $BAM \
    -O ${OUT}.g.vcf \
    --native-pair-hmm-threads 10 \
    -ERC GVCF \

echo BAM file $BAM was processed by the GATK Haplotype Caller by Slurm Task ID $SLURM_ARRAY_TASK_ID on $(date)
```
```
sbatch run_haplotypecaller_array.qsh
```

9. To ensure the HaplotypeCaller jobs are completing, you can manually check the `.err` files with a text editor, which shows the entire task log, or use the following shortcut and ensure that your actual number of completed jobs matches with your expected number of completed jobs.
```
grep 'Traversal complete' *err | wc -l
```

10. Move array `.err` and `.out` files to a new directory to reduce clutter.
```
mkdir gatk_haplotypecaller_outs
mv *err *out gatk_haplotypecaller_outs
```

11. Make a new directory and move samples that should be excluded from the study due to low read counts into it.
```
mkdir excluded_samples
mv FL1_004* excluded_samples/
mv MS-L4_006* excluded_samples/
mv NJ-CAM_016* excluded_samples/
```

</details>
</p>

***

<details><summary> 06. GATK: GenomicsDBImport & GenotypeGVCFs </summary>
<p>

### 06. GATK:  GenomicsDBImport & GenotypeGVCFs

1. Make a new directory in the analysis directory to house this analysis.
```
cd ..
mkdir 03_gatk_genotype_gvcfs
cd 03_gatk_genotype_gvcfs
```

2. Link the `g.vcf` files, and the reference genome into the new directory.
```
ln -s ../02_gatk/*g.vcf* .
ln -s ../02_gatk/Tplicata_572_v3.* .
```

3. Create a list of all `g.vcf` files to feed into the CombineGVCFs tool.
```
ls *g.vcf > gvcfs.list
```

4. Next, I ran the CombineGVCFs tool (which failed--see below).
```
nano run_combine_gvcfs.qsh
```
```
#!/bin/bash
#SBATCH --job-name=CT_CombineGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks=50
#SBATCH --mem=150G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus-bigmem
#SBATCH --qos=campus-bigmem
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load gatk

gatk CombineGVCFs \
        --java-options "-Xmx140g" \
        -R Tplicata_572_v3.fa \
        -variant gvcfs.list \
        -O Cthyoides_ref2.g.vcf.gz
```
```
sbatch run_combine_gvcfs.qsh
```

5. Unfortunately, the CombineGVCFs tool did not perform well on this dataset, so I attempted to use GenomicsDBImport, another GATK module. First I prepared an intervals file to feed into the tool. The following script will print the name of each scaffold in the assembly, remove the `>` in front of the scaffold name, as well as remove any genomic coordinate data, such as ">29382484 **1871870 0 29346470+,...,29353346+**". Failure to remove these features will result in failure of the GenomicsDBImport tool as GATK's GenotypeGVCFs tool does not include any genomic coordinate data in the resulting `g.vcf` files--only the scaffold name, in the above case: **29382484**.
* Note: This file must end in a `.list` suffix to be properly read by GenomicsDBImport. 

```
nano create_intervals_file.qsh
```
```
# Create the intervals.list file.
grep '>' Tplicata_572_v3.fa | sed 's/>//' >> intervals.list

# Replace 'input.txt' with the path to your input text file.
input_file="intervals.list"

while IFS= read -r line; do
  # Surround each line with backticks and print the result
  echo "$line" | cut -d' ' -f1 >> trimmed_intervals.list
done < "$input_file"

rm intervals.list
```
```
bash create_intervals_file.qsh
```

6. Now, run the GenomicsImportDB tool using the following script. The following tool requires quite a bit of nuance to run appropriately, and it does not produce a simple combined GVCF, instead it creates a database file that includes all of joint-call data.  
```
nano run_genomicsdb_import.qsh
```
```
#!/bin/bash
#SBATCH --job-name=CT_GenomicsDBImport
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --mem=500G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus-bigmem
#SBATCH --qos=campus-bigmem
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load gatk

gatk \
        --java-options "-Xmx460g -Xms460g" \
        GenomicsDBImport \
        --genomicsdb-shared-posixfs-optimizations \
        --reference Tplicata_572_v3.fa \
        --variant gvcfs.list \
        --genomicsdb-workspace-path CT_database \
        --overwrite-existing-genomicsdb-workspace \
        --intervals trimmed_intervals.list \
        --batch-size 50 \
        --merge-contigs-into-num-partitions 25 \
        --reader-threads 56
```
```
sbatch run_genomicsdb_import.qsh
```
* `--java-options "-Xmx460g -Xms460g"` Sets the heap memory for java to 460 Gb (~10% less than the total requested on ISAAC).
* `GenomicsDBImport` specifies the GATK tool.
* `--genomicsdb-shared-posixfs-optimizations` toggles optimizations for HPC systems, such as those run on lustre (like ISAAC).
* `--reference Tplicata_572_v3.fa` designates the reference sequence.
* `--variant gvcfs.list` designates the list of gvcf files.
* `--genomicsdb-workspace-path CT_database` designates the workspace output, aka the genomicsdb database directory. This will be your input file for GenotypeGVCFs.
* `--overwrite-existing-genomicsdb-workspace`. This tool fails if the database already exists. I enabled this due to several failed trial runs due to forgetting to remove the previous database file. Use at your own risk.
* `--intervals trimmed_intervals.list` designates the genomic intervals (scaffolds, in my case) in the reference to read over.
* `--batch-size 50` sets the batch size of files to process. 50 is the recommended value used in-house at the Broad Institute (GATK's developer). Using the default number (all files) caused me to exceed the SLURM walltime. More info here: https://gatk.broadinstitute.org/hc/en-us/articles/360056138571-GenomicsDBImport-usage-and-performance-guidelines
* `--merge-contigs-into-num-partitions 25`. This flag merges smaller contigs (but supposedly not larger ones) to improve speed and efficiency of the tool. This setting was recommended by the Broad Institute.
* `--reader-threads 56`. This designates the number of parallel threads for GenomicsDBImport to run across.

6.5. As a bonus, here's the slurm resource utilization from the above script. It appears that less resources were required than anticipated!
```
Job ID: 734083
Cluster: isaac
User/Group: zsmith10/tug2106
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 56
CPU Utilized: 07:48:26
CPU Efficiency: 1.87% of 17-09:26:00 core-walltime
Job Wall-clock time: 07:27:15
Memory Utilized: 157.33 GB
Memory Efficiency: 31.47% of 500.00 GB
```

7. Run the GenotypeGVCFs tool again on the combined GVCF file.
```
nano run_gatk_genotypegvcfs.qsh
```
```
#!/bin/bash
#SBATCH --job-name=CT_GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --mem=192G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

module load gatk

gatk \
   --java-options "-Xmx180g" \
    GenotypeGVCFs \
   -R Tplicata_572_v3.fa  \
   -V gendb://CT_database \
   -O Cthyoides_refaligned.vcf.gz
```
```
sbatch run_gatk_genotypegvcfs.qsh
```

</details>
</p>

***


</details>


***

<details><summary> 01. HDplot: Paralog Pruning </summary>
<p>

Documentation: https://github.com/edgardomortiz/paralog-finder

1. First, I downloaded the master.zip from the paralog-finder Github. This method is based on McKinney et al., 2017, but generally improved. After unzipping the file, it returns a paralog-finder-master directory with all required scripts.
```
wget https://github.com/edgardomortiz/paralog-finder/archive/refs/heads/master.zip
unzip master.zip
```

2. Next, I linked my VCFs with GATK filters applied. My naming wasn't perfectly consistent in the previous SNP calling analyses, so I've renamed them to be consistent moving forward here.
```
nano link_files.sh
```
```
ln -s ../../03_C.obtusa_reference-aligned/analysis/07_bcftools/Cthyoides_Cobtusa-refaligned_intervals_GATKfilters.vcf.gz .
ln -s ../../02_Ct_denovo_dDocent2GATK/02_analysis_C0.96_k3_k8_current/09_bcftools/Cthyoides_denovo-refaligned_GATKfilters.vcf.gz .
ln -s ../../01_T.plicata_reference-aligned/analysis/04_bcftools/Cthyoides_refaligned.GATKfilters.vcf.gz .

mv Cthyoides_Cobtusa-refaligned_intervals_GATKfilters.vcf.gz Cthyoides_Cobtusa-refaligned_GATKfilters.vcf.gz
mv Cthyoides_refaligned.GATKfilters.vcf.gz Cthyoides_Tplicata-refaligned_GATKfilters.vcf.gz
```
```
bash link_files.sh
```

3. I then created a conda environment for HDplot.
```
conda create -n paralog-finder -c conda-forge \
        python=2.7 \
        scipy \
        pandas \
        numpy \
        statsmodels
conda install r-ggplot2
conda install r-ggExtra
conda install r-argparse
```

4. I also need to add a few packages to my R environment, so I did that like so.
```
conda install r-ggExtra
conda install r-argparse
conda install r-ggplot2
```

5. The original paralog-finder python implementation of HDplot was unfortunately not compatible with GATK VCF files--it is compatible with STACKS and ipyrad only. Thus, I adapted these scripts to be compatible with GATK and downstream BCFtools SNP filtering. A copy of the modified paralog-finder scripts are included below.

* Modified paralog-finder scripts: [paralog-finder-modified.tar.gz](https://github.com/statonlab/StatonLabDocs/files/14216569/paralog-finder-modified.tar.gz)

6. Before running paralog-finder, I realized I needed to annotate my VCFs with loci names, so I did this with BCFtools. I named each locus id as the scaffolds name and its corresponding bp position.
```
nano run_bcftools_annotate.qsh
```
```
#!/bin/bash
#SBATCH --job-name=BCFtools_annotate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu
#SBATCH --array=0-2
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *vcf.gz))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# define the vcf.gz file from the array file
IN=${ARRAY_FILE}

# create an output file name
OUT=$(basename $IN | sed 's/.vcf.gz/_annotated.vcf.gz/' )

echo "IN: $IN"
echo "OUT: $OUT"

eval "$(conda shell.bash hook)"
conda activate bcftools

bcftools annotate --set-id '%CHROM\_%POS' -o $OUT -O z $IN
```
```
sbatch run_bcftools_annotate.qsh
```

7. Paralog-finder is comprised of three scripts. The first identifies paralogs and prints them into a `.depthsBias` file.
```
nano run_paralog-finder_identify-paralogs.qsh
```
```
#!/bin/bash
#SBATCH --job-name=paralog-finder
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=40G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edux

# Load R environment
eval "$(conda shell.bash hook)"
conda activate paralog-finder

# Step 1: Run modified HDplot script for GATK-formatted VCF
for file in *_annotated*vcf*
do
        # create an out directory
        out_dir=$(basename $file | sed 's/_GATKfilters_annotated.vcf.gz/.paralog-finder.out/')
        out_file=$(basename $file | sed 's/.vcf.gz/.depthsBias/')
        mkdir $out_dir

        python paralog-finder-modified/HDplot_process_vcf_modified.py \
                -i $file
        mv $out_file $out_dir
done
```
```
sbatch run_paralog-finder_identify-paralogs.qsh
```

8. Once the depth of each locus has been calculated, we can plot them with the second script. 
```
nano run_paralog-finder_plot-paralogs.qsh
```
```
#!/bin/bash
#SBATCH --job-name=paralog-finder
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=40G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edux

# Step 2: Create HDplot plots
eval "$(conda shell.bash hook)"
conda activate R

Cobtusa=Cthyoides_Cobtusa-refaligned_GATKfilters_annotated.depthsBias
Tplicata=Cthyoides_Tplicata-refaligned_GATKfilters_annotated.depthsBias
denovo=Cthyoides_denovo-refaligned_GATKfilters_annotated.depthsBias

# Step 3: Create SNP blacklists
# create an out directory
Cobtusa_out_dir=$(basename $Cobtusa | sed 's/_GATKfilters_annotated.depthsBias/.paralog-finder.out/')
# create an out directory
Tplicata_out_dir=$(basename $Tplicata | sed 's/_GATKfilters_annotated.depthsBias/.paralog-finder.out/')
# create an out directory
denovo_out_dir=$(basename $denovo | sed 's/_GATKfilters_annotated.depthsBias/.paralog-finder.out/')


# Run HDplot graphing
Rscript paralog-finder-modified/HDplot_graphs_modified.R \
        -i $Cobtusa_out_dir/$Cobtusa \
        --transp 0.5
#move pngs to out directory
mv *png $Cobtusa_out_dir


# Run HDplot graphing
Rscript paralog-finder-modified/HDplot_graphs_modified.R \
        -i $Tplicata_out_dir/$Tplicata \
        --transp 0.5
#move pngs to out directory
mv *png $Tplicata_out_dir


# Run HDplot graphing
Rscript paralog-finder-modified/HDplot_graphs_modified.R \
        -i $denovo_out_dir/$denovo \
        --transp 0.5
#move pngs to out directory
mv *png $denovo_out_dir
```
```
sbatch run_paralog-finder_plot-paralogs.qsh
```

9. This script generates four plots:
* Plot 1: 
* Plot 2:
* Plot 3:
* Plot 4:

10. Finally, after selecting cut-off thresholds based on the heterozygosity and read deviation (primarily from Plot 1), I created the blacklists and whitelists of loci for each file with the third and final script. While paralog-finder can help identify both diverged and undiverged paralogs (ancient and recent paralogs, respectively). I did not observe clear signs of undiverged paralogs, thus, I selected a threshold of maxH 0.99 to filter only those paralogs based on heterozygosity. Read deviation thresholds were selected for each of the three datasets by excluding the upper and lower 2.5th percentile of loci (i.e., outliers) denoted by thickened black lines/points on graph 1. 

```
nano run_paralog-finder_blacklist-paralogs.qsh
```
```
#!/bin/bash
#SBATCH --job-name=paralog-finder
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edux

# Load R environment
eval "$(conda shell.bash hook)"
conda activate paralog-finder

# depthBias files
Cobtusa=Cthyoides_Cobtusa-refaligned_GATKfilters_annotated.depthsBias
Tplicata=Cthyoides_Tplicata-refaligned_GATKfilters_annotated.depthsBias
denovo=Cthyoides_denovo-refaligned_GATKfilters_annotated.depthsBias

# Step 3: Create SNP blacklists

# create an out directory
Cobtusa_out_dir=$(basename $Cobtusa | sed 's/_GATKfilters_annotated.depthsBias/.paralog-finder.out/')

python paralog-finder-modified/blacklist_paralogs_modified.py \
                -i $Cobtusa_out_dir/$Cobtusa \
                --maxH 0.99 \
                --minN 1 \
                --minD -29 \
                --maxD 43

#move blacklists/whitelists to out directory
mv *.*list $Cobtusa_out_dir

# create an out directory
Tplicata_out_dir=$(basename $Tplicata | sed 's/_GATKfilters_annotated.depthsBias/.paralog-finder.out/')

python paralog-finder-modified/blacklist_paralogs_modified.py \
                -i $Tplicata_out_dir/$Tplicata \
                --maxH 0.99 \
                --minN 1 \
                --minD -20 \
                --maxD 31

#move blacklists/whitelists to out directory
mv *.*list $Tplicata_out_dir

# create an out directory
denovo_out_dir=$(basename $denovo | sed 's/_GATKfilters_annotated.depthsBias/.paralog-finder.out/')

python paralog-finder-modified/blacklist_paralogs_modified.py \
                -i $denovo_out_dir/$denovo \
                --maxH 0.99 \
                --minN 1 \
                --minD -18 \
                --maxD 200

#move blacklists/whitelists to out directory
mv *.*list $denovo_out_dir
```
```
bash run_paralog-finder_blacklist-paralogs.qsh
```

</p>
</details> 

***

<details><summary> 02. BCFtools: Filtering </summary>
<p>

Directory: Directory: /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/02_bcftools/01_first_pass

1. First, I linked the relevant GATK-marked VCFs and their corresponding blacklists from my hdplot directory. 
```
ln -s ../01_hdplot/*_GATKfilters_annotated.vcf.gz .
cp ../01_hdplot/Cthyoides_*-refaligned.paralog-finder.out/*blacklist .
```

2. Next, I ran several incremental filters on my VCFs, including one that retains all high-quality loci per individual to generate summary statistics.

* There is a mistake here that is fixed in my second pass at this. The "unfiltered" option should actually just be filtered for everything except MAF. Thus, I would recommend screening multiple levels of missing data without minor allele frequency filtering also. I only did the single threshold of 90% on my second pass as I was already quite sure the missing data percentages don't affect the dimensionality of the data. 

```
nano run_bcftools_filter.qsh
```
```
#!/bin/bash
#SBATCH --job-name=BCFtools_filter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edux
#SBATCH --array=0-2
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate bcftools

# create an array variable containing the file names
FILES=($(ls -1 *_GATKfilters_annotated.vcf.gz))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}
OUTPUT_FILE_UNFILTERED=$(basename $ARRAY_FILE | sed 's/.vcf.gz/_minDP3_PASSonly_removed-indv_removed-paralogs.vcf/')
OUTPUT_FILE_FILTERED=$(basename $ARRAY_FILE | sed 's/.vcf.gz/_biallelic_maf05_minDP3_PASSonly_removed-indv_removed-paralogs/')
BLACKLIST=$(basename $ARRAY_FILE | sed 's/_GATKfilters_annotated.vcf.gz/_paralogs.blacklist/')
WHITELIST=$(basename $ARRAY_FILE | sed 's/_GATKfilters_annotated.vcf.gz/_singletons.whitelist/')

echo "IN: $ARRAY_FILE"
echo "UNFILTERED_OUT: $OUTPUT_FILE_UNFILTERED"
echo "FILTERED_OUT: $OUTPUT_FILE_FILTERED"
echo "BLACKLIST: $BLACKLIST"

# UNFILTERED for summary statistics
bcftools view -i "ID!=@$BLACKLIST" $ARRAY_FILE | \
bcftools view -S ^excluded_individuals.txt | \
bcftools filter -i 'TYPE="snp" & FILTER="PASS" & FORMAT/DP>=3' -o $OUTPUT_FILE_UNFILTERED \
&& bgzip $OUTPUT_FILE_UNFILTERED \
&& bcftools index $OUTPUT_FILE_UNFILTERED\.gz

# FILTERED for population analyses - 80%
bcftools view -i "ID!=@$BLACKLIST" $ARRAY_FILE | \
bcftools view -m2 -M2 -S ^excluded_individuals.txt | \
bcftools filter -i 'TYPE="snp" & MAF>=0.05 & F_MISSING<=0.20 & FILTER="PASS" & FORMAT/DP>=3' -o ${OUTPUT_FILE_FILTERED}_mm80.vcf \
&& bgzip ${OUTPUT_FILE_FILTERED}_mm80.vcf \
&& bcftools index ${OUTPUT_FILE_FILTERED}_mm80.vcf.gz

# FILTERED for population analyses - 85%
bcftools view -i "ID!=@$BLACKLIST" $ARRAY_FILE | \
bcftools view -m2 -M2 -S ^excluded_individuals.txt | \
bcftools filter -i 'TYPE="snp" & MAF>=0.05 & F_MISSING<=0.15 & FILTER="PASS" & FORMAT/DP>=3' -o ${OUTPUT_FILE_FILTERED}_mm85.vcf \
&& bgzip ${OUTPUT_FILE_FILTERED}_mm85.vcf \
&& bcftools index ${OUTPUT_FILE_FILTERED}_mm85.vcf.gz

# FILTERED for population analyses - 90%
bcftools view -i "ID!=@$BLACKLIST" $ARRAY_FILE | \
bcftools view -m2 -M2 -S ^excluded_individuals.txt | \
bcftools filter -i 'TYPE="snp" & MAF>=0.05 & F_MISSING<=0.10 & FILTER="PASS" & FORMAT/DP>=3' -o ${OUTPUT_FILE_FILTERED}_mm90.vcf \
&& bgzip ${OUTPUT_FILE_FILTERED}_mm90.vcf \
&& bcftools index ${OUTPUT_FILE_FILTERED}_mm90.vcf.gz

# FILTERED for population analyses - 95%
bcftools view -i "ID!=@$BLACKLIST" $ARRAY_FILE | \
bcftools view -m2 -M2 -S ^excluded_individuals.txt | \
bcftools filter -i 'TYPE="snp" & MAF>=0.05 & F_MISSING<=0.05 & FILTER="PASS" & FORMAT/DP>=3' -o ${OUTPUT_FILE_FILTERED}_mm95.vcf \
&& bgzip ${OUTPUT_FILE_FILTERED}_mm95.vcf \
&& bcftools index ${OUTPUT_FILE_FILTERED}_mm95.vcf.gz

```
```
sbatch run_bcftools_filter.qsh
```

3. Just to check that I have the expected number of individuals in each file, I ran the following script.
```
for file in *vcf.gz; do echo $file; bcftools query -l $file | wc -l; done
```

4. Finally, to ensure that the blacklisted loci were removed, I ran a loop to grep for each line of the mm80 files filtered blacklist on its corresponding filtered maximum of 20% missing data (mm80) VCF file.
```
nano check_paralogs.sh
```
```
#!/bin/bash
#SBATCH --job-name=count-paralogs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=1:00:00

# Denovo
while read -r line; do
    if grep -q "$line" Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_maf05_mm80_minDP3_PASSonly_removed-i$
        echo "Found: $line" >> denovo_paralog_check.txt
    else
        echo "Not found: $line" >> denovo_paralogs.txt
    fi
done < Cthyoides_denovo-refaligned_paralogs.blacklist


#Tplicata
while read -r line; do
    if grep -q "$line" Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_maf05_mm80_minDP3_PASSonly_removed$
        echo "Found: $line" >> Tplicata_paralogs.txt
    else
        echo "Not found: $line" >> Tplicata_paralog_check.txt
    fi
done < Cthyoides_Tplicata-refaligned_paralogs.blacklist


#Cobtusa
while read -r line; do
    if grep -q "$line" Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf05_mm80_minDP3_PASSonly_removed-$
        echo "Found: $line" >> Cobtusa_paralogs.txt
    else
        echo "Not found: $line" >> Cobtusa_paralog_check.txt
    fi
done < Cthyoides_Cobtusa-refaligned_paralogs.blacklist

for file in *_paralog_check.txt
do
        echo $file >> paralog_validation_report.txt
        grep -c 'Found:' $file >> paralog_validation_report.txt
done
```
```
sbatch check_paralogs.sh
```

</p>
</details> 

***

<details><summary> 03. Plink: Linkage Disequilibrium Pruning </summary>
<p>

Directory: `/lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/03_LD_pruning`

1. Link the filtered VCFs from each respective SNP calling pipeline. Note: I've copied these VCFs to make it more simplistic for linking through the following analyses.
```
ln -s ../02_bcftools/*vcf.gz .
```

2. To prune for linkage disequilibrium, I used plink. I pruned in windows of 50 SNPs, shifting the window forward 5 SNPs after each iteration, and removed any SNPs with an r-squared higher than 0.4 (see, https://onlinelibrary.wiley.com/doi/full/10.1111/mec.15407?saml_referrer). It is very important that you have annotated your loci with unique names (which I did with BCFtools), otherwise GATK assigns your loci as ".", which will result in all loci being filtered out.
```
nano run_plink.qsh
```
```
#!/bin/bash
#SBATCH --job-name=plink_LD_prune
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

eval "$(conda shell.bash hook)"
conda activate plink

for input_file in *vcf.gz
do
        base_file=$(basename $input_file | sed 's/.vcf.gz//')
        output_file=$(basename $input_file | sed 's/.vcf.gz/_LD-50-5-0.4/')

        plink \
                --vcf $input_file \
                --make-founders \
                --indep-pairwise  50 5 0.4 \
                --allow-extra-chr \
                --out $base_file \
        && \
        plink \
                --vcf $input_file \
                --exclude $base_file\.prune.out \
                --recode vcf bgz \
                --allow-extra-chr \
                --out $output_file
done
```
```
sbatch run_plink.qsh
```

3. I counted the SNPs in each file to get a rough estimate of how many SNPs were lost through plink's sliding-window-based LD pruning.
```
for file in *LD-50-5-0.4.vcf.gz
do
    echo $file >> SNPcounts.txt
    zgrep -c -v '^#' $file >> SNPcounts.txt
done
```

4. I included all of the SNP counts for each of my VCF filtering steps here: 

* https://docs.google.com/spreadsheets/d/1LyT2Z5pgxM3Lfv-SsWOXsxO3wW77mOvVqHmW6vG5NbE/edit#gid=1861693885

</p>
</details> 

***

<details><summary> 04. King: Evaluate Pairwise Kinship </summary>
<p>

1. Before using King, it's necessary to filter and prune to ensure a high-quality basis for inferring kinship. Thus, I linked the VCF files I intended to carry forward, which used a missing data cut-off of 10%. These files have paralogs removed and are linkage pruned.
```
ln -s ../03_LD_pruning/*_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz .
```

2. Next, I ran king as implemented in plink2. Notably, this analysis requires the VCF be converted to BED format, which is included below in the first command.
```
nano run_plink2_king.qsh
```
```
#!/bin/bash
#SBATCH --job-name=plink_king
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

eval "$(conda shell.bash hook)"
conda activate plink2

for input_file in *vcf.gz
do
        bed_file=$(basename $input_file | sed 's/.vcf.gz//')
        output_file=$(basename $input_file | sed 's/.vcf.gz/.king/')

        plink2 \
                --vcf $input_file \
                --make-bed \
                --allow-extra-chr \
                --out $bed_file \
        && \
        plink2 \
                --bfile $bed_file \
                --make-king-table \
                #--king-cutoff 0.30 \
                --allow-extra-chr \
                --out $output_file
done
```
```
sbatch run_plink2_king.qsh
```

3. This runs lightning fast, so results are pretty quick to assess. To easily evaluate the output pairwise kinship files (`.kin0` files), I used visidata to sort for the maximum kinship values. I'm including the output for the C. obtusa aligned SNPs below, but I repeated the visualization for all 3 of the kinship files. Plink's documentation suggests the following kinship cut-offs:

* Duplicate samples/Monozygotic Twins: ~0.354 (the geometric mean of 0.5 and 0.25)
* 1st Degree Relations (Parent-Child or Full Siblings): ~0.25
* 2nd Degree Relations (Half-siblings): ~0.125

I suspect the kinship values are slightly inflated in this dataset, due to observed heterozygosity being low. The algorithm also evaluates Isolation-by-Descent (IBD) and Isolation-by-State (IBS). **Of note, the column IBS0 is actually the number of loci that have an IBS of 0, thus, lower number indicate higher IBS.**
* Isolation-by-Descent (IBD): Identical alleles among individuals inherited from a common ancestor.
* Isolation-by-State (IBS): Identical alleles among individuals that may be identical by chance.

Here's a few relative values according to ChatGPT:
IBD (Identical By Descent):
1. 0% IBD: Individuals are unrelated or very distantly related.
2. 25% IBD: Suggests half-sibling relationship or grandparent-grandchild relationship.
3. 50% IBD: Suggests full-sibling relationship or parent-offspring relationship.
4. 100% IBD: Individuals are identical twins or clones.
IBS (Identical By State):
1. 0% IBS: Individuals do not share any alleles at the analyzed loci.
2. 100% IBS: Individuals are identical at all analyzed loci.

```
vd Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf05_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.king.kin0
```
* File excerpt:
![image](https://github.com/statonlab/StatonLabDocs/assets/115577973/66abe3bd-dbf1-49a7-b895-a3fbf8763a76)

4. Based on these results, earlier PCA/structure evaluation of this dataset was seems to be validated in that the the samples from the DE (now VA) population and MS4_1, MS4_2, MS4_3, MS4_4, MS4_6, and MS4_7 appear to be more related than half-siblings.

5. The previous command also produces pruning files, so I reran the command with the `--king-cutoff 0.30` to select individuals with more than 30% kinship. The individuals listed in the output file `Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf05_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.king.king.cutoff.out.id` will be added to excluded individuals when refiltering with BCFtools.

</p>
</details> 

***

<details><summary> 02.1. BCFtools: Removing Related Individuals </summary>
<p>

Directory: /lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/02_bcftools/02_second_pass

1. I linked my VCFs for summary stats, as well as my 10% maximum missing genotypes VCFs into a new directory.
```
ln -s ../../01_hdplot/*_GATKfilters_annotated.vcf.gz .
cp ../../01_hdplot/Cthyoides_*-refaligned.paralog-finder.out/*blacklist .
cp ../01_first_pass/excluded_individuals.txt
cp ../../04_king/Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf05_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.king.king.cutoff.out.id .
```

2. Next, I added the individuals from `Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf05_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.king.king.cutoff.out.id` to `excluded_individuals.txt`.
* New `excluded_individuals.txt`:

3. Then, I reran BCFtools filtering.

* There was a mistake in my first pass that is fixed here. The "unfiltered" option should actually just be filtered for everything except MAF. Thus, I would recommend screening multiple levels of missing data without minor allele frequency filtering also. I only did the single threshold of 90% on my second pass as I was already quite sure the missing data percentages don't affect the dimensionality of the data. 

```
nano run_bcftools_filter.qsh
```
```

```
```
sbatch run_bcftools_filter.qsh
```

</p>
</details> 

***

<details><summary> 03.1 Linkage Pruning - Second Pass </summary>
<p>

Directory: `/lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/03_LD_pruning/02_second_pass`

1. I linked filtered files from my second pass at pruning (after removing closely related individuals).
```
ln -s ../../02_bcftools/02_second_pass/*_removed-paralogs*.vcf.gz .
```

2. Next, I redid linkage pruning with plink.
```
nano run_plink.qsh
```
```
#!/bin/bash
#SBATCH --job-name=plink_LD_prune
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=3:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

eval "$(conda shell.bash hook)"
conda activate plink

for input_file in *vcf.gz
do
        base_file=$(basename $input_file | sed 's/.vcf.gz//')
        output_file=$(basename $input_file | sed 's/.vcf.gz/_LD-50-5-0.4/')

        plink \
                --vcf $input_file \
                --make-founders \
                --indep-pairwise  50 5 0.4 \
                --allow-extra-chr \
                --out $base_file \
        && \
        plink \
                --vcf $input_file \
                --exclude $base_file\.prune.out \
                --recode vcf bgz \
                --allow-extra-chr \
                --out $output_file
done
```
```
sbatch run_plink.qsh
```

3. Just to get an idea of how the number of SNPs is changing, I counted the SNPs of all the files (pruned and unpruned).
```
nano count_snps.sh
```
```
#!/bin/bash
#SBATCH --job-name=count-paralogs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=5G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=1:00:00

for file in *_mm*vcf.gz
do
        echo $file >> snp_counts.txt
        zgrep -c -v '^#' $file >> snp_counts.txt
done
```
```
sbatch count_snps.sh
```

</p>
</details> 

***


<details><summary> 05. R: Interactive PCAs </summary>
<p>

Directory: `/lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/05_pca/01_first_pass`

01. First, I loaded my conda environment for base R. The following is the recipe I used to create my environment. 
```
conda create -n R -c bioconda R
conda install r-markdown
conda install -c conda-forge pandoc
conda install r-tidyverse
```

2. Next, I linked all of my filtered VCF files to evaluate how filtering affected the data.
```
ln -s ../02_bcftools/*maf05* .
```

3. PCAs were performed using an Rscript and constructed directly from a VCF. The following methodology is robust to mixed ploidy VCFs. The raw PCA script is included in section 5.1.
```
nano run_pca.sh
```
```
#!/bin/bash
#SBATCH --job-name=PCAs
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

# Load R environment
eval "$(conda shell.bash hook)"
conda activate R

# Launch PCA Script
/nfs/home/zsmith10/.conda/envs/vcf2ploidy/lib/R/bin/R --no-restore --file=pcas.R

# Remove temporary GDS files
rm *.gds
```
```
sbatch run_pca.sh
```

</p>
</details> 

***

<details><summary> 05.1 R: Interactive PCAs </summary>
<p>

Directory: `/lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/05_pca/02_second_pass`

1. Once again, I linked over the linkage pruned and unpruned VCFs to visualize on PCAs to check overall data spread.
```
ln -s ../../03_LD_pruning/02_second_pass/*vcf.gz .
```

2. I recreated a popfile with the following file. I manually changed the DE samples to VA as their true provenance is from Virginia.
```
nano make_popfile.sh
```
```
eval "$(conda shell.bash hook)"
conda activate bcftools

input_file=Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_minDP3_PASSonly_removed-indv_removed-paralogs_mm90.vcf.gz

# Create popfile
bcftools query -l $input_file > temp_file
awk '{print $0, substr($1, 1, 2)}' temp_file > popfile.txt && rm temp_file

# Create popfile.csv
cat popfile.txt | sed 's/\s/,/g' > popfile.csv
```
```
sbatch make_popfile.sh
```

3. Next, I reran PCAs.
```
nano run_pca.sh
```
```
#!/bin/bash
#SBATCH --job-name=PCAs
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

# Load R environment
eval "$(conda shell.bash hook)"
conda activate R

# Launch PCA Script
/nfs/home/zsmith10/.conda/envs/vcf2ploidy/lib/R/bin/R --no-restore --file=pcas.R

# Remove temporary GDS files
rm *.gds
```
```
sbatch run_pca.sh
```

* Rscript: `pcas.R`
```
#!/usr/bin/env Rscript

##########################################
#### Installing and loading libraries ####
##########################################

# Define a CRAN mirror to download packages from.
cran_mirrors <- c("https://cloud.r-project.org/")

# Install BiocManager (if required) and load SNPRelate.
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = cran_mirrors)
BiocManager::install("SNPRelate")

# Specify other libraries
packages_to_load <- c("data.table", "plotly", "rgl", "scales", "htmlwidgets", "htmltools", "SNPRelate", "rmarkdown")

# Check if each package is installed, and if not, install it
for (package in packages_to_load) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = cran_mirrors)
  }
}

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

############################
#### Performing the PCA ####
############################

## Loop through each VCF and perform PCAs
# Set the path to your directory containing VCF files
directory_path <- "."
# Get a list of all files in the directory with a .vcf extension
vcf_files <- list.files(directory_path, pattern = "\\.vcf.gz$", full.names = TRUE)

# Loop through each VCF file
for (vcf_file in vcf_files) {

        # 1. Convert VCF to gds
        ## 2a. Create a timestamp
        timestamp <- format(Sys.time(), "%H%M%S")
        # Dynamically generate the output filename with a timestamp
        gds_name <- paste0("temp_", timestamp, ".gds")
        ## 2b. Write to gds
        snpgdsVCF2GDS(vcf_file, gds_name, method="biallelic.only")

        # 2. Assign population file & name of output file to objects
        pops_file <- "popfile.txt"
        input_filename <- basename(vcf_file)
        output_name <- paste(input_filename, "_PCA", sep = "")

        # 3. Open GDS file
        genofile <- snpgdsOpen(gds_name)
    
        # 4. Run PCA
        pca <- snpgdsPCA(genofile, num.thread=1, autosome.only=F)

        # 5. Extract eigenvalues
        eigenval_1 <- pca$eigenval[1]
        eigenval_2 <- pca$eigenval[2]
        eigenval_3 <- pca$eigenval[3]
        #eigenval_4 <- pca$eigenval[4]
        #eigenval_5 <- pca$eigenval[5]
        pc.percent <- pca$varprop * 100
        print(round(pc.percent, 2))

        # 6. Open figure driver
        pdf(paste(output_name, ".pdf", sep=""))

        # 7.Create colors for PCA
        ## 7a. Identify samples and population
        sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
        pop_code <- read.table(pops_file, sep=" ")
        sorted_pops <- pop_code$V2[order(match(pop_code$V1, sample.id))]

        print(length(pca$sample.id))
        print(length(sorted_pops))
        print(nrow(pca$eigenvect))

        ## 7b. Create color vector
        site_colors <- c("#88CCEE", "#332288", "#882255", "#CC6677", "#AA4499", "#661100", "#DDCC77", "#44AA99", "#6699CC", "#117733", "#999933", "#888888")
        ## 7c. Assign site colors to a hue palatte
        site_colors <- hue_pal()(12)
        ## 7d. Format populations from a population file as factors
        site_to_color <- factor(sorted_pops, levels = unique(sorted_pops))
        ## 7e. Assign color levels for each site to the site color column
        site_col <- site_colors[site_to_color]

        # 8. Plots PCA
        ## IF MISSING POPFILE
        if (!is.na(pops_file)) {
                #sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
                #pop_code <- read.table(pops_file, sep=",")
                #sorted_pops <- pop_code$V2#[order(match(pop_code$V1, sample.id))]
                tab <- data.frame(sample.id = pca$sample.id, pop = sorted_pops, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3 = pca$eigenvect[,3], stringsAsFactors=F)
                #save(tab, file=paste(output_name, ".Rdata", sep=""))
        p <- plot_ly(tab, x=tab$EV1, y=tab$EV2, text=tab$sample.id, color=tab$pop, colors=site_colors, type = "scatter", mode = "markers")
        p <- layout(p, title="PCA", xaxis=list(title=paste("PC 1(", round(eigenval_1, d=2) , "%)")), yaxis=list(title=paste("PC 1(", round(eigenval_2, d=2) , "%)")))
        htmlwidgets::saveWidget(as_widget(p), paste(output_name, ".html", sep=""))

        ## IF POPFILE PRESENT
        } else {
                tab <- data.frame(sample.id = pca$sample.id, EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], EV3 = pca$eigenvect[,3], stringsAsFactors=F)
                #print(pca$sample.id)
                #print(tab$EV1)
                #print(tab$EV2)
                p <- plot_ly(tab, x=tab$EV1, y=tab$EV2,
                text=tab$sample.id, type = "scatter", mode = "markers")
                p <- layout(p, title="PCA", xaxis=list(title=paste("PC 1(", round(eigenval_1, d=2) , "%)"), yaxis=list(title=paste("PC 2(", round(eigenval_2, d=2) , "%)"))))
        }
        # 9. Save widget to object & write interactive PCA widget to html
        p <- htmlwidgets::appendContent(p, htmltools::tags$input(id='inputText', value='', ''), htmltools::tags$button(id='buttonSearch', 'Search'))
        p <- htmlwidgets::appendContent(p, htmltools::tags$script(HTML('document.getElementById("buttonSearch").addEventListener("click", function()
        { var i = 0; var j = 0; var found = []; var myDiv = document.getElementsByClassName("js-plotly-plot")[0] var data = JSON.parse(document.querySelectorAll("script[type=\'application/json\']")[0].innerHTML); console.log(data.x.data)
        for (i = 0 ;i < data.x.data.length; i += 1) {
        for (j = 0; j < data.x.data[i].text.length; j += 1) { if (data.x.data[i].text[j].indexOf(document.getElementById("inputText").value) !== -1) { found.push({curveNumber: i, pointNumber: j});
                                        }
                                }
                        }
                        Plotly.Fx.hover(myDiv, found);
                }
        );')))
        htmlwidgets::saveWidget(as_widget(p), paste(output_name, ".html", sep=""))

        # 10. remove temporary gds file
        #file.remove("/lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/02_pca/temp.gds")
}

```

4. All of the PCAs look remarkably consistent among all three SNP calling methods. The single individual left from the related cluster now appropriately clusters with other individuals from Mississippi.

</p>
</details> 

***

<details><summary> 06. Summary Statistics </summary>
<p>

Directory: `/lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/06_summary_stats`

1. To begin, I created a conda environment for base R.
```
conda create -n R -c bioconda R
conda install r-markdown
conda install -c conda-forge pandoc
conda install r-tidyverse # tentative
conda install -c conda-forge gdal # tentative
```

2. First, I linked the VCF files that were NOT linkage pruned as we want to retain all of that information for calculating summary statistics.
```
ln -s ../02_bcftools/02_second_pass/*_biallelic_minDP3_PASSonly_removed-indv_removed-paralogs_mm90.vcf.gz .
```

3. Next, I created a population file.
```
nano make_popfile.qsh
```
```
eval "$(conda shell.bash hook)"
conda activate bcftools

# Check to make sure all files have the same number of individuals
for file in *vcf.gz
do
        bcftools query -l Cthyoides_Cobtusa-refaligned_intervals_GATKfilters_biallelic_maf05_mm95_minDP3_PASSonly_removed-indv_annotated.vcf.gz | wc -l
done

# Create popfile
bcftools query -l Cthyoides_Cobtusa-refaligned_intervals_GATKfilters_biallelic_maf05_mm95_minDP3_PASSonly_removed-indv_annotated.vcf.gz > temp_file
awk '{print $0, substr($1, 1, 2)}' temp_file > popfile.txt && rm temp_file
#awk '{print $0, substr($1, 1, 2)}' bcftools_popfile.txt > temp_file && mv temp_file popfile.txt

# Create popfile.csv
cat popfile.txt | sed 's/\s/,/g' > popfile.csv
```
```
bash make_popfile.qsh
```

3. Next, I ran basic summary stats using hierfstat.
```
nano run_summary_stats.qsh
```
```
#!/bin/bash
#SBATCH --job-name=summary_stats
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=500G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus-bigmem
#SBATCH --qos=campus-bigmem
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

# Load R environment
eval "$(conda shell.bash hook)"
conda activate R

# Launch PCA Script
/nfs/home/zsmith10/.conda/envs/R/lib/R/bin/R --no-restore --file=summary_stats.R
```
```
sbatch run_summary_stats.qsh
```

* Rscript: `summary_stats.R`
```
<insert after troubleshooting>
```

</p>
</details> 

***

<details><summary> 07. PGDspider: VCF to STRUCTURE </summary>
<p>

1. I linked the 10% missing, LD-pruned VCFs in this directory.
```
ln -s ../03_LD_pruning/02_second_pass/*maf05*mm90*LD*.vcf.gz .
```

2. I converted the files from VCF to STRUCTURE format using PGDspider.
```
nano run_pgdspider.sh
```
```
#!/bin/bash
#SBATCH --job-name=pgdspider
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

eval "$(conda shell.bash hook)"
conda activate pgdspider

for file in *vcf.gz
do
        # Make temp uncompressed VCF file
        input_vcf_temp=$(mktemp)
        gunzip -c $file > $input_vcf_temp

        output_faststr=$(basename $file | sed 's/.vcf.gz/.faststr/')

        # Convert VCF to fastSTR format.
        python /nfs/home/zsmith10/.conda/envs/pgdspider/bin/PGDSpider2-cli \
                -inputfile $input_vcf_temp \
                -inputformat VCF \
                -outputfile $output_faststr \
                -outputformat STRUCTURE \
                -spid VCF2fastSTR.spid

        # Remove temp file
        rm $input_vcf_temp
done
```
```
sbatch run_pgdspider.sh
```

* Spid file:
```
# spid-file generated: Tue Aug 22 17:48:13 EDT 2023

# VCF Parser questions
PARSER_FORMAT=VCF

# Only output SNPs with a phred-scaled quality of at least:
VCF_PARSER_QUAL_QUESTION=0
# Select population definition file:
VCF_PARSER_POP_FILE_QUESTION=
# What is the ploidy of the data?
VCF_PARSER_PLOIDY_QUESTION=DIPLOID
# Do you want to include a file with population definitions?
VCF_PARSER_POP_QUESTION=false
# Output genotypes as missing if the phred-scale genotype quality is below:
VCF_PARSER_GTQUAL_QUESTION=
# Do you want to include non-polymorphic SNPs?
VCF_PARSER_MONOMORPHIC_QUESTION=true
# Only output following individuals (ind1, ind2, ind4, ...):
VCF_PARSER_IND_QUESTION=
# Only input following regions (refSeqName:start:end, multiple regions: whitespace separated):
VCF_PARSER_REGION_QUESTION=
# Output genotypes as missing if the read depth of a position for the sample is below:
VCF_PARSER_READ_QUESTION=
# Take most likely genotype if "PL" or "GL" is given in the genotype field?
VCF_PARSER_PL_QUESTION=false
# Do you want to exclude loci with only missing data?
VCF_PARSER_EXC_MISSING_LOCI_QUESTION=false

# STRUCTURE Writer questions
WRITER_FORMAT=STRUCTURE

# Specify the locus/locus combination you want to write to the STRUCTURE file:
STRUCTURE_WRITER_LOCUS_COMBINATION_QUESTION=
# Do you want to include inter-marker distances?
STRUCTURE_WRITER_LOCI_DISTANCE_QUESTION=false
# Specify which data type should be included in the STRUCTURE file  (STRUCTURE can only analyze one data type per file):
STRUCTURE_WRITER_DATA_TYPE_QUESTION=SNP
# Save more specific fastSTRUCTURE format?
STRUCTURE_WRITER_FAST_FORMAT_QUESTION=true
```

</p>
</details> 

***

<details><summary> 08. fastStructure </summary>
<p>

Directory: `/lustre/isaac/proj/UTK0032/zsmith10/chamaecyparis_thyoides/04_popgen_analysis/08_fastStructure`

1. First, I linked my fastStructure files into my fastStructure directory.
```
ln -s ../07_pgdspider/*.faststr .
```

2. Structure_threader requires a file with all individual names for plot labeling, which I made like so.
```
nano create_indfile.sh
```
```
eval "$(conda shell.bash hook)"
conda activate bcftools

bcftools query -l ../07_pgdspider/Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf05_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz > indfile.txt
```
```
bash create_indfile.sh
```

3. Next, I ran the simple model of fastStructure as implemented in structure_threader. I don't suspect weak population structure, so I did not run the logistic model (which assigns individuals more rigidly to populations and takes exponentially longer).
```
nano run_structure_threader_run.qsh
```
```
#!/bin/bash
#SBATCH --job-name=STRUCTURE
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus
#SBATCH --qos=campus
#SBATCH --time=6:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edux
#SBATCH --array=0-2
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

echo "host name : " `hostname`
echo This is array task number $SLURM_ARRAY_TASK_ID

# create an array variable containing the file names
FILES=($(ls -1 *.faststr))

# get specific file name, assign it to the array function
        # note that FILE variable is 0-indexed so
        # for convenience we also began the task IDs with 0
ARRAY_FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

mkdir fastStructure_simple_out

structure_threader \
        run \
        -fs /nfs/home/zsmith10/.local/bin/fastStructure \
        -K 10 \
        -R 10 \
        -i $ARRAY_FILE \
        -o ${ARRAY_FILE}_simple_out \
        -t 20 \
        --ind indfile.txt \
        --seed 1234 \
```
```
sbatch run_structure_threader_run.qsh
```

</p>
</details> 

***

<details><summary> 09. DAPC </summary>
<p>

1. DAPC is essentially used to evaluate between-population differences, while ignoring within-population differences by performing discriminant analysis on identified principal components. Thus, this method takes advantage of variation as synthetic variables segregate populations. I 

2. As usual, I linked over my VCF 3 files, specifically the VCFs filtered for MAF.
```
ln -s ../03_LD_pruning/04_fourth_pass/*_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz .
```

3. In addition, I used the following population file to assign populations for the DAPC. Here, I used the respective state of each sample to verify that there were indeed two metapopulations (North & South), including two subpopulations among the southern metapopulation.

```
head -n 10 popfile.csv
```
```
AL-L1_001,AL
AL-L1_002,AL
AL-L1_003,AL
AL-L1_004,AL
AL-L1_006,AL
AL-L1_007,AL
AL-L1_008,AL
AL-L1_009,AL
AL-L1_010,AL
AL-L1_011,AL
```

4. Next, I constructed an R script run a DAPC using adegenet.
```
nano dapc.R
```
```
# Define a CRAN mirror to download packages from.
cran_mirrors <- c('https://cloud.r-project.org/', 'http://cran.us.r-project.org')

# Deep dependency of poppr - libicui18n C library - hard-install the first time only
#install.packages("stringi", configure.args="--disable-pkg-config", repos = cran_mirrors)

# Specify libraries
packages_to_load <- c('adegenet', 'tidyverse', 'vcfR', 'poppr')

# Install packages, if required.
#for (package in packages_to_load) {
#    if (!require(package, quietly = TRUE))
#        install.packages(packages_to_load, dependencies = TRUE, repos = cran_mirrors)
#}

# Check if each package is installed, and if not, install it
for (package in packages_to_load) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = cran_mirrors)
  }
}

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

############################
#### Performing a DAPC ####
############################

# References: https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html

## Loop through each VCF and perform PCAs
# Set the path to your directory containing VCF files
directory_path <- "."

# Get a list of all files in the directory with a .vcf extension
vcf_files <- list.files(directory_path, pattern = "\\.vcf.gz$", full.names = TRUE)
print(vcf_files)

popfile <- read.csv("popfile.csv", header = FALSE)

# Loop through each VCF file
for (vcf_file in vcf_files) {

        # 1. Convert VCF to genind
        awc_vcf <- vcfR::read.vcfR(vcf_file)
        awc_genind <- vcfR::vcfR2genind(awc_vcf)
        head(awc_genind, 5)

        # 2. Assign pops
        pop(awc_genind) <- popfile$V2
        print(popfile$V2)

        # 3. Create DAPC object
        dapc.awc <- dapc(awc_genind, var.contrib = TRUE, scale = FALSE, n.pca = 2, n.da = nPop(awc_genind) - 1)
        print("dapc.awc")
        print(dapc.awc)

        # Extract sample coordinates & populations
        dapc_scores <- as.data.frame(dapc.awc$ind.coord, header = FALSE)
        print("dapc_scores")
        print(dapc_scores)

        # Combine into a single df
        dapc.awc_df <- cbind(dapc_scores, group = popfile$V2)
        glimpse(dapc.awc_df)

        # Plot - with 95% confidence intervals
        dapc.awc_plot <- ggplot(dapc.awc_df, aes(x = LD1, y = LD2, color = group)) +
                geom_point(size = 1) +
                stat_ellipse(type = "norm", linetype = 3) +
                labs(title = "DAPC among Chamaecyparis thyoides populations",
                        x = "LD1",
                        y = "LD2") +
                theme_minimal() +
                theme(plot.title = element_text(hjust = 0.5))

        # Save the file
        ggsave(filename = paste0(vcf_file, "_DAPC_plot.png"), plot = dapc.awc_plot, width = 8, height = 6)
}
```

5. To run this on ISAAC and still produce informative outfiles, I used a slurm script to run it.
```
nano run_dapc.qsh
```
```
#!/bin/bash
#SBATCH --job-name=DAPC_adegenet
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=20G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=short
#SBATCH --qos=short
#SBATCH --time=1:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edu

# Load R environment
eval "$(conda shell.bash hook)"
conda activate R

# Launch PCA Script
/nfs/home/zsmith10/.conda/envs/R/lib/R/bin/R --no-restore --file=dapc.R
```
```
sbatch run_dapc.qsh
```
6. I am attaching the results for each below. Of note, the results are virtually the same, supporting the previous fastStructure results. Specifically, there is a split in the southern metapopulation somewhere in Alabama, and a singular metapopulation in the northern Atlantic Coastal Plain (New Jersey/Maryland/Delaware).

De novo genome:
![Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0 4 vcf gz_DAPC_plot](https://github.com/user-attachments/assets/9f295fd1-0acb-430a-be8d-20123c6ac255)

T. plicata genome:
![Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0 4 vcf gz_DAPC_plot](https://github.com/user-attachments/assets/f8498656-f13e-4df3-a15b-a525096a17c1)

C. obtusa genome:
![Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0 4 vcf gz_DAPC_plot](https://github.com/user-attachments/assets/7830c7a1-adf9-4379-9e9e-66ba582960bb)

</p>
</details> 

***

<details><summary> 10. IQTree </summary>
<p>

</p>
</details> 

***

<details><summary> 11. dartR: Effective Population Size </summary>
<p>

1. First, I copied over my VCF files that were NOT MAF filtered (to capture as much variation throughout the population as possible.) I then unzipped these files for dartR.
```
cp ../02_bcftools/04_fourth_pass_HDplot-MaxH_0.50/*_biallelic_minDP3_PASSonly_removed-indv_removed-paralogs_mm90.vcf.gz .
```

2. Next, I reformatted my population file to pertain to putative subspecies (subsp. henryae in the South and subsp. thyoides in the North). The file looks like so:
```
head -n 10 popfile_subspecies.tsv
```
```
AL-L1_001       henryae
AL-L1_002       henryae
AL-L1_003       henryae
AL-L1_004       henryae
AL-L1_006       henryae
AL-L1_007       henryae
AL-L1_008       henryae
AL-L1_009       henryae
AL-L1_010       henryae
AL-L1_011       henryae
```

3. Since dartR is an R package, I constructed an R script to calculate the effective population size of each putative subspecies.

```
nano dartR.R
```
```
#!/usr/bin/env Rscript

##########################################
#### Installing and loading libraries ####
##########################################

# Define a CRAN mirror to download packages from.
cran_mirrors <- c("https://cran.r-project.org/","https://cloud.r-project.org/")

# Specify other libraries - rgdal is essential for dartR and will give install issues without directly downloading it.
packages_to_load <- c("dartRverse", "vcfR", "dartR.popgen", "poppr")

# Check if each package is installed, and if not, install it
for (package in packages_to_load) {
  if (!requireNamespace(package, quietly = TRUE)){
        install.packages(package, repos = cran_mirrors, dependencies = TRUE)
  }
}

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

# Load all dartR modules
dartRverse_install("all")

#############################
####### Data Loading #######
#############################

# Data Load-in & Curation - Read in VCFs with vcfR - gl.read.vcf() from dartR didn't work for me.
denovo_vcf <- "Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_minDP3_PASSonly_removed-indv_removed-paralogs_mm90.vcf"
Cobtusa_vcf <- "Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_minDP3_PASSonly_removed-indv_removed-paralogs_mm90.vcf"
Tplicata_vcf <- "Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_minDP3_PASSonly_removed-indv_removed-paralogs_mm90.vcf"

# Subsetting to northern and southern populations
#northern_states <- c("DE", "NJ", "MD")
#southern_states <- c("AL", "FL", "MS")

# Assign populations & read in genlight
pop_csv <- read.csv("popfile_subspecies.tsv", sep = '\t', header = FALSE)
pops <- pop_csv$V2

######################################
############ De novo ################
######################################

# Assign populations & read in genlight
denovo_vcf_R <- read.vcfR(denovo_vcf)
denovo_gl <- vcfR2genlight(denovo_vcf_R)
pop(denovo_gl) <- pops

# Save and restore genlight objects if rerunning analysis
# Save an object to a file
saveRDS(denovo_gl, file = "denovo_gl.rds")
# Restore the object
#denovo_gl <- readRDS(file = "denovo_gl.rds")

# Subsets
#denovo_gl_north <- popsub(denovo_gl, sublist = northern_states)
#denovo_gl_south <- popsub(denovo_gl, sublist = southern_states)

# Run dartR
dartR.popgen::gl.LDNe(denovo_gl,
        outfile = "denovo_LDNe.txt",
        outpath = "denovo_out",
        neest.path = "/nfs/home/zsmith10/software/NeEstimator")

###################################
########### C. obtusa #############
###################################

# Assign populations & read in genlight
Cobtusa_vcf_R <- read.vcfR(Cobtusa_vcf)
Cobtusa_gl <- vcfR2genlight(Cobtusa_vcf_R)
pop(Cobtusa_gl) <- pops

# Save and restore genlight objects if rerunning analysis
# Save an object to a file
saveRDS(Cobtusa_gl, file = "Cobtusa_gl.rds")
# Restore the object
#Cobtusa_gl <- readRDS(file = "Cobtusa_gl.rds")

# Subsets
#Cobtusa_gl_north <- popsub(Cobtusa_gl, sublist = northern_states)
#Cobtusa_gl_south <- popsub(Cobtusa_gl, sublist = southern_states)

# Run dartR
dartR.popgen::gl.LDNe(Cobtusa_gl,
        outfile = "Cobtusa_LDNe.txt",
        outpath = "Cobtusa_out",
        neest.path = "/nfs/home/zsmith10/software/NeEstimator")

######################################
########### T. plicata ###############
####################a#################

# Assign populations & read in genlight
Tplicata_vcf_R <- read.vcfR(Tplicata_vcf)
Tplicata_gl <- vcfR2genlight(Tplicata_vcf_R)
pop(Tplicata_gl) <- pops

# Save and restore genlight objects if rerunning analysis
# Save an object to a file
saveRDS(Tplicata_gl, file = "Tplicata_gl.rds")
# Restore the object
#Tplicata_gl <- readRDS(file = "Tplicata_gl.rds")

# Subsets
#Tplicata_gl_north <- popsub(Tplicata_gl, sublist = northern_states)
#Tplicata_gl_south <- popsub(Tplicata_gl, sublist = southern_states)

# Run dartR
dartR.popgen::gl.LDNe(Tplicata_gl,
        outfile = "Tplicata_LDNe.txt",
        outpath = "Tplicata_out",
        neest.path = "/nfs/home/zsmith10/software/NeEstimator")

```

4. To run this on ISAAC and still produce informative outfiles, I used a slurm script to run it. Note that some of my VCF files are quite large, so I ran this on the bigmem node with 200G of RAM--just to be safe.

```
run_dartR.sh
```
```
#!/bin/bash
#SBATCH --job-name=dartR
#SBATCH --nodes=1
#SBATCH --ntasks=25
#SBATCH --mem=200G
#SBATCH -A ACF-UTK0032
#SBATCH --partition=campus-bigmem
#SBATCH --qos=campus-bigmem
#SBATCH --time=24:00:00
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=zsmith10@vols.utk.edux

# Load R environment
eval "$(conda shell.bash hook)"
conda activate R

mkdir denovo_out Cobtusa_out Tplicata_out

# Launch PCA Script
/nfs/home/zsmith10/.conda/envs/R/lib/R/bin/R --no-restore --file=dartR.R

```
```
sbatch run_dartR.sh
```

5. These are the preliminary results from the de novo SNP calls.
```
$`henryae_AL-L1`
                    Statistic Frequency 1 Frequency 2
 Lowest Allele Frequency Used          0+       No S*
    Harmonic Mean Sample Size       102.5       102.9
      Independent Comparisons    41982490    15109397
                  OverAll r^2    0.010584    0.011008
          Expected r^2 Sample    0.010009    0.009978
                Estimated Ne^         578       321.3
            CI low Parametric       573.5       318.9
           CI high Parametric       582.6       323.8
             CI low JackKnife       468.3       266.2
            CI high JackKnife       750.3       402.1

$thyoides_DE1
                    Statistic Frequency 1 Frequency 2
 Lowest Allele Frequency Used          0+       No S*
    Harmonic Mean Sample Size         129       129.3
      Independent Comparisons    27086066     9484569
                  OverAll r^2     0.00855    0.008393
          Expected r^2 Sample    0.007934    0.007921
                Estimated Ne^       539.8         704
            CI low Parametric       535.8       692.9
           CI high Parametric       543.8       715.5
             CI low JackKnife       361.5         596
            CI high JackKnife      1011.9       856.8
```

</p>
</details> 

***

<details><summary> 12. Isolation-by-Distance: Mantel Test </summary>
<p>

1. Mantel test:
```
  GNU nano 2.9.8                                                                                                                                                                                                                                                                        dartR_IBD_loop.R                                                                                                                                                                                                                                                                                   

#!/usr/bin/env Rscript

## Installing and loading libraries

# Install the remotes package if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Define a CRAN mirror to download packages from.
cran_mirrors <- c("https://cran.r-project.org/","https://cloud.r-project.org/")

# Specify other libraries
packages_to_load <- c("vcfR", "poppr", "dartRverse", "raster")

# Check if each package is installed, and if not, install it
#for (package in packages_to_load) {
#  if (!requireNamespace(package, quietly = TRUE)){
#    install.packages(package, repos = cran_mirrors, dependencies = TRUE)
#  }
#}

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

## Data Loading

# Load population and coordinate data
pop_csv <- read.csv("popfile.tsv", sep = '\t', header = FALSE)
pops <- pop_csv$V2

coord_csv <- read.csv("coord_file.tsv", sep = '\t', header = TRUE)
coords <- as.data.frame(coord_csv[, c(3, 2)])

## Perform Mantel Test on VCF files

# List all VCF files in the directory
vcf_files <- list.files(pattern = "*.vcf.gz")

for (vcf_file in vcf_files) {
  # Read in the VCF file
  vcf_R <- read.vcfR(vcf_file)
  all_gl <- vcfR2genlight(vcf_R)

  # Assign populations and lat/long
  pop(all_gl) <- pops
  all_gl@other$latlon <- coords

  # Create genlight subsets
  ACP_gl <- all_gl[pop(all_gl) %in% c("ACP"), ]
  WGCP_gl <- all_gl[pop(all_gl) %in% c("WGCP"), ]
  EGCP_gl <- all_gl[pop(all_gl) %in% c("EGCP"), ]
  ACP_WGCP_gl <- all_gl[pop(all_gl) %in% c("ACP", "WGCP"), ]
  ACP_EGCP_gl <- all_gl[pop(all_gl) %in% c("ACP", "EGCP"), ]
  WGCP_EGCP_gl <- all_gl[pop(all_gl) %in% c("WGCP", "EGCP"), ]

  # Save the genlight object
  output_rds <- sub(".vcf.gz", "_gl.rds", vcf_file)
  saveRDS(all_gl, file = output_rds)

  # Perform population-based Mantel tests
  print("RUNNING all_ibd:")
  all_ibd <- dartR::gl.ibd(all_gl, distance = "euclidean", paircols = 'pop', permutations = 10000)
  all_ibd

  print("RUNNING ACP_ibd:")
  ACP_ibd <- dartR::gl.ibd(ACP_gl, distance = "euclidean", paircols = 'pop', permutations = 10000)
  ACP_ibd

  print("RUNNING WGCP_ibd:")
  WGCP_ibd <- dartR::gl.ibd(WGCP_gl, distance = "euclidean", paircols = 'pop', permutations = 10000)
  WGCP_ibd

  print("RUNNING EGCP_ibd:")
  EGCP_ibd <- dartR::gl.ibd(EGCP_gl, distance = "euclidean", paircols = 'pop', permutations = 10000)
  EGCP_ibd

  print("RUNNING ACP_WGCP_ibd:")
  ACP_WGCP_ibd <- dartR::gl.ibd(ACP_WGCP_gl, distance = "euclidean", paircols = 'pop', permutations = 10000)
  ACP_WGCP_ibd

  print("RUNNING ACP_EGCP_ibd:")
  ACP_EGCP_ibd <- dartR::gl.ibd(ACP_EGCP_gl, distance = "euclidean", paircols = 'pop', permutations = 10000)
  ACP_EGCP_ibd

  print("RUNNING WGCP_EGCP_ibd:")
  WGCP_EGCP_ibd <- dartR::gl.ibd(WGCP_EGCP_gl, distance = "euclidean", paircols = 'pop', permutations = 10000)
  WGCP_EGCP_ibd
```


</p>
</details> 

***

<details><summary> 13. RDA </summary>
<p>

Directory: `/pickett_sphinx/projects/zsmith10/chamaecyparis_thyoides/12_RDA`

1. RDA is essentially a dimensionality reduction analysis that examine the correlation of the genotype matrix to scaled environmental matrices to evaluate candidate environmental drivers of genetic diversity, aka adaptation. 

Unfortunately, I had substantial issues installing the updated geospatial packages on ISAAC, so I did this analysis on Sphinx. These are the files I coped over from ISAAC & my metadata.
* VCFs
```
Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz
```
```
Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz
```
```
Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz
```

* Coordinate data
```
head -n 5 coord_file.tsv
```
```
sample  lat     long
AL-L1_001       30.64124446     -88.39755033
AL-L1_002       30.6413274      -88.39754653
AL-L1_003       30.64131946     -88.3975295
AL-L1_004       30.64134089     -88.3976078

```

* Population file
```
head -n 5 popfile_subspecies.tsv
```
```
AL-L1_001       henryae
AL-L1_002       henryae
AL-L1_003       henryae
AL-L1_004       henryae
AL-L1_006       henryae
```

2. I then downloaded environmental data from the 30-year WorldClim normals (the 19 BioClim variables; 1970-2000).
* WorldClim link: https://www.worldclim.org/data/worldclim21.html

Additionally, I used the WorldClim 30-year normals for monthly vapor pressure deficit, solar radiation, and wind speed in addition to the digital elevational model.

I am currently using the 5 arcminute layers as there is some uncertainty around some of the coordinate information. However, I downloaded all of the data, in case I wanted to use it.

I downloaded this data in the following directory: `/pickett_sphinx/projects/zsmith10/chamaecyparis_thyoides/WorldClim_8.4.24`

```
nano download_bioclim.sh
```
```
#30s
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_wind.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_vapr.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_srad.zip >> download.log &

#2.5m
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_vapr.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_wind.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_srad.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_bio.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_2.5m_elev.zip >> download.log &

#5m
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_5m_elev.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_5m_bio.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_5m_vapr.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_5m_wind.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_5m_srad.zip >> download.log &

#10m
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_elev.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_vapr.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_wind.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_srad.zip >> download.log &
wget https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip >> download.log &
```
```
bash download_bioclim/sh
```

3. Next, I used the coordinate file to extract environmental variables at each coordinate point. 
```
nano run_step1_get-WorldClim.qsh
```
```
eval "$(conda shell.bash hook)"
conda activate R

# Launch PCA Script
Rscript WorldClim.R
```
```
bash run_step1_get-WorldClim.qsh
```

* WorldClim.R
```
# Define a CRAN mirror to download packages from.
cran_mirrors <- c('http://cran.us.r-project.org', 'https://cloud.r-project.org/')

# Specify libraries
packages_to_load <- c('raster', 'sp', 'geodata', 'terra')

# Check if each package is installed, and if not, install it
for (package in packages_to_load) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = cran_mirrors)
  }
}

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

# Load coordinates from CSV
coords <- read.csv("coord_file.tsv", sep = '\t', header = TRUE)

# Convert specific columns to numeric
coords$lat <- as.numeric(coords$lat)
coords$long <- as.numeric(coords$long)

# Convert coordinates to SpatialPoints
sample_points <- SpatialPoints(coords[, c("long", "lat")], proj4string = CRS("+proj=longlat +datum=WGS84"))
head(sample_points)

# Convert from SpatPoints to SpatVector for value extraction
sample_points_vector <- vect(sample_points)
head(sample_points_vector)

# Download the dataset using geodata function Replace 'path_to_files' with the actual path to your GeoTIFF files
tiff_files <- list.files(path = "/pickett_sphinx/projects/zsmith10/chamaecyparis_thyoides/WorldClim_8.4.24/5m", pattern = "\\.tif$", full.names = TRUE)

# Read each file as a SpatRaster object
raster_list <- lapply(tiff_files, function(file) rast(file))

# Combine into a raster stack
bioclim <- c(raster_list)
bioclim

# Reproject the WorldClim raster layers - DEPRECATED
bioclim_reproj <- project(bioclim, y = "+proj=longlat +datum=WGS84") Extract values from bioclim at sample_points_vector
bioclim_data <- lapply(bioclim, function(r)
        { extract(r, vect(sample_points))
})

head(bioclim_data, 10)

# Combine coordinates and bioclimatic data
result <- cbind(coords, bioclim_data)
#View(result)

# Write the result to a new CSV
write.csv(result, "AWC_coords_bioclim-vars.csv", row.names = FALSE)

```

4. Next, I wrote a script to (manually) awk the appropriate data columns into a new spreadsheet. Oddly, the previous script includes ID columns between each data column.
```
nano run_step2_extract_bioclim-vars.sh
```
```
awk -F','  '{OFS=","; print $1, $5, $7, $9, $11, $13, $15, $17, $19, $21, $23, $25, $27, $29, $31, $33, $35, $37, $39, $41, $43, $45, $47, $49, $51, $53, $55, $57, $59, $61, $63, $65, $67, $69, $71, $73, $75, $77, $79, $81, $83, $85, $87, $89, $91, $93, $95, $97, $99, $101, $103, $105, $107, $109, $111, $113, $115}' AWC_coords_bioclim-vars.csv | sed 's/"//g' > AWC_bioclim-vars.csv
```
```
bash run_step2_extract_bioclim-vars.sh
```

5. Time for the RDA! As always, this was ran in triplicate for all three reference VCFs. I did include the `nice` command and `cpulimit -l 4000` function to attempt to limit my job to 40 cores. It didn't really work as it still took all of sphinx's cores.
```
nano run_step3_vegan-RDA.qsh
```
```
 Load R environment
eval "$(conda shell.bash hook)"
conda activate R

# Launch PCA Script
nice cpulimit -l 4000 Rscript vegan-RDA_5.R > vegan-RDA.log
```
```
bash run_step3_vegan-RDA.qsh
```

* `vegan-RDA_5.R`
```
# Define a CRAN mirror to download packages from.
cran_mirrors <- c('http://cran.us.r-project.org', 'https://cloud.r-project.org/')

# Specify libraries
packages_to_load <- c('vegan', 'mice', 'vcfR', 'ggplot2', 'GGally', 'dplyr', 'tidyr', 'ggrepel')

# Check if each package is installed, and if not, install it
for (package in packages_to_load) {
  if (!requireNamespace(package, quietly = TRUE)){
        install.packages(package, repos = cran_mirrors, dependencies = TRUE)
  }
}

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

## Loop through each VCF and perform PCAs
# Set the path to your directory containing VCF files
directory_path <- "."

# Assign VCF names to objects.
Cobtusa_vcf <- read.vcfR("Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz")
Tplicata_vcf <- read.vcfR("Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz")
denovo_vcf <- read.vcfR("Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf.gz")

popfile <- read.csv("popfile_subspecies.tsv", sep = '\t', header = FALSE)
pops <- popfile$V2

# Extract genotype data
Cobtusa_genotypes <- extract.gt(Cobtusa_vcf, element = "GT")
Tplicata_genotypes <- extract.gt(Tplicata_vcf, element = "GT")
denovo_genotypes <- extract.gt(denovo_vcf, element = "GT")

Cobtusa_genotypes_df <- as.data.frame(Cobtusa_genotypes)
Tplicata_genotypes_df <- as.data.frame(Tplicata_genotypes)
denovo_genotypes_df <- as.data.frame(denovo_genotypes)

#head(Cobtusa_genotypes_df)
#head(Tplicata_genotypes_df)
#head(denovo_genotypes_df)

# Format genotype data into 0,1,2,NA format Example conversion function
# Create conversion function:
convert_genotypes <- function(genotype) {
  genotype <- as.character(genotype)
  genotype[genotype == "0/0"] <- 0
  genotype[genotype == "0/1"] <- 1
  genotype[genotype == "1/1"] <- 2
  genotype[genotype == "./."] <- NA
  return(as.numeric(genotype))
}

# Now convert each genotype for all genotype dataframes
Cobtusa_converted_genotypes_df <- as.data.frame(lapply(Cobtusa_genotypes_df, convert_genotypes))
Tplicata_converted_genotypes_df <- as.data.frame(lapply(Tplicata_genotypes_df, convert_genotypes))
denovo_converted_genotypes_df <- as.data.frame(lapply(denovo_genotypes_df, convert_genotypes))

Cobtusa_converted_cleaned_genotypes_df <- Cobtusa_converted_genotypes_df %>%
  drop_na(everything())
#nrow(Cobtusa_converted_genotypes_df)
#nrow(Cobtusa_converted_cleaned_genotypes_df)

Tplicata_converted_cleaned_genotypes_df <- Tplicata_converted_genotypes_df %>%
  drop_na(everything())
#nrow(Tplicata_converted_genotypes_df)
#nrow(Tplicata_converted_cleaned_genotypes_df)

denovo_converted_cleaned_genotypes_df <- denovo_converted_genotypes_df %>%
  drop_na(everything())
#nrow(denovo_converted_genotypes_df)
#nrow(denovo_converted_cleaned_genotypes_df)

##### BioClim Variables ######
AWC_bioclim_vars <- read.csv("AWC_bioclim-vars.csv", header = TRUE)
ncol(AWC_bioclim_vars)

#Here, I scaled the variables to make sure values are the most comparably for comparisons.
AWC_bioclim_vars[2:57] <- scale(AWC_bioclim_vars[2:57])
#AWC_bioclim_vars <- as.numeric(AWC_bioclim_vars[2:20])
AWC_bioclim_vars <- as.data.frame(AWC_bioclim_vars)

# Reassign names
names(AWC_bioclim_vars) <- c("Sample",
                        "Annual_Mean_Temperature",
                        "Mean_Diurnal_Range",
                        "Isothermality",
                        "Temperature_Seasonality",
                        "Max_Temperature_of_Warmest_Month",
                        "Min_Temperature_of_Coldest_Month",
                        "Temperature_Annual_Range",
                        "Mean_Temperature_of_Wettest_Quarter",
                        "Mean_Temperature_of_Driest_Quarter",
                        "Mean_Temperature_of_Warmest_Quarter",
                        "Mean_Temperature_of_Coldest_Quarter",
                        "Annual_Precipitation",
                        "Precipitation_of_Wettest_Month",
                        "Precipitation_of_Driest_Month",
                        "Precipitation_Seasonality",
                        "Precipitation_of_Wettest_Quarter",
                        "Precipitation_of_Driest_Quarter",
                        "Precipitation_of_Warmest_Quarter",
                        "Precipitation_of_Coldest_Quarter",
                        "Elevation",
                        "January_Solar_Radiation",
                        "February_Solar_Radiation",
                        "March_Solar_Radiation",
                        "April_Solar_Radiation",
                        "May_Solar_Radiation",
                        "June_Solar_Radiation",
                        "July_Solar_Radiation",
                        "August_Solar_Radiation",
                        "September_Solar_Radiation",
                        "October_Solar_Radiation",
                        "November_Solar_Radiation",
                        "December_Solar_Radiation",
                        "January_Vapor_Pressure_Deficit",
                        "February_Vapor_Pressure_Deficit",
                        "March_Vapor_Pressure_Deficit",
                        "April_Vapor_Pressure_Deficit",
                        "May_Vapor_Pressure_Deficit",
                        "June_Vapor_Pressure_Deficit",
                        "July_Vapor_Pressure_Deficit",
                        "August_Vapor_Pressure_Deficit",
                        "September_Vapor_Pressure_Deficit",
                        "October_Vapor_Pressure_Deficit",
                        "November_Vapor_Pressure_Deficit",
                        "December_Vapor_Pressure_Deficit",
                        "January_Wind_Speed",
                        "February_Wind_Speed",
                        "March_Wind_Speed",
                        "April_Wind_Speed",
                        "May_Wind_Speed",
                        "June_Wind_Speed",
                        "July_Wind_Speed",
                        "August_Wind_Speed",
                        "September_Wind_Speed",
                        "October_Wind_Speed",
                        "November_Wind_Speed",
                        "December_Wind_Speed"
                        )

length(names(AWC_bioclim_vars))
#summary(AWC_bioclim_vars)

##### Variable Collinearity #####
#Evaluate collinearity among predictors
# Create the pairs plot using ggpairs
bioclim_collinearity_plot <- ggpairs(AWC_bioclim_vars, columns = c(2:57), color = "cyl", size = 10)
ggsave("bioclim_collinearity_plot.png", plot = bioclim_collinearity_plot, width = 40, height = 30, dpi = 300)

#bioclim_collinearity_subset_plot <- ggpairs(AWC_bioclim_vars, columns = c(6,7,13,15,19))
#ggsave("bioclim_collinearity_subset_plot.png", plot = bioclim_collinearity_subset_plot, width = 20, height = 15, dpi = 300)

# Subset down the bioclim variables
#AWC_bioclim_vars <- AWC_bioclim_vars[, c(1,6,7,13,15,19)]
nrow(AWC_bioclim_vars)
ncol(AWC_bioclim_vars)

##### RDA  #####
# Perform RDA
transposed_Cobtusa_converted_cleaned_genotypes_df <- t(Cobtusa_converted_cleaned_genotypes_df)
transposed_Tplicata_converted_cleaned_genotypes_df <- t(Tplicata_converted_cleaned_genotypes_df)
transposed_denovo_converted_cleaned_genotypes_df <- t(denovo_converted_cleaned_genotypes_df)
#nrow(transposed_denovo_converted_cleaned_genotypes_df)

Cobtusa_rda <- rda(transposed_Cobtusa_converted_cleaned_genotypes_df ~
                                AWC_bioclim_vars$Annual_Mean_Temperature +
                                AWC_bioclim_vars$Mean_Diurnal_Range +
                                AWC_bioclim_vars$Isothermality +
                                AWC_bioclim_vars$Temperature_Seasonality +
                                AWC_bioclim_vars$Max_Temperature_of_Warmest_Month +
                                AWC_bioclim_vars$Min_Temperature_of_Coldest_Month +
                                AWC_bioclim_vars$Temperature_Annual_Range +
                                AWC_bioclim_vars$Mean_Temperature_of_Wettest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Driest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Warmest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Coldest_Quarter +
                                AWC_bioclim_vars$Annual_Precipitation +
                                AWC_bioclim_vars$Precipitation_of_Wettest_Month +
                                AWC_bioclim_vars$Precipitation_of_Driest_Month +
                                AWC_bioclim_vars$Precipitation_Seasonality +
                                AWC_bioclim_vars$Precipitation_of_Wettest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Driest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Warmest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Coldest_Quarter +
                                AWC_bioclim_vars$Elevation +
                                AWC_bioclim_vars$January_Solar_Radiation +
                                AWC_bioclim_vars$February_Solar_Radiation +
                                AWC_bioclim_vars$March_Solar_Radiation +
                                AWC_bioclim_vars$April_Solar_Radiation +
                                AWC_bioclim_vars$May_Solar_Radiation +
                                AWC_bioclim_vars$June_Solar_Radiation +
                                AWC_bioclim_vars$July_Solar_Radiation +
                                AWC_bioclim_vars$August_Solar_Radiation +
                                AWC_bioclim_vars$September_Solar_Radiation +
                                AWC_bioclim_vars$October_Solar_Radiation +
                                AWC_bioclim_vars$November_Solar_Radiation +
                                AWC_bioclim_vars$December_Solar_Radiation +
                                AWC_bioclim_vars$January_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$February_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$March_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$April_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$May_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$June_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$July_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$August_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$September_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$October_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$November_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$December_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$January_Wind_Speed +
                                AWC_bioclim_vars$February_Wind_Speed +
                                AWC_bioclim_vars$March_Wind_Speed +
                                AWC_bioclim_vars$April_Wind_Speed +
                                AWC_bioclim_vars$May_Wind_Speed +
                                AWC_bioclim_vars$June_Wind_Speed +
                                AWC_bioclim_vars$July_Wind_Speed +
                                AWC_bioclim_vars$August_Wind_Speed +
                                AWC_bioclim_vars$September_Wind_Speed +
                                AWC_bioclim_vars$October_Wind_Speed +
                                AWC_bioclim_vars$November_Wind_Speed +
                                AWC_bioclim_vars$December_Wind_Speed,
                                data = AWC_bioclim_vars)

Tplicata_rda <- rda(transposed_Tplicata_converted_cleaned_genotypes_df ~
                                AWC_bioclim_vars$Annual_Mean_Temperature +
                                AWC_bioclim_vars$Mean_Diurnal_Range +
                                AWC_bioclim_vars$Isothermality +
                                AWC_bioclim_vars$Temperature_Seasonality +
                                AWC_bioclim_vars$Max_Temperature_of_Warmest_Month +
                                AWC_bioclim_vars$Min_Temperature_of_Coldest_Month +
                                AWC_bioclim_vars$Temperature_Annual_Range +
                                AWC_bioclim_vars$Mean_Temperature_of_Wettest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Driest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Warmest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Coldest_Quarter +
                                AWC_bioclim_vars$Annual_Precipitation +
                                AWC_bioclim_vars$Precipitation_of_Wettest_Month +
                                AWC_bioclim_vars$Precipitation_of_Driest_Month +
                                AWC_bioclim_vars$Precipitation_Seasonality +
                                AWC_bioclim_vars$Precipitation_of_Wettest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Driest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Warmest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Coldest_Quarter +
                                AWC_bioclim_vars$Elevation +
                                AWC_bioclim_vars$January_Solar_Radiation +
                                AWC_bioclim_vars$February_Solar_Radiation +
                                AWC_bioclim_vars$March_Solar_Radiation +
                                AWC_bioclim_vars$April_Solar_Radiation +
                                AWC_bioclim_vars$May_Solar_Radiation +
                                AWC_bioclim_vars$June_Solar_Radiation +
                                AWC_bioclim_vars$July_Solar_Radiation +
                                AWC_bioclim_vars$August_Solar_Radiation +
                                AWC_bioclim_vars$September_Solar_Radiation +
                                AWC_bioclim_vars$October_Solar_Radiation +
                                AWC_bioclim_vars$November_Solar_Radiation +
                                AWC_bioclim_vars$December_Solar_Radiation +
                                AWC_bioclim_vars$January_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$February_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$March_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$April_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$May_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$June_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$July_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$August_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$September_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$October_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$November_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$December_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$January_Wind_Speed +
                                AWC_bioclim_vars$February_Wind_Speed +
                                AWC_bioclim_vars$March_Wind_Speed +
                                AWC_bioclim_vars$April_Wind_Speed +
                                AWC_bioclim_vars$May_Wind_Speed +
                                AWC_bioclim_vars$June_Wind_Speed +
                                AWC_bioclim_vars$July_Wind_Speed +
                                AWC_bioclim_vars$August_Wind_Speed +
                                AWC_bioclim_vars$September_Wind_Speed +
                                AWC_bioclim_vars$October_Wind_Speed +
                                AWC_bioclim_vars$November_Wind_Speed +
                                AWC_bioclim_vars$December_Wind_Speed,
                                data = AWC_bioclim_vars)

denovo_rda <- rda(transposed_denovo_converted_cleaned_genotypes_df ~
                                AWC_bioclim_vars$Annual_Mean_Temperature +
                                AWC_bioclim_vars$Mean_Diurnal_Range +
                                AWC_bioclim_vars$Isothermality +
                                AWC_bioclim_vars$Temperature_Seasonality +
                                AWC_bioclim_vars$Max_Temperature_of_Warmest_Month +
                                AWC_bioclim_vars$Min_Temperature_of_Coldest_Month +
                                AWC_bioclim_vars$Temperature_Annual_Range +
                                AWC_bioclim_vars$Mean_Temperature_of_Wettest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Driest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Warmest_Quarter +
                                AWC_bioclim_vars$Mean_Temperature_of_Coldest_Quarter +
                                AWC_bioclim_vars$Annual_Precipitation +
                                AWC_bioclim_vars$Precipitation_of_Wettest_Month +
                                AWC_bioclim_vars$Precipitation_of_Driest_Month +
                                AWC_bioclim_vars$Precipitation_Seasonality +
                                AWC_bioclim_vars$Precipitation_of_Wettest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Driest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Warmest_Quarter +
                                AWC_bioclim_vars$Precipitation_of_Coldest_Quarter +
                                AWC_bioclim_vars$Elevation +
                                AWC_bioclim_vars$January_Solar_Radiation +
                                AWC_bioclim_vars$February_Solar_Radiation +
                                AWC_bioclim_vars$March_Solar_Radiation +
                                AWC_bioclim_vars$April_Solar_Radiation +
                                AWC_bioclim_vars$May_Solar_Radiation +
                                AWC_bioclim_vars$June_Solar_Radiation +
                                AWC_bioclim_vars$July_Solar_Radiation +
                                AWC_bioclim_vars$August_Solar_Radiation +
                                AWC_bioclim_vars$September_Solar_Radiation +
                                AWC_bioclim_vars$October_Solar_Radiation +
                                AWC_bioclim_vars$November_Solar_Radiation +
                                AWC_bioclim_vars$December_Solar_Radiation +
                                AWC_bioclim_vars$January_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$February_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$March_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$April_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$May_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$June_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$July_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$August_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$September_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$October_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$November_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$December_Vapor_Pressure_Deficit +
                                AWC_bioclim_vars$January_Wind_Speed +
                                AWC_bioclim_vars$February_Wind_Speed +
                                AWC_bioclim_vars$March_Wind_Speed +
                                AWC_bioclim_vars$April_Wind_Speed +
                                AWC_bioclim_vars$May_Wind_Speed +
                                AWC_bioclim_vars$June_Wind_Speed +
                                AWC_bioclim_vars$July_Wind_Speed +
                                AWC_bioclim_vars$August_Wind_Speed +
                                AWC_bioclim_vars$September_Wind_Speed +
                                AWC_bioclim_vars$October_Wind_Speed +
                                AWC_bioclim_vars$November_Wind_Speed +
                                AWC_bioclim_vars$December_Wind_Speed,
                                data = AWC_bioclim_vars)

##### Re-evaluating model expectations #####

forward.selection_denovo <- ordiR2step(rda(transposed_denovo_converted_cleaned_genotypes_df ~ 1),
                                scope = formula(denovo_rda),
                                direction = "forward",
                                R2scope = TRUE,
                                pstep = 1000,
                                trace = FALSE)
print("Printing result for denovo forward selection:")
forward.selection_denovo$call

forward.selection_Tplicata <- ordiR2step(rda(transposed_Tplicata_converted_cleaned_genotypes_df ~ 1),
                                scope = formula(Tplicata_rda),
                                direction = "forward",
                                R2scope = TRUE,
                                pstep = 1000,
                                trace = FALSE)
print("Printing result for Tplicata forward selection:")
forward.selection_Tplicata$call

forward.selection_Cobtusa <- ordiR2step(rda(transposed_Cobtusa_converted_cleaned_genotypes_df ~ 1),
                                scope = formula(Cobtusa_rda),
                                direction = "forward",
                                R2scope = TRUE,
                                pstep = 1000,
                                trace = FALSE)
print("Printing result for Cobtusa forward selection:")
forward.selection_Cobtusa$call

# Check for 0 variance variables
#print("Checking for zero variance variables:")
#apply(AWC_bioclim_vars[, -1], 2, var)

# Print RDA result
## Cobtusa
# Extract site scores (samples)
Cobtusa_site_scores <- as.data.frame(scores(Cobtusa_rda, display = "sites"))
Cobtusa_site_scores$Sample <- as.data.frame(scores(Cobtusa_rda, display = "sites"))
Cobtusa_site_scores$pops <- factor(pops)
#head(Cobtusa_site_scores, 5)

# Extract biplot scores (explanatory variables) - bp arg = biplot
Cobtusa_biplot_scores <- as.data.frame(scores(Cobtusa_rda, display = "bp"))

# Determine if any variables were dropped from the model
rownames(Cobtusa_biplot_scores)
Cobtusa_biplot_scores$variable <- rownames(Cobtusa_biplot_scores)
# Yep, Precipitation_of_Driest_Quarter was dropped, so let's remove it

# Plot RDA
Cobtusa_rda_plot <- ggplot() +
  # Plot site scores (samples)
  geom_point(data = Cobtusa_site_scores, aes(x = RDA1, y = RDA2, color = pops), size = 2) +
  # Add arrows for biplot scores (explanatory variables)
  geom_segment(data = Cobtusa_biplot_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                                                arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "darkblue") +
  # Add labels for biplot scores
  geom_text_repel(data = Cobtusa_biplot_scores, aes(x = RDA1, y = RDA2, label = rownames(Cobtusa_biplot_scores)), size = 1, color = "blue") +
  # Customize the plot
  labs(title = "Chamaecyparis obtusa WorldClim RDA Biplot", x = "RDA1", y = "RDA2", color = "Subspecies") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("Cobtusa_RDA.png", plot = Cobtusa_rda_plot, width = 8, height = 6, dpi = 300)

## Tplicata
# Extract site scores (samples)
Tplicata_site_scores <- as.data.frame(scores(Tplicata_rda, display = "sites"))
Tplicata_site_scores$Sample <- as.data.frame(scores(Tplicata_rda, display = "sites"))
Tplicata_site_scores$pops <- factor(pops)

# Extract biplot scores (explanatory variables) - bp arg = biplot
Tplicata_biplot_scores <- as.data.frame(scores(Tplicata_rda, display = "bp"))

# Determine if any variables were dropped from the model
rownames(Tplicata_biplot_scores)
Tplicata_biplot_scores$variable <- rownames(Tplicata_biplot_scores)
# Yep, Precipitation_of_Driest_Quarter was dropped, so let's remove it -- Deprecated & no longer necessary rownames() removes it by default.

# Plot RDA
Tplicata_rda_plot <- ggplot() +
  # Plot site scores (samples)
  geom_point(data = Tplicata_site_scores, aes(x = RDA1, y = RDA2, color = pops), size = 2) +
  # Add arrows for biplot scores (explanatory variables)
  geom_segment(data = Tplicata_biplot_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                                                arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "darkblue") +
  # Add labels for biplot scores
  geom_text_repel(data = Tplicata_biplot_scores, aes(x = RDA1, y = RDA2, label = rownames(Tplicata_biplot_scores)), size = 1, color = "blue") +
  # Customize the plot
  labs(title = "Thuja plicata WorldClim RDA Biplot", x = "RDA1", y = "RDA2", color = "Subspecies") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("Tplicata_RDA.png", plot = Tplicata_rda_plot, width = 8, height = 6, dpi = 300)

## Denovo
# Extract site scores (samples)
denovo_site_scores <- as.data.frame(scores(denovo_rda, display = "sites"))
denovo_site_scores$Sample <- as.data.frame(scores(denovo_rda, display = "sites"))
denovo_site_scores$pops <- factor(pops)

# Extract biplot scores (explanatory variables) - bp arg = biplot
denovo_biplot_scores <- as.data.frame(scores(denovo_rda, display = "bp"))

# Determine if any variables were dropped from the model
rownames(denovo_biplot_scores)
denovo_biplot_scores$variable <- rownames(denovo_biplot_scores)
# Yep, Precipitation_of_Driest_Quarter was dropped, so let's remove it

# Plot RDA
denovo_rda_plot <- ggplot() +
  # Plot site scores (samples)
  geom_point(data = denovo_site_scores, aes(x = RDA1, y = RDA2, color = pops), size = 2) +
  # Add arrows for biplot scores (explanatory variables)
  geom_segment(data = denovo_biplot_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
                                                arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "darkblue") +
  # Add labels for biplot scores
  geom_text_repel(data = denovo_biplot_scores, aes(x = RDA1, y = RDA2, label = rownames(denovo_biplot_scores)), size = 1, color = "blue") +
  # Customize the plot
  labs(title = "De novo WorldClim RDA Biplot", x = "RDA1", y = "RDA2", color = "Subspecies") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("denovo_RDA.png", plot = denovo_rda_plot, width = 8, height = 6, dpi = 300)

##### Print RDA results #####
print("Cobtusa RDA Results")
summary(Cobtusa_rda)
RsquareAdj(Cobtusa_rda)

print("Tplicata RDA Results")
summary(Tplicata_rda)
RsquareAdj(Tplicata_rda)

print("denovo RDA Results")
summary(denovo_rda)
RsquareAdj(denovo_rda)

```

</p>
</details> 

***

<details><summary> 14. LEA: LFMM </summary>
<p>

1. I first copied over my VCFs and unzipped them (compressed VCFs aren't compatible with LEA).

2. Next, I mined climate data from WorldClim for each coordinate point.

3. I ran a GWAS-style latent factor mixed model to evaluate Genotype Environment Associations (GEA) for each reference approach. Furthermore, I ran this on the overall dataset, the henryae subset, and the thyoides subset with the following three scripts:

* All samples
```
# Define a CRAN mirror to download packages from.
cran_mirrors <- c('http://cran.us.r-project.org', 'https://cloud.r-project.org/')

#library(devtools)
#library(remotes)
#devtools::find_rtools()
#devtools::install_github('bcm-uga/lfmm')

#install.packages('RSpectra')

# Specify libraries
packages_to_load <- c('BiocManager', 'tidyr', 'vcfR', 'LEA', 'ggplot2', 'dplyr', 'lfmm')

setwd("C:/Users/zaerr/Downloads/test_LEA")

# Install packages, if required.
#for (package in packages_to_load) {
#    if (!require(package, quietly = TRUE))
#        install.packages(packages_to_load, dependencies = FALSE, repos = cran_mirrors)
#}

### Add a function to read in env files.
# Function to read the environmental variable files
read_env_file <- function(filename) {
  read.csv(filename, header = FALSE)
}

#BiocManager::install("vcfR")

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

Cobtusa_all_vcf <- "Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"
Tplicata_all_vcf <- "Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"
denovo_all_vcf <- "Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"

vcf2geno(Cobtusa_all_vcf, "Cobtusa_all.geno")
vcf2geno(Tplicata_all_vcf, "Tplicata_all.geno")
vcf2geno(denovo_all_vcf, "denovo_all.geno")

vcf2lfmm(Cobtusa_all_vcf, output.file = "Cobtusa_all.lfmm", force = TRUE)
vcf2lfmm(Tplicata_all_vcf, output.file = "Tplicata_all.lfmm", force = TRUE)
vcf2lfmm(denovo_all_vcf, output.file = "denovo_all.lfmm", force = TRUE)

Cobtusa_all_geno <- LEA::read.geno("Cobtusa_all.geno")
Cobtusa_all_geno
Tplicata_all_geno <- LEA::read.geno("Tplicata_all.geno")
Tplicata_all_geno
denovo_all_geno <- LEA::read.geno("denovo_all.geno")
denovo_all_geno

#Y <- Cobtusa_all_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

#Y <- Tplicata_all_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

#Y <- denovo_all_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

denovo_all.lfmm <- "Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
#denovo_all.lfmm
Tplicata_all.lfmm <- "Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
#Tplicata_all.lfmm
Cobtusa_all.lfmm <- "Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
#Cobtusa_all.lfmm

denovo_all.snmf <- snmf("denovo_all.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#denovo_all.snmf
Tplicata_all.snmf <- snmf("Tplicata_all.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#Tplicata_all.snmf
Cobtusa_all.snmf <- snmf("Cobtusa_all.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#Cobtusa_all.snmf

denovo_all_best <- which.min(cross.entropy(denovo_all.snmf, K = 3))
#denovo_all_best
Tplicata_all_best <- which.min(cross.entropy(Tplicata_all.snmf, K = 3))
#Tplicata_all_best
Cobtusa_all_best <- which.min(cross.entropy(Cobtusa_all.snmf, K = 3))
#Cobtusa_all_best

impute(denovo_all.snmf, denovo_all.lfmm, method = "random", K = 3, run = denovo_all_best)
impute(Tplicata_all.snmf, Tplicata_all.lfmm, method = "random", K = 3, run = Tplicata_all_best)
impute(Cobtusa_all.snmf, Cobtusa_all.lfmm, method = "random", K = 3, run = Cobtusa_all_best)

all_bioclim_vars <- read.csv("all_bioclim-vars.csv", header = TRUE)
all_bioclim_vars[2:20] <- scale(all_bioclim_vars[2:20])
#all_bioclim_vars <- as.data.frame(all_bioclim_vars) %>% drop_na()

### Subsp. thyoides
all_Annual_Mean_Temperature <- all_bioclim_vars[2]
all_Mean_Diurnal_Range <- all_bioclim_vars[3]
all_Isothermality <- all_bioclim_vars[4]
all_Temperature_Seasonality <- all_bioclim_vars[5]
all_Max_Temperature_of_Warmest_Month <- all_bioclim_vars[6]
all_Min_Temperature_of_Coldest_Month <- all_bioclim_vars[7]
all_Temperature_Annual_Range <- all_bioclim_vars[8]
all_Mean_Temperature_of_Wettest_Quarter <- all_bioclim_vars[9]
all_Mean_Temperature_of_Driest_Quarter <- all_bioclim_vars[10]
all_Mean_Temperature_of_Warmest_Quarter <- all_bioclim_vars[11]
all_Mean_Temperature_of_Coldest_Quarter <- all_bioclim_vars[12]
all_Annual_Precipitation <- all_bioclim_vars[13]
all_Precipitation_of_Wettest_Month <- all_bioclim_vars[14]
all_Precipitation_of_Driest_Month <- all_bioclim_vars[15]
all_Precipitation_Seasonality <- all_bioclim_vars[16]
all_Precipitation_of_Wettest_Quarter <- all_bioclim_vars[17]
all_Precipitation_of_Driest_Quarter <- all_bioclim_vars[18]
all_Precipitation_of_Warmest_Quarter <- all_bioclim_vars[19]
all_Precipitation_of_Coldest_Quarter <- all_bioclim_vars[20]

LEA::write.env(all_Annual_Mean_Temperature, "all_Annual_Mean_Temperature.env")
LEA::write.env(all_Mean_Diurnal_Range, "all_Mean_Diurnal_Range.env")
LEA::write.env(all_Isothermality, "all_Isothermality.env")
LEA::write.env(all_Temperature_Seasonality, "all_Temperature_Seasonality.env")
LEA::write.env(all_Max_Temperature_of_Warmest_Month, "all_Max_Temperature_of_Warmest_Month.env")
LEA::write.env(all_Min_Temperature_of_Coldest_Month, "all_Min_Temperature_of_Coldest_Month.env")
LEA::write.env(all_Temperature_Annual_Range, "all_Temperature_Annual_Range.env")
LEA::write.env(all_Mean_Temperature_of_Wettest_Quarter, "all_Mean_Temperature_of_Wettest_Quarter.env")
LEA::write.env(all_Mean_Temperature_of_Driest_Quarter, "all_Mean_Temperature_of_Driest_Quarter.env")
LEA::write.env(all_Mean_Temperature_of_Warmest_Quarter, "all_Mean_Temperature_of_Warmest_Quarter.env")
LEA::write.env(all_Mean_Temperature_of_Coldest_Quarter, "all_Mean_Temperature_of_Coldest_Quarter.env")
LEA::write.env(all_Annual_Precipitation, "all_Annual_Precipitation.env")
LEA::write.env(all_Precipitation_of_Wettest_Month, "all_Precipitation_of_Wettest_Month.env")
LEA::write.env(all_Precipitation_of_Driest_Month, "all_Precipitation_of_Driest_Month.env")
LEA::write.env(all_Precipitation_Seasonality, "all_Precipitation_Seasonality.env")
LEA::write.env(all_Precipitation_of_Wettest_Quarter, "all_Precipitation_of_Wettest_Quarter.env")
LEA::write.env(all_Precipitation_of_Driest_Quarter, "all_Precipitation_of_Driest_Quarter.env")
LEA::write.env(all_Precipitation_of_Warmest_Quarter, "all_Precipitation_of_Warmest_Quarter.env")
LEA::write.env(all_Precipitation_of_Coldest_Quarter, "all_Precipitation_of_Coldest_Quarter.env")

all_Annual_Mean_Temperature.env <- "all_Annual_Mean_Temperature.env"
all_Mean_Diurnal_Range.env <- "all_Mean_Diurnal_Range.env"
all_Isothermality.env <- "all_Isothermality.env"
all_Temperature_Seasonality.env <- "all_Temperature_Seasonality.env"
all_Max_Temperature_of_Warmest_Month.env <- "all_Max_Temperature_of_Warmest_Month.env"
all_Min_Temperature_of_Coldest_Month.env <- "all_Min_Temperature_of_Coldest_Month.env"
all_Temperature_Annual_Range.env <- "all_Temperature_Annual_Range.env"
all_Mean_Temperature_of_Wettest_Quarter.env <- "all_Mean_Temperature_of_Wettest_Quarter.env"
all_Mean_Temperature_of_Driest_Quarter.env <- "all_Mean_Temperature_of_Driest_Quarter.env"
all_Mean_Temperature_of_Warmest_Quarter.env <- "all_Mean_Temperature_of_Warmest_Quarter.env"
all_Mean_Temperature_of_Coldest_Quarter.env <- "all_Mean_Temperature_of_Coldest_Quarter.env"
all_Annual_Precipitation.env <- "all_Annual_Precipitation.env"
all_Precipitation_of_Wettest_Month.env <- "all_Precipitation_of_Wettest_Month.env"
all_Precipitation_of_Driest_Month.env <- "all_Precipitation_of_Driest_Month.env"
all_Precipitation_Seasonality.env <- "all_Precipitation_Seasonality.env"
all_Precipitation_of_Wettest_Quarter.env <- "all_Precipitation_of_Wettest_Quarter.env"
all_Precipitation_of_Driest_Quarter.env <- "all_Precipitation_of_Driest_Quarter.env"
all_Precipitation_of_Warmest_Quarter.env <- "all_Precipitation_of_Warmest_Quarter.env"
all_Precipitation_of_Coldest_Quarter.env <- "all_Precipitation_of_Coldest_Quarter.env"

all_bioclim_vars_list <- list(
  "all_Annual_Mean_Temperature.env",
  "all_Mean_Diurnal_Range.env",
  "all_Isothermality.env",
  "all_Temperature_Seasonality.env",
  "all_Max_Temperature_of_Warmest_Month.env",
  "all_Min_Temperature_of_Coldest_Month.env",
  "all_Temperature_Annual_Range.env",
  "all_Mean_Temperature_of_Wettest_Quarter.env",
  "all_Mean_Temperature_of_Driest_Quarter.env",
  "all_Mean_Temperature_of_Warmest_Quarter.env",
  "all_Mean_Temperature_of_Coldest_Quarter.env",
  "all_Annual_Precipitation.env",
  "all_Precipitation_of_Wettest_Month.env",
  "all_Precipitation_of_Driest_Month.env",
  "all_Precipitation_Seasonality.env",
  "all_Precipitation_of_Wettest_Quarter.env",
  "all_Precipitation_of_Driest_Quarter.env",
  "all_Precipitation_of_Warmest_Quarter.env",
  "all_Precipitation_of_Coldest_Quarter.env"
)

#lfmm2geno(input.file, output.file = NULL, force = TRUE)
denovo_all_imputed.lfmm <- "Cthyoides_denovo-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
denovo_all_imputed.lfmm
lfmm2geno(denovo_all_imputed.lfmm, output.file = "denovo_all_imputed.geno", force = TRUE)
denovo_all_imputed.geno <- read.geno("denovo_all_imputed.geno")

Tplicata_all_imputed.lfmm <- "Cthyoides_Tplicata-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
Tplicata_all_imputed.lfmm
lfmm2geno(Tplicata_all_imputed.lfmm, output.file = "Tplicata_all_imputed.geno", force = TRUE)
Tplicata_all_imputed.geno <- read.geno("Tplicata_all_imputed.geno")

Cobtusa_all_imputed.lfmm <- "Cthyoides_Cobtusa-refaligned_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
Cobtusa_all_imputed.lfmm
lfmm2geno(Cobtusa_all_imputed.lfmm, output.file = "Cobtusa_all_imputed.geno", force = TRUE)
Cobtusa_all_imputed.geno <- read.geno("Cobtusa_all_imputed.geno")

###############################
########## thyoides ###########
###############################
options(digits = 22)

# Create/initialize lists
input_files <- list(denovo_all_imputed.geno, Tplicata_all_imputed.geno, Cobtusa_all_imputed.geno)
outfiles <- list("denovo_all_LFMM", "Tplicata_all_LFMM", "Cobtusa_all_LFMM")
all_results_list <- list()

for (i in seq_along(input_files)) {  
  
  for (var in all_bioclim_vars_list) {
    
    outfile <- outfiles[[i]]
    geno <- input_files[[i]]
    
    #Print variable being processed
    print(var)
    
    # Read the environmental variable file
    env_data <- read_env_file(var)
    
    # Run LFMM
    test.mod.k3 <- lfmm::lfmm_ridge(Y = geno, X = env_data, K = 3)
    
    # Run the LFMM test for each bioclimatic variable
    test_result <- lfmm::lfmm_test(lfmm = test.mod.k3, Y = geno, X = env_data, calibrate = "gif")
    
    # Ensure p-values are properly transformed for plotting
    calibrated_pvalues <- data.frame(test_result$calibrated.pvalue)
    
    # Ensure p-values are properly transformed for plotting
    neg_log10_pvalues <- data.frame(-log10(test_result$calibrated.pvalue))
    
    # Add significant column based on a p-value threshold
    threshold <- 5E-08
    significant <- ifelse(calibrated_pvalues < threshold, "Significant", "Not Significant")
    
    # Convert significant to a factor to ensure proper handling in ggplot
    significant <- factor(significant, levels = c("Not Significant", "Significant"))
    
    # Create index for SNPs
    index <- 1:length(test_result$calibrated.pvalue)
    
    # Prepare data frame for plotting
    results_df <- data.frame(index = index, calibrated_pvalues = calibrated_pvalues, neg_log10_pvalues = neg_log10_pvalues, significant = significant, variable = var)
    colnames(results_df) <- c("index", "calibrated_pvalues", "neg_log10_pvalues", "significant", "variable")
    results_df
    # Store results in list
    all_results_list[[var]] <- results_df
    
  }
  
  # Combine all results into one data frame
  all_results_df <- bind_rows(all_results_list)
  all_results_df
  
  # Define color palette for significant points
  significant_colors <- c(
    "all_Annual_Mean_Temperature.env" = "red",
    "all_Mean_Diurnal_Range.env" = "blue",
    "all_Isothermality.env" = "green",
    "all_Temperature_Seasonality.env" = "purple",
    "all_Max_Temperature_of_Warmest_Month.env" = "orange",
    "all_Min_Temperature_of_Coldest_Month.env" = "brown",
    "all_Temperature_Annual_Range.env" = "pink",
    "all_Mean_Temperature_of_Wettest_Quarter.env" = "yellow",
    "all_Mean_Temperature_of_Driest_Quarter.env" = "cyan",
    "all_Mean_Temperature_of_Warmest_Quarter.env" = "magenta",
    "all_Mean_Temperature_of_Coldest_Quarter.env" = "darkgreen",
    "all_Annual_Precipitation.env" = "darkblue",
    "all_Precipitation_of_Wettest_Month.env" = "darkred",
    "all_Precipitation_of_Driest_Month.env" = "darkorange",
    "all_Precipitation_Seasonality.env" = "hotpink",
    "all_Precipitation_of_Wettest_Quarter.env" = "darkcyan",
    "all_Precipitation_of_Driest_Quarter.env" = "darkmagenta",
    "all_Precipitation_of_Warmest_Quarter.env" = "lightgreen",
    "all_Precipitation_of_Coldest_Quarter.env" = "lightblue"
  )
  
  # Create names for significant and non-significant combined palette
  significant_names <- paste0(names(significant_colors), ".Significant")
  combined_palette <- c(setNames(significant_colors, significant_names), "Not Significant" = "gray")
  
  all_results_df_noNA <- na.omit(all_results_df)
  glimpse(all_results_df_noNA)
  
  # Plot the results with ggplot2
  lfmm_faceted <- ggplot(all_results_df_noNA, aes(x = index, y = neg_log10_pvalues)) +
    geom_point(aes(color = ifelse(significant == "Significant", paste0(variable, ".Significant"), "Not Significant")), size = 1.5) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = combined_palette) +
    labs(title = "Significant SNPs by environment", x = "SNP Index", y = "-log10(p-value)") +
    theme_minimal()
  
  # Extract the base name of the LFMM file (remove directory and extension)
  #lfmm_base_name <- tools::file_path_sans_ext(basename("geno"))
  
  # Generate the output file name
  output_file <- paste0(outfile, ".png")
  
  # Save the plot
  ggsave(output_file, plot = lfmm_faceted, width = 32, height = 24, dpi = 300)
  
}

```

* Subsp. henryae
```
# Define a CRAN mirror to download packages from.
cran_mirrors <- c('http://cran.us.r-project.org', 'https://cloud.r-project.org/')

#library(devtools)
#library(remotes)
#devtools::find_rtools()
#devtools::install_github('bcm-uga/lfmm')

#install.packages('RSpectra')

# Specify libraries
packages_to_load <- c('BiocManager', 'tidyr', 'vcfR', 'LEA', 'ggplot2', 'dplyr', 'lfmm')

setwd("C:/Users/zaerr/Downloads/test_LEA")

# Install packages, if required.
#for (package in packages_to_load) {
#    if (!require(package, quietly = TRUE))
#        install.packages(packages_to_load, dependencies = FALSE, repos = cran_mirrors)
#}

### Add a function to read in env files.
# Function to read the environmental variable files
read_env_file <- function(filename) {
  read.csv(filename, header = FALSE)
}

#BiocManager::install("vcfR")

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

Cobtusa_henryae_vcf <- "Cthyoides_Cobtusa-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"
Tplicata_henryae_vcf <- "Cthyoides_Tplicata-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"
denovo_henryae_vcf <- "Cthyoides_denovo-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"

vcf2geno(Cobtusa_henryae_vcf, "Cobtusa_henryae.geno")
vcf2geno(Tplicata_henryae_vcf, "Tplicata_henryae.geno")
vcf2geno(denovo_henryae_vcf, "denovo_henryae.geno")

vcf2lfmm(Cobtusa_henryae_vcf, output.file = "Cobtusa_henryae.lfmm", force = TRUE)
vcf2lfmm(Tplicata_henryae_vcf, output.file = "Tplicata_henryae.lfmm", force = TRUE)
vcf2lfmm(denovo_henryae_vcf, output.file = "denovo_henryae.lfmm", force = TRUE)

Cobtusa_henryae_geno <- LEA::read.geno("Cobtusa_henryae.geno")
Cobtusa_henryae_geno
Tplicata_henryae_geno <- LEA::read.geno("Tplicata_henryae.geno")
Tplicata_henryae_geno
denovo_henryae_geno <- LEA::read.geno("denovo_henryae.geno")
denovo_henryae_geno

#Y <- Cobtusa_henryae_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

#Y <- Tplicata_henryae_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

#Y <- denovo_henryae_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

denovo_henryae.lfmm <- "Cthyoides_denovo-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
denovo_henryae.lfmm
Tplicata_henryae.lfmm <- "Cthyoides_Tplicata-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
Tplicata_henryae.lfmm
Cobtusa_henryae.lfmm <- "Cthyoides_Cobtusa-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
Cobtusa_henryae.lfmm

denovo_henryae.snmf <- snmf("denovo_henryae.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#denovo_henryae.snmf
Tplicata_henryae.snmf <- snmf("Tplicata_henryae.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#Tplicata_henryae.snmf
Cobtusa_henryae.snmf <- snmf("Cobtusa_henryae.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#Cobtusa_henryae.snmf

denovo_henryae_best <- which.min(cross.entropy(denovo_henryae.snmf, K = 3))
#denovo_henryae_best
Tplicata_henryae_best <- which.min(cross.entropy(Tplicata_henryae.snmf, K = 3))
#Tplicata_henryae_best
Cobtusa_henryae_best <- which.min(cross.entropy(Cobtusa_henryae.snmf, K = 3))
#Cobtusa_henryae_best

impute(denovo_henryae.snmf, denovo_henryae.lfmm, method = "random", K = 3, run = denovo_henryae_best)
impute(Tplicata_henryae.snmf, Tplicata_henryae.lfmm, method = "random", K = 3, run = Tplicata_henryae_best)
impute(Cobtusa_henryae.snmf, Cobtusa_henryae.lfmm, method = "random", K = 3, run = Cobtusa_henryae_best)

henryae_bioclim_vars <- read.csv("henryae_bioclim-vars.csv", header = TRUE)
henryae_bioclim_vars[2:20] <- scale(henryae_bioclim_vars[2:20])
#henryae_bioclim_vars <- as.data.frame(henryae_bioclim_vars) %>% drop_na()

### Subsp. thyoides
henryae_Annual_Mean_Temperature <- henryae_bioclim_vars[2]
henryae_Mean_Diurnal_Range <- henryae_bioclim_vars[3]
henryae_Isothermality <- henryae_bioclim_vars[4]
henryae_Temperature_Seasonality <- henryae_bioclim_vars[5]
henryae_Max_Temperature_of_Warmest_Month <- henryae_bioclim_vars[6]
henryae_Min_Temperature_of_Coldest_Month <- henryae_bioclim_vars[7]
henryae_Temperature_Annual_Range <- henryae_bioclim_vars[8]
henryae_Mean_Temperature_of_Wettest_Quarter <- henryae_bioclim_vars[9]
henryae_Mean_Temperature_of_Driest_Quarter <- henryae_bioclim_vars[10]
henryae_Mean_Temperature_of_Warmest_Quarter <- henryae_bioclim_vars[11]
henryae_Mean_Temperature_of_Coldest_Quarter <- henryae_bioclim_vars[12]
henryae_Annual_Precipitation <- henryae_bioclim_vars[13]
henryae_Precipitation_of_Wettest_Month <- henryae_bioclim_vars[14]
henryae_Precipitation_of_Driest_Month <- henryae_bioclim_vars[15]
henryae_Precipitation_Seasonality <- henryae_bioclim_vars[16]
henryae_Precipitation_of_Wettest_Quarter <- henryae_bioclim_vars[17]
henryae_Precipitation_of_Driest_Quarter <- henryae_bioclim_vars[18]
henryae_Precipitation_of_Warmest_Quarter <- henryae_bioclim_vars[19]
henryae_Precipitation_of_Coldest_Quarter <- henryae_bioclim_vars[20]

LEA::write.env(henryae_Annual_Mean_Temperature, "henryae_Annual_Mean_Temperature.env")
LEA::write.env(henryae_Mean_Diurnal_Range, "henryae_Mean_Diurnal_Range.env")
LEA::write.env(henryae_Isothermality, "henryae_Isothermality.env")
LEA::write.env(henryae_Temperature_Seasonality, "henryae_Temperature_Seasonality.env")
LEA::write.env(henryae_Max_Temperature_of_Warmest_Month, "henryae_Max_Temperature_of_Warmest_Month.env")
LEA::write.env(henryae_Min_Temperature_of_Coldest_Month, "henryae_Min_Temperature_of_Coldest_Month.env")
LEA::write.env(henryae_Temperature_Annual_Range, "henryae_Temperature_Annual_Range.env")
LEA::write.env(henryae_Mean_Temperature_of_Wettest_Quarter, "henryae_Mean_Temperature_of_Wettest_Quarter.env")
LEA::write.env(henryae_Mean_Temperature_of_Driest_Quarter, "henryae_Mean_Temperature_of_Driest_Quarter.env")
LEA::write.env(henryae_Mean_Temperature_of_Warmest_Quarter, "henryae_Mean_Temperature_of_Warmest_Quarter.env")
LEA::write.env(henryae_Mean_Temperature_of_Coldest_Quarter, "henryae_Mean_Temperature_of_Coldest_Quarter.env")
LEA::write.env(henryae_Annual_Precipitation, "henryae_Annual_Precipitation.env")
LEA::write.env(henryae_Precipitation_of_Wettest_Month, "henryae_Precipitation_of_Wettest_Month.env")
LEA::write.env(henryae_Precipitation_of_Driest_Month, "henryae_Precipitation_of_Driest_Month.env")
LEA::write.env(henryae_Precipitation_Seasonality, "henryae_Precipitation_Seasonality.env")
LEA::write.env(henryae_Precipitation_of_Wettest_Quarter, "henryae_Precipitation_of_Wettest_Quarter.env")
LEA::write.env(henryae_Precipitation_of_Driest_Quarter, "henryae_Precipitation_of_Driest_Quarter.env")
LEA::write.env(henryae_Precipitation_of_Warmest_Quarter, "henryae_Precipitation_of_Warmest_Quarter.env")
LEA::write.env(henryae_Precipitation_of_Coldest_Quarter, "henryae_Precipitation_of_Coldest_Quarter.env")

henryae_Annual_Mean_Temperature.env <- "henryae_Annual_Mean_Temperature.env"
henryae_Mean_Diurnal_Range.env <- "henryae_Mean_Diurnal_Range.env"
henryae_Isothermality.env <- "henryae_Isothermality.env"
henryae_Temperature_Seasonality.env <- "henryae_Temperature_Seasonality.env"
henryae_Max_Temperature_of_Warmest_Month.env <- "henryae_Max_Temperature_of_Warmest_Month.env"
henryae_Min_Temperature_of_Coldest_Month.env <- "henryae_Min_Temperature_of_Coldest_Month.env"
henryae_Temperature_Annual_Range.env <- "henryae_Temperature_Annual_Range.env"
henryae_Mean_Temperature_of_Wettest_Quarter.env <- "henryae_Mean_Temperature_of_Wettest_Quarter.env"
henryae_Mean_Temperature_of_Driest_Quarter.env <- "henryae_Mean_Temperature_of_Driest_Quarter.env"
henryae_Mean_Temperature_of_Warmest_Quarter.env <- "henryae_Mean_Temperature_of_Warmest_Quarter.env"
henryae_Mean_Temperature_of_Coldest_Quarter.env <- "henryae_Mean_Temperature_of_Coldest_Quarter.env"
henryae_Annual_Precipitation.env <- "henryae_Annual_Precipitation.env"
henryae_Precipitation_of_Wettest_Month.env <- "henryae_Precipitation_of_Wettest_Month.env"
henryae_Precipitation_of_Driest_Month.env <- "henryae_Precipitation_of_Driest_Month.env"
henryae_Precipitation_Seasonality.env <- "henryae_Precipitation_Seasonality.env"
henryae_Precipitation_of_Wettest_Quarter.env <- "henryae_Precipitation_of_Wettest_Quarter.env"
henryae_Precipitation_of_Driest_Quarter.env <- "henryae_Precipitation_of_Driest_Quarter.env"
henryae_Precipitation_of_Warmest_Quarter.env <- "henryae_Precipitation_of_Warmest_Quarter.env"
henryae_Precipitation_of_Coldest_Quarter.env <- "henryae_Precipitation_of_Coldest_Quarter.env"

henryae_bioclim_vars_list <- list(
  "henryae_Annual_Mean_Temperature.env",
  "henryae_Mean_Diurnal_Range.env",
  "henryae_Isothermality.env",
  "henryae_Temperature_Seasonality.env",
  "henryae_Max_Temperature_of_Warmest_Month.env",
  "henryae_Min_Temperature_of_Coldest_Month.env",
  "henryae_Temperature_Annual_Range.env",
  "henryae_Mean_Temperature_of_Wettest_Quarter.env",
  "henryae_Mean_Temperature_of_Driest_Quarter.env",
  "henryae_Mean_Temperature_of_Warmest_Quarter.env",
  "henryae_Mean_Temperature_of_Coldest_Quarter.env",
  "henryae_Annual_Precipitation.env",
  "henryae_Precipitation_of_Wettest_Month.env",
  "henryae_Precipitation_of_Driest_Month.env",
  "henryae_Precipitation_Seasonality.env",
  "henryae_Precipitation_of_Wettest_Quarter.env",
  "henryae_Precipitation_of_Driest_Quarter.env",
  "henryae_Precipitation_of_Warmest_Quarter.env",
  "henryae_Precipitation_of_Coldest_Quarter.env"
)

#lfmm2geno(input.file, output.file = NULL, force = TRUE)
denovo_henryae_imputed.lfmm <- "Cthyoides_denovo-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
denovo_henryae_imputed.lfmm
lfmm2geno(denovo_henryae_imputed.lfmm, output.file = "denovo_henryae_imputed.geno", force = TRUE)
denovo_henryae_imputed.geno <- read.geno("denovo_henryae_imputed.geno")

Tplicata_henryae_imputed.lfmm <- "Cthyoides_Tplicata-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
Tplicata_henryae_imputed.lfmm
lfmm2geno(Tplicata_henryae_imputed.lfmm, output.file = "Tplicata_henryae_imputed.geno", force = TRUE)
Tplicata_henryae_imputed.geno <- read.geno("Tplicata_henryae_imputed.geno")

Cobtusa_henryae_imputed.lfmm <- "Cthyoides_Cobtusa-refaligned_henryae_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
Cobtusa_henryae_imputed.lfmm
lfmm2geno(Cobtusa_henryae_imputed.lfmm, output.file = "Cobtusa_henryae_imputed.geno", force = TRUE)
Cobtusa_henryae_imputed.geno <- read.geno("Cobtusa_henryae_imputed.geno")

###############################
########## thyoides ###########
###############################
options(digits = 22)

# Create/initialize lists
input_files <- list(denovo_henryae_imputed.geno, Tplicata_henryae_imputed.geno, Cobtusa_henryae_imputed.geno)
outfiles <- list("denovo_henryae_LFMM", "Tplicata_henryae_LFMM", "Cobtusa_henryae_LFMM")
henryae_results_list <- list()

for (i in seq_along(input_files)) {  
  
  for (var in henryae_bioclim_vars_list) {
    
    outfile <- outfiles[[i]]
    geno <- input_files[[i]]
    
    #Print variable being processed
    print(var)
    
    # Read the environmental variable file
    env_data <- read_env_file(var)
    
    # Run LFMM
    test.mod.k3 <- lfmm::lfmm_ridge(Y = geno, X = env_data, K = 3)
    
    # Run the LFMM test for each bioclimatic variable
    test_result <- lfmm::lfmm_test(lfmm = test.mod.k3, Y = geno, X = env_data, calibrate = "gif")
    
    # Ensure p-values are properly transformed for plotting
    calibrated_pvalues <- data.frame(test_result$calibrated.pvalue)
    
    # Ensure p-values are properly transformed for plotting
    neg_log10_pvalues <- data.frame(-log10(test_result$calibrated.pvalue))
    
    # Add significant column based on a p-value threshold
    threshold <- 5E-08
    significant <- ifelse(calibrated_pvalues < threshold, "Significant", "Not Significant")
    
    # Convert significant to a factor to ensure proper handling in ggplot
    significant <- factor(significant, levels = c("Not Significant", "Significant"))
    
    # Create index for SNPs
    index <- 1:length(test_result$calibrated.pvalue)
    
    # Prepare data frame for plotting
    results_df <- data.frame(index = index, calibrated_pvalues = calibrated_pvalues, neg_log10_pvalues = neg_log10_pvalues, significant = significant, variable = var)
    colnames(results_df) <- c("index", "calibrated_pvalues", "neg_log10_pvalues", "significant", "variable")
    results_df
    # Store results in list
    henryae_results_list[[var]] <- results_df
    
  }
  
  # Combine all results into one data frame
  henryae_results_df <- bind_rows(henryae_results_list)
  henryae_results_df
  
  # Define color palette for significant points
  significant_colors <- c(
    "henryae_Annual_Mean_Temperature.env" = "red",
    "henryae_Mean_Diurnal_Range.env" = "blue",
    "henryae_Isothermality.env" = "green",
    "henryae_Temperature_Seasonality.env" = "purple",
    "henryae_Max_Temperature_of_Warmest_Month.env" = "orange",
    "henryae_Min_Temperature_of_Coldest_Month.env" = "brown",
    "henryae_Temperature_Annual_Range.env" = "pink",
    "henryae_Mean_Temperature_of_Wettest_Quarter.env" = "yellow",
    "henryae_Mean_Temperature_of_Driest_Quarter.env" = "cyan",
    "henryae_Mean_Temperature_of_Warmest_Quarter.env" = "magenta",
    "henryae_Mean_Temperature_of_Coldest_Quarter.env" = "darkgreen",
    "henryae_Annual_Precipitation.env" = "darkblue",
    "henryae_Precipitation_of_Wettest_Month.env" = "darkred",
    "henryae_Precipitation_of_Driest_Month.env" = "darkorange",
    "henryae_Precipitation_Seasonality.env" = "hotpink",
    "henryae_Precipitation_of_Wettest_Quarter.env" = "darkcyan",
    "henryae_Precipitation_of_Driest_Quarter.env" = "darkmagenta",
    "henryae_Precipitation_of_Warmest_Quarter.env" = "lightgreen",
    "henryae_Precipitation_of_Coldest_Quarter.env" = "lightblue"
  )
  
  # Create names for significant and non-significant combined palette
  significant_names <- paste0(names(significant_colors), ".Significant")
  combined_palette <- c(setNames(significant_colors, significant_names), "Not Significant" = "gray")
  
  henryae_results_df_noNA <- na.omit(henryae_results_df)
  glimpse(henryae_results_df_noNA)
  
  # Plot the results with ggplot2
  lfmm_faceted <- ggplot(henryae_results_df_noNA, aes(x = index, y = neg_log10_pvalues)) +
    geom_point(aes(color = ifelse(significant == "Significant", paste0(variable, ".Significant"), "Not Significant")), size = 1.5) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = combined_palette) +
    labs(title = "Significant SNPs by environment", x = "SNP Index", y = "-log10(p-value)") +
    theme_minimal()
  
  # Extract the base name of the LFMM file (remove directory and extension)
  #lfmm_base_name <- tools::file_path_sans_ext(basename("geno"))
  
  # Generate the output file name
  output_file <- paste0(outfile, ".png")
  
  # Save the plot
  ggsave(output_file, plot = lfmm_faceted, width = 32, height = 24, dpi = 300)
  
}

```

* Subsp. thyoides
```
# Define a CRAN mirror to download packages from.
cran_mirrors <- c('http://cran.us.r-project.org', 'https://cloud.r-project.org/')

#library(devtools)
#library(remotes)
#devtools::find_rtools()
#devtools::install_github('bcm-uga/lfmm')

#install.packages('RSpectra')

# Specify libraries
packages_to_load <- c('BiocManager', 'tidyr', 'vcfR', 'LEA', 'ggplot2', 'dplyr', 'lfmm')

setwd("C:/Users/zaerr/Downloads/test_LEA")

# Install packages, if required.
#for (package in packages_to_load) {
#    if (!require(package, quietly = TRUE))
#        install.packages(packages_to_load, dependencies = FALSE, repos = cran_mirrors)
#}

### Add a function to read in env files.
# Function to read the environmental variable files
read_env_file <- function(filename) {
    read.csv(filename, header = FALSE)
}

#BiocManager::install("vcfR")

# Load all the packages
sapply(packages_to_load, library, character.only = TRUE)

Cobtusa_thyoides_vcf <- "Cthyoides_Cobtusa-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"
Tplicata_thyoides_vcf <- "Cthyoides_Tplicata-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"
denovo_thyoides_vcf <- "Cthyoides_denovo-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.vcf"

vcf2geno(Cobtusa_thyoides_vcf, "Cobtusa_thyoides.geno")
vcf2geno(Tplicata_thyoides_vcf, "Tplicata_thyoides.geno")
vcf2geno(denovo_thyoides_vcf, "denovo_thyoides.geno")

vcf2lfmm(Cobtusa_thyoides_vcf, output.file = "Cobtusa_thyoides.lfmm", force = TRUE)
vcf2lfmm(Tplicata_thyoides_vcf, output.file = "Tplicata_thyoides.lfmm", force = TRUE)
vcf2lfmm(denovo_thyoides_vcf, output.file = "denovo_thyoides.lfmm", force = TRUE)

Cobtusa_thyoides_geno <- LEA::read.geno("Cobtusa_thyoides.geno")
Cobtusa_thyoides_geno
Tplicata_thyoides_geno <- LEA::read.geno("Tplicata_thyoides.geno")
Tplicata_thyoides_geno
denovo_thyoides_geno <- LEA::read.geno("denovo_thyoides.geno")
denovo_thyoides_geno

#Y <- Cobtusa_thyoides_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

#Y <- Tplicata_thyoides_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

#Y <- denovo_thyoides_geno
#pc <- prcomp(Y)
#plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
#points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")

denovo_thyoides.lfmm <- "Cthyoides_denovo-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
#denovo_thyoides.lfmm
Tplicata_thyoides.lfmm <- "Cthyoides_Tplicata-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
#Tplicata_thyoides.lfmm
Cobtusa_thyoides.lfmm <- "Cthyoides_Cobtusa-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm"
#Cobtusa_thyoides.lfmm

denovo_thyoides.snmf <- snmf("denovo_thyoides.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#denovo_thyoides.snmf
Tplicata_thyoides.snmf <- snmf("Tplicata_thyoides.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#Tplicata_thyoides.snmf
Cobtusa_thyoides.snmf <- snmf("Cobtusa_thyoides.geno", K = 3, entropy = TRUE, repetitions = 10, project = "new")
#Cobtusa_thyoides.snmf

denovo_thyoides_best <- which.min(cross.entropy(denovo_thyoides.snmf, K = 3))
#denovo_thyoides_best
Tplicata_thyoides_best <- which.min(cross.entropy(Tplicata_thyoides.snmf, K = 3))
#Tplicata_thyoides_best
Cobtusa_thyoides_best <- which.min(cross.entropy(Cobtusa_thyoides.snmf, K = 3))
#Cobtusa_thyoides_best

impute(denovo_thyoides.snmf, denovo_thyoides.lfmm, method = "random", K = 3, run = denovo_thyoides_best)
impute(Tplicata_thyoides.snmf, Tplicata_thyoides.lfmm, method = "random", K = 3, run = Tplicata_thyoides_best)
impute(Cobtusa_thyoides.snmf, Cobtusa_thyoides.lfmm, method = "random", K = 3, run = Cobtusa_thyoides_best)

thyoides_bioclim_vars <- read.csv("thyoides_bioclim-vars.csv", header = TRUE)
thyoides_bioclim_vars[2:20] <- scale(thyoides_bioclim_vars[2:20])
#thyoides_bioclim_vars <- as.data.frame(thyoides_bioclim_vars) %>% drop_na()

### Subsp. thyoides
thyoides_Annual_Mean_Temperature <- thyoides_bioclim_vars[2]
thyoides_Mean_Diurnal_Range <- thyoides_bioclim_vars[3]
thyoides_Isothermality <- thyoides_bioclim_vars[4]
thyoides_Temperature_Seasonality <- thyoides_bioclim_vars[5]
thyoides_Max_Temperature_of_Warmest_Month <- thyoides_bioclim_vars[6]
thyoides_Min_Temperature_of_Coldest_Month <- thyoides_bioclim_vars[7]
thyoides_Temperature_Annual_Range <- thyoides_bioclim_vars[8]
thyoides_Mean_Temperature_of_Wettest_Quarter <- thyoides_bioclim_vars[9]
thyoides_Mean_Temperature_of_Driest_Quarter <- thyoides_bioclim_vars[10]
thyoides_Mean_Temperature_of_Warmest_Quarter <- thyoides_bioclim_vars[11]
thyoides_Mean_Temperature_of_Coldest_Quarter <- thyoides_bioclim_vars[12]
thyoides_Annual_Precipitation <- thyoides_bioclim_vars[13]
thyoides_Precipitation_of_Wettest_Month <- thyoides_bioclim_vars[14]
thyoides_Precipitation_of_Driest_Month <- thyoides_bioclim_vars[15]
thyoides_Precipitation_Seasonality <- thyoides_bioclim_vars[16]
thyoides_Precipitation_of_Wettest_Quarter <- thyoides_bioclim_vars[17]
thyoides_Precipitation_of_Driest_Quarter <- thyoides_bioclim_vars[18]
thyoides_Precipitation_of_Warmest_Quarter <- thyoides_bioclim_vars[19]
thyoides_Precipitation_of_Coldest_Quarter <- thyoides_bioclim_vars[20]

LEA::write.env(thyoides_Annual_Mean_Temperature, "thyoides_Annual_Mean_Temperature.env")
LEA::write.env(thyoides_Mean_Diurnal_Range, "thyoides_Mean_Diurnal_Range.env")
LEA::write.env(thyoides_Isothermality, "thyoides_Isothermality.env")
LEA::write.env(thyoides_Temperature_Seasonality, "thyoides_Temperature_Seasonality.env")
LEA::write.env(thyoides_Max_Temperature_of_Warmest_Month, "thyoides_Max_Temperature_of_Warmest_Month.env")
LEA::write.env(thyoides_Min_Temperature_of_Coldest_Month, "thyoides_Min_Temperature_of_Coldest_Month.env")
LEA::write.env(thyoides_Temperature_Annual_Range, "thyoides_Temperature_Annual_Range.env")
LEA::write.env(thyoides_Mean_Temperature_of_Wettest_Quarter, "thyoides_Mean_Temperature_of_Wettest_Quarter.env")
LEA::write.env(thyoides_Mean_Temperature_of_Driest_Quarter, "thyoides_Mean_Temperature_of_Driest_Quarter.env")
LEA::write.env(thyoides_Mean_Temperature_of_Warmest_Quarter, "thyoides_Mean_Temperature_of_Warmest_Quarter.env")
LEA::write.env(thyoides_Mean_Temperature_of_Coldest_Quarter, "thyoides_Mean_Temperature_of_Coldest_Quarter.env")
LEA::write.env(thyoides_Annual_Precipitation, "thyoides_Annual_Precipitation.env")
LEA::write.env(thyoides_Precipitation_of_Wettest_Month, "thyoides_Precipitation_of_Wettest_Month.env")
LEA::write.env(thyoides_Precipitation_of_Driest_Month, "thyoides_Precipitation_of_Driest_Month.env")
LEA::write.env(thyoides_Precipitation_Seasonality, "thyoides_Precipitation_Seasonality.env")
LEA::write.env(thyoides_Precipitation_of_Wettest_Quarter, "thyoides_Precipitation_of_Wettest_Quarter.env")
LEA::write.env(thyoides_Precipitation_of_Driest_Quarter, "thyoides_Precipitation_of_Driest_Quarter.env")
LEA::write.env(thyoides_Precipitation_of_Warmest_Quarter, "thyoides_Precipitation_of_Warmest_Quarter.env")
LEA::write.env(thyoides_Precipitation_of_Coldest_Quarter, "thyoides_Precipitation_of_Coldest_Quarter.env")

thyoides_Annual_Mean_Temperature.env <- "thyoides_Annual_Mean_Temperature.env"
thyoides_Mean_Diurnal_Range.env <- "thyoides_Mean_Diurnal_Range.env"
thyoides_Isothermality.env <- "thyoides_Isothermality.env"
thyoides_Temperature_Seasonality.env <- "thyoides_Temperature_Seasonality.env"
thyoides_Max_Temperature_of_Warmest_Month.env <- "thyoides_Max_Temperature_of_Warmest_Month.env"
thyoides_Min_Temperature_of_Coldest_Month.env <- "thyoides_Min_Temperature_of_Coldest_Month.env"
thyoides_Temperature_Annual_Range.env <- "thyoides_Temperature_Annual_Range.env"
thyoides_Mean_Temperature_of_Wettest_Quarter.env <- "thyoides_Mean_Temperature_of_Wettest_Quarter.env"
thyoides_Mean_Temperature_of_Driest_Quarter.env <- "thyoides_Mean_Temperature_of_Driest_Quarter.env"
thyoides_Mean_Temperature_of_Warmest_Quarter.env <- "thyoides_Mean_Temperature_of_Warmest_Quarter.env"
thyoides_Mean_Temperature_of_Coldest_Quarter.env <- "thyoides_Mean_Temperature_of_Coldest_Quarter.env"
thyoides_Annual_Precipitation.env <- "thyoides_Annual_Precipitation.env"
thyoides_Precipitation_of_Wettest_Month.env <- "thyoides_Precipitation_of_Wettest_Month.env"
thyoides_Precipitation_of_Driest_Month.env <- "thyoides_Precipitation_of_Driest_Month.env"
thyoides_Precipitation_Seasonality.env <- "thyoides_Precipitation_Seasonality.env"
thyoides_Precipitation_of_Wettest_Quarter.env <- "thyoides_Precipitation_of_Wettest_Quarter.env"
thyoides_Precipitation_of_Driest_Quarter.env <- "thyoides_Precipitation_of_Driest_Quarter.env"
thyoides_Precipitation_of_Warmest_Quarter.env <- "thyoides_Precipitation_of_Warmest_Quarter.env"
thyoides_Precipitation_of_Coldest_Quarter.env <- "thyoides_Precipitation_of_Coldest_Quarter.env"

thyoides_bioclim_vars_list <- list(
    "thyoides_Annual_Mean_Temperature.env",
    "thyoides_Mean_Diurnal_Range.env",
    "thyoides_Isothermality.env",
    "thyoides_Temperature_Seasonality.env",
    "thyoides_Max_Temperature_of_Warmest_Month.env",
    "thyoides_Min_Temperature_of_Coldest_Month.env",
    "thyoides_Temperature_Annual_Range.env",
    "thyoides_Mean_Temperature_of_Wettest_Quarter.env",
    "thyoides_Mean_Temperature_of_Driest_Quarter.env",
    "thyoides_Mean_Temperature_of_Warmest_Quarter.env",
    "thyoides_Mean_Temperature_of_Coldest_Quarter.env",
    "thyoides_Annual_Precipitation.env",
    "thyoides_Precipitation_of_Wettest_Month.env",
    "thyoides_Precipitation_of_Driest_Month.env",
    "thyoides_Precipitation_Seasonality.env",
    "thyoides_Precipitation_of_Wettest_Quarter.env",
    "thyoides_Precipitation_of_Driest_Quarter.env",
    "thyoides_Precipitation_of_Warmest_Quarter.env",
    "thyoides_Precipitation_of_Coldest_Quarter.env"
)

#lfmm2geno(input.file, output.file = NULL, force = TRUE)
denovo_thyoides_imputed.lfmm <- "Cthyoides_denovo-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
denovo_thyoides_imputed.lfmm
lfmm2geno(denovo_thyoides_imputed.lfmm, output.file = "denovo_thyoides_imputed.geno", force = TRUE)
denovo_thyoides_imputed.geno <- read.geno("denovo_thyoides_imputed.geno")
  
Tplicata_thyoides_imputed.lfmm <- "Cthyoides_Tplicata-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
Tplicata_thyoides_imputed.lfmm
lfmm2geno(Tplicata_thyoides_imputed.lfmm, output.file = "Tplicata_thyoides_imputed.geno", force = TRUE)
Tplicata_thyoides_imputed.geno <- read.geno("Tplicata_thyoides_imputed.geno")

Cobtusa_thyoides_imputed.lfmm <- "Cthyoides_Cobtusa-refaligned_thyoides_GATKfilters_annotated_biallelic_maf01_minDP3_PASSonly_removed-indv_removed-paralogs_mm90_LD-50-5-0.4.lfmm_imputed.lfmm"
Cobtusa_thyoides_imputed.lfmm
lfmm2geno(Cobtusa_thyoides_imputed.lfmm, output.file = "Cobtusa_thyoides_imputed.geno", force = TRUE)
Cobtusa_thyoides_imputed.geno <- read.geno("Cobtusa_thyoides_imputed.geno")

###############################
########## thyoides ###########
###############################
options(digits = 22)

# Create/initialize lists
input_files <- list(denovo_thyoides_imputed.geno, Tplicata_thyoides_imputed.geno, Cobtusa_thyoides_imputed.geno)
outfiles <- list("denovo_thyoides_LFMM", "Tplicata_thyoides_LFMM", "Cobtusa_thyoides_LFMM")
thyoides_results_list <- list()

for (i in seq_along(input_files)) {  

  for (var in thyoides_bioclim_vars_list) {
      
    outfile <- outfiles[[i]]
    geno <- input_files[[i]]
      
    #Print variable being processed
    print(var)
    
    # Read the environmental variable file
    env_data <- read_env_file(var)
    
    # Run LFMM
    test.mod.k3 <- lfmm::lfmm_ridge(Y = geno, X = env_data, K = 3)
    
    # Run the LFMM test for each bioclimatic variable
    test_result <- lfmm::lfmm_test(lfmm = test.mod.k3, Y = geno, X = env_data, calibrate = "gif")
    
    # Ensure p-values are properly transformed for plotting
    calibrated_pvalues <- data.frame(test_result$calibrated.pvalue)
    
    # Ensure p-values are properly transformed for plotting
    neg_log10_pvalues <- data.frame(-log10(test_result$calibrated.pvalue))
    
    # Add significant column based on a p-value threshold
    threshold <- 5E-08
    significant <- ifelse(calibrated_pvalues < threshold, "Significant", "Not Significant")
    
    # Convert significant to a factor to ensure proper handling in ggplot
    significant <- factor(significant, levels = c("Not Significant", "Significant"))
    
    # Create index for SNPs
    index <- 1:length(test_result$calibrated.pvalue)
    
    # Prepare data frame for plotting
    results_df <- data.frame(index = index, calibrated_pvalues = calibrated_pvalues, neg_log10_pvalues = neg_log10_pvalues, significant = significant, variable = var)
    colnames(results_df) <- c("index", "calibrated_pvalues", "neg_log10_pvalues", "significant", "variable")
    results_df
    # Store results in list
    thyoides_results_list[[var]] <- results_df

  }

  # Combine all results into one data frame
  thyoides_results_df <- bind_rows(thyoides_results_list)
  thyoides_results_df

  # Define color palette for significant points
  significant_colors <- c(
    "thyoides_Annual_Mean_Temperature.env" = "red",
    "thyoides_Mean_Diurnal_Range.env" = "blue",
    "thyoides_Isothermality.env" = "green",
    "thyoides_Temperature_Seasonality.env" = "purple",
    "thyoides_Max_Temperature_of_Warmest_Month.env" = "orange",
    "thyoides_Min_Temperature_of_Coldest_Month.env" = "brown",
    "thyoides_Temperature_Annual_Range.env" = "pink",
    "thyoides_Mean_Temperature_of_Wettest_Quarter.env" = "yellow",
    "thyoides_Mean_Temperature_of_Driest_Quarter.env" = "cyan",
    "thyoides_Mean_Temperature_of_Warmest_Quarter.env" = "magenta",
    "thyoides_Mean_Temperature_of_Coldest_Quarter.env" = "darkgreen",
    "thyoides_Annual_Precipitation.env" = "darkblue",
    "thyoides_Precipitation_of_Wettest_Month.env" = "darkred",
    "thyoides_Precipitation_of_Driest_Month.env" = "darkorange",
    "thyoides_Precipitation_Seasonality.env" = "hotpink",
    "thyoides_Precipitation_of_Wettest_Quarter.env" = "darkcyan",
    "thyoides_Precipitation_of_Driest_Quarter.env" = "darkmagenta",
    "thyoides_Precipitation_of_Warmest_Quarter.env" = "lightgreen",
    "thyoides_Precipitation_of_Coldest_Quarter.env" = "lightblue"
  )

  # Create names for significant and non-significant combined palette
  significant_names <- paste0(names(significant_colors), ".Significant")
  combined_palette <- c(setNames(significant_colors, significant_names), "Not Significant" = "gray")
  
  thyoides_results_df_noNA <- na.omit(thyoides_results_df)
  glimpse(thyoides_results_df_noNA)

  # Plot the results with ggplot2
  lfmm_faceted <- ggplot(thyoides_results_df_noNA, aes(x = index, y = neg_log10_pvalues)) +
    geom_point(aes(color = ifelse(significant == "Significant", paste0(variable, ".Significant"), "Not Significant")), size = 1.5) +
    facet_wrap(~ variable, scales = "free_y") +
    scale_color_manual(values = combined_palette) +
    labs(title = "Significant SNPs by environment", x = "SNP Index", y = "-log10(p-value)") +
    theme_minimal()

  # Extract the base name of the LFMM file (remove directory and extension)
  #lfmm_base_name <- tools::file_path_sans_ext(basename("geno"))

  # Generate the output file name
  output_file <- paste0(outfile, ".png")

  # Save the plot
  ggsave(output_file, plot = lfmm_faceted, width = 32, height = 24, dpi = 300)

}

```




</p>
</details> 

***

