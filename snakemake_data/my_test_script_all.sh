
#!/bin/zsh
#based on https://carpentries-incubator.github.io/snakemake-novice-bioinformatics
########################################################################################################################
########################################################################################################################
############ Episode 1 Running commands with Snakemake
############
echo "Episode 1 "

if [ -d yeast ]; then #better to run under yeast directory so go to yeast directory...
    cd yeast
    ls reads
fi
if [ -d original_reads ]; then #revert to original reads
    rm -rf reads
    mv original_reads reads
fi


#some basic shell commands
head -n8 reads/ref1_1.fq
wc -l reads/ref1_1.fq
wc -l reads/ref1_1.fq > ref1_1.fq.count
head -v *.count

#The sample dataset represents a transcriptomics experiment in brewer’s yeast (Saccharomyces cerevisiae) under three conditions.
#Three replicates.


cat << 'EOF' > Snakefile
rule countlines:
    output: "ref1_1.fq.count"
    input:  "reads/ref1_1.fq"
    shell:
        "wc -l reads/ref1_1.fq > ref1_1.fq.count"
EOF

snakemake -j1 -F -p ref1_1.fq.count

#quotation, indendation, and colons are important in Snakefile
#snakemake --help | less #to see all options
#count fastq reads rather than lines
echo $(( $(wc -l <reads/ref1_1.fq) / 4 ))



cat << 'EOF' > Snakefile
rule countreads:
    output: "ref1_1.fq.count"
    input:  "reads/ref1_1.fq"
    shell:
        "echo $(( $(wc -l <reads/ref1_1.fq) / 4 )) > ref1_1.fq.count"

rule countreads2:
  output: "etoh60_1_1.fq.count"
  input:  "reads/etoh60_1_1.fq"
  shell:
    "echo $(( $(wc -l <reads/etoh60_1_1.fq) / 4 )) > etoh60_1_1.fq.count"
EOF

snakemake -j1 -F -p ref1_1.fq.count etoh60_1_1.fq.count


#Before running Snakemake you need to write a Snakefile
#A Snakefile is a text file which defines a list of rules
#Rules have inputs, outputs, and shell commands to be run
#You tell Snakemake what file to make and it will run the shell command defined in the appropriate rule

########################################################################################################################
########################################################################################################################
############ Episode 2 Placeholders and wildcards
############
echo "Episode 2 "

# New generic read counter
cat << 'EOF' > Snakefile
rule countreads:
    output: "{myfile}.fq.count"
    input:  "reads/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"
EOF

#input output are placeholders
#myfile is a wildcard

#you may check https://docs.google.com/presentation/d/1GrrW0lbFre4NzRgS3UlgV_Zr3eBQIvMiQa3s-acNT6g/edit?usp=sharing

snakemake -j1 -F -p temp33_1_1.fq.count

cat temp33_1_1.fq.count


echo "this will fail due to missing input file"
snakemake -j1 -F -p wibble_1.fq.count

#dry run, dry run, dry run
snakemake -n -F -p temp33_1_1.fq.count

#snakemake -f --dag temp33_1_1.fq.count | dot -Tpng > dag.png 
#snakemake --rulegraph -j1 -p temp33_1_1.fq.count | dot -Tpng > dag.png


#lets another rule to trim reads BUT they are not connected, yet...
cat << 'EOF' > Snakefile
rule countreads:
    output: "{myfile}.fq.count"
    input:  "reads/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# Trim any FASTQ reads for base quality
rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
EOF

snakemake -j1 -F -p trimmed/ref1_1.fq trimmed/ref1_2.fq



########################################################################################################################
########################################################################################################################
############ Episode 3 Chaining rules
############
echo "Episode 3 "

cat << 'EOF' > Snakefile
# New even-more-generic read counter
rule countreads:
    output: "{indir}.{myfile}.fq.count"
    input:  "{indir}/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
EOF

snakemake -j1 -F -p trimmed.ref1_1.fq.count
cat trimmed.ref1_1.fq.count
#you may run for reads.ref1_1.fq.count


# now we can add another rule to run kallisto and chain them
# kallisto requires trimmed reads
cat << 'EOF' > Snakefile
rule countreads:
    output: "{indir}.{myfile}.fq.count"
    input:  "{indir}/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"

rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "kallisto quant -i {input.index} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"

EOF

echo "this will fail due to missing index"
snakemake -j1 -F -p kallisto.ref1/abundance.h5


#now we should add kallisto index rule and chain it
cat << 'EOF' > Snakefile
rule countreads:
    output: "{indir}.{myfile}.fq.count"
    input:  "{indir}/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"

rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "kallisto quant -i {input.index} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"

rule kallisto_index:
    output:
        idx = "{strain}.kallisto_index",
        log = "{strain}.kallisto_log",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    shell:
        "kallisto index -i {output.idx} {input.fasta} >& {output.log}"
EOF


#multiple samples to process, it is possible to run parallel
snakemake -j1 -F -p kallisto.ref1/abundance.h5  kallisto.temp33_1/abundance.h5

#visualize the DAG
#snakemake -f --dag kallisto.ref1/abundance.h5 | dot -Tpng > dag.png 
#snakemake --rulegraph -j1 -p kallisto.ref1/abundance.h5 | dot -Tpng > dag.png


#beware of typos
#snakemake -j1 -F -p kallisto.ref1/abundances.h5

########################################################################################################################
########################################################################################################################
############ Episode 4 How Snakemake plans what jobs to run
############
echo "Episode 4 "

snakemake -j1 -F -p kallisto.temp33_1/abundance.h5 #just in case create again
snakemake -j1 -p kallisto.temp33_1/abundance.h5 #check what happens

rm kallisto.temp33_1/* #remove all created files in kallisto directory
snakemake -j1 -p kallisto.temp33_1/abundance.h5 #now it should only run kallisto for temp33_1

rm trimmed/temp33_*.fq #remove all trimmed files
snakemake -j1 -p kallisto.temp33_1/abundance.h5 #now what happens???

touch transcriptome/*.fa.gz #touch fasta file under transcriptome
snakemake -j1 -p kallisto.temp33_1/abundance.h5 #now what happens???



#The -R flag allows you to explicitly tell Snakemake that a rule has changed and that all outputs from that rule need to be re-evaluated.

snakemake -j1 -R trimreads -p kallisto.temp33_1/abundance.h5 #forcing to run trimreads rule


#forcing multiple rules
snakemake -j1 -R trimreads kallisto_index -p kallisto.temp33_1/abundance.h5 kallisto.temp33_2/abundance.h5


snakemake -j1 -f -p kallisto.temp33_1/abundance.h5 #forcing only the target

snakemake -f --dag kallisto.etoh60_1/abundance.h5 | dot -Tpng > dag.png #create a dag plot for etoh60_1



########################################################################################################################
########################################################################################################################
############ Episode 5 Processing lists of inputs
############
echo "Episode 5 "

#we have to do some processing to the reads due to inconsistency in the naming
#Snakemake prefers to work with consistent file names and extensions
if [ ! -d original_reads ]; then
    mv reads original_reads
    mkdir reads
    cd reads
    ln -s ../original_reads/* .
    rename -v -s ref ref_ * #rename command is required, available in Bioconda
    cd ..
fi

# Input conditions and replicates to process
# CONDITIONS = ["ref", "etoh60", "temp33"]
# REPLICATES = ["1", "2", "3"]




cat << 'EOF' > Snakefile
# Input conditions and replicates to process
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]

rule all_counts:
    input: expand("trimmed.{cond}_{rep}_1.fq.count", cond=CONDITIONS, rep=REPLICATES)

rule countreads:
    output: "{indir}.{myfile}.fq.count"
    input:  "{indir}/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"

rule kallisto_quant:
    output:
        h5   = "kallisto.{sample}/abundance.h5",
        tsv  = "kallisto.{sample}/abundance.tsv",
        json = "kallisto.{sample}/run_info.json",
    input:
        index = "Saccharomyces_cerevisiae.R64-1-1.kallisto_index",
        fq1   = "trimmed/{sample}_1.fq",
        fq2   = "trimmed/{sample}_2.fq",
    shell:
        "kallisto quant -i {input.index} -o kallisto.{wildcards.sample} {input.fq1} {input.fq2}"

rule kallisto_index:
    output:
        idx = "{strain}.kallisto_index",
        log = "{strain}.kallisto_log",
    input:
        fasta = "transcriptome/{strain}.cdna.all.fa.gz"
    shell:
        "kallisto index -i {output.idx} {input.fasta} >& {output.log}"
EOF



# snakemake -j1 -p all_counts
# snakemake -j1 -p -n #dry run to check
# snakemake  -j1 --dag | dot -Tpng > dag.png  #generate dag plot and check which files are generated
# are the any missing??

snakemake -j1 -p # without telling the target, it will run the first rule in the Snakefile


# Input conditions and replicates to process
# CONDITIONS = ["ref", "etoh60", "temp33"]
# REPLICATES = ["1", "2", "3"]
# READ_ENDS  = ["1", "2"]
# COUNT_DIR  = ["reads", "trimmed"]

# rule all_counts:
#     input: expand("{indir}.{cond}_{rep}_{end}.fq.count", indir=COUNT_DIR, cond=CONDITIONS, rep=REPLICATES, end=READ_ENDS)


cat << 'EOF' > Snakefile
CONDITIONS = ["ref", "etoh60", "temp33"]
REPLICATES = ["1", "2", "3"]
READ_ENDS  = ["1", "2"]
COUNT_DIR  = ["reads", "trimmed"]

rule all_counts:
    input:
        untrimmed = expand( "reads.{cond}_{rep}_{end}.fq.count",   cond  = CONDITIONS,
                                                                   rep   = REPLICATES,
                                                                   end   = READ_ENDS ),
        trimmed   = expand( "trimmed.{cond}_{rep}_{end}.fq.count", cond  = CONDITIONS,
                                                                   rep   = REPLICATES,
                                                                   end   = READ_ENDS ),
    output:
        untrimmed = "untrimmed_counts_concatenated.txt",
        trimmed   = "trimmed_counts_concatenated.txt",
    shell:
        "cat {input.untrimmed} > {output.untrimmed} ; cat {input.trimmed} > {output.trimmed}"

rule countreads:
    output: "{indir}.{myfile}.fq.count"
    input:  "{indir}/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
EOF

######################################################################
#Dynamically determine the input files for the all_counts rule
#  CONDITIONS = glob_wildcards("reads/{condition}_1_1.fq").condition
#  print("Conditions are: ", CONDITIONS)
######################################################################


########################################################################################################################
########################################################################################################################
############ Episode 6 Handling awkward programs
############
echo "Episode 6 "



fastqc reads/ref_1_1.fq
# placeholders vs wildcards

#see sample fastqc.Snakefile: ../../sample_snakefiles/fastqc.Snakefile

########################################################################################################################
########################################################################################################################
############ Episode 7 Finishing the basic workflow
############
echo "Episode 7 "


#see sample ep07.Snakefile: ../../reference_snakefiles/ep07.Snakefile

cat ../../reference_snakefiles/ep07.Snakefile > Snakefile

#to generate the dag plot
snakemake --rulegraph -j1 -p multiqc | dot -Tpng > dag.png

#run the pipeline
snakemake -j5 -p multiqc

#multiqc rule can be moved to up trigger auto run
#it is possible to use CONDITIONS = glob_wildcards("reads/{condition}_1_1.fq").condition

#test with different number of cores
#time snakemake -j5 -p multiqc
#snakemake -j5 multiqc  
#snakemake -j1 multiqc





#add rule all counts for trimmed back
#rule all_counts:
#    input:
#        expand("trimmed.{cond}_{rep}_{end}.fq.count", cond=CONDITIONS, rep=REPLICATES, end=["1","2"])

