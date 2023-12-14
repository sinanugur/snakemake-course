
#!/bin/zsh

########################################################################################################################
########################################################################################################################
############ Episode 1
############
echo "Episode 1 "

if [ -d yeast ]; then #better to run under yeast directory
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

cat << 'EOF' > Snakefile
rule countlines:
    output: "ref1_1.fq.count"
    input:  "reads/ref1_1.fq"
    shell:
        "wc -l reads/ref1_1.fq > ref1_1.fq.count"
EOF

snakemake -j1 -F -p ref1_1.fq.count

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


########################################################################################################################
########################################################################################################################
############ Episode 2
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

snakemake -j1 -F -p temp33_1_1.fq.count

cat temp33_1_1.fq.count


echo "this will fail due to missing input file"
snakemake -j1 -F -p wibble_1.fq.count

#dry run
snakemake -n -F -p temp33_1_1.fq.count

cat << 'EOF' > Snakefile
# Trim any FASTQ reads for base quality
rule countreads:
    output: "{myfile}.fq.count"
    input:  "reads/{myfile}.fq"
    shell:
        "echo $(( $(wc -l <{input}) / 4 )) > {output}"

rule trimreads:
    output: "trimmed/{myfile}.fq"
    input:  "reads/{myfile}.fq"
    shell:
        "fastq_quality_trimmer -t 20 -l 100 -o {output} <{input}"
EOF

snakemake -j1 -F -p trimmed/ref1_1.fq trimmed/ref1_2.fq


########################################################################################################################
########################################################################################################################
############ Episode 3
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


echo "this will fail due to typo" #abundances.h5
snakemake -j1 -F -p kallisto.ref1/abundances.h5



#multiple samples to process
snakemake -j1 -F -p kallisto.ref1/abundance.h5  kallisto.temp33_1/abundance.h5


########################################################################################################################
########################################################################################################################
############ Episode 4
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
############ Episode 5
############
echo "Episode 5 "


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
############ Episode 6
############
echo "Episode 6 "



fastqc reads/ref_1_1.fq


