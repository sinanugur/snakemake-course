
#!/bin/zsh

############
############
############ Episode 1
############
echo "Episode 1 "
cd yeast
ls reads

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


############
############
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


echo "this will fail"
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


############
############
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
EOF

snakemake -j1 -F -p trimmed.ref1_1.fq.count
cat trimmed.ref1_1.fq.count

cat << 'EOF' > Snakefile
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


snakemake -j1 -F -p kallisto.ref1/abundance.h5


cat << 'EOF' > Snakefile
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

snakemake -j1 -F -p kallisto.ref1/abundance.h5