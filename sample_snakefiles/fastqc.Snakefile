rule fastqc:
    output:
        html = "{indir}/{myfile}_fastqc.html",
        zip  = "{indir}/{myfile}_fastqc.zip",
    input:  "{indir}/{myfile}.fq"
    shell:
        "fastqc {input}"

rule fastqc:
    output:
        html = "{indir}.{myfile}_fastqc.html",
        zip  = "{indir}.{myfile}_fastqc.zip"
    input:  "{indir}/{myfile}.fq"
    shell:
        """fastqc -o . {input}
           mv {wildcards.myfile}_fastqc.html {output.html}
           mv {wildcards.myfile}_fastqc.zip  {output.zip}
        """

#snakemake -s ../../sample_snakefiles/fastqc.Snakefile -j1 -p reads/etoh60_1_1_fastqc.html reads/etoh60_1_2_fastqc.html