configfile: "config/config.yaml"

rule download_vs2_db:
    output:
        "resources/vs2_db/done"
    params:
        db_dir = "resources/vs2_db"
    threads: 5
    conda:
        "../envs/virsorter2.yaml"
    shell:
        '''
        virsorter setup -d {params.db_dir} -j 4
        touch {output}
        '''


rule run_virsorter2:
    input:
        rules.download_vs2_db.output,
        fasta = "results/3_bacterial_contigs/{cl}.fasta"
    output:
        "results/4_prophages/{cl}/final-viral-boundary.tsv"
    params:
        outdir = "results/4_prophages/{cl}",
        db = "resources/vs2_db"
    threads: 20
    conda:
        "../envs/virsorter2.yaml"
    shell:
        '''
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            mkdir -p {params.outdir}
            touch {output}
        else
            virsorter run -w {params.outdir} -i {input.fasta} --min-length 1500 -j {threads} -d {params.db} all
        fi
        '''
