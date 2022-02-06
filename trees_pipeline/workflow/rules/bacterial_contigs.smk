configfile: "config/config.yaml"


rule download_bacterial_contigs:
    # input:
    #     "results/3_bacterial/{cl}_mrcas_sister.txt"
    output:
        ipg = "results/3_bacterial/{cl}.ipg",
        fasta = "results/3_bacterial/{cl}.fasta",
        to_download = "results/3_bacterial/{cl}.todownl"
    threads: 10
    conda:
        "../envs/ncbi_download_biopython.yaml"
    script:
        "../scripts/download_bacterial_contigs.py"

rule filter_bacterial_contigs:
    input:
        ipg = "results/3_bacterial/{cl}.ipg",
        fasta = "results/3_bacterial/{cl}.fasta"
    output:
        "results/4_prophages/{cl}_bacterial.fasta"
    conda:
        "../envs/ncbi_download_biopython.yaml"
    script:
        "../scripts/filter_bacterial_contigs.py"

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
        fasta = rules.filter_bacterial_contigs.output
    output:
        "results/4_prophages/{cl}/final-viral-boundary.tsv"
    params:
        outdir = "results/4_prophages/{cl}",
        db = "resources/vs2_db"
    threads: 23
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

rule parse_virsorter2_results:
    input:
        vs2 = rules.run_virsorter2.output,
        all_contigs_fasta = rules.download_bacterial_contigs.output.fasta,
        filt_contigs_fasta = rules.filter_bacterial_contigs.output,
        ipg = rules.download_bacterial_contigs.output.ipg
    output:
        "results/4_prophages/{cl}_bacterial.table"
    params:
        score_cutoff = 0.5
    conda:
        "../envs/ncbi_download_biopython.yaml"
    script:
        "../scripts/parse_vs2_output.py"
