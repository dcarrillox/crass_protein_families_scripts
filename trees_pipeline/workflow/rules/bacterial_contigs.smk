
rule download_bacterial_contigs:
    input:
        get_input_faa
    output:
        ipg = "results/3_bacterial_contigs/{cl}.ipg",
        fasta = "results/3_bacterial_contigs/{cl}.fasta"
    conda:
        "../envs/ncbi_download_biopython.yaml"
    script:
        "../scripts/download_bacterial_contigs.py"
