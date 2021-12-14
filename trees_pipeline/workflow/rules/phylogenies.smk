
rule multiple_sequence_alignment:
    input:
        get_input_faa
    output:
        "results/1_msas/{cl}.mafft"
    threads: 25
    conda:
        "../envs/phylogenies.yaml"
    log:
        "logs/1_msas/{cl}.log"
    shell:
        '''
        nseqs=$(grep -c ">" {input})
        if [[ $nseqs -lt 5000 ]]
        then
            if [[ $nseqs -gt 1000 ]]
            then
                mafft-fftnsi --maxiterate 200 --thread {threads} {input} > {output} 2> {log}
            else
                mafft-einsi --thread {threads} {input} > {output} 2> {log}
            fi
        else
            touch {output}
        fi
        '''

rule trimming:
    input:
        rules.multiple_sequence_alignment.output
    output:
        "results/1_msas/{cl}_trimmed.mafft"
    threads: 1
    conda:
        "../envs/phylogenies.yaml"
    shell:
        '''
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            touch {output}
        else
            trimal -in {input} -out {output} -gt 0.5
        fi
        '''

rule construct_phylogeny:
    input:
        rules.trimming.output
    output:
        "results/2_trees/{cl}_trimmed.nw"
    threads: 5
    log:
        "logs/2_trees/{cl}_trimmed.log"
    conda:
        "../envs/phylogenies.yaml"
    shell:
        '''
        size=$(stat --printf="%s" {input})
        if [[ $size -eq 0 ]]
        then
            touch {output}
        else
            FastTreeMP {input} > {output} 2> {log}
        fi
        '''

rule get_crassvirales_mrcas_sister:
    input:
        rules.construct_phylogeny.output
    output:
        "results/3_bacterial_sister_clades/{cl}_mrcas_sister.txt"
    conda:
        "../envs/parse_trees.yaml"
    script:
        "../scripts/get_sister_bacteria.py"
