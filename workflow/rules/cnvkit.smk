rule cnvkit_access:
    input:
        config["ref_fasta"],
    output:
        "results/access-mappable.bed",
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/acces/access.log"
    params:
        access_param=config["params"]["cnvkit"]["access_param"]
    shell:
        "(cnvkit.py access {input} {params.access_param} -o {output}) 2>{log}"




rule cnvkit_batch:
    input:
        tumor=config["bam_directory"]+"{sample}_T.sorted.bam",
        normal=config["bam_directory"]+"{sample}_N.sorted.bam",
        bai_T=config["bam_directory"]+"{sample}_T.sorted.bai",
        bai_N=config["bam_directory"]+"{sample}_N.sorted.bai",
        fasta=config["ref_fasta"],
        targets=get_targets,
        access=rules.cnvkit_access.output,
    output:
        folder=directory("results/cnvkit/{sample}/batch/"),
        cns="results/cnvkit/{sample}/batch/{sample}_T.sorted.cns",
        cnr="results/cnvkit/{sample}/batch/{sample}_T.sorted.cnr"
    params:
        reference_out="results/cnvkit/{sample}/batch/{sample}_reference.cnn",
        batch=config["params"]["cnvkit"]["batch"],
        chr_sex=get_chr_sex
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/batch/{sample}.log"
    threads: workflow.cores
    shell:
        "(cnvkit.py batch {input.tumor} --normal {input.normal} --targets {input.targets} --fasta {input.fasta} "
        "--output-reference {params.reference_out} --access {input.access} --output-dir {output.folder} "
        "--diagram --scatter -p {threads} {params.chr_sex} {params.batch}) 2>{log}"


ruleorder: cnvkit_call_cns > filter_tumor_suppressor > filter_oncogene > filter_vgp


rule cnvkit_call_cns:
    input:
        "results/cnvkit/{sample}/batch/{sample}_T.sorted.cns"
    output:
        "results/cnvkit/cnv_call/{sample}.cns"
    params:
        call_param=config["params"]["cnvkit"]["call_param"],
        tumor_purity=get_tumor_purity,
        chr_sex=get_chr_sex,
        sample_sex=get_sample_sex
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/call/{sample}.cns.log"
    threads: 1
    shell:
        "(cnvkit.py call {params.chr_sex} -x {params.sample_sex} --drop-low-coverage -m clonal --purity {params.tumor_purity} {input} -o {output}) 2>{log}"




rule cnvkit_call_cnr:
    input:
        "results/cnvkit/{sample}/batch/{sample}_T.sorted.cnr"
    output:
        "results/cnvkit/cnv_call/{sample}.cnr"
    params:
        call_param=config["params"]["cnvkit"]["call_param"],
        tumor_purity=get_tumor_purity,
        chr_sex=get_chr_sex,
        sample_sex=get_sample_sex
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/call/{sample}.cnr.log"
    threads: 1
    shell:
        "(cnvkit.py call {params.chr_sex} -x {params.sample_sex} --drop-low-coverage -m clonal --purity {params.tumor_purity} {input} -o {output}) 2>{log}"



rule cnvkit_export_seg:
    input:
        cns=expand("results/cnvkit_call/{sample}.cns", sample=TUMOR_SAMPLES),
    output:
        "results/cnvkit/segments.seg"
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/segments.log"
    threads: 1
    shell:
        "(cnvkit.py export seg {input.cns} -o {output}) 2>{log}"




rule gistic_get_intput:
    input:
        "results/cnvkit/segments.seg"
    output:
        "results/cnvkit/gistic.input.seg"
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/cnvkit/pandas.log"
    threads: 1
    script:
        "../scripts/gistic2_input.py"