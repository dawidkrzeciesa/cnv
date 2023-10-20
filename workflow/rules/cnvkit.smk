rule cnvkit_access:
    input:
        fasta=get_reference, 
    output:
        "results/access-mappable.bed",
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/acces/access.log",
    params:
        access_param=config["cnvkit"]["access_param"],
    shell:
        "(cnvkit.py access {input} {params.access_param} -o {output}) 2>{log}"


rule cnvkit_batch:
    input:
        tumor=get_cnvkit_batch_input,
        normal=lambda w: get_cnvkit_batch_input(w, sample_type="normal"),
        bai_T=lambda w: get_cnvkit_batch_input(w, ext="bai"),
        bai_N=lambda w: get_cnvkit_batch_input(w, sample_type="normal", ext="bai"),
        fasta=get_reference,
        targets=get_targets,
        access=rules.cnvkit_access.output,
    output:
        cns="results/cnvkit_batch/{sample}.{group}.cns",
        cnr="results/cnvkit_batch/{sample}.{group}.cnr",
        cnn="results/cnvkit_batch/{sample}.{group}.cnn",
    params:
        folder=lambda wc, output: path.dirname(output.cns),
        basename=lambda wc, input: path.splitext(path.basename(input.tumor))[0],
        normal=lambda wc, input: input.normal if (input.normal != input.tumor) else " ",
        batch=config["cnvkit"]["batch"],
        chr_sex=get_chr_sex,
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit_batch/{sample}.{group}.log",
    threads: 64
    shell:
        "(cnvkit.py batch {input.tumor} "
        "  --normal {params.normal} "
        "  --targets {input.targets} "
        "  --fasta {input.fasta} "
        "  --output-reference {output.cnn} "
        "  --access {input.access} "
        "  --output-dir {params.folder} "
        "  --diagram "
        "  --scatter "
        "  -p {threads} "
        "  {params.chr_sex} "
        "  {params.batch}; "
        " mv {params.folder}/{params.basename}.cns {output.cns}; "
        " mv {params.folder}/{params.basename}.cnr {output.cnr}; "
        ") 2>{log}"


rule present_bcf_to_vcf:
    input:
        get_varlociraptor_present_bcf,
    output:
        "results/theta2/{sample}.{group}.vcf",
    log:
        "logs/theta2/{sample}.{group}.vcf.log",
    params:
        extra="--types snps",
    wrapper:
        "v2.6.0/bio/bcftools/view"


rule cnvkit_to_theta2:
    input:
        cns="results/cnvkit_batch/{sample}.{group}.cns",
        cnn="results/cnvkit_batch/{sample}.{group}.cnn",
        vcf="results/theta2/{sample}.{group}.vcf",
    output:
        interval_count="results/cnvkit_batch/{sample}.{group}.interval_count",
        tumor_snps="results/cnvkit_batch/{sample}.{group}.tumor.snp_formatted.txt",
        normal_snps="results/cnvkit_batch/{sample}.{group}.normal.snp_formatted.txt",
    log:
        "logs/theta2/{sample}.{group}.input.log",
    conda:
        "../envs/cnvkit.yaml"
    params:
        tumor_alias=lambda wc: samples.loc[samples["sample_name"] == wc.sample, "alias"].squeeze(),
        normal_alias=lambda wc: get_normal_alias_of_group(wc.group),
    shell:
        "(cnvkit.py export theta "
        "  --reference {input.cnn} "
        "  --vcf {input.vcf} "
        "  --sample-id {params.tumor_alias} "
        "  --normal-id {params.normal_alias} "
        "  {input.cns} "
        ") 2>{log}"


rule theta2_purity_estimation:
    input:
        interval_count="results/cnvkit_batch/{sample}.{group}.interval_count",
        tumor_snps="results/cnvkit_batch/{sample}.{group}.tumor.snp_formatted.txt",
        normal_snps="results/cnvkit_batch/{sample}.{group}.normal.snp_formatted.txt",
    output:
        res="results/theta2/{sample}.{group}.BEST.results",
    log:
        "logs/theta2/{sample}.{group}.BEST.log",
    conda:
        "../envs/theta2.yaml"
    params:
        out_dir=lambda wc, output: path.dirname(output.res),
    threads: 8
    shell:
        "(RunTHetA.py {input.interval_count} "
        "  --TUMOR_FILE {input.tumor_snps} "
        "  --NORMAL_FILE {input.normal_snps} "
        "  --DIR {params.out_dir} "
        "  --BAF "
        "  --NUM_PROCESSES {threads} "
        "  --FORCE "
        ") 2>{log}"


rule theta2_to_cnvkit:
    input:
        res="results/theta2/{sample}.{group}.BEST.results",
    output:
        cns="results/cnvkit_batch/{sample}.{group}.purity_adjusted.cns",
    log:
        "logs/cnvkit_batch/{sample}.{group}.purity_adjusted.log",
    conda:
        "../envs/cnvkit.yaml"
    params:
        out_dir=lambda wc, output: path.dirname(output.cns),
    shell:
        "(cnvkit.py import-theta "
        "  --output-dir {params.out_dir} "
        "  {output.cns} "
        "  {input.res} "
        ") 2>{log}"


rule cnvkit_call:
    input:
        cns=get_cnvkit_call_input,
    output:
        cns="results/cnvkit_call/{sample}.{group}.cns",
    params:
        call_param=config["cnvkit"]["call_param"],
        tumor_purity=get_tumor_purity_setting,
        chr_sex=get_chr_sex,
        sample_sex=get_sample_sex,
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/call/{sample}.{group}.cns.log",
    threads: 1
    shell:
        "(cnvkit.py call "
        "  {params.chr_sex} "
        "  -x {params.sample_sex} "
        "  --drop-low-coverage "
        "  -m clonal "
        "  {params.tumor_purity} "
        "  {input.cns} "
        "  -o {output.cns} "
        ") 2>{log}"


rule cnvkit_export_seg:
    input:
        cns=expand("results/cnvkit_call/{sample_group}.cns",
            sample_group=get_tumor_sample_group_pairs()
        ),
    output:
        "results/cnvkit/segments.seg",
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/segments.log",
    threads: 1
    shell:
        "(cnvkit.py export seg {input.cns} -o {output}) 2>{log}"


rule gistic_get_input:
    input:
        "results/cnvkit/segments.seg",
    output:
        "results/cnvkit/gistic.input.seg",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/cnvkit/pandas.log",
    threads: 1
    script:
        "../scripts/gistic2_input.py"
