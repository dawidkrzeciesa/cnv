rule cnvkit_access:
    input:
        config["ref_fasta"],
    output:
        "results/refs/access-mappable.bed",
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/acces/access.log"
    params:
        access_param=config["params"]["cnvkit"]["access_param"]
    shell:
        "(cnvkit.py access {input} {params.access_param} -o {output}) 2>{log}"


rule cnvkit_autobin:
    input:
        ref=config["ref_fasta"],
        bams=get_bams_autobin,
        targets=get_targets_autobin,
        access=rules.cnvkit_access.output,
    output:
        antitarget_bed=temp("antitarget.{target}.bed"),
        target_bed=temp("target.{target}.bed"),
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/autobin/pooled.{target}.log"
    params:
        access_param=config["params"]["cnvkit"]["access_param"],
        antitarget_bed_name="antitarget.{target}.bed",
        target_bed_name="target.{target}.bed",
    shell:
        "(cnvkit.py autobin {input.bams} -t {input.targets} -g {input.access} "
            "--antitarget-output-bed {params.antitarget_bed_name} --target-output-bed {params.target_bed_name}) 2>{log}"


if config["CNVkit"]["mode"] == "pool":
    print("mode: pool normals")


    rule make_cnvkit_pool_ref:
        input:
            ref=config["ref_fasta"],
            bams=get_normal_bams,
            targets="target.{target}.bed",
            antitargets="antitarget.{target}.bed",
            access=rules.cnvkit_access.output,
        output:
            "results/refs/pooled_{target}.cnn",
        conda:
            "../envs/cnvkit.yaml"
        log:
            "logs/cnvkit/refs/{target}.log"
        params:
            output_dir="results/refs/pooled_{target}/"
        threads: workflow.cores 
        shell:
            "(cnvkit.py batch -n {input.bams} --output-reference {output} -t {input.targets} -a {input.antitargets} -f {input.ref} "
            "-g {input.access} --output-dir {params.output_dir} -p {threads}) 2>{log}"


    rule run_cnvkit_pool_ref:
        input:
            ref=get_ref,
            bam=get_tumor_bam,
        output:
            cns="results/cnvkit/{group}/{sample}.cns",
            cnr="results/cnvkit/{group}/{sample}.cnr",
        conda:
            "../envs/cnvkit.yaml"
        log:
            "logs/cnvkit/batch/{group}.{sample}.ref.log"
        params:
            run_cnvkit=config["params"]["cnvkit"]["run_cnvkit"],
            out_dir="results/cnvkit/{group}/",
        threads: workflow.cores 
        shell:
            "(cnvkit.py batch {input.bam} -r {input.ref} -d {params.out_dir} -p {threads} --diagram --scatter {params.run_cnvkit}) 2> {log}"



elif config["CNVkit"]["mode"] == "tumor_only":
    print("mode: tumor only")

    rule make_cnvkit_flat_ref:
        input:
            ref=config["ref_fasta"],
            targets="target.{target}.bed",
            antitargets="antitarget.{target}.bed",
        output:
            "results/refs/flatref_{target}.cnn",
        conda:
            "../envs/cnvkit.yaml"
        log:
            "logs/cnvkit/refs/flat_ref.{target}.log"
        shell:
            "(cnvkit.py reference -o {output} -t {input.targets} -a {input.antitargets} -f {input.ref}) 2>{log}"


    rule run_cnvkit_tumor_only:
        input:
            ref=get_ref,
            bam=get_tumor_bam,
        output:
            cns="results/cnvkit/{group}/{sample}.cns",
            cnr="results/cnvkit/{group}/{sample}.cnr",
        conda:
            "../envs/cnvkit.yaml"
        log:
            "logs/cnvkit/batch/{group}.{sample}.ref.log"
        params:
            run_cnvkit=config["params"]["cnvkit"]["run_cnvkit"],
            out_dir="results/cnvkit/{group}/",
        threads: workflow.cores 
        shell:
            "(cnvkit.py batch {input.bam} -r {input.ref} -d {params.out_dir} -p {threads} --diagram --scatter {params.run_cnvkit}) 2> {log}"



elif config["CNVkit"]["mode"] == "batch":
    print("batch tumor normal")
    rule cnvkit_batch_tumor_normal:
        input:
            tumor=get_tumor_bam,
            normal=get_matched_normal_bam,
            fasta=config["ref_fasta"],
            targets=get_targets,
            access=rules.cnvkit_access.output,
        output:
            cns="results/cnvkit/{group}/{sample}.cns",
            cnr="results/cnvkit/{group}/{sample}.cnr",
        params:
            run_cnvkit=config["params"]["cnvkit"]["run_cnvkit"],
            out_dir="results/cnvkit/{group}/",
        conda:
            "../envs/cnvkit.yaml"
        log:
            "logs/cnvkit/batch/{group}.{sample}.log"
        threads: workflow.cores
        shell:
            "(cnvkit.py batch {input.tumor} --normal {input.normal} --targets {input.targets} --fasta {input.fasta} "
            "--access {input.access} --output-dir {params.out_dir} --diagram --scatter -p {threads} {params.run_cnvkit}) 2>{log}"


#ruleorder: cnvkit_call_cns > filter_tumor_suppressor > filter_oncogene > filter_vgp


rule cnvkit_call_cns:
    input:
        get_cns
    output:
        "results/cnvkit/cnv_call/{group}.cns"
    params:
        call_param=config["params"]["cnvkit"]["call_param"],
        tumor_purity=get_tumor_purity,
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/call/{group}.cns.log"
    threads: 1
    shell:
        "(cnvkit.py call --drop-low-coverage -m clonal --purity {params.tumor_purity} {input} -o {output}) 2>{log}"




rule cnvkit_call_cnr:
    input:
        get_cnr
    output:
        "results/cnvkit/cnv_call/{group}.cnr"
    params:
        call_param=config["params"]["cnvkit"]["call_param"],
        tumor_purity=get_tumor_purity,
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/call/{group}.cnr.log"
    threads: 1
    shell:
        "(cnvkit.py call --drop-low-coverage -m clonal --purity {params.tumor_purity} {input} -o {output}) 2>{log}"



rule cnvkit_export_seg:
    input:
        cns=expand("results/cnvkit/cnv_call/{group}.cns", group=GROUP),
        cnr=expand("results/cnvkit/cnv_call/{group}.cnr", group=GROUP)
    output:
        "results/cnvkit/segments.seg"
    conda:
        "../envs/cnvkit.yaml"
    log:
        "logs/cnvkit/segments.log"
    threads: 1
    shell:
        "(cnvkit.py export seg {input.cns} -o {output}) 2>{log}"


