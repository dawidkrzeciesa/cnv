rule filter_oncogene:
    input:
        cns="results/cnvkit_call/{sample}.cns",
    output:
        cns_oncogene="results/cnvkit_call/filtered/oncogene/{sample}_oncogene_only.cns",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/filter_oncogene/{sample}.log"
    threads: 1
    script:
        "../scripts/filter_oncogene.py"




rule filter_tumor_suppressor:
    input:
        cns="results/cnvkit_call/{sample}.cns",
    output:
        cns_tsg="results/cnvkit_call/filtered/tumor_supressor/{sample}_tumor_supressor_only.cns",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/filter_tumor_supressor/{sample}.log"
    threads: 1
    script:
        "../scripts/filter_tumor_supressor.py"



rule filter_vgp:
    input:
        cns="results/cnvkit_call/{sample}.cns",
    output:
        cns_vgp="results/cnvkit_call/filtered/vgp/{sample}_vgp.cns",
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/filter_vgp/{sample}.log"
    threads: 1
    script:
        "../scripts/filter_vgp.py"



#ruleorder: cnvkit_call_cns > filter_oncogene 



rule build_matrix_tumor_suppressor:
    input:
        expand("results/cnvkit/cnv_call/filtered/tumor_supressor/{sample}_tumor_supressor_only.cns", sample=SAMPLES)
    output:
        matrix_tsg = "results/oncoprint/matrix/tumor_supressor_matrix.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_matrix_oncoprint_tumor_suppressor.py"




rule oncoprint_tumor_suppressors:
    input:
        matrix = "results/oncoprint/matrix/tumor_supressor_matrix.tsv",
    output:
        report("results/oncoprint/tumor_supressor_deletions_only.pdf", category="oncoprints", caption="../report/oncoprint_fusions.rst")
    conda:
        "../envs/oncoprint.yaml"
    params:
        title = config["oncoprint"]["title"]
    script:
        "../scripts/oncoprint_tumor_suppressor.R"



rule build_matrix_oncogene:
    input:
        expand("results/cnvkit/cnv_call/filtered/oncogene/{sample}_oncogene_only.cns", sample=SAMPLES)
    output:
        matrix_tsg = "results/oncoprint/matrix/oncogene_matrix.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_matrix_oncoprint_oncogene.py"




rule oncoprint_oncogene:
    input:
        matrix = "results/oncoprint/matrix/oncogene_matrix.tsv",
    output:
        report("results/oncoprint/oncogene_amplifications_only.pdf", category="oncoprints", caption="../report/oncoprint_fusions.rst")
    conda:
        "../envs/oncoprint.yaml"
    params:
        title = config["oncoprint"]["title"]
    script:
        "../scripts/oncoprint_oncogene.R"



rule build_matrix_vgp:
    input:
        expand("results/cnvkit/cnv_call/filtered/vgp/{sample}_vgp.cns", sample=SAMPLES)
    output:
        matrix_vgp = "results/oncoprint/matrix/vgp_matrix.tsv"
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/build_matrix_oncoprint_vgp.py"




rule oncoprint_vgp:
    input:
        matrix = "results/oncoprint/matrix/vgp_matrix.tsv",
    output:
        report("results/oncoprint/oncogene_vgp.pdf", category="oncoprints", caption="../report/oncoprint_fusions.rst")
    conda:
        "../envs/oncoprint.yaml"
    params:
        title = config["oncoprint"]["title"]
    script:
        "../scripts/oncoprint_vgp.R"