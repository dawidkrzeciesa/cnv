import glob
from os import path

import pandas as pd



def get_tumor_purity(wildcards):
    purity=samples.loc[(samples["group"] == wildcards.group) & (samples["alias"] == config["aliases"]["tumor"])]
    purity=float(purity["purity"][0])
    return purity



def get_bams_autobin(wildcards):
    bam_names=samples.loc[samples["target_name"] == wildcards.target].index.to_list()
    bam_names= [config["bam_directory"] + i + ".bam" for i in bam_names]
    return bam_names


def get_targets_autobin(wildcards):
    bed_path=samples.loc[samples["target_name"] == wildcards.target]
    bed_path=bed_path["path_to_target"].unique()
    return bed_path

def get_normal_bams(wildcards):
    bam_names_N=samples.loc[(samples["target_name"] == wildcards.target) & (samples["alias"] == config["aliases"]["normal"])].index.to_list()
    bam_names_N= [config["bam_directory"] + i + ".bam" for i in bam_names_N]
    return bam_names_N

def get_tumor_bam(wildcards):
    bam_name_T=samples.loc[(samples["alias"] == config["aliases"]["tumor"]) & (samples["group"] == wildcards.group)].index.to_list()
    bam_name_T= [config["bam_directory"] + i + ".bam" for i in bam_name_T]
    return bam_name_T


def get_ref(wildcards):
    target=samples.loc[(samples["alias"] == config["aliases"]["tumor"]) & (samples["group"] == wildcards.group)]
    target=target.iloc[0]["target_name"]
    
    if config["CNVkit"]["mode"] == "pool":
        return f"results/refs/pooled_{target}.cnn"
    
    elif config["CNVkit"]["mode"] == "tumor_only":
        return f"results/refs/flatref_{target}.cnn"
        
    
def get_targets(wildcards):
    targets=samples.loc[(samples["alias"] == config["aliases"]["normal"]) & (samples["group"] == wildcards.group)]
    targets=targets.iloc[0]["path_to_target"]
    return targets
    
    
def get_matched_normal_bam(wildcards):
    bam_name_N=samples.loc[(samples["alias"] == config["aliases"]["normal"]) & (samples["group"] == wildcards.group)].index.to_list()
    bam_name_N= [config["bam_directory"] + i + ".bam" for i in bam_name_N]
    return bam_name_N


def get_cns(wildcards):
    in_name=samples.loc[(samples["alias"] == config["aliases"]["tumor"]) & (samples["group"] == wildcards.group)].index.to_list()
    in_name=in_name[0]
    return f"results/cnvkit/{wildcards.group}/{in_name}.cns"

def get_cnr(wildcards):
    in_name=samples.loc[(samples["alias"] == config["aliases"]["tumor"]) & (samples["group"] == wildcards.group)].index.to_list()
    in_name=in_name[0]
    return f"results/cnvkit/{wildcards.group}/{in_name}.cnr"