import glob
from os import path

import pandas as pd



def get_targets(wildcards):
    return samples.loc[wildcards.sample, "target_bed"]


def get_tumor_purity_setting(wildcards):
    if "tumor_purity" in samples.columns:
        return f"--purity {samples.loc[wildcards.sample, 'tumor_purity']}"
    else:
        return ""


def get_chr_sex(wildcards):
    sex = samples.loc[wildcards.sample, "sex"]
    if sex == "male":
        return "-y"
    else:
        return ""


def get_sample_sex(wildcards):
    sex = samples.loc[wildcards.sample, "sex"]
    if sex == "male":
        return "male"
    else:
        return "female"


def get_group_of_sample(sample):
    return samples.loc[samples["sample_name"] == sample, "group"]




def get_cnvkit_call_input(wildcards):
    # no purity specified for this sample
    if len(get_tumor_purity_setting(wildcards.sample)) == 0:
        return f"results/cnvkit_batch/{wildcards.sample}.purity_adjusted.cns"
    else:
        return f"results/cnvkit_batch/{wildcards.sample}.cns"


def get_varlociraptor_present_bcf(wildcards):
    group = get_group_of_sample(wildcards.sample)
    event = config["cnvkit"]["joint_event"]
    return f"results/final-calls/{group}.{event}.fdr-controlled.bcf"
