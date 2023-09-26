import glob
from os import path

import pandas as pd

# read in sample sheet

samples = (
    pd.read_csv(
        config["samples"], sep="\t", dtype={"sample_name": str, "target_bed": str}
    )
    .set_index(["sample_name"], drop=False)
    .sort_index()
)

### set global variables

TUMOR_SAMPLES = set(
    samples.loc[
        samples["alias"].str.startswith(config["alias_prefixes"]["tumor"]), "sample_name"
    ]
)

### wildcard constraints

wildcard_constraints:
    sample="|".join(samples["sample_name"].drop_duplicates()),
    group="|".join(samples["group"].drop_duplicates()),

### helper functions


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


def get_group_sample_type(group, sample_type):
    return samples.loc[
        (samples["group"] == group)
        & (samples["alias"].str.startswith(config["alias_prefixes"][sample_type])),
        "sample_name",
    ]


def get_normal_alias_of_group(group):
    return samples.loc[
        (samples["group"] == group)
        & (samples["alias"].str.startswith(config["alias_prefixes"]["normal"])),
        "alias",
    ]


### input functions


def get_cnvkit_batch_input(wildcards, sample_type="tumor", ext="bam"):
    group = get_group_of_sample(wildcards.sample)
    sample_name = get_group_sample_type(group, sample_type)
    if (len(sample_name) == 0) and (sample_type == "normal"):
        # fall back to a `panel_of_normals` alias, if no matched normal sample is present
        sample_name = get_group_sample_type(group, "panel_of_normals")
        if len(sample_name) == 0:
            # if not even a `panel_of_normal` sample is available, stick the tumor sample
            # here to be able to automatically trigger the rule with an empty `--normal`
            # specification via the params.normal lambda function
            sample_name = get_group_sample_type(group, "tumor")
    if len(sample_name) == 0:
        raise AssertionError(
            f"Group {group} does not have a sample with whose alias has the tumor prefix '{config['alias_prefixes']['tumor']}'"
        )
    return f"results/recal/{sample_name}.{ext}"


def get_cnvkit_call_input(wildcards):
    # no purity specified for this sample
    if len(get_tumor_purity_setting(wildcards.sample)) == 0:
        return f"results/cnvkit_batch/{wildcards.sample}.purity_adjusted.cns"
    else:

def get_reference(wildcards):
    return config["ref"]




def get_varlociraptor_present_bcf(wildcards):
    group = get_group_of_sample(wildcards.sample)
    event = config["cnvkit"]["joint_event"]
    return f"results/final-calls/{group}.{event}.fdr-controlled.bcf"
