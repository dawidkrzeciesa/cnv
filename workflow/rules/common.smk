import glob
from os import path

import pandas as pd



def get_targets(wildcards):
    targets=samples.loc[wildcards.sample, "target_bed"]
    return targets



def get_tumor_purity(wildcards):
    purity=samples.loc[wildcards.sample, "tumor_purity"]
    return purity



def get_chr_sex(wildcards):
    sex=samples.loc[wildcards.sample, "sex"]
    if sex == "male":
        return "-y"
    else:
        return "" 


def get_sample_sex(wildcards):
    sex=samples.loc[wildcards.sample, "sex"]
    if sex == "male":
        return "male"
    else:
        return "female" 