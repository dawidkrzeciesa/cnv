import pandas as pd
import glob,os


####### OncoKB #######

inport_columns=["Hugo Symbol", "Is Oncogene", "Is Tumor Suppressor Gene"]

OncoKB = pd.read_csv(snakemake.input["oncokb"], sep="\t", usecols=inport_columns)

OncoKB = OncoKB.rename(columns={"Hugo Symbol": "gene", "Is Oncogene": "oncogene", "Is Tumor Suppressor Gene": "tsg"})

OncoKB = OncoKB[OncoKB.tsg != "No"]

tsg_lst = OncoKB["gene"].tolist()

####### filter cns #######
cns = pd.read_csv(snakemake.input["cns"], sep="\t")

cns["gene"]=cns.gene.str.split(",")
cns = cns.explode("gene")

cns_tsg = cns[cns["gene"].isin(tsg_lst)]

cns_tsg.to_csv(snakemake.output["cns_tsg"], sep="\t", index=False)

####### filter cnr #######

cnr = pd.read_csv(snakemake.input["cnr"], sep="\t", low_memory=False)

cnr_tsg = cnr[cnr["gene"].isin(tsg_lst)]

cnr_tsg.to_csv(snakemake.output["cnr_tsg"], sep="\t", index=False)