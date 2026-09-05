#!/usr/bin/env python3
# stl
# 2022.10.11
import os
import sys
import argparse
import pandas as pd
import numpy as np

def getargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dnafile", action="store", dest="dnafile", help="")
    parser.add_argument("--outdir", action="store", dest="outdir", help="", default="")
    p = parser.parse_args()
    return p

def filter(data):
    ## filter depth < 100x & af >1% in GnomAD/1KGEAS(1,000 Genomes project phase 3 East Asians)
    #data = data[(~data["GnomAD"].isin(["-","NA"])) | (~data["1KGEAS"].isin(["-","NA"])) |(~data["dbSNP"].isin(["-","NA"])) | (~data["COSMIC"].isin(["-","NA"]))]
    data["GnomAD"] = data["GnomAD"].fillna(-1) ## NA
    data["1KGEAS"] = data["1KGEAS"].fillna(-1)
    index1 = data["GnomAD"].isin(['-'])
    data.loc[index1,"GnomAD"] = -99
    index2 = data["1KGEAS"].isin(['-'])
    data.loc[index2,"1KGEAS"] = -99
    data = data[(~data["Depth"].isin(["-"])) & (~data["Freq"].isin(["-"]))]
    if data[~data["GnomAD"].isin(["-"])].empty and data[~data["1KGEAS"].isin(["-"])].empty:  ## no GnomAD and 1KGEAS annotation
        data = data[(data["Depth"] > 100)]
    elif not data[~data["GnomAD"].isin(["-"])].empty and data[~data["1KGEAS"].isin(["-"])].empty: ## have GnomAD anno but 1KGEAS no
        data = data[~data["GnomAD"].isin(["-"])]
        data = data[(data["Depth"] > 100) & (data["GnomAD"] < 0.01)]
    elif data[~data["GnomAD"].isin(["-"])].empty and not data[~data["1KGEAS"].isin(["-"])].empty: ## have 1KGEAS anno but GnomAD no
        data = data[~data["1KGEAS"].isin(["-"])]
        data = data[(data["Depth"] > 100) & (data["1KGEAS"] < 0.01)]
    elif not (data[~data["GnomAD"].isin(["-"])].empty and data[~data["1KGEAS"].isin(["-"])].empty): ## have GnomAD and 1KGEAS annotation
        data = data[(~data["GnomAD"].isin(["-"])) & (~data["1KGEAS"].isin(["-"]))]
        data = data[(data["Depth"] > 100) & (data["GnomAD"] < 0.01) & (data["1KGEAS"] < 0.01)]
    data1 = []
    data2 = []
    data1 = pd.DataFrame(data1)
    data2 = pd.DataFrame(data2)
    ## hotspot loci
    if not data[data["Tag"].str.contains("HotSpot")].empty:
        data1 = data[data["Tag"].str.contains("HotSpot")]
        data1 = data1[data1["Freq"] >= 0.0002]
    ## non-hotspot loci, select 0.2% <= freq <= 10%
    if not data[~data["Tag"].str.contains("HotSpot")].empty:
        data2 = data[~data["Tag"].str.contains("HotSpot")]
        # data2 = data2[(data2["Freq"] >= 0.002) & (data2["Freq"] <= 0.1)]
        data2 = data2[data2["Freq"] >= 0.002]
        data2 = data2[data2["Tag"].map(lambda x: "WBC" not in x.split(";"))]
        data2 = data2[data2["Tag"].map(lambda x: "Germline" not in x.split(";"))]
        data2 = data2[data2["Tag"].map(lambda x: "Polymorphism" not in x.split(";"))]
        data2 = data2[data2["Tag"].map(lambda x: "BackGround" not in x.split(";"))]
        # data2 = data2[data2["Tag"].map(lambda x: "PW" not in x.split(";"))]
        data2 = data2[
            ~((data2["Tag"].str.contains("Polymer|STR")) & (data2["Type"].str.contains("Intronic|3'UTR|5'UTR")))]
        data2 = data2[~data2["Tag"].str.contains("Bad|OutOfReg|HighClip|HighEN|PoorSB|LowRatio|Noise")]
        # data = data[~data["Tag"].str.contains("Bad|OutOfReg|Polym|HighClip|HighEN|PoorSB|BackGround|LowRatio|Noise|STR")]
        data2 = data2[~((data2["Gene"].str.contains("MUC16")) & (data2["Type"].str.contains("Intronic|3'UTR|5'UTR")))]
    if not (data1.empty and data2.empty):
        data = pd.concat([data1, data2])
    elif data1.empty and not data2.empty:
        data = data2
    elif not data1.empty and data2.empty:
        data = data1
    else:
        print("Not data output")
        exit(0)

    return data


def transform(data):
    dic = {}
    samples = data["Sample"]
    chrs = data["Chr"]
    starts = data["Start"]
    ends = data["End"]
    refs = data["Ref"]
    alts = data["Alt"]
    deps = data["Depth"]
    freqs = data["Freq"]
    genes = data["Gene"]
    geneids = data["Gene"]
    types = data["Type"]

    uniq_samples = set(list(samples))
    samples = np.array(samples)
    title = "##fileformat=VCFv4.1\n" + "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])+"\n"

    for sample in uniq_samples:
        indx_old = list(np.where(samples==sample))
        ## filter index that start != end(i.e indel locus but not SNP locus )
        indx = []
        start_old = list(np.array(starts)[indx_old])
        end_old = list(np.array(ends)[indx_old])
        check_pos_list = list(zip(start_old,end_old))
        for i,loc in enumerate(check_pos_list):
            if loc[0] == loc[1]:  ##select start equal end locus
                indx.append(indx_old[0][i])
        indx = np.array(indx)
        if indx.any():        
            chr = [x.replace("chr","").replace("Chr","") for x in list(np.array(chrs)[indx])]
            start = list(np.array(starts)[indx])
        #end = list(np.array(ends)[indx])
            ref = list(np.array(refs)[indx])
            alt = list(np.array(alts)[indx])
            length = len(indx)
            id = ["."]*length
            qual = ["."]*length
            filter = ["."]*length
            dep = [ "Depth="+str(x) for x in list(np.array(deps)[indx])]
            freq = [ "Freq="+str(x) for x in list(np.array(freqs)[indx])]
            gene = ["Gene="+str(x) for x in list(np.array(genes)[indx])]
            geneid = ["GeneID="+str(x) for x in list(np.array(geneids)[indx])]
            type = ["Variant_Classification="+str(x) for x in list(np.array(types)[indx])]

            vcf_line = ["\t".join([str(y) for y in x]) for x in list(zip(chr, start, id, ref, alt,qual,filter))]
            info_line = [";".join([str(y) for y in x]) for x in list(zip(gene, geneid, dep, freq,type))]
            vcf_and_info_line = title+ "\n".join(["\t".join([y for y in x ]) for x in list(zip(vcf_line,info_line))])
            dic[sample] = vcf_and_info_line

    return dic

if __name__ == "__main__":
    args = getargs()
    data = pd.read_excel(args.dnafile, sheet_name=["Sheet1"])
    #hotsomatic = data['SNVIndelHotSomatic']
    #nonhotsomatic = data['SNVIndelSomatic']
    #all_somatic = pd.concat([hotsomatic,nonhotsomatic])
    all_somatic = data["Sheet1"]
    #dic = transform(filter(all_somatic))
    dic = transform(all_somatic)
    print(dic)
    for sample in dic:
        #print(dic[sample])
        fw = open(sample+'.vcf','w').write(dic[sample])
