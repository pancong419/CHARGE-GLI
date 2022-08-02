#!/usr/bin/python

###################################
### Copyright (C) 2022 Cong Pan ###
### version 0.2, April 21, 2022   ###
###################################

# NOTE: this script is used to convert GEM (v1.4.1 or later) output to the meta-analysis-ready format for Phase2 CHARGE Gene-Lifestyle Interactions projects
# Please run using python3
# You can check the python version using
# python --version
# Example command:
# python GEM2CHARGEGLI.py Phase2.ARIC.AA.HDL.LTST.M1.COMBINED.chr22.out Phase2.ARIC.AA.HDL.LTST.M2.COMBINED.chr22.out chr22.info.gz PHASE2.ARIC.AA.HDL.LTST.COMBINED.20220216.txt LTST SNP Rsq Genotyped tab
# For questions or comments, please contact Cong.Pan.2 AT uth.tmc.edu
# v0.2: added formatting quantitative exposure 


import sys
import gzip
from scipy import stats

infile1 = str(sys.argv[1]) # input file name(s) of GEM Model 1 output
infile2 = str(sys.argv[2]) # input file name(s) of GEM Model 2 output
snpinfofile = str(sys.argv[3]) # SNP info file should at least contain the following information: SNPID, INFO, IMPUTED
# NOTE: if GEM was run separately by each chromosome, please only use the file names for chr22 output and SNP info file here: the script assumes results from all other 21 autosomes are located in the same directories respectively, with the same naming convention and will loop over starting from chr1, so there is no need to specify 66 different input files
# e.g., Phase2.ARIC.AA.HDL.LTST.M1.COMBINED.chr22.out
# and Phase2.ARIC.AA.HDL.LTST.M2.COMBINED.chr22.out
# and chr22.info.gz
# If you have concatenated results from all 22 autosomes into a single file, just provide three single file names for Models 1 and 2 here (without "chr22")
# e.g., Phase2.ARIC.AA.HDL.LTST.M1.COMBINED.all.out
# and Phase2.ARIC.AA.HDL.LTST.M2.COMBINED.all.out
# and snpinfo.txt
outfile = str(sys.argv[4]) # your output file name (without .gz), the results will be gzipped automatically
# e.g., PHASE2.ARIC.AA.HDL.LTST.COMBINED.20220216.txt
# the final gzipped file would be PHASE2.ARIC.AA.HDL.LTST.COMBINED.20220216.txt.gz
Ename = str(sys.argv[5]) # variable name for the environmental exposure used in the gene-environment interaction analyses, as shown in your GEM output (NOTE: this is CaSE-senSItiVe)
# e.g., LTST, STST, or ltst, stst
SNPIDname = str(sys.argv[6]) # SNPID column name in snpinfofile, e.g., SNP
INFOname = str(sys.argv[7]) # INFO column name in snpinfofile, e.g., Rsq
IMPUTEDname = str(sys.argv[8]) # IMPUTED column name in snpinfofile, e.g., Genotyped
# NOTE: if the IMPUTED column is coded in your snpinfofile (instead of using "Imputed" and "Genotyped" as values), please make sure that 1 = Imputed and 0 = Genotyped
snpinfofile_delim = str(sys.argv[9]) # delimiter for the snpinfofile, e.g., tab, space or ,
na_string = "."
maf_cutoff = 0.001
info_cutoff = 0.3

print ("Output file: %s.gz"  % outfile)
print ("Environmental variable name in GEM output: %s" % Ename)
print ("SNPID column name in %s: %s" %(snpinfofile, SNPIDname))
print ("INFO column name in %s: %s" % (snpinfofile, INFOname))
print ("IMPUTED column name in %s: %s" % (snpinfofile, IMPUTEDname))
print ("Delimiter in %s: %s" % (snpinfofile, snpinfofile_delim))
if snpinfofile_delim == "tab":
    snpinfofile_delim = "\t"
elif snpinfofile_delim == "space":
    snpinfofile_delim = " "
print ("Missing value string in the output file: %s" % na_string)
print ("Minor allele frequency cut-off: %f" % maf_cutoff)
print ("INFO cut-off: %f" % info_cutoff)

outfile_handle = gzip.open(outfile + ".gz", "wb")
if not "chr22" in infile1 and not "chr22" in infile2 and not "chr22" in snpinfofile:
    print ("Reading a single SNP info file: %s" % snpinfofile)
    if snpinfofile.endswith(".gz"):
        snpinfofile_handle = gzip.open(snpinfofile, "rt", encoding='utf-8')
    else:
        snpinfofile_handle = open(snpinfofile, "r")
        
    line = snpinfofile_handle.readline()
    line = line.strip()
    fields = line.split(snpinfofile_delim)
    columns = {}
    imputed = {}
    info = {}
    for i in range(len(fields)):
        columns[fields[i]] = i
    while line:
        line = snpinfofile_handle.readline()
        line = line.strip()
        fields = line.split(snpinfofile_delim)
        if line == "":
            break
        if fields[columns[INFOname]].upper() == "NA" or fields[columns[INFOname]].upper() == "NAN" or fields[columns[INFOname]].upper() == "." or float(fields[columns[INFOname]]) < info_cutoff:
            continue
        if fields[columns[IMPUTEDname]].upper() == "IMPUTED":
            imputed[fields[columns[SNPIDname]]] = "1"
        elif fields[columns[IMPUTEDname]].upper() == "GENOTYPED":
            imputed[fields[columns[SNPIDname]]] = "0"
        else:
            imputed[fields[columns[SNPIDname]]] = fields[columns[IMPUTEDname]]
        info[fields[columns[SNPIDname]]] = fields[columns[INFOname]]
    snpinfofile_handle.close()

    print ("Reading GEM output...")
    print ("Model 1 input file: %s" % infile1)
    print ("Model 2 input file: %s" % infile2)
    if infile1.endswith(".gz"):
        infile1_handle = gzip.open(infile1, "rt", encoding='utf-8')
    else:
        infile1_handle = open(infile1, "r")
    line1 = infile1_handle.readline()
    line1 = line1.strip()
    fields1 = line1.split("\t")
    columns1 = {}
    for i in range(len(fields1)):
        columns1[fields1[i]] = i
    if not "AF_" + Ename + "_0" in columns1 or not "AF_" + Ename + "_1" in columns1 or not "N_" + Ename + "_1" in columns1:
        quantE = True
        print("quantE is True")
    else:
        quantE = False
        print("quantE is False")
    if infile2.endswith(".gz"):
        infile2_handle = gzip.open(infile2, "rt", encoding='utf-8')
    else:
        infile2_handle = open(infile2, "r")
    line2 = infile2_handle.readline()
    line2 = line2.strip()
    fields2 = line2.split("\t")
    columns2 = {}
    for i in range(len(fields2)):
        columns2[fields2[i]] = i
    if quantE:
        outfile_handle.write(b"SNPID\tCHR\tPOS\tINFO\tIMPUTED\tEFFECT_ALLELE\tNON_EFFECT_ALLELE\tEAF_ALL\tN\tBETA_SNP_M2\tSE_SNP_M2\tP_SNP_M2\tBETA_SNP_M1\tSE_SNP_M1_MB\tP_SNP_M1_MB\tSE_SNP_M1_ROBUST\tP_SNP_M1_ROBUST\tBETA_INT\tSE_INT_MB\tP_INT_MB\tSE_INT_ROBUST\tP_INT_ROBUST\tP_JOINT_MB\tCOV_SNP_INT_MB\tP_JOINT_ROBUST\tCOV_SNP_INT_ROBUST\n")
    else:
        outfile_handle.write(b"SNPID\tCHR\tPOS\tINFO\tIMPUTED\tEFFECT_ALLELE\tNON_EFFECT_ALLELE\tEAF_ALL\tEAF_E0\tEAF_E1\tN\tN_EXP\tBETA_SNP_M2\tSE_SNP_M2\tP_SNP_M2\tBETA_SNP_M1\tSE_SNP_M1_MB\tP_SNP_M1_MB\tSE_SNP_M1_ROBUST\tP_SNP_M1_ROBUST\tBETA_INT\tSE_INT_MB\tP_INT_MB\tSE_INT_ROBUST\tP_INT_ROBUST\tP_JOINT_MB\tCOV_SNP_INT_MB\tP_JOINT_ROBUST\tCOV_SNP_INT_ROBUST\n")
    line_ct = 1
    while line1 and line2:
        line1 = infile1_handle.readline()
        line1 = line1.strip()
        fields1 = line1.split("\t")
        if line1 == "":
            break
        for i in range(len(fields1)):
            if fields1[i].upper() == "NA" or fields1[i].upper() == "NAN" or fields1[i].upper() == "-NAN" or fields1[i].upper() == "INF":
                fields1[i] = na_string
        line2 = infile2_handle.readline()
        line2 = line2.strip()
        fields2 = line2.split("\t")
        if line2 == "":
            break
        for i in range(len(fields2)):
            if fields2[i].upper() == "NA" or fields2[i].upper() == "NAN" or fields2[i].upper() == "-NAN" or fields2[i].upper() == "INF":
                fields2[i] = na_string
        line_ct = line_ct + 1   
        if fields1[columns1["SNPID"]] != fields2[columns2["SNPID"]]:
            print ("Error occurred at line %d" % line_ct)
            print ("SNPID in Model 1 GEM output: %s" % fields1[columns1["SNPID"]])
            print ("SNPID in Model 2 GEM output: %s" % fields2[columns2["SNPID"]])
            sys.exit(1)
        if fields1[columns1["CHR"]] != fields2[columns2["CHR"]]:
            print ("Error occurred at line %d" % line_ct)
            print ("CHR in Model 1 GEM output: %s" % fields1[columns1["CHR"]])
            print ("CHR in Model 2 GEM output: %s" % fields2[columns2["CHR"]])
            sys.exit(1)
        if fields1[columns1["POS"]] != fields2[columns2["POS"]]:
            print ("Error occurred at line %d" % line_ct)
            print ("POS in Model 1 GEM output: %s" % fields1[columns1["POS"]])
            print ("POS in Model 2 GEM output: %s" % fields2[columns2["POS"]])
            sys.exit(1)
        if fields1[columns1["Non_Effect_Allele"]] != fields2[columns2["Non_Effect_Allele"]]:
            print ("Error occurred at line %d" % line_ct)
            print ("Non_Effect_Allele in Model 1 GEM output: %s" % fields1[columns1["Non_Effect_Allele"]])
            print ("Non_Effect_Allele in Model 2 GEM output: %s" % fields2[columns2["Non_Effect_Allele"]])
            sys.exit(1)
        if fields1[columns1["Effect_Allele"]] != fields2[columns2["Effect_Allele"]]:
            print ("Error occurred at line %d" % line_ct)
            print ("Effect_Allele in Model 1 GEM output: %s" % fields1[columns1["Effect_Allele"]])
            print ("Effect_Allele in Model 2 GEM output: %s" % fields2[columns2["Effect_Allele"]])
            sys.exit(1)
        if fields1[columns1["N_Samples"]] != fields2[columns2["N_Samples"]]:
            print ("Error occurred at line %d" % line_ct)
            print ("N_Samples in Model 1 GEM output: %s" % fields1[columns1["N_Samples"]])
            print ("N_Samples in Model 2 GEM output: %s" % fields2[columns2["N_Samples"]])
            sys.exit(1)
        if fields1[columns1["AF"]] != fields2[columns2["AF"]]:
            print ("Error occurred at line %d" % line_ct)
            print ("AF in Model 1 GEM output: %s" % fields1[columns1["AF"]])
            print ("AF in Model 2 GEM output: %s" % fields2[columns2["AF"]])
            sys.exit(1)
        if not fields1[columns1["SNPID"]] in info or not fields1[columns1["SNPID"]] in imputed:
            continue
        if fields1[columns1["AF"]] == na_string or float(fields1[columns1["AF"]]) < maf_cutoff or float(fields1[columns1["AF"]]) > 1 - maf_cutoff:
            continue
        if fields1[columns1["Beta_G"]] == na_string or fields1[columns1["SE_Beta_G"]] == na_string or float(fields1[columns1["SE_Beta_G"]]) <= 0:
            p_snp_m1_mb = na_string
        else:
            p_snp_m1_mb = "{:.6}".format(stats.chi2.sf((float(fields1[columns1["Beta_G"]])/float(fields1[columns1["SE_Beta_G"]]))**2, 1))
            
        if fields1[columns1["Beta_G"]] == na_string or fields1[columns1["robust_SE_Beta_G"]] == na_string or float(fields1[columns1["robust_SE_Beta_G"]]) <= 0:    
            p_snp_m1_robust = na_string
        else:
            p_snp_m1_robust = "{:.6}".format(stats.chi2.sf((float(fields1[columns1["Beta_G"]])/float(fields1[columns1["robust_SE_Beta_G"]]))**2, 1))
        if quantE:
            outfile_handle.write(str.encode(fields1[columns1["SNPID"]]) + b"\t" + str.encode(fields1[columns1["CHR"]]) + b"\t" + str.encode(fields1[columns1["POS"]]) + b"\t" + str.encode(info[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(imputed[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(fields1[columns1["Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["Non_Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["AF"]]) + b"\t"  + str.encode(fields1[columns1["N_Samples"]]) + b"\t"  + str.encode(fields2[columns2["Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["SE_Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["P_Value_Marginal"]]) + b"\t" + str.encode(fields1[columns1["Beta_G"]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_mb) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_robust) + b"\t" + str.encode(fields1[columns1["Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["Cov_Beta_G_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["robust_Cov_Beta_G_G-" + Ename]]) + b"\n")
        else:
            outfile_handle.write(str.encode(fields1[columns1["SNPID"]]) + b"\t" + str.encode(fields1[columns1["CHR"]]) + b"\t" + str.encode(fields1[columns1["POS"]]) + b"\t" + str.encode(info[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(imputed[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(fields1[columns1["Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["Non_Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["AF"]]) + b"\t" + str.encode(fields1[columns1["AF_" + Ename + "_0"]]) + b"\t" + str.encode(fields1[columns1["AF_" + Ename + "_1"]]) + b"\t" + str.encode(fields1[columns1["N_Samples"]]) + b"\t" + str.encode(fields1[columns1["N_" + Ename + "_1"]]) + b"\t" + str.encode(fields2[columns2["Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["SE_Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["P_Value_Marginal"]]) + b"\t" + str.encode(fields1[columns1["Beta_G"]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_mb) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_robust) + b"\t" + str.encode(fields1[columns1["Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["Cov_Beta_G_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["robust_Cov_Beta_G_G-" + Ename]]) + b"\n")
    infile1_handle.close()
    infile2_handle.close()

elif "chr22" in infile1 and "chr22" in infile2 and "chr22" in snpinfofile:
    snpinfofile_parts = snpinfofile.split("chr22")
    infile1_parts = infile1.split("chr22")
    infile2_parts = infile2.split("chr22")
    for chr in range(1, 23):
        current_snpinfofile = snpinfofile_parts[0] + "chr" + str(chr) + snpinfofile_parts[1]
        print ("Reading SNP info file: %s" % current_snpinfofile)    
        if current_snpinfofile.endswith(".gz"):
            snpinfofile_handle = gzip.open(current_snpinfofile, "rt", encoding='utf-8')
        else:
            snpinfofile_handle = open(current_snpinfofile, "r")
        line = snpinfofile_handle.readline()
        line = line.strip()
        
        fields = line.split(snpinfofile_delim)
        columns = {}
        imputed = {}
        info = {}
        for i in range(len(fields)):
            columns[fields[i]] = i
        while line:
            line = snpinfofile_handle.readline()
            line = line.strip()
            fields = line.split(snpinfofile_delim)
            if line == "":
                break
            if fields[columns[INFOname]].upper() == "NA" or fields[columns[INFOname]].upper() == "NAN" or fields[columns[INFOname]].upper() == "." or float(fields[columns[INFOname]]) < info_cutoff:
                continue
            if fields[columns[IMPUTEDname]].upper() == "IMPUTED":
                imputed[fields[columns[SNPIDname]]] = "1"
            elif fields[columns[IMPUTEDname]].upper() == "GENOTYPED":
                imputed[fields[columns[SNPIDname]]] = "0"
            else:
                imputed[fields[columns[SNPIDname]]] = fields[columns[IMPUTEDname]]
            info[fields[columns[SNPIDname]]] = fields[columns[INFOname]]
        snpinfofile_handle.close()

        current_infile1 = infile1_parts[0] + "chr" + str(chr) + infile1_parts[1]
        current_infile2 = infile2_parts[0] + "chr" + str(chr) + infile2_parts[1]
        print ("Reading GEM output...")
        print ("Model 1 input file: %s" % current_infile1)
        print ("Model 2 input file: %s" % current_infile2)
        if current_infile1.endswith(".gz"):
            infile1_handle = gzip.open(current_infile1, "rt", encoding='utf-8')
        else:
            infile1_handle = open(current_infile1, "r")
        line1 = infile1_handle.readline()
        line1 = line1.strip()
        fields1 = line1.split("\t")
        columns1 = {}
        for i in range(len(fields1)):
            columns1[fields1[i]] = i
        if current_infile2.endswith(".gz"):
            infile2_handle = gzip.open(current_infile2, "rt", encoding='utf-8')
        else:
            infile2_handle = open(current_infile2, "r")
        line2 = infile2_handle.readline()
        line2 = line2.strip()
        fields2 = line2.split("\t")
        columns2 = {}
        for i in range(len(fields2)):
            columns2[fields2[i]] = i
        if chr == 1:
            if not "AF_" + Ename + "_0" in columns1 or not "AF_" + Ename + "_1" in columns1 or not "N_" + Ename + "_1" in columns1:
                quantE = True
            else:
                quantE = False
            if quantE:
                outfile_handle.write(b"SNPID\tCHR\tPOS\tINFO\tIMPUTED\tEFFECT_ALLELE\tNON_EFFECT_ALLELE\tEAF_ALL\tN\tBETA_SNP_M2\tSE_SNP_M2\tP_SNP_M2\tBETA_SNP_M1\tSE_SNP_M1_MB\tP_SNP_M1_MB\tSE_SNP_M1_ROBUST\tP_SNP_M1_ROBUST\tBETA_INT\tSE_INT_MB\tP_INT_MB\tSE_INT_ROBUST\tP_INT_ROBUST\tP_JOINT_MB\tCOV_SNP_INT_MB\tP_JOINT_ROBUST\tCOV_SNP_INT_ROBUST\n")
            else:
                outfile_handle.write(b"SNPID\tCHR\tPOS\tINFO\tIMPUTED\tEFFECT_ALLELE\tNON_EFFECT_ALLELE\tEAF_ALL\tEAF_E0\tEAF_E1\tN\tN_EXP\tBETA_SNP_M2\tSE_SNP_M2\tP_SNP_M2\tBETA_SNP_M1\tSE_SNP_M1_MB\tP_SNP_M1_MB\tSE_SNP_M1_ROBUST\tP_SNP_M1_ROBUST\tBETA_INT\tSE_INT_MB\tP_INT_MB\tSE_INT_ROBUST\tP_INT_ROBUST\tP_JOINT_MB\tCOV_SNP_INT_MB\tP_JOINT_ROBUST\tCOV_SNP_INT_ROBUST\n")
  
        line_ct = 1
        while line1 and line2:
            line1 = infile1_handle.readline()
            line1 = line1.strip()
            fields1 = line1.split("\t")
            if line1 == "":
                break
            for i in range(len(fields1)):
                if fields1[i].upper() == "NA" or fields1[i].upper() == "NAN" or fields1[i].upper() == "-NAN" or fields1[i].upper() == "INF":
                    fields1[i] = na_string
            line2 = infile2_handle.readline()
            line2 = line2.strip()
            fields2 = line2.split("\t")
            if line2 == "":
                break
            for i in range(len(fields2)):
                if fields2[i].upper() == "NA" or fields2[i].upper() == "NAN" or fields2[i].upper() == "-NAN" or fields2[i].upper() == "INF":
                    fields2[i] = na_string
            line_ct = line_ct + 1
            if fields1[columns1["SNPID"]] != fields2[columns2["SNPID"]]:
                print ("Error occurred at line %d" % line_ct)
                print ("SNPID in Model 1 GEM output: %s" % fields1[columns1["SNPID"]])
                print ("SNPID in Model 2 GEM output: %s" % fields2[columns2["SNPID"]])
                sys.exit(1)
            if fields1[columns1["CHR"]] != fields2[columns2["CHR"]]:
                print ("Error occurred at line %d" % line_ct)
                print ("CHR in Model 1 GEM output: %s" % fields1[columns1["CHR"]])
                print ("CHR in Model 2 GEM output: %s" % fields2[columns2["CHR"]])
                sys.exit(1)
            if fields1[columns1["POS"]] != fields2[columns2["POS"]]:
                print ("Error occurred at line %d" % line_ct)
                print ("POS in Model 1 GEM output: %s" % fields1[columns1["POS"]])
                print ("POS in Model 2 GEM output: %s" % fields2[columns2["POS"]])
                sys.exit(1)
            if fields1[columns1["Non_Effect_Allele"]] != fields2[columns2["Non_Effect_Allele"]]:
                print ("Error occurred at line %d" % line_ct)
                print ("Non_Effect_Allele in Model 1 GEM output: %s" % fields1[columns1["Non_Effect_Allele"]])
                print ("Non_Effect_Allele in Model 2 GEM output: %s" % fields2[columns2["Non_Effect_Allele"]])
                sys.exit(1)
            if fields1[columns1["Effect_Allele"]] != fields2[columns2["Effect_Allele"]]:
                print ("Error occurred at line %d" % line_ct)
                print ("Effect_Allele in Model 1 GEM output: %s" % fields1[columns1["Effect_Allele"]])
                print ("Effect_Allele in Model 2 GEM output: %s" % fields2[columns2["Effect_Allele"]])
                sys.exit(1)
            if fields1[columns1["N_Samples"]] != fields2[columns2["N_Samples"]]:
                print ("Error occurred at line %d" % line_ct)
                print ("N_Samples in Model 1 GEM output: %s" % fields1[columns1["N_Samples"]])
                print ("N_Samples in Model 2 GEM output: %s" % fields2[columns2["N_Samples"]])
                sys.exit(1)
            if fields1[columns1["AF"]] != fields2[columns2["AF"]]:
                print ("Error occurred at line %d" % line_ct)
                print ("AF in Model 1 GEM output: %s" % fields1[columns1["AF"]])
                print ("AF in Model 2 GEM output: %s" % fields2[columns2["AF"]])
                sys.exit(1)
            if not fields1[columns1["SNPID"]] in info or not fields1[columns1["SNPID"]] in imputed:
                continue
            if fields1[columns1["AF"]] == na_string or float(fields1[columns1["AF"]]) < maf_cutoff or float(fields1[columns1["AF"]]) > 1 - maf_cutoff:
                continue
            if fields1[columns1["Beta_G"]] == na_string or fields1[columns1["SE_Beta_G"]] == na_string or float(fields1[columns1["SE_Beta_G"]]) <= 0:
                p_snp_m1_mb = na_string
            else:
                p_snp_m1_mb = "{:.6}".format(stats.chi2.sf((float(fields1[columns1["Beta_G"]])/float(fields1[columns1["SE_Beta_G"]]))**2, 1))            
            if fields1[columns1["Beta_G"]] == na_string or fields1[columns1["robust_SE_Beta_G"]] == na_string or float(fields1[columns1["robust_SE_Beta_G"]]) <= 0:    
                p_snp_m1_robust = na_string
            else:
                p_snp_m1_robust = "{:.6}".format(stats.chi2.sf((float(fields1[columns1["Beta_G"]])/float(fields1[columns1["robust_SE_Beta_G"]]))**2, 1))
            if quantE:
                outfile_handle.write(str.encode(fields1[columns1["SNPID"]]) + b"\t" + str.encode(fields1[columns1["CHR"]]) + b"\t" + str.encode(fields1[columns1["POS"]]) + b"\t" + str.encode(info[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(imputed[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(fields1[columns1["Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["Non_Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["AF"]]) + b"\t"  + str.encode(fields1[columns1["N_Samples"]]) + b"\t"  + str.encode(fields2[columns2["Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["SE_Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["P_Value_Marginal"]]) + b"\t" + str.encode(fields1[columns1["Beta_G"]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_mb) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_robust) + b"\t" + str.encode(fields1[columns1["Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["Cov_Beta_G_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["robust_Cov_Beta_G_G-" + Ename]]) + b"\n")
            else:
                outfile_handle.write(str.encode(fields1[columns1["SNPID"]]) + b"\t" + str.encode(fields1[columns1["CHR"]]) + b"\t" + str.encode(fields1[columns1["POS"]]) + b"\t" + str.encode(info[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(imputed[fields1[columns1["SNPID"]]]) + b"\t" + str.encode(fields1[columns1["Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["Non_Effect_Allele"]]) + b"\t" + str.encode(fields1[columns1["AF"]]) + b"\t" + str.encode(fields1[columns1["AF_" + Ename + "_0"]]) + b"\t" + str.encode(fields1[columns1["AF_" + Ename + "_1"]]) + b"\t" + str.encode(fields1[columns1["N_Samples"]]) + b"\t" + str.encode(fields1[columns1["N_" + Ename + "_1"]]) + b"\t" + str.encode(fields2[columns2["Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["SE_Beta_Marginal"]]) + b"\t" + str.encode(fields2[columns2["P_Value_Marginal"]]) + b"\t" + str.encode(fields1[columns1["Beta_G"]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_mb) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G"]]) + b"\t" + str.encode(p_snp_m1_robust) + b"\t" + str.encode(fields1[columns1["Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["robust_SE_Beta_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Interaction"]]) + b"\t" + str.encode(fields1[columns1["P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["Cov_Beta_G_G-" + Ename]]) + b"\t" + str.encode(fields1[columns1["robust_P_Value_Joint"]]) + b"\t" + str.encode(fields1[columns1["robust_Cov_Beta_G_G-" + Ename]]) + b"\n")
        infile1_handle.close()
        infile2_handle.close()

else:
    print ("Error: check your input file and SNP info file names! They must all contain \"chr22\" or all NOT contain \"chr22\"")
    print ("Model 1 input file: %s" % infile1)
    print ("Model 2 input file: %s" % infile2)
    print ("SNP info file: %s" % snpinfofile)
    sys.exit(1)

outfile_handle.close()
    
