import random
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

random.seed(0)

components = ["Virus","ACE2","PKC","ANG_2_T1R","ANG_2","ADAM_17","SIL6R","TLR4","RIG1","NFKB","TNF","TNF_alpha", "IRF3","STAT1","ISG","C_FLIP","INF_A_B","NRLP3","CASP1","FOXO3A","IFNR","ProCASP8","DISC","BCL_2","tBID","Bax_Bak","CASP9","ROS","TNFR","FADD","Gasdermin","Pyroptosis","IL1","IL1R","IL6", "MLKL","Necroptosis","RIPK1_3","CASP8","ProCASP3_7","CASP3_7","Apoptosis"]
x_dict = {k : 0 for k in components}

external = ["TNF_e", "Virus_e", "ROS_e", "IFN_e"]
external_dict = {k : 0 for k in external}
external_dict["Virus_e"] = 1

x_dict["Virus"] = 0
#x_dict["ROS"] = 1

#x_begin = x_dict.copy()

n_iter = 25
mat = np.zeros((n_iter + 1, len(x_dict)))
for j in range(0, n_iter):
    r = random.sample(list(x_dict), len(x_dict))
    for i in r:
        if (i == "Virus"):
            x_dict[i] = (not x_dict["STAT1"] and x_dict["Virus"]) or external_dict["Virus_e"] # or x_dict["RIG1"]
            #x_dict[i] = 1
            # x_dict[i] = 0 # possibly make virus = 0 
        if (i == "ACE2"):
            # PKC mediates ACE2 shedding from tubular cells
            x_dict[i] = x_dict["PKC"] and not (x_dict["ADAM_17"] or x_dict["RIG1"])
            # x_dict[i] = not x_dict["Virus"] 
            # not sure about if there is a relation since Virus just relies on ACE2 to enter cells, the presence of ACE2 promotes disease prog
        if (i == "PKC"):
            x_dict[i] = x_dict["ANG_2_T1R"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9000463/
        if (i == "ANG_2_T1R"):
            x_dict[i] = x_dict["ANG_2"] 
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8193025/
            # Ang II acts on angiotensin type 1 (AT1") receptors and activates the NADPH-oxidase complex producing superoxide and promoting cell pro-oxidative and pro-inflammatory responses
        if (i == "ANG_2"):
            # ACE2 converts Ang II into Ang-(1–7")
            x_dict[i] = not x_dict["ACE2"]
        if (i == "ADAM_17"):
            # A2_t1r activates ADAM 17, which promotes ACE2 shedding (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8185693/")
            x_dict[i] = x_dict["ANG_2_T1R"]
        if (i == "SIL6R"):
            x_dict[i] = x_dict["ADAM_17"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7224649/
        if (i == "TLR4"): 
            x_dict[i] = x_dict["Virus"]
            # Spike glycoprotein, the major infective surface protein of SARS-CoV-2 has been found as a ligand for human TLR4
            # https://www.futuremedicine.com/doi/full/10.2217/fvl-2021-0249
        if (i == "RIG1"):
            x_dict[i] = x_dict["Virus"]
            # antiviral activity of RIG-1 may comprise inhibition of viral entry into the host cell by preventing the expression of its receptor, ACE2
            # https://www.news-medical.net/news/20210215/RIG-1-like-receptors-may-play-dominant-role-in-suppressing-SARS-CoV-2-infection.aspx
        if (i == "NFKB"): 
            x_dict[i] = x_dict["TNFR"] or x_dict["C_FLIP"] or x_dict["IL1R"] #x_dict["ANG_2_T1R"] or x_dict["PKC"] or x_dict["RIG1"] or
            # should be good, may need to find what exactly inhibits NFKB
            # common drug therapy is inhibiting NFKB (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7648206/")
            
            ############################## NEW ##############################
            # not sure if "or" or "and"
        
        if (i == "TNF"):
            x_dict[i] = x_dict["NFKB"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7114322/
        if (i == "IRF3"):
            x_dict[i] = x_dict["RIG1"] and not x_dict["Virus"]
            # https://journals.asm.org/doi/10.1128/CMR.00299-20
            # SARS-CoV-2 membrane protein binds to importin karyopherin subunit alpha-6 (KPNA6) to inhibit interferon regulatory factor 3(IRF3) nuclear translocation
            # https://www.frontiersin.org/articles/10.3389/fcimb.2021.766922/full
        if (i == "STAT1"):
            x_dict[i] = (x_dict["ISG"] or x_dict["IFNR"]) and not x_dict["Virus"]
            # After the infection, STAT1 activity is inhibited by the SARS-CoV-2 proteins, NSP1, and ORF6 
            # https://www.nature.com/articles/s41418-020-00633-7
        if (i == "ISG"):
            x_dict[i] = x_dict["STAT1"] # x_dict["Virus"]
            # https://www.nature.com/articles/s41586-021-03234-7
        if (i == "C_FLIP"):
            x_dict[i] = x_dict["NFKB"] and not x_dict["FOXO3A"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == "INF_A_B"):
            x_dict[i] = x_dict["IRF3"] and not x_dict["Virus"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7995242/
        if (i == "NRLP3"):
            x_dict[i] = x_dict["NFKB"] and x_dict["Virus"]
            # https://www.nature.com/articles/ni.3772
        if (i == "CASP1"):
            x_dict[i] = x_dict["NRLP3"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6651423/
        if (i == "FOXO3A"):
            x_dict[i] = x_dict["Virus"] 
            ## this should be right if "Virus" causes stress
            #x_dict[i] = x_dict[i]
            # no direct relation (drug targetting in place - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8187014/")
        if (i == "IFNR"):
            x_dict[i] = x_dict["INF_A_B"] or external_dict["IFN_e"]
            # need to confirm more
            # https://www.frontiersin.org/articles/10.3389/fimmu.2020.606456/full
        if (i == "ProCASP8"):
            x_dict[i] = x_dict["DISC"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8899608/
        if (i == "DISC"):
            x_dict[i] = (x_dict["ProCASP8"] or x_dict["FADD"]) and not x_dict["C_FLIP"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == "BCL_2"):
            x_dict[i] = x_dict["NFKB"]
            # https://www.nature.com/articles/1204926
        if (i == "tBID"):
            x_dict[i] = x_dict["CASP8"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4882451/
        if (i == "Bax_Bak"):
            x_dict[i] = x_dict["BCL_2"] or x_dict["tBID"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == "CASP9"):
            x_dict[i] = x_dict["Bax_Bak"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4219646/
        if (i == "ROS"):
            x_dict[i] = external_dict["ROS_e"] or x_dict["SIL6R"] ## new
            # x_dict["IL6"] or x_dict["ANG_2_T1R"]

            # apparently IL6 and ANG2R (by NOX) promotes ROS -- NEW
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7768996/ -- look at paper for more details
        if (i == "TNFR"):
            x_dict[i] = x_dict["TNF"] or external_dict["TNF_e"]
            # Do more research later
            # https://www.frontiersin.org/articles/10.3389/fimmu.2020.585880/full
        if (i == "TNF_alpha"):
            x_dict[i] = x_dict["NFKB"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7114322/ -- NEW
            # https://www.frontiersin.org/articles/10.3389/fpubh.2022.833967/full#B22
            # ^^ vague info
            ## info mostly about how drugs target TNF alpha but no reason stated
        if (i == "FADD"):
            x_dict[i] = x_dict["TNFR"]
            # Do more research on this later
            # https://pubmed.ncbi.nlm.nih.gov/9430227/
        if (i == "Gasdermin"):
            x_dict[i] = x_dict["CASP1"]
            # Do more research on this later
            # https://www.cell.com/immunity/fulltext/S1074-7613(20")30237-5?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1074761320302375%3Fshowall%3Dtrue
        if (i == "Pyroptosis"):
            x_dict[i] = x_dict["Gasdermin"]
            # https://www.nature.com/articles/s41467-019-09753-2
        if (i == "IL1"):
            x_dict[i] = x_dict["MLKL"] or x_dict["NFKB"] # changed from "or" - NEW
            # Look into this more later
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
        if (i == "IL1R"):
            x_dict[i] = x_dict["IL1"]
        if (i == "IL6"):
            x_dict[i] = x_dict["NFKB"] and x_dict["ROS"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7768996/ ---- NEW
        if (i == "MLKL"):
            x_dict[i] = x_dict["RIPK1_3"] and not x_dict["CASP1"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
            # read more of this later
        if (i == "Necroptosis"):
            x_dict[i] = x_dict["MLKL"]
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5321867/
        if (i == "RIPK1_3"):
            x_dict[i] = x_dict["RIG1"] and not x_dict["CASP8"]
            # read more later
            # https://www.sciencedirect.com/science/article/pii/S1097276514008661
        if (i == "CASP8"):
            x_dict[i] = x_dict["FADD"] and x_dict["ROS"] and x_dict["DISC"] 
            #(x_dict["FADD"] and x_dict["Virus"]) or x_dict["ROS"] or x_dict["DISC"]
            # too many, will assume all is true, only certain drugs inhibit casp 8
        if (i == "ProCASP3_7"):
            x_dict[i] = x_dict["CASP8"]
            # https://pubmed.ncbi.nlm.nih.gov/9765224/
        if (i == "CASP3_7"):
            x_dict[i] = x_dict["ProCASP3_7"]
            # obvious (trivial") relationship based on naming
        if (i == "Apoptosis"):
            x_dict[i] = x_dict["CASP3_7"]
            # https://pubmed.ncbi.nlm.nih.gov/10200555/
    mat[j,:] = list(x_dict.values())
mat[n_iter, :] = np.average(mat[0:n_iter-1,:],axis=0)
print("Complete")
yticklabels = [str(x) for x in range(1,n_iter + 1)]
yticklabels.append("Average")
ax = sns.heatmap(mat, cmap = "viridis", linewidths = .05, xticklabels = components, yticklabels=yticklabels)
ax.tick_params(axis='y', which='major', labelsize= 10)

#diff = np.subtract(list(x_dict.values()),list(x_begin.values()))
#print("Nothing changed") if sum(diff) == 0 else print(diff)
print(mat[n_iter, :])

fig = plt.gcf()
fig.set_size_inches(10, 7, forward=True)
plt.tight_layout()
plt.show()