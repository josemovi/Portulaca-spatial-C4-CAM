

########################################################################
#     	   Supplementary Script for Moreno-Villena et al. (2021)       #
# This script used the newly-built flux balance models of C4-CAM to    #
# perform the analyses for the paper,                                  #
#                                                                      #
#########################################################################


from scipy import stats
from cobra.core import Metabolite, Reaction
from studyFunctions import *
import scobra

# import libraries
from libsbml import readSBML
from cobra import io, flux_analysis
import re
import cobra.test

from cobra.flux_analysis import flux_variability_analysis


################################ MAIN ######################################

# pFBA for the drought condition--"Drought"
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0


C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Drought.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")
  
# FVA for the drought condition--"Drought_FVA"
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=flux_variability_analysis(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Drought_FVA.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tminimum\tmaximum\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.minimum.get(rxn.id))+"\t"+str(C4CAMSolution.maximum.get(rxn.id))+"\n")

# pFBA for the well-watered condition--"Wet"
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Wet.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# FVA for the well-watered condition--"Wet_FVA"  
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=flux_variability_analysis(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Wet_FVA.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tminimum\tmaximum\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.minimum.get(rxn.id))+"\t"+str(C4CAMSolution.maximum.get(rxn.id))+"\n")

# pFBA for the C3 anatomy condition--"C3+CAM_ C3" 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -40
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 40
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -171.4
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 171.4
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = -30
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 30
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = -128.6
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 128.6

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("C3+CAM_C3.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for CAM condition under the C3 anatomy condition--"CAM_C3"
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:5.15})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:5.15})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:5.15})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:5.15})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -171.4
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 171.4
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = -128.6
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 128.6

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("CAM_C3.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for the CAM condition under the C4 anatomy condition--"CAM_C4" 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("CAM_C4.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for blocking malate transfer under the drought condition--"Drought_bMal" 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Drought_bMal.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for blocking malate transfer under the well-water condition--"Wet_bMal" 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_c_1_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Wet_bMal.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for blocking malate storage in mesophyll under the drought condition--"Drought_bM" 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_M").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Drought_bM.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")


# pFBA for blocking malate storage in bundle sheath under the drought condition--"Drought_bBS"
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_BS").upper_bound = 0


C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Drought_bBS.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for blocking malate storage in mesophyll and bundle sheath under the drought condition--"Drought_bMBS" 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_p2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("MAL_v_dielTransfer_M").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("Drought_bMBS.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for sensitivity analysis of Vc/Vo=10 in the bundle sheath--"VcVo10": 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:10})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("VcVo10.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")


# pFBA for sensitivity analysis of Vc/Vo=20 in the bundle sheath--"VcVo20": 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:20})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:20})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("VcVo20.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for sensitivity analysis of Vc/Vo=40 in the bundle sheath--"VcVo40": 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:40})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:40})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("VcVo40.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for sensitivity analysis of Vc/Vo=80 in the bundle sheath--"VcVo80": 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:80})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:80})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("VcVo80.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")

# pFBA for sensitivity analysis of Vc/Vo=200 in the bundle sheath--"VcVo200": 
C4CAM_model = io.read_sbml_model("C4_CAM.xml")
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_SUCROSE_e2_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e1_boundary_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("EX_GLC_e2_boundary_M").upper_bound = 0

Rubisco_balance_M = Metabolite("rubisco_bal_p1_M", name = "Weights to balance RuBP carboxygenase oxygenase balance_M", compartment = "p1_M")
C4CAM_model.reactions.get_by_id("RXN_961_p1_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_M").add_metabolites({Rubisco_balance_M:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_M").add_metabolites({Rubisco_balance_M:3})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_M").add_metabolites({Rubisco_balance_M:-1})

Rubisco_balance_BS = Metabolite("rubisco_bal_p1_BS", name = "Weights to balance RuBP carboxygenase oxygenase balance_BS", compartment = "p1_BS")
C4CAM_model.reactions.get_by_id("RXN_961_p1_BS").add_metabolites({Rubisco_balance_BS:200})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p1_BS").add_metabolites({Rubisco_balance_BS:-1})
C4CAM_model.reactions.get_by_id("RXN_961_p2_BS").add_metabolites({Rubisco_balance_BS:200})
C4CAM_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p2_BS").add_metabolites({Rubisco_balance_BS:-1})

Maintenance_constraint11 = Metabolite("ATPase_NADPHoxidase_constraint_c1_M",name =  "ATPase_NADPHoxidase_constraint_c1_M", compartment = "c1")
Maintenance_constraint12 = Metabolite("ATPase_NADPHoxidase_constraint_c2_M",name =  "ATPase_NADPHoxidase_constraint_c2_M", compartment = "c2")
Maintenance_constraint13 = Metabolite("Light_dark_maintainence_constraint_M",name =  "Light_dark_maintainence_constraint_M", compartment = "c1")

Maintenance_constraint21 = Metabolite("ATPase_NADPHoxidase_constraint_c1_BS",name =  "ATPase_NADPHoxidase_constraint_c1_BS", compartment = "c1")
Maintenance_constraint22 = Metabolite("ATPase_NADPHoxidase_constraint_c2_BS",name =  "ATPase_NADPHoxidase_constraint_c2_BS", compartment = "c2")
Maintenance_constraint23 = Metabolite("Light_dark_maintainence_constraint_BS",name =  "Light_dark_maintainence_constraint_BS", compartment = "c1")

C4CAM_model.reactions.get_by_id("ATPase_tx1_M").add_metabolites({Maintenance_constraint11:1,Maintenance_constraint13:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").add_metabolites({Maintenance_constraint12:1,Maintenance_constraint13:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_M").add_metabolites({Maintenance_constraint12:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_M").add_metabolites({Maintenance_constraint11:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_M").add_metabolites({Maintenance_constraint12:-3})     

C4CAM_model.reactions.get_by_id("ATPase_tx1_BS").add_metabolites({Maintenance_constraint21:1,Maintenance_constraint23:1})
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").add_metabolites({Maintenance_constraint22:1,Maintenance_constraint23:-1})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxc_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxm_tx2_BS").add_metabolites({Maintenance_constraint22:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx1_BS").add_metabolites({Maintenance_constraint21:-3})
C4CAM_model.reactions.get_by_id("NADPHoxp_tx2_BS").add_metabolites({Maintenance_constraint22:-3})

rxn = C4CAM_model.reactions.get_by_id("diel_biomass_BS")
rxn.add_metabolites({C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t1_t_BS"):1.75, C4CAM_model.metabolites.get_by_id("X_Phloem_contribution_t2_t_BS"):-1.25})
 
C4CAM_model.reactions.get_by_id("Photon_tx2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("Photon_tx1_M").lower_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_M").upper_bound = 1120
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").lower_bound = 480
C4CAM_model.reactions.get_by_id("Photon_tx1_BS").upper_bound = 480

C4CAM_model.reactions.get_by_id("CO2_tx1_M").lower_bound = -70
C4CAM_model.reactions.get_by_id("CO2_tx1_M").upper_bound = 70
C4CAM_model.reactions.get_by_id("CO2_tx2_M").lower_bound = -300
C4CAM_model.reactions.get_by_id("CO2_tx2_M").upper_bound = 300
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("CO2_tx2_BS").upper_bound = 0

C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_BS").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c1_M").upper_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").lower_bound = 0
C4CAM_model.reactions.get_by_id("PEPCARBOXYKIN_RXN_c2_M").upper_bound = 0

C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("CARBON_DIOXIDE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").lower_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_BS").upper_bound=23.73
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").lower_bound=10.17
C4CAM_model.reactions.get_by_id("ATPase_tx2_M").upper_bound=10.17

C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("PHOSPHO_ENOL_PYRUVATE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("2_PG_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("GAP_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("OXYGEN_MOLECULE_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_1_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("WATER_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("DIHYDROXY_ACETONE_PHOSPHATE_c_2_BSMTransfer").upper_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").lower_bound = 0
C4CAM_model.reactions.get_by_id("G3P_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").lower_bound = -1000
C4CAM_model.reactions.get_by_id("MAL_c_2_BSMTransfer").upper_bound = 0

C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").lower_bound = 0
C4CAM_model.reactions.get_by_id("PYRUVATEORTHOPHOSPHATE_DIKINASE_RXN_c2_BS").upper_bound = 0

C4CAMSolution=cobra.flux_analysis.pfba(C4CAM_model)

fin = open("rxnpathwaydict.csv","r")
pathclass = dict()
pathname=dict()
for line in fin:
  line=line.replace("\n","")
  lineparts=line.split(",")
  pathname[lineparts[0]]=lineparts[1]


fin = open("pathwayprioritizer.csv","r")
pathindex = dict()

for line in fin:
  line=line.replace("\n","")
  lineparts=line.split("\t")
  pathindex[lineparts[0]]=lineparts[1]


fout = open("VcVo200.csv","w")
fout.write("Reaction ID\tReaction Name\tEquation\tCompartment\tPathway\tPathway index\tE.C number\tC3 optimal\tC3 starch-hydrolyzing & maltose-exporting\tC3 starch-hydrolyzing & glucose-exporting\tC3 no hydrolysis\tCAM optimal\tCAM starch-hydrolyzing & maltose-exporting\tCAM starch-hydrolyzing & glucose-exporting\tCAM no hydrolysis\n")
for rxn in C4CAM_model.reactions:
  compSet=set()
  for met in rxn.metabolites.keys():
    compSet.add(met.compartment)
  comp="";
  for C in compSet:
    comp=comp+C
  path=""
  index=""
  EC=""
  RXN=rxn.id
  if(rxn.id[len(rxn.id)-1]=="M"):
    RXN = rxn.id[0:len(rxn.id)-3]
  if(rxn.id[len(rxn.id)-1]=="S"):
    RXN = rxn.id[0:len(rxn.id)-4]
  if(pathname.keys().__contains__(RXN)):
    path=pathname.get(RXN)
  if pathindex.keys().__contains__(path):
    index=pathindex.get(path)
  if(not str(rxn.notes)=="{}"):
    EC=rxn.notes.get("PROTEIN CLASS")[0]
  fout.write(rxn.id+"\t"+rxn.name+"\t"+rxn.reaction+"\t"+comp+"\t"+path+"\t"+index+"\t"+EC+"\t"+str(C4CAMSolution.fluxes.get(rxn.id))+"\n")
