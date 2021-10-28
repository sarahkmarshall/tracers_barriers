# -*- coding: utf-8 -*-
"""
This is a script that takes the output from the models and produces the figures
for the paper.

"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy
import flopy.utils.binaryfile as bf
import pandas as pd
import os
import sys
from matplotlib import ticker
import matplotlib.cm as cm
from matplotlib import gridspec
from scipy.ndimage.filters import gaussian_filter


print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('pandas version: {}'.format(pd.__version__))
print('flopy version: {}'.format(flopy.__version__))

#==== DIRECTORIES =============================================================

modelName           = "TB_107_c" # 
modelName_mt        = "TB_107_c_MT3D" # 

proj_2_folder       = r'C:\workspace\Proj2_TracersBarriers/'
MT3D_USGS_folder    = r'C:\Program Files\mt3d-usgs_Distribution\bin\MT3D-USGS_64'
exe_name            = r'C:\workspace\Proj2_TracersBarriers\MT3D_Reaction_sngl' 
mt3d_version        = 'mt3d-usgs' 

this_model_folder = proj_2_folder + modelName

# Set working directory
os.chdir(this_model_folder)

print("Current working directory is: " + str(os.getcwd()))

figureDirectory = os.path.join(this_model_folder, "Figures")

figureDirectory2 = r"C:\SMarshall_PhD\Papers\2_Paper_Two\Figures\Maybe" 
if not os.path.exists(figureDirectory2):
    os.makedirs(figureDirectory2)

dataDirectory = os.path.join(this_model_folder, "Data") 

#==== MODEL PROPERTIES REQUIRED FOR PLOTS =====================================

nlay    = 12            # Number of layers.
Lx      = 10000.        # Length of the model sides in metres.
Ly      = 5000          # Rectangular domain.
   
delr    = 10.           # Spacing length across rows in metres.

delc    = delr          # The spacing along columns and rows are equal.
ncol    = int(Lx/delr)  # Number of columns
nrow    = int(Ly/delc)  # Number of rows
ncells  = nlay*ncol*nrow # Total number of cells in the model.
ncells_per_layer = ncol*nrow # Number of cells per layer.
surface_area = Lx*Ly # m^2
    
ztop = 300. # Top of model.
zbot = 0. # Base of model.
delv = (ztop - zbot) / nlay # Length of cells in z direction.
botm = np.linspace(ztop, zbot, nlay + 1) 

column_of_barrier = [(int(ncol/2)-2), (int(ncol/2)-1),
                     (int(ncol/2)), (int(ncol/2)+1), (int(ncol/2)+2)] 
                     
middle_of_barrier = ncol/2
vol_barr_fullypen = 50*4000*300  
vol_barr_buried   = 50*5000*150

sa_barr_fullypen = 50*4000 # m2
sa_barr_buried   = 50*5000 # m2

xsectn_area_fp   = 4000*300 # m2
xsectn_area_bu   = 5000*150 # m2

len_barr_fullypen = 4000
len_barr_buried   = 5000 

#==== HYDRAULIC PROPERTIES ====================================================

sHead = 290 # Starting head across aquifer.

sy = 0.1 ## Specific yield
ss = sy/ztop ## Specific storage

hk_aquifer = 1.  # Hydraulic conductvity along rows (m/day) It worked ar 0.09
vka = hk_aquifer/10. # Vertical hydraulic conductivity

prsity = sy

#==== DATABASE SET UP =========================================================

number_of_recharge_scenarios = 2
number_of_casestudies = 3

# ------------------------------
# DICTIONARIES
headDict = {}
ageDict = {}
massDict = {}

diffHeadAbsDict = {} # Absolute difference in head barrier/non-barrier
diffHeadRelDict = {} # Relative difference in head barrier/non-barrier

diffAgeAbsDict = {} # Absolute difference in age barrier/non-barrier
diffAgeRelDict = {} # Relative difference in age barrier/non-barrier

gammaHeadDict = {}
gammaAgeDict = {}

prctDiffHeadDict = {} # Percent change in head with barrier compared to without 
prctDiffAgeDict = {} # Percent change in age with barrier compared to without 

statsDict = {}

modelObjectDict = {}

cellBudgetDict = {}
frfDict = {} # Flow right face
fffDict = {} # Flow front face
flfDict = {} # Flow left face

nodeCoordsDict = {}

# ------------------------------

rechargeScenarioNames = [] # Create the keys for the dictionary
for i in range(number_of_recharge_scenarios):
    rechargeScenarioNames.append('RECHARGE_'+str(i+1))

caseStudyNames = []
for i in range(number_of_casestudies):
    caseStudyNames.append('CS_'+str(i+1))


# Making dictionary with head files for each recharge scenario and case study
for rechargeScenario in range(number_of_recharge_scenarios): # number_of_recharge_scenarios
    headDict[rechargeScenarioNames[rechargeScenario]] = {}
    ageDict[rechargeScenarioNames[rechargeScenario]] = {}
    massDict[rechargeScenarioNames[rechargeScenario]] = {}
    modelObjectDict[rechargeScenarioNames[rechargeScenario]] = {}
    cellBudgetDict[rechargeScenarioNames[rechargeScenario]] = {}
    frfDict[rechargeScenarioNames[rechargeScenario]] = {}
    fffDict[rechargeScenarioNames[rechargeScenario]] = {}
    flfDict[rechargeScenarioNames[rechargeScenario]] = {}
    nodeCoordsDict[rechargeScenarioNames[rechargeScenario]] = {}

    for caseStudy in range(number_of_casestudies): #number_of_casestudies

        modelname = modelName
        modelname_mt = modelName_mt
        modelname = (str(modelname) + "_R" + str(rechargeScenario) 
        + "_CS" + str(caseStudy))
        modelname_mt = (modelname_mt + "_R" + str(rechargeScenario) 
        + "_CS" + str(caseStudy))
       
        # HEADFILES - FROM FLOW MODEL
        headobj = bf.HeadFile(modelname + '.hds')
        head_times = headobj.get_times()
        head_raster = headobj.get_data(totim=head_times[0])
        headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = head_raster

        # CONCENTRATION/AGE FILES - FROM TRANSPORT MODEL
        directage_file = modelname_mt + str('.ucn') 
        ucnobjc_direct_age = flopy.utils.binaryfile.UcnFile(directage_file, precision="double")
        age_times = ucnobjc_direct_age.get_times()
        age_raster = ucnobjc_direct_age.get_data(totim=age_times[0])
        ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = age_raster
        
        # MASS BUDGET FILES - FROM TRANSPORT MODEL
        mass_file = modelname_mt + str('.mas') 
        mas = flopy.mt3d.Mt3dms.load_mas(mass_file)
        massDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = mas

        # MODEL OBJECTS (USE FOR X-SECTN PLOTS) - FROM FLOW MODEL
        ml = flopy.modflow.Modflow.load(str(modelname) + ".nam")
        modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = ml
        node_coordinates = ml.dis.get_node_coordinates()
        nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = node_coordinates

        # CELL BUDGET FILES WHICH I CAN USE TO DRAW VELOCITY ARROWS ON MY PLOTS
        cbb = bf.CellBudgetFile(modelname+'.cbc', precision='double')
        cellBudgetDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = cbb

        frf = cbb.get_data(text='FLOW RIGHT FACE', totim=head_times[0])[0]
        frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = frf
        fff = cbb.get_data(text='FLOW FRONT FACE', totim=head_times[0])[0]
        fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = fff
        flf = cbb.get_data(text='FLOW LOWER FACE', totim=head_times[0])[0]
        flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = flf


#------------------------------------------------------------------------------

# MAKING ARRAYS FOR AGE AND HEAD: DIFFERENCE BETWEEN CASE STUDIES
# Automate this so that it doesn't matter how many scenarios there are
# Add the data straight away to a dictionary.

for rechargeScenario in range(number_of_recharge_scenarios):
    diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]] = {}
    diffHeadRelDict[rechargeScenarioNames[rechargeScenario]] = {}
    diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]] = {}
    diffAgeRelDict[rechargeScenarioNames[rechargeScenario]] = {}
    
    for caseStudy in range(number_of_casestudies): 
        
        (diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]]   
         [caseStudyNames[caseStudy]]) = ((headDict[rechargeScenarioNames[rechargeScenario]]
         [caseStudyNames[caseStudy]]) - (headDict[rechargeScenarioNames[rechargeScenario]]
         [caseStudyNames[0]]))

        diffHeadRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = ((headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] - headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]])/headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]])

        (diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]]
         [caseStudyNames[caseStudy]]) = ((ageDict[rechargeScenarioNames[rechargeScenario]]
         [caseStudyNames[caseStudy]]) - (ageDict[rechargeScenarioNames[rechargeScenario]]
         [caseStudyNames[0]]))

        diffAgeRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = ((ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] - ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]])/ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]])
        
#------------------------------------------------------------------------------                          
#------------- RELATIVE/FRACTIONAL RESULTS ------------------------------------
#------------------------------------------------------------------------------

# Zinn and Konikow normalisation: call this "Gamma"

for rechargeScenario in range(number_of_recharge_scenarios):
    gammaAgeDict[rechargeScenarioNames[rechargeScenario]] = {}
    gammaHeadDict[rechargeScenarioNames[rechargeScenario]] = {}
    for caseStudy in range(number_of_casestudies): 
        gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = np.ones((nlay,nrow,ncol),dtype=np.float32)
        gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = np.ones((nlay,nrow,ncol),dtype=np.float32)
for rechargeScenario in range(number_of_recharge_scenarios):
    for caseStudy in range(number_of_casestudies): 
        for il in range(nlay): #nlay
            for ir in range(nrow): #nrow
                for ic in range(ncol): #ncol                 
                    age_with_barr = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][il][ir, ic]
                    age_no_barr = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]][il][ir, ic]                       
                    
                    head_with_barr = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][il][ir, ic]
                    head_no_barr = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]][il][ir, ic]                       
                    
                    if diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][il][ir, ic] >= 0:
                        gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][il][ir, ic] = (age_with_barr-age_no_barr)/age_no_barr
                        gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][il][ir, ic] = (head_with_barr-head_no_barr)/head_no_barr
                    else:
                        gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][il][ir, ic] = (age_with_barr-age_no_barr)/age_with_barr
                        gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][il][ir, ic] = (head_with_barr-head_no_barr)/head_with_barr

#------------------------------------------------------------------------------

# SPATIAL AREA FOR THE DIFFERENCE IN HEAD AND AGE
   # SPATIAL AREA FOR THE DIFFERENCE IN HEAD AND AGE
def myround(x, base=5):
    return int(base * np.ceil(x/base))
    
def myround_down(x, base=5):
    return int(base * np.floor(x/base))

total_no_cells = nlay*nrow*ncol
cell_volume = delr*delc*delv # in cubic metres

# ABSOLUTE AGE --> i.e. just difference, not greater than or less than
spatAbsAgeDict = {}

for rechargeScenario in range(number_of_recharge_scenarios): # number_of_recharge_scenarios
    spatAbsAgeDict[rechargeScenarioNames[rechargeScenario]] = {}
    for caseStudy in range(number_of_casestudies): # 
        if caseStudy == 0: # No barrier
            vol_barr = 0
            sa_barr = 0
            len_barr = 0
            x_sectn_area = 0
        elif caseStudy == 1: # Fully-penetrating
            vol_barr = vol_barr_fullypen
            sa_barr = sa_barr_fullypen
            len_barr = len_barr_fullypen
            x_sectn_area = xsectn_area_fp
        elif caseStudy == 2: # Buried
            vol_barr = vol_barr_buried
            sa_barr = sa_barr_buried
            len_barr = len_barr_buried
            x_sectn_area = xsectn_area_bu
        else:
            vol_barr = vol_barr_fullypen
            sa_barr = sa_barr_fullypen
            len_barr = len_barr_fullypen
            x_sectn_area = xsectn_area_fp
    
        print("Recharge: " + str(rechargeScenario) + ". CS: " + str(caseStudy))

        y = gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
        abs_y = np.absolute(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]])
        
        abs_max_y = np.amax(np.absolute(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]))
        print("The absolute maximum value of relative difference is: " + str(abs_max_y))
      
        max_to_nearest_5 = myround(abs_max_y*10000)
        max_to_nearest_5 = max_to_nearest_5/10000
        print("rounded up to: " + str(max_to_nearest_5))
        
        min_value = 0
        
        brackets = np.arange(min_value, max_to_nearest_5, 0.0005)
        
        number_cells_in_each_bracket = []
        area_of_aquifer_per_bracket = []
        volume_aquifer = []
        volume_water = []
        ratio_volumes = []
        ratio_surfacearea = []
        ratio_length = []
        ratio_xsectn = []
        if len(brackets) > 0:
            for i in brackets:
                print(i)
                a = (abs_y>i).sum()
                number_cells_in_each_bracket.append(a)
                area_of_aquifer_per_bracket.append((a/total_no_cells)*100)
                volume_aquifer.append(a*cell_volume) # in cubic metres
                volume_water.append((a*cell_volume)*prsity) # Where 0.1 is the porosity
                ratio_volumes.append((a*cell_volume)/vol_barr)
                ratio_surfacearea.append((a*cell_volume)/(sa_barr*1000)) # so it's in km
                ratio_length.append((a*cell_volume)/(len_barr*1000000)) # In km^2
                ratio_xsectn.append((a*cell_volume)/(x_sectn_area*1000)) # In km

            spatAbsAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = (pd.DataFrame({
                                            'percent_area_aquifer': area_of_aquifer_per_bracket, 
                                            'cells_aquifer': number_cells_in_each_bracket,
                                            'brackets_change_in_age': brackets,
                                            'volume_aquifer': volume_aquifer,
                                            'volume_water': volume_water,
                                            'ratio_volumes': ratio_volumes,
                                            'ratio_sa': ratio_surfacearea,
                                            'ratio_length': ratio_length,
                                            'ratio_xsectn': ratio_xsectn})) 
        
        else:
            pass
            
# ABSOLUTE HEAD --> i.e. just difference, not greater than or less than
spatAbsHeadDict = {}

for rechargeScenario in range(number_of_recharge_scenarios): # number_of_recharge_scenarios
    spatAbsHeadDict[rechargeScenarioNames[rechargeScenario]] = {}
    for caseStudy in range(number_of_casestudies): # 
        if caseStudy == 0:
            vol_barr = 0
            sa_barr = 0
            len_barr = 0
            x_sectn_area = 0
        elif caseStudy == 1:
            vol_barr = vol_barr_fullypen
            sa_barr = sa_barr_fullypen
            len_barr = len_barr_fullypen
            x_sectn_area = xsectn_area_fp
        elif caseStudy == 2:
            vol_barr = vol_barr_buried
            sa_barr = sa_barr_buried
            len_barr = len_barr_buried
            x_sectn_area = xsectn_area_bu
        else:
            vol_barr = vol_barr_fullypen
            sa_barr = sa_barr_fullypen
            len_barr = len_barr_fullypen
            x_sectn_area = xsectn_area_fp
        print("Recharge: " + str(rechargeScenario) + ". CS: " + str(caseStudy))

        y = diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
        abs_y = np.absolute(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]])
        
        abs_max_y = np.amax(np.absolute(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]))
        print("The absolute maximum value of relative difference is: " + str(abs_max_y))
      
        max_to_nearest_5 = myround(abs_max_y*1000)
        max_to_nearest_5 = max_to_nearest_5/1000
        print("rounded up to: " + str(max_to_nearest_5))
        
        min_value = 0
        
        brackets = np.arange(min_value, max_to_nearest_5, 0.005)
        
        number_cells_in_each_bracket = []
        area_of_aquifer_per_bracket = []
        volume_aquifer = []
        volume_water = []
        ratio_volumes = []
        ratio_surfacearea = []
        ratio_length = []
        ratio_xsectn = []
        if len(brackets) > 0:
            for i in brackets:
                print(i)
                a = (abs_y>i).sum()
                number_cells_in_each_bracket.append(a)
                area_of_aquifer_per_bracket.append((a/total_no_cells)*100)
                volume_aquifer.append(a*cell_volume) # in cubic metres
                volume_water.append((a*cell_volume)*0.1) # Where 0.1 is the porosity
                ratio_volumes.append((a*cell_volume)/vol_barr)
                ratio_surfacearea.append((a*cell_volume)/(sa_barr*1000)) # This is in km, standard in m
                ratio_length.append((a*cell_volume)/(len_barr*1000000)) # In km^2
                ratio_xsectn.append((a*cell_volume)/(x_sectn_area*1000)) # In km

            spatAbsHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = (pd.DataFrame({
                                            'percent_area_aquifer': area_of_aquifer_per_bracket, 
                                            'cells_aquifer': number_cells_in_each_bracket,
                                            'brackets_change_in_age': brackets,
                                            'volume_aquifer': volume_aquifer,
                                            'volume_water': volume_water,
                                            'ratio_volumes': ratio_volumes,
                                            'ratio_sa': ratio_surfacearea,
                                            'ratio_length': ratio_length,
                                            'ratio_xsectn': ratio_xsectn})) 
        else:
            pass
            
#------------------------------------------------------------------------------
    
# HEAD AND AGE CHANGE ACROSS BARRIER
# NO LONGER A GRADIENT - INSTEAD JUST A CHANGE 

age_type = ageDict # gammaAgeDict # ageDict

length_gradient = 100
row_start_for_grad = 0
row_end_for_grad = nrow

middle_of_barrier = ncol/2
left_grad_sample_col = middle_of_barrier - (length_gradient/delr)/2 
right_grad_sample_col = middle_of_barrier + (length_gradient/delr)/2

list_grad_h_dict = []
list_grad_a_dict = []

for sample_layer in range(nlay):
    headLeftDict = {}
    headRightDict = {}
    ageLeftDict = {}
    ageRightDict = {}
    headGradDict = {}
    ageGradDict = {}
    for rechargeScenario in range(number_of_recharge_scenarios): # number_of_recharge_scenarios
        headLeftDict[rechargeScenarioNames[rechargeScenario]] = {}
        headRightDict[rechargeScenarioNames[rechargeScenario]] = {}
        ageLeftDict[rechargeScenarioNames[rechargeScenario]] = {}
        ageRightDict[rechargeScenarioNames[rechargeScenario]] = {}
        headGradDict[rechargeScenarioNames[rechargeScenario]] = {}
        ageGradDict[rechargeScenarioNames[rechargeScenario]] = {}

        for caseStudy in range(number_of_casestudies):
            list_head_left = []
            list_head_right = []
            list_age_left = []
            list_age_right = []
            for i in range(row_start_for_grad, row_end_for_grad):
                list_head_left.append(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][sample_layer][i, int(left_grad_sample_col)])
                list_head_right.append(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][sample_layer][i, int(right_grad_sample_col)])
                list_age_left.append(age_type[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][sample_layer][i, int(left_grad_sample_col)])
                list_age_right.append(age_type[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][sample_layer][i, int(right_grad_sample_col)])                
                
            headLeftDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_head_left
            headRightDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_head_right
            ageLeftDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_age_left
            ageRightDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_age_right    

            list_age = []
            list_head = []
            for i in range(row_end_for_grad-row_start_for_grad):
                list_head.append((headRightDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i]-
                                 headLeftDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i]))
                list_age.append((ageRightDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i]-
                                 ageLeftDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i])/
                                ((ageRightDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i]+
                                 ageLeftDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i])/2))
            headGradDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_head
            ageGradDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_age

    list_grad_h_dict.append(headGradDict)
    list_grad_a_dict.append(ageGradDict)

#------------------------------------------------------------------------------
    
# AGE CHANGE ACROSS BARRIER --> NOT NORMALISED BY AVERAGE --> USE FOR DISCUSSION
# AND CALCULATION OF VELOCITIES

age_type = ageDict # gammaAgeDict # ageDict

list_grad_a_dict_nonorm = []

for sample_layer in range(nlay):
    ageLeftDict = {}
    ageRightDict = {}
    ageGradDict = {}
    for rechargeScenario in range(number_of_recharge_scenarios): # number_of_recharge_scenarios
        ageLeftDict[rechargeScenarioNames[rechargeScenario]] = {}
        ageRightDict[rechargeScenarioNames[rechargeScenario]] = {}
        ageGradDict[rechargeScenarioNames[rechargeScenario]] = {}

        for caseStudy in range(number_of_casestudies):
            list_age_left = []
            list_age_right = []
            for i in range(row_start_for_grad, row_end_for_grad):
                list_age_left.append(age_type[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][sample_layer][i, int(left_grad_sample_col)])
                list_age_right.append(age_type[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][sample_layer][i, int(right_grad_sample_col)])                
                
            ageLeftDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_age_left
            ageRightDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_age_right    

            list_age = []
            for i in range(row_end_for_grad-row_start_for_grad):
                
                list_age.append(ageRightDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i]-
                                 ageLeftDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][i])
            ageGradDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]] = list_age

    list_grad_a_dict_nonorm.append(ageGradDict)

# # # =========================================================================
        
# KEY STATISTICS

for rechargeScenario in range(number_of_recharge_scenarios):
    
    statsDict[rechargeScenarioNames[rechargeScenario]] = pd.DataFrame({

# Head
    'max_head': [
    np.amax(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

    'min_head': [
    np.amin(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amin(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amin(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

    'mean_head': [
    np.mean(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.mean(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.mean(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],
            
    'std_head': [
    np.std(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.std(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.std(headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

# Age
    'max_age': [
    np.amax(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

    'min_age': [
    np.amin(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amin(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amin(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

    'mean_age': [
    np.mean(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.mean(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.mean(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],
            
    'std_age': [
    np.std(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.std(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.std(ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

# Abs difference head                        
    'max_diff_head_abs': [
    np.amax(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])], 

    'mean_diff_head_abs': [
    np.mean(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.mean(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.mean(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],  

    'std_diff_head_abs': [
    np.std(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.std(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.std(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],    

    'min_diff_head_abs': [
    np.min(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.min(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.min(diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],                  

# Relative difference head                        
    'max_diff_head_rel': [
    np.amax(diffHeadRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(diffHeadRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(diffHeadRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],                       

# Abs difference age                        
    'max_diff_age_abs': [
    np.amax(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],        
    
    'min_diff_age_abs': [
    np.amin(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amin(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amin(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],   
            
    'mean_diff_age_abs': [
    np.mean(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.mean(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.mean(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],
            
    'std_diff_age_abs': [
    np.std(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.std(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.std(diffAgeAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

# Relative difference age                        
    'max_diff_age_rel': [
    np.amax(diffAgeRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(diffAgeRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(diffAgeRelDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

# Gamma age
    'max_diff_age_gamma': [
    np.amax(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],
            
    'mean_diff_age_gamma': [
    np.mean(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.mean(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.mean(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],
            
    'std_diff_age_gamma': [
    np.std(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.std(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.std(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

    'min_diff_age_gamma': [
    np.amin(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amin(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amin(gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

# Gamma head
    'max_diff_head_gamma': [
    np.amax(gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amax(gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amax(gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])],

    'min_diff_head_gamma': [
    np.amin(gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[0]]),
    np.amin(gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[1]]),
    np.amin(gammaHeadDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[2]])]
    }, index = ["CS_1", "CS_2", "CS_3"])
    
# Before I run 114 - Need to save certain data
spatAbsHeadDict_Standard = spatAbsHeadDict
spatAbsAgeDict_Standard = spatAbsAgeDict
    
#==============================================================================
# PLOT SETTINGS 
# Universal settings for the plot (width of axes, font size etc)

fontSize = 14

def cm2inch(value):
    return value/2.54
    
mpl.rc("axes", lw=0.5, edgecolor='grey', labelsize=fontSize)
mpl.rc('lines', lw=1, c='k')
#mpl.rc("grid", c='gray', ls='-', lw='0.5')
#plt.rc('xtick', direction='out', color='0.1')
#plt.rc('ytick', direction='out', color='0.1')

# A column in journal is ~8.2 cm wide

mpl.rc('font', size=fontSize)
mpl.rc('legend', labelspacing=0.1, frameon=False, fontsize=fontSize)
mpl.rc("figure", figsize=(cm2inch(8.25), cm2inch(8.25*1.25)), titlesize='medium')
#mpl.rc('font', **{'sans-serif' : 'Helvetica','family' : 'sans-serif'}) 
mpl.rc('font',family='Arial')
# # # Settings

colormap_options = ["gist_ncar", "RdGy_r", "gist_stern", "hsv", 
                    "nipy_spectral", 'gist_rainbow', "seismic",
                    "PiYG_r", "PuOr", "Spectral"]

contour_colors = "k"
velocity_arrows_color = "white"

head_cmap = 'coolwarm'
age_cmap = 'jet'

diff_head_cmap = colormap_options[7]
diff_age_cmap = colormap_options[6]

gamma_cmap = colormap_options[6] 

row_cross_section = int(nrow/2)

horizontal_colour = "darkcyan"

buried_colour = "firebrick"

### Velocity arrow parameters

iskipCols = 75
iskipLayers = 2
iskipRows = 10
layerDisMultiplyer = 20
colDisMultiplyer = .8
rowDisMultiplyer = .5

## Max and min head and age for plots --> because I want them all to be 
# directly comparable, i.e. to have exactly the same colourbars

age_min = 0
age_max = 21000

head_min = 290.0
head_max = 296.0

gamma_age_min = -1.6 
gamma_age_max = 1.6

# Contour labels 
gamma_cs_labels = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]
levels = np.arange(0, 20000, 5000) # Plotting contour set
levels_2 = np.arange(0, 20000, 2500) # Plotting contour set

# Some of the contours I want to have manual plotting, but this is annoying when I
# 'm playing around with them so I can change some of them here
manual_plot = True# True/False

dpi_value = 600

#==== PLOTTING ================================================================

#------------------------------------------------------------------------------
### FIGURE 3
#------------------------------------------------------------------------------

# These values are to display on the plots
mean_a = (statsDict[rechargeScenarioNames[0]].loc["CS_1", "mean_head"])
mean_b = int(statsDict[rechargeScenarioNames[0]].loc["CS_1", "mean_age"])
mean_c = (statsDict[rechargeScenarioNames[1]].loc["CS_1", "mean_head"])
mean_d = int(statsDict[rechargeScenarioNames[1]].loc["CS_1", "mean_age"])
std_a = (statsDict[rechargeScenarioNames[0]].loc["CS_1", "std_head"])
std_b = int(statsDict[rechargeScenarioNames[0]].loc["CS_1", "std_age"])
std_c = (statsDict[rechargeScenarioNames[1]].loc["CS_1", "std_head"])
std_d = int(statsDict[rechargeScenarioNames[1]].loc["CS_1", "std_age"])

# These values are for colour-bar evaluation purposes only
max_head_plot1 = [(statsDict[rechargeScenarioNames[0]].loc["CS_1", "max_head"]), 
                  (statsDict[rechargeScenarioNames[1]].loc["CS_1", "max_head"])]
max_age_plot1 = [int(statsDict[rechargeScenarioNames[0]].loc["CS_1", "max_age"]), 
                 int(statsDict[rechargeScenarioNames[1]].loc["CS_1", "max_age"])]   
               
min_head_plot1 = [(statsDict[rechargeScenarioNames[0]].loc["CS_1", "min_head"]), 
                  (statsDict[rechargeScenarioNames[1]].loc["CS_1", "min_head"])]
min_age_plot1 = [int(statsDict[rechargeScenarioNames[0]].loc["CS_1", "min_age"]), 
                 int(statsDict[rechargeScenarioNames[1]].loc["CS_1", "min_age"])]  

#==================================
fig = plt.figure(figsize = (10, 7))
#===================================

# # # Recharge Scenario 1: Head
caseStudy = 0
rechargeScenario = 0
ax1 = fig.add_subplot(4, 1, 1)

qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qz_avg = np.empty(flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qz_avg[1:, :, :] = 0.5 * (flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0:nlay-1, :, :] + 
                    flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][1:nlay, :, :])
qz_avg[0, :, :] = 0.5 * flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0, :, :]

ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img1 = ax1.imshow(np.flipud(raster), cmap=head_cmap, extent=extent, aspect=5, vmin = head_min, vmax = head_max) # 
axes1 = plt.gca()
axes1.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes1.get_xaxis().set_visible(False)


y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Z = np.meshgrid(x, z[:, 0, 0])

plt.quiver(X[::iskipLayers, ::iskipCols], Z[::iskipLayers, ::iskipCols],
           colDisMultiplyer*qx_avg[::iskipLayers, row_cross_section, ::iskipCols], 
            layerDisMultiplyer*(-qz_avg[::iskipLayers, row_cross_section, ::iskipCols]),
           color='w', scale=5, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# # # Recharge Scenario 1: Age

ax2 = fig.add_subplot(4, 1, 2)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img2 = ax2.imshow(np.flipud(raster), cmap=age_cmap, extent=extent, aspect=5, vmin = age_min, vmax = age_max) # 

cs = plt.contour(raster, levels=levels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes2 = plt.gca()
axes2.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes2.get_xaxis().set_visible(False)

# # # Recharge Scenario 2: Head

caseStudy = 0
rechargeScenario = 1
ax3 = fig.add_subplot(4, 1, 3)

qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qz_avg = np.empty(flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qz_avg[1:, :, :] = 0.5 * (flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0:nlay-1, :, :] + 
                    flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][1:nlay, :, :])
qz_avg[0, :, :] = 0.5 * flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0, :, :]

ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img3 = ax3.imshow(np.flipud(raster), cmap=head_cmap, extent=extent, aspect=5) # , vmin = head_min, vmax = head_max
axes3 = plt.gca()
axes3.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes3.get_xaxis().set_visible(False)

y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Z = np.meshgrid(x, z[:, 0, 0])

plt.quiver(X[::iskipLayers, ::iskipCols], Z[::iskipLayers, ::iskipCols],
           colDisMultiplyer*qx_avg[::iskipLayers, row_cross_section, ::iskipCols], 
            layerDisMultiplyer*(-qz_avg[::iskipLayers, row_cross_section, ::iskipCols]),
           color='w', scale=5, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# # # Recharge Scenario 2: Age

ax4 = fig.add_subplot(4, 1, 4)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img4 = ax4.imshow(np.flipud(raster), cmap=age_cmap, extent=extent, aspect=5, vmin = age_min, vmax = age_max)# 


#print(get ready for manual contour labels)
cs = plt.contour(raster, levels=levels_2, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, inline=1, manual = manual_plot, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes4 = plt.gca()
axes4.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
plt.xlabel('$x$ [m]')

# Adding colour-bars

plt.subplots_adjust(bottom=0.2)

cbaxes1 = fig.add_axes([0.133, 0.1, 0.35, 0.03]) 
cb1 = plt.colorbar(img1, cax = cbaxes1, orientation="horizontal")  
cb1.set_label('Hydraulic head [m]')
tick_locator = ticker.MaxNLocator(nbins = 5)
cb1.locator = tick_locator
cb1.update_ticks()

cbaxes2 = fig.add_axes([0.533, 0.1, 0.35, 0.03]) 
cb2 = plt.colorbar(img2, cax = cbaxes2, orientation="horizontal")  
cb2.set_label('Age [y]')
cb2.locator = tick_locator
cb2.update_ticks()

# Adding text to indicate which recharge scenario
#plt.gcf().text(0.92, (0.8), "Recharge 1", fontsize=14, rotation=90*3)        
#plt.gcf().text(0.92, (0.45), "Recharge 2", fontsize=14, rotation=90*3) 
bbox_props = {'facecolor':'white', 'alpha':0.75, 'lw':0}
axes1.text(7500, 200, (r"$\mu$ = %.0f, $\sigma$ = %.2f" % (mean_a, std_a)),
           bbox=bbox_props, zorder=15)
axes2.text(7500, 200, (r"$\mu$ = " + str(mean_b) + r", $\sigma$ = " +str(std_b)),
           bbox=bbox_props, zorder=15)
axes3.text(7500, 200, (r"$\mu$ = %.0f, $\sigma$ = %.2f" % (mean_c, std_c)),
           bbox=bbox_props, zorder=15)
axes4.text(7500, 200, (r"$\mu$ = " + str(mean_d) + r", $\sigma$ = " +str(std_d)),
           bbox=bbox_props, zorder=15)

# Plotting text on the plots
pos1 = axes1.get_position()
pos2 = axes2.get_position()
pos3 = axes3.get_position()
pos4 = axes4.get_position()

gap_amount = .146
# Use position values for the text, i.e. (x, (y + height), ...
plt.gcf().text(0.125, (0.7544685976449275+gap_amount), "(a) Hydraulic head, Uniform recharge")
plt.gcf().text(0.125, (0.5718599019927537+gap_amount), "(b) Groundwater age, Uniform recharge")
plt.gcf().text(0.125, (0.3892512063405797+gap_amount), "(c) Hydraulic head, Upgradient recharge")
plt.gcf().text(0.125, (0.20664251068840572+gap_amount), "(d) Groundwater age, Upgradient recharge")

# Saving the plot
name = "fig_3"
plt.savefig((figureDirectory + "\_" + str(name)), format="pdf", dpi = dpi_value)


#------------------------------------------------------------------------------
### FIGURE 4
#------------------------------------------------------------------------------

### CASE STUDY horizontal barrier TRIFECTA PLOT

mean_a = (statsDict[rechargeScenarioNames[0]].loc["CS_2", "mean_head"])
mean_b = int(statsDict[rechargeScenarioNames[0]].loc["CS_2", "mean_age"])
mean_c = statsDict[rechargeScenarioNames[0]].loc["CS_2", "mean_diff_age_gamma"]
mean_d = (statsDict[rechargeScenarioNames[1]].loc["CS_2", "mean_head"])
mean_e = int(statsDict[rechargeScenarioNames[1]].loc["CS_2", "mean_age"])
mean_f = statsDict[rechargeScenarioNames[1]].loc["CS_2", "mean_diff_age_gamma"]

std_a = (statsDict[rechargeScenarioNames[0]].loc["CS_2", "std_head"])
std_b = int(statsDict[rechargeScenarioNames[0]].loc["CS_2", "std_age"])
std_c = statsDict[rechargeScenarioNames[0]].loc["CS_2", "std_diff_age_gamma"]
std_d = (statsDict[rechargeScenarioNames[1]].loc["CS_2", "std_head"])
std_e = int(statsDict[rechargeScenarioNames[1]].loc["CS_2", "std_age"])
std_f = statsDict[rechargeScenarioNames[1]].loc["CS_2", "std_diff_age_gamma"]

# These values are for colour-bar evaluation purposes only
max_head_plot1 = [(statsDict[rechargeScenarioNames[0]].loc["CS_2", "max_head"]), 
                  (statsDict[rechargeScenarioNames[1]].loc["CS_2", "max_head"])]
max_age_plot1 = [int(statsDict[rechargeScenarioNames[0]].loc["CS_2", "max_age"]), 
                 int(statsDict[rechargeScenarioNames[1]].loc["CS_2", "max_age"])]   
               
min_head_plot1 = [(statsDict[rechargeScenarioNames[0]].loc["CS_2", "min_head"]), 
                  (statsDict[rechargeScenarioNames[1]].loc["CS_2", "min_head"])]
min_age_plot1 = [int(statsDict[rechargeScenarioNames[0]].loc["CS_2", "min_age"]), 
                 int(statsDict[rechargeScenarioNames[1]].loc["CS_2", "min_age"])]  
                     
#max_diff_age_plot = [np.amax(gammaAgeDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][:, row_cross_section, :]),
#       np.amax(gammaAgeDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][:, row_cross_section, :])]
#
#min_diff_age_plot = [np.amin(gammaAgeDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][:, row_cross_section, :]),
#       np.amin(gammaAgeDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][:, row_cross_section, :])]

#==================================
fig = plt.figure(figsize = (10, 12))
#===================================

# # # Recharge Scenario 1: Head
caseStudy = 1
rechargeScenario = 0
ax1 = fig.add_subplot(6, 1, 1)

qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qz_avg = np.empty(flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qz_avg[1:, :, :] = 0.5 * (flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0:nlay-1, :, :] + 
                    flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][1:nlay, :, :])
qz_avg[0, :, :] = 0.5 * flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0, :, :]

ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img1 = ax1.imshow(np.flipud(raster), cmap=head_cmap, extent=extent, aspect=5, vmin = head_min, vmax = head_max) # 
axes1 = plt.gca()
axes1.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes1.get_xaxis().set_visible(False)


y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Z = np.meshgrid(x, z[:, 0, 0])

plt.quiver(X[::iskipLayers, ::iskipCols], Z[::iskipLayers, ::iskipCols],
           colDisMultiplyer*qx_avg[::iskipLayers, row_cross_section, ::iskipCols], 
            layerDisMultiplyer*(-qz_avg[::iskipLayers, row_cross_section, ::iskipCols]),
           color='w', scale=6, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# # # Recharge Scenario 1: Age

ax2 = fig.add_subplot(6, 1, 2)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img2 = ax2.imshow(np.flipud(raster), cmap=age_cmap, extent=extent, aspect=5, vmin = age_min, vmax = age_max) # 

cs = plt.contour(raster, levels=levels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes2 = plt.gca()
axes2.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes2.get_xaxis().set_visible(False)

# # # Recharge Scenario 1: Difference in Age (gamma)

# ---
data_3 = gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
sigma = 4.0 # this depends on how noisy your data is, play with it!
filter_3 = gaussian_filter(data_3, sigma)
# ---

gamma_cs_labels_2 = np.linspace(-0.1, 0.1, 21)

ax3 = fig.add_subplot(6, 1, 3)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = filter_3
img3 = ax3.imshow(np.flipud(raster), cmap=diff_age_cmap, extent=extent, aspect=5, vmin = gamma_age_min, vmax = gamma_age_max) # 

cs = plt.contour(raster, levels=gamma_cs_labels_2, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, gamma_cs_labels_2, inline=1, fontsize=10, fmt='%1.2f', zorder=11, colors=contour_colors)

axes3 = plt.gca()
axes3.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes3.get_xaxis().set_visible(False)

# # # Recharge Scenario 2: Head

rechargeScenario = 1
ax4 = fig.add_subplot(6, 1, 4)

qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qz_avg = np.empty(flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qz_avg[1:, :, :] = 0.5 * (flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0:nlay-1, :, :] + 
                    flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][1:nlay, :, :])
qz_avg[0, :, :] = 0.5 * flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0, :, :]

ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img4 = ax4.imshow(np.flipud(raster), cmap=head_cmap, extent=extent, aspect=5, vmin = head_min, vmax = head_max) # 
axes4 = plt.gca()
axes4.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes4.get_xaxis().set_visible(False)

y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Z = np.meshgrid(x, z[:, 0, 0])

plt.quiver(X[::iskipLayers, ::iskipCols], Z[::iskipLayers, ::iskipCols],
           colDisMultiplyer*qx_avg[::iskipLayers, row_cross_section, ::iskipCols], 
            layerDisMultiplyer*(-qz_avg[::iskipLayers, row_cross_section, ::iskipCols]),
           color='w', scale=5, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# # # Recharge Scenario 2: Age

ax5 = fig.add_subplot(6, 1, 5)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img5 = ax5.imshow(np.flipud(raster), cmap=age_cmap, extent=extent, aspect=5, vmin = age_min, vmax = age_max) 

#print(get ready for manual contour labels)
cs = plt.contour(raster, levels=levels_2, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, inline=1, manual=manual_plot, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes5 = plt.gca()
axes5.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes5.get_xaxis().set_visible(False)

# # # Recharge Scenario 2: Gamma difference in age

ax6 = fig.add_subplot(6, 1, 6)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img6 = ax6.imshow(np.flipud(raster), cmap=diff_age_cmap, extent=extent, aspect=5, vmin = gamma_age_min, vmax = gamma_age_max) 

#print(get ready for manual contour labels)
cs = plt.contour(raster, levels=gamma_cs_labels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, gamma_cs_labels, inline=1, manual=manual_plot, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes6 = plt.gca()
axes6.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
plt.xlabel('$x$ [m]')

# Adding colour-bars

plt.subplots_adjust(bottom=0.2)

cbaxes1 = fig.add_axes([.1, 0.1, 0.24, 0.03]) 
cb1 = plt.colorbar(img1, cax = cbaxes1, orientation="horizontal")  
cb1.set_label('Hydraulic head [m]')
tick_locator = ticker.MaxNLocator(nbins = 5)
cb1.locator = tick_locator
cb1.update_ticks()

cbaxes2 = fig.add_axes([.38, 0.1, 0.24, 0.03]) 
cb2 = plt.colorbar(img2, cax = cbaxes2, orientation="horizontal")  
cb2.set_label('Age [y]')
cb2.locator = tick_locator
cb2.update_ticks()

#plt.rcParams['font.family'] = "sans-serif"
#plt.rcParams['font.sans-serif'] = "Helvetica"
#plt.rcParams['text.usetex'] = False

cbaxes3 = fig.add_axes([0.66, 0.1, 0.24, 0.03]) 
cb3 = plt.colorbar(img3, cax = cbaxes3, orientation="horizontal")  
cb3.set_label(r'$\Gamma\rm{^a}$ [-]') #r'$\Gamma^a$ [-]')  u'${\u0413}\\rm{^a}$ [-]'
cb3.locator = tick_locator
cb3.update_ticks()

# Adding text to indicate which recharge scenario
#plt.gcf().text(0.92, (0.8), "Recharge 1", fontsize=14, rotation=90*3)        
#plt.gcf().text(0.92, (0.45), "Recharge 2", fontsize=14, rotation=90*3) 
bbox_props = {'facecolor':'white', 'alpha':0.75, 'lw':0}
axes1.text(7500, 200, (r"$\mu$ = %.0f, $\sigma$ = %.2f" % (mean_a, std_a)),
           bbox=bbox_props, zorder=15)
axes2.text(7500, 200, (r"$\mu$ = " + str(mean_b) + r", $\sigma$ = " +str(std_b)),
           bbox=bbox_props, zorder=15)
axes3.text(7500, 200, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_c, std_c)),
           bbox=bbox_props, zorder=15)
axes4.text(7500, 200, (r"$\mu$ = %.0f, $\sigma$ = %.2f" % (mean_d, std_d)),
           bbox=bbox_props, zorder=15)
axes5.text(7500, 200, (r"$\mu$ = " + str(mean_e) + r", $\sigma$ = " +str(std_e)),
           bbox=bbox_props, zorder=15)
axes6.text(7500, 200, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_f, std_f)),
           bbox=bbox_props, zorder=15)

# Add text to the plot to label the plots
pos1 = axes1.get_position() # shows x, y, width, height
pos2 = axes2.get_position()
pos3 = axes3.get_position()
pos4 = axes4.get_position()
pos5 = axes5.get_position()
pos6 = axes6.get_position()

# Use position values for the text, i.e. (x, (y + height), ...
gap_amount = .086
plt.gcf().text(0.125, (0.8095511474164927+gap_amount), "(a) Hydraulic head, Uniform recharge")
plt.gcf().text(0.125, (0.6895511474164927+gap_amount), "(b) Groundwater age, Uniform recharge")
plt.gcf().text(0.125, (0.5695511474164927+gap_amount), "(c) Groundwater age difference, Uniform recharge")
plt.gcf().text(0.125, (0.4495511474164927+gap_amount), "(d) Hydraulic head, Upgradient recharge")
plt.gcf().text(0.125, (0.3295511474164927+gap_amount), "(e) Groundwater age, Upgradient recharge")
plt.gcf().text(0.125, (0.20955114741649272+gap_amount), "(f) Groundwater age difference, Upgradient recharge")

# Saving the plot
name = "fig_4"
plt.savefig((figureDirectory + "\_" + str(name)), format="pdf", dpi = dpi_value)

#------------------------------------------------------------------------------
### FIGURE 5
#------------------------------------------------------------------------------

## ABSOLUTE DIFFERENCE IN HEAD/AGE - Case Study 2
# Uniform recharge 
min_abs_head_diff = -1.2
max_abs_head_diff = 1.2
min_abs_age_diff = -100
max_abs_age_diff = 100

caseStudy = 1
rechargeScenario = 0
nLayer = 1 #2 Middle of aquifer

mean_a = (statsDict[rechargeScenarioNames[0]].loc["CS_2", "mean_diff_head_abs"])
mean_b = (statsDict[rechargeScenarioNames[0]].loc["CS_2", "mean_diff_age_gamma"])
std_a = (statsDict[rechargeScenarioNames[0]].loc["CS_2", "std_diff_head_abs"])
std_b = (statsDict[rechargeScenarioNames[0]].loc["CS_2", "std_diff_age_gamma"])

# These values are for colour-bar evaluation purposes only

max_head_plot = [np.amax(diffHeadAbsDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amax(diffHeadAbsDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]

max_age_plot = [np.amax(gammaAgeDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amax(gammaAgeDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]

min_head_plot = [np.amin(diffHeadAbsDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amin(diffHeadAbsDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]

min_age_plot = [np.amin(gammaAgeDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amin(gammaAgeDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]


qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qy_avg = np.empty(fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)

qy_avg[:, 1:, :] = 0.5 * (fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, 0:nrow-1, :] + 
                    fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, 1:nrow, :])

qy_avg[:, 0, :] = 0.5 * fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, 0, :]               
               
               
              
levels_1 = np.arange(-1.2, 1.2, 0.2)
#levels_2 = np.arange(-100, 100, 20)
diff_head_contours = np.linspace(-2, 2, 41)

#==================================
fig = plt.figure(figsize = (7, 7))
#==================================

ax1 = fig.add_subplot(2, 1, 1)
extent = (delr/2., Lx - delr/2., Ly - delc/2., delc/2.)
raster1 = diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]
img1 = ax1.imshow(np.flipud(raster1), cmap=diff_head_cmap, extent=extent, vmin = min_abs_head_diff, vmax = max_abs_head_diff) # 
axes1 = plt.gca()
plt.ylabel("$y$ [m]")
axes1.axes.get_xaxis().set_visible(False)
#plt.xlabel("x [m]")

#print(manual contours ahead)
cs = plt.contour(raster1, levels=diff_head_contours, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, diff_head_contours, inline=1, manual=manual_plot, fontsize=10, fmt='%1.1f', zorder=11, colors=contour_colors) #   manual = True, 

# -------- AGE PLOT --------- # 

# ---
data_2 = gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]
sigma = 4.0 # this depends on how noisy your data is, play with it!
filter_2 = gaussian_filter(data_2, sigma)
# ---

#gamma_cs_labels_fig5 = [-0.05, 0.0, 0.05]

ax2 = fig.add_subplot(2, 1, 2)
extent = (delr/2., Lx - delr/2., Ly - delc/2., delc/2.)
raster2 = filter_2
img2 = ax2.imshow(np.flipud(raster2), cmap=diff_age_cmap, extent=extent, vmin = -0.03, vmax =0.03) # , vmin = gamma_age_min, vmax = gamma_age_max
axes2 = plt.gca()
plt.xlabel("$x$ [m]")
plt.ylabel("$y$ [m]")

#print(warning manual contours ahead)
cs = plt.contour(raster2, levels=gamma_cs_labels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, gamma_cs_labels, manual=manual_plot, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors) #  manual=True, 

y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Y = np.meshgrid(x, y)

iskipRows = 37
iskipCols = 50

plt.quiver(X[::iskipRows, ::iskipCols], Y[::iskipRows, ::iskipCols],
           colDisMultiplyer*qx_avg[nLayer, ::iskipRows, ::iskipCols], 
            rowDisMultiplyer*(qy_avg[nLayer, ::iskipRows, ::iskipCols]),
           color='k', scale=5, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# ---------------------------- #
plt.subplots_adjust(bottom=0.21)

cbaxes1 = fig.add_axes([0.133, 0.1, 0.32, 0.03]) 
cb1 = plt.colorbar(img1, cax = cbaxes1, orientation="horizontal", format='%1.1f')  
cb1.set_label(r'$\Delta h$ [m]')
tick_locator = ticker.MaxNLocator(nbins = 5)
cb1.locator = tick_locator
cb1.update_ticks()

cbaxes2 = fig.add_axes([0.547, 0.1, 0.32, 0.03]) 
cb2 = plt.colorbar(img2, cax = cbaxes2, orientation="horizontal", format='%1.2f')  
cb2.set_label(r'$\Gamma\mathrm{^a}$ [-]')
cb2.locator = tick_locator
cb2.update_ticks()

axes1.text(500, 1000, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_a, std_a)),
           bbox=bbox_props, zorder=15)
axes2.text(500, 1000, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_b, std_b)),
           bbox=bbox_props, zorder=15)

# Add text to the plot to label the plots
pos1 = axes1.get_position() # shows x, y, width, height
pos2 = axes2.get_position()

# Use position values for the text, i.e. (x, (y + height), ...
gap_amount = 0.33
plt.gcf().text(0.19513686463836777, (0.5818181818181819+gap_amount), r"(a) Hydraulic head difference $\Delta h$")
plt.gcf().text(0.19513686463836777, (0.20000000000000018+gap_amount), r"(b) Normalised age difference $\Gamma\mathrm{^a}$")

# Saving the plot
name = "fig_5"
plt.savefig((figureDirectory + "\_" + str(name)), format='pdf', dpi = dpi_value)

#------------------------------------------------------------------------------
### FIGURE 6  
#------------------------------------------------------------------------------

## ABSOLUTE DIFFERENCE IN HEAD/AGE - Case Study 2
# Uniform recharge 
min_abs_head_diff = -2.1
max_abs_head_diff = 2.1


caseStudy = 1
rechargeScenario = 1
nLayer = 2 # Middle of aquifer

mean_a = (statsDict[rechargeScenarioNames[1]].loc["CS_2", "mean_diff_head_abs"])
mean_b = (statsDict[rechargeScenarioNames[1]].loc["CS_2", "mean_diff_age_gamma"])
std_a = (statsDict[rechargeScenarioNames[1]].loc["CS_2", "std_diff_head_abs"])
std_b = (statsDict[rechargeScenarioNames[1]].loc["CS_2", "std_diff_age_gamma"])

# These values are for colour-bar evaluation purposes only

max_head_plot = [np.amax(diffHeadAbsDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amax(diffHeadAbsDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]

max_age_plot = [np.amax(gammaAgeDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amax(gammaAgeDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]

min_head_plot = [np.amin(diffHeadAbsDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amin(diffHeadAbsDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]

min_age_plot = [np.amin(gammaAgeDict[rechargeScenarioNames[0]][caseStudyNames[caseStudy]][nLayer, :, :]),
       np.amin(gammaAgeDict[rechargeScenarioNames[1]][caseStudyNames[caseStudy]][nLayer, :, :])]

levels_1 = np.arange(-1.2, 1.2, 0.2)

#==================================
fig = plt.figure(figsize = (7, 7))
#==================================

diff_head_contours = np.linspace(-2, 2, 41)


qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qy_avg = np.empty(fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)

qy_avg[:, 1:, :] = 0.5 * (fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, 0:nrow-1, :] + 
                    fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, 1:nrow, :])

qy_avg[:, 0, :] = 0.5 * fffDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, 0, :]


ax1 = fig.add_subplot(2, 1, 1)
extent = (delr/2., Lx - delr/2., Ly - delc/2., delc/2.)
raster1 = diffHeadAbsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]
img1 = ax1.imshow(np.flipud(raster1), cmap=diff_head_cmap, extent=extent, vmin = min_abs_head_diff, vmax = max_abs_head_diff) # 
axes1 = plt.gca()
plt.ylabel("$y$ [m]")
axes1.axes.get_xaxis().set_visible(False)
#plt.xlabel("x [m]")

#print(manual contours ahead)
cs = plt.contour(raster1, levels=diff_head_contours, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, diff_head_contours, inline=1, manual=manual_plot, fontsize=10, fmt='%1.1f', zorder=11, colors=contour_colors) #   manual = True, 

# -------- AGE PLOT --------- # 

#gamma_cs_labels_fig5 = [-0.05, 0.0, 0.05]

ax2 = fig.add_subplot(2, 1, 2)
extent = (delr/2., Lx - delr/2., Ly - delc/2., delc/2.)
raster2 = gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]
img2 = ax2.imshow(np.flipud(raster2), cmap=diff_age_cmap, extent=extent, vmin = gamma_age_min, vmax = gamma_age_max) # 
axes2 = plt.gca()
plt.xlabel("$x$ [m]")
plt.ylabel("$y$ [m]")

cs = plt.contour(raster2, levels=gamma_cs_labels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, gamma_cs_labels, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors) #  manual=True, 

y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Y = np.meshgrid(x, y)

iskipRows = 37
iskipCols = 50

plt.quiver(X[::iskipRows, ::iskipCols], Y[::iskipRows, ::iskipCols],
           colDisMultiplyer*qx_avg[nLayer, ::iskipRows, ::iskipCols], 
            rowDisMultiplyer*(qy_avg[nLayer, ::iskipRows, ::iskipCols]),
           color='k', scale=5, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# ------------------------------ #
plt.subplots_adjust(bottom=0.21)

cbaxes1 = fig.add_axes([0.133, 0.1, 0.32, 0.03]) 
cb1 = plt.colorbar(img1, cax = cbaxes1, orientation="horizontal", format='%1.1f')  
cb1.set_label(r'$\Delta h$ [m]')
tick_locator = ticker.MaxNLocator(nbins = 5)
cb1.locator = tick_locator
cb1.update_ticks()

cbaxes2 = fig.add_axes([0.547, 0.1, 0.32, 0.03]) 
cb2 = plt.colorbar(img2, cax = cbaxes2, orientation="horizontal", format='%1.2f')  
cb2.set_label(r'$\Gamma\mathrm{^a}$ [-]')
cb2.locator = tick_locator
cb2.update_ticks()


axes1.text(500, 1000, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_a, std_a)),
           bbox=bbox_props, zorder=15)
axes2.text(500, 1000, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_b, std_b)),
           bbox=bbox_props, zorder=15)


# Add text to the plot to label the plots
pos1 = axes1.get_position() # shows x, y, width, height
pos2 = axes2.get_position()

# Use position values for the text, i.e. (x, (y + height), ...
gap_amount = 0.33
plt.gcf().text(0.19513686463836777, (0.5818181818181819+gap_amount), r"(a) Hydraulic head difference $\Delta h$")
plt.gcf().text(0.19513686463836777, (0.20000000000000018+gap_amount), r"(b) Normalised age difference $\Gamma\mathrm{^a}$")

# Saving the plot
name = "fig_6"
plt.savefig((figureDirectory + "\_" + str(name)), format="pdf", dpi = dpi_value)


#------------------------------------------------------------------------------
### FIGURE 7
#------------------------------------------------------------------------------

mean_a = (statsDict[rechargeScenarioNames[0]].loc["CS_3", "mean_head"])
mean_b = int(statsDict[rechargeScenarioNames[0]].loc["CS_3", "mean_age"])
mean_c = statsDict[rechargeScenarioNames[0]].loc["CS_3", "mean_diff_age_gamma"]
mean_d = (statsDict[rechargeScenarioNames[1]].loc["CS_3", "mean_head"])
mean_e = int(statsDict[rechargeScenarioNames[1]].loc["CS_3", "mean_age"])
mean_f = statsDict[rechargeScenarioNames[1]].loc["CS_3", "mean_diff_age_gamma"]

std_a = (statsDict[rechargeScenarioNames[0]].loc["CS_3", "std_head"])
std_b = int(statsDict[rechargeScenarioNames[0]].loc["CS_3", "std_age"])
std_c = statsDict[rechargeScenarioNames[0]].loc["CS_3", "std_diff_age_gamma"]
std_d = (statsDict[rechargeScenarioNames[1]].loc["CS_3", "std_head"])
std_e = int(statsDict[rechargeScenarioNames[1]].loc["CS_3", "std_age"])
std_f = statsDict[rechargeScenarioNames[1]].loc["CS_3", "std_diff_age_gamma"]

                    
#==================================
fig = plt.figure(figsize = (10, 12))
#===================================

# # # Recharge Scenario 1: Head
caseStudy = 2
rechargeScenario = 0
ax1 = fig.add_subplot(6, 1, 1)

qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qz_avg = np.empty(flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qz_avg[1:, :, :] = 0.5 * (flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0:nlay-1, :, :] + 
                    flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][1:nlay, :, :])
qz_avg[0, :, :] = 0.5 * flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0, :, :]

ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img1 = ax1.imshow(np.flipud(raster), cmap=head_cmap, extent=extent, aspect=5, vmin = head_min, vmax = head_max) # 
axes1 = plt.gca()
axes1.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes1.get_xaxis().set_visible(False)

y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Z = np.meshgrid(x, z[:, 0, 0])

plt.quiver(X[::iskipLayers, ::iskipCols], Z[::iskipLayers, ::iskipCols],
           colDisMultiplyer*qx_avg[::iskipLayers, row_cross_section, ::iskipCols], 
            layerDisMultiplyer*(-qz_avg[::iskipLayers, row_cross_section, ::iskipCols]),
           color='w', scale=5, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# # # Recharge Scenario 1: Age

ax2 = fig.add_subplot(6, 1, 2)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img2 = ax2.imshow(np.flipud(raster), cmap=age_cmap, extent=extent, aspect=5, vmin = age_min, vmax = age_max) # 

cs = plt.contour(raster, levels=levels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes2 = plt.gca()
axes2.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes2.get_xaxis().set_visible(False)

# # # Recharge Scenario 1: Difference in Age (gamma)

ax3 = fig.add_subplot(6, 1, 3)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img3 = ax3.imshow(np.flipud(raster), cmap=diff_age_cmap, extent=extent, aspect=5, vmin = gamma_age_min, vmax = gamma_age_max) # 

cs = plt.contour(raster, levels=gamma_cs_labels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, gamma_cs_labels, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes3 = plt.gca()
axes3.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes3.get_xaxis().set_visible(False)

# # # Recharge Scenario 2: Head

rechargeScenario = 1
ax4 = fig.add_subplot(6, 1, 4)

qx_avg = np.empty(frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qx_avg[:, :, 1:] = 0.5 * (frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0:ncol-1] + 
                            frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 1:ncol])
qx_avg[:, :, 0] = 0.5 * frfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, :, 0]

qz_avg = np.empty(flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].shape, 
                  dtype=flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]].dtype)
qz_avg[1:, :, :] = 0.5 * (flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0:nlay-1, :, :] + 
                    flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][1:nlay, :, :])
qz_avg[0, :, :] = 0.5 * flfDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][0, :, :]

ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = headDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img4 = ax4.imshow(np.flipud(raster), cmap=head_cmap, extent=extent, aspect=5, vmin = head_min, vmax = head_max) # 
axes4 = plt.gca()
axes4.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes4.get_xaxis().set_visible(False)

y, x, z = nodeCoordsDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
X, Z = np.meshgrid(x, z[:, 0, 0])

plt.quiver(X[::iskipLayers, ::iskipCols], Z[::iskipLayers, ::iskipCols],
           colDisMultiplyer*qx_avg[::iskipLayers, row_cross_section, ::iskipCols], 
            layerDisMultiplyer*(-qz_avg[::iskipLayers, row_cross_section, ::iskipCols]),
           color='w', scale=5, headwidth=3, headlength=2,
           headaxislength=2, width=0.0025)

# # # Recharge Scenario 2: Age

ax5 = fig.add_subplot(6, 1, 5)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = ageDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img5 = ax5.imshow(np.flipud(raster), cmap=age_cmap, extent=extent, aspect=5, vmin = age_min, vmax = age_max) 

#print(get ready for manual contour labels)
cs = plt.contour(raster, levels=levels_2, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, inline=1, manual=manual_plot, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes5 = plt.gca()
axes5.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
axes5.get_xaxis().set_visible(False)

# # # Recharge Scenario 2: Gamma difference in age

ax6 = fig.add_subplot(6, 1, 6)
ml = modelObjectDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
mm = flopy.plot.crosssection.ModelCrossSection(model=ml, line={'row':row_cross_section})                
extent = (delr/2., Lx - delr/2., ztop - delc/2., delc/2.)  # values from model
raster = gammaAgeDict[rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][:, row_cross_section, :]
img6 = ax6.imshow(np.flipud(raster), cmap=diff_age_cmap, extent=extent, aspect=5, vmin = gamma_age_min, vmax = gamma_age_max) 

cs = plt.contour(raster, levels=gamma_cs_labels, extent=extent, zorder=10, colors=contour_colors) 
plt.clabel(cs, gamma_cs_labels, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)

axes6 = plt.gca()
axes6.set_ylim([25, 275])
plt.ylabel('$z$ [m]')
plt.xlabel('$x$ [m]')

# Adding colour-bars

plt.subplots_adjust(bottom=0.2)

cbaxes1 = fig.add_axes([.1, 0.1, 0.24, 0.03]) 
cb1 = plt.colorbar(img1, cax = cbaxes1, orientation="horizontal")  
cb1.set_label('Hydraulic head [m]')
tick_locator = ticker.MaxNLocator(nbins = 5)
cb1.locator = tick_locator
cb1.update_ticks()

cbaxes2 = fig.add_axes([.38, 0.1, 0.24, 0.03]) 
cb2 = plt.colorbar(img2, cax = cbaxes2, orientation="horizontal")  
cb2.set_label('Age [y]')
cb2.locator = tick_locator
cb2.update_ticks()

cbaxes3 = fig.add_axes([0.66, 0.1, 0.24, 0.03]) 
cb3 = plt.colorbar(img3, cax = cbaxes3, orientation="horizontal")  
cb3.set_label(r'$\Gamma\mathrm{^a}$ [-]')
cb3.locator = tick_locator
cb3.update_ticks()

# Adding text to indicate which recharge scenario
#plt.gcf().text(0.92, (0.8), "Recharge 1", fontsize=14, rotation=90*3)        
#plt.gcf().text(0.92, (0.45), "Recharge 2", fontsize=14, rotation=90*3) 
bbox_props = {'facecolor':'white', 'alpha':0.75, 'lw':0}
axes1.text(7500, 200, (r"$\mu$ = %.0f, $\sigma$ = %.2f" % (mean_a, std_a)),
           bbox=bbox_props, zorder=15)
axes2.text(7500, 200, (r"$\mu$ = " + str(mean_b) + r", $\sigma$ = " +str(std_b)),
           bbox=bbox_props, zorder=15)
axes3.text(7500, 200, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_c, std_c)),
           bbox=bbox_props, zorder=15)
axes4.text(7500, 200, (r"$\mu$ = %.0f, $\sigma$ = %.2f" % (mean_d, std_d)),
           bbox=bbox_props, zorder=15)
axes5.text(7500, 200, (r"$\mu$ = " + str(mean_e) + r", $\sigma$ = " +str(std_e)),
           bbox=bbox_props, zorder=15)
axes6.text(7500, 200, (r"$\mu$ = %.2f, $\sigma$ = %.2f" % (mean_f, std_f)),
           bbox=bbox_props, zorder=15)

# Add text to the plot to label the plots
pos1 = axes1.get_position() # shows x, y, width, height
pos2 = axes2.get_position()
pos3 = axes3.get_position()
pos4 = axes4.get_position()
pos5 = axes5.get_position()
pos6 = axes6.get_position()

# Use position values for the text, i.e. (x, (y + height), ...
gap_amount = .086
plt.gcf().text(0.125, (0.8095511474164927+gap_amount), "(a) Hydraulic head, Uniform recharge")
plt.gcf().text(0.125, (0.6895511474164927+gap_amount), "(b) Groundwater age, Uniform recharge")
plt.gcf().text(0.125, (0.5695511474164927+gap_amount), "(c) Groundwater age difference, Uniform recharge")
plt.gcf().text(0.125, (0.4495511474164927+gap_amount), "(d) Hydraulic head, Upgradient recharge")
plt.gcf().text(0.125, (0.3295511474164927+gap_amount), "(e) Groundwater age, Upgradient recharge")
plt.gcf().text(0.125, (0.20955114741649272+gap_amount), "(f) Groundwater age difference, Upgradient recharge")

# Saving the plot
name = "fig_7"
plt.savefig((figureDirectory + "\_" + str(name)), format="pdf", dpi = dpi_value)


#------------------------------------------------------------------------------
### FIGURE 8 - normalised to xsection area of barrier
#------------------------------------------------------------------------------
lw_spatial_plot = 2

##################################
plt.figure(figsize=[6,8])
##################################


## PLOTS OF THE SPATIAL CHANGE IN HEAD THROUGH THE AQUIFER

plt.subplot(2, 1, 1)

plt.plot(spatAbsHeadDict[rechargeScenarioNames[1]][caseStudyNames[1]]["ratio_xsectn"],
         spatAbsHeadDict[rechargeScenarioNames[1]][caseStudyNames[1]]["brackets_change_in_age"],
         color = horizontal_colour, lw=lw_spatial_plot, label = "$Fully-penetrating$; $Upgradient$")
plt.plot(spatAbsHeadDict[rechargeScenarioNames[1]][caseStudyNames[2]]["ratio_xsectn"],
         spatAbsHeadDict[rechargeScenarioNames[1]][caseStudyNames[2]]["brackets_change_in_age"],
         color = buried_colour, lw=lw_spatial_plot, label = "$Buried$; $Upgradient$")

plt.plot(spatAbsHeadDict[rechargeScenarioNames[0]][caseStudyNames[1]]["ratio_xsectn"],
         spatAbsHeadDict[rechargeScenarioNames[0]][caseStudyNames[1]]["brackets_change_in_age"],
         ls="--", color = horizontal_colour, lw=lw_spatial_plot, label = "$Fully-penetrating$; $Uniform$")
plt.plot(spatAbsHeadDict[rechargeScenarioNames[0]][caseStudyNames[2]]["ratio_xsectn"],
         spatAbsHeadDict[rechargeScenarioNames[0]][caseStudyNames[2]]["brackets_change_in_age"], 
         ls="--", color = buried_colour, lw=lw_spatial_plot, label = "$Buried$; $Uniform$")

plt.ylabel(r"$|\Delta h|$ [m]")
axes = plt.gca()
axes.set_ylim([0.1, 2.5])
axes.set_xlim([0, 15])
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
plt.axvspan(0, 15, ymin=0, ymax=(0.4/(2.4)), color="0.5", alpha=0.3)

axes.set_xticklabels([])
axes.text(0.25, 2.2, ("(a) Hydraulic head difference"), zorder=15)

################################
# SPATIAL CHANGE IN AGE
plt.subplot(2, 1, 2)

plt.plot(spatAbsAgeDict[rechargeScenarioNames[1]][caseStudyNames[1]]["ratio_xsectn"],
         spatAbsAgeDict[rechargeScenarioNames[1]][caseStudyNames[1]]["brackets_change_in_age"],
         color = horizontal_colour, lw=lw_spatial_plot, label = "Fully-penetrating; Upgradient")
plt.plot(spatAbsAgeDict[rechargeScenarioNames[1]][caseStudyNames[2]]["ratio_xsectn"],
         spatAbsAgeDict[rechargeScenarioNames[1]][caseStudyNames[2]]["brackets_change_in_age"],
         color = buried_colour, lw=lw_spatial_plot, label = "Buried; Upgradient")

plt.plot(spatAbsAgeDict[rechargeScenarioNames[0]][caseStudyNames[1]]["ratio_xsectn"],
         spatAbsAgeDict[rechargeScenarioNames[0]][caseStudyNames[1]]["brackets_change_in_age"],
         ls="--", color = horizontal_colour, lw=lw_spatial_plot, label = "Fully-penetrating; Uniform")
plt.plot(spatAbsAgeDict[rechargeScenarioNames[0]][caseStudyNames[2]]["ratio_xsectn"],
         spatAbsAgeDict[rechargeScenarioNames[0]][caseStudyNames[2]]["brackets_change_in_age"], 
         ls="--", color = buried_colour, lw=lw_spatial_plot, label = "Buried; Uniform")


plt.ylabel(r"$|\Gamma\mathrm{^a}|$ [-]")

axes = plt.gca()
axes.set_yscale("log")
plt.xlabel("Effective distance from the barrier [km]")
axes.set_xlim([0, 15])
axes.set_ylim([.01, 5])
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
plt.axvspan(0, 15, ymin=0, ymax=((np.log(0.2/0.01))/np.log(5/0.01)), color="0.5", alpha=0.3)

axes.text(0.25, 2, ("(b) Normalised age difference"), zorder=15)

# LEGEND
mpl.rc('legend', labelspacing=0.1,  frameon=False, facecolor="white", edgecolor = "white")  # , framealpha=0.5
plt.subplots_adjust(bottom=0.2, wspace = 0.1, hspace = 0.1)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.02, -0.2))

name = "fig_8"
plt.savefig((figureDirectory + "\_" + str(name)), format="pdf", dpi = dpi_value)

# This is to use for figure 12 below
hk_1_head_x = spatAbsHeadDict[rechargeScenarioNames[1]][caseStudyNames[1]]["ratio_xsectn"]
hk_1_head_y = spatAbsHeadDict[rechargeScenarioNames[1]][caseStudyNames[1]]["brackets_change_in_age"]
hk_1_age_x = spatAbsAgeDict[rechargeScenarioNames[1]][caseStudyNames[1]]["ratio_xsectn"]
hk_1_age_y = spatAbsAgeDict[rechargeScenarioNames[1]][caseStudyNames[1]]["brackets_change_in_age"]

#------------------------------------------------------------------------------
### FIGURE 9
#------------------------------------------------------------------------------

# Setting up the colour settings

alpha_value = 1.
lw_value = 2
zero_lw = 5
c1 = "k" 
c2 = horizontal_colour 
c3 = buried_colour 
c4 = "k" 
c5 = horizontal_colour 
c6 = buried_colour 
layer_f9 = 2

# Setting up y axis

y_y = []
for i in range(row_start_for_grad, row_end_for_grad):
    y_y.append(i*delr)

x_solid_y = [0]*len(y_y)

# Setting  up data: HEAD, FULLY-PENETRATING (a)

# ----- Setting up the data for Recharge Scenario 1 (Unifrom recharge)

# Start with an empty array
x1_h1_sum_array = np.zeros(len(list_grad_h_dict[0][rechargeScenarioNames[0]][caseStudyNames[1]]))
x2_h1_sum_array = np.zeros(len(list_grad_h_dict[0][rechargeScenarioNames[0]][caseStudyNames[1]]))

# Now add the array for each layer so that you get a sum of all of the layers
for il in range(nlay):
    x1_h1 = []
    x2_h1 = [] 
    for ir in range(len(list_grad_h_dict[il][rechargeScenarioNames[0]][caseStudyNames[1]])):
        x1_h1.append((list_grad_h_dict[il][rechargeScenarioNames[0]][caseStudyNames[0]][ir]))        
        x2_h1.append((list_grad_h_dict[il][rechargeScenarioNames[0]][caseStudyNames[1]][ir]))
    x1_h1_sum_array = x1_h1_sum_array + x1_h1
    x2_h1_sum_array = x2_h1_sum_array + x2_h1

# Now take the mean of all layers
x1_h1_mean_array = x1_h1_sum_array/nlay
x2_h1_mean_array = x2_h1_sum_array/nlay
   
# ----- Setting up the data for Recharge Scenario 2 (Upgradient recharge)

# Start with an empty array 
x1_h2_sum_array = np.zeros(len(list_grad_h_dict[0][rechargeScenarioNames[1]][caseStudyNames[1]]))
x2_h2_sum_array = np.zeros(len(list_grad_h_dict[0][rechargeScenarioNames[1]][caseStudyNames[1]]))

# Now add the array for each layer so that you get a sum of all of the layers
for il in range(nlay):
    x1_h2 = []
    x2_h2 = [] 
    for ir in range(len(list_grad_h_dict[il][rechargeScenarioNames[1]][caseStudyNames[1]])):
        x1_h2.append((list_grad_h_dict[il][rechargeScenarioNames[1]][caseStudyNames[0]][ir]))  
        x2_h2.append((list_grad_h_dict[il][rechargeScenarioNames[1]][caseStudyNames[1]][ir]))
    x1_h2_sum_array = x1_h2_sum_array + x1_h2
    x2_h2_sum_array = x2_h2_sum_array + x2_h2

# Now take the mean of all layers
x1_h2_mean_array = x1_h2_sum_array/nlay
x2_h2_mean_array = x2_h2_sum_array/nlay
    
# Setting up data: AGE, FULLY-PENETRATING (b)

# ----- Setting up the data for Recharge Scenario 1 (Unifrom recharge)

# Start with an empty array
x1_a1_sum_array = np.zeros(len(list_grad_a_dict[0][rechargeScenarioNames[0]][caseStudyNames[1]]))
x2_a1_sum_array = np.zeros(len(list_grad_a_dict[0][rechargeScenarioNames[0]][caseStudyNames[1]]))

# Now add the array for each layer so that you get a sum of all of the layers
for il in range(nlay):
    x1_a1 = []
    x2_a1 = [] 
    for ir in range(len(list_grad_a_dict[il][rechargeScenarioNames[0]][caseStudyNames[1]])):
        x1_a1.append((list_grad_a_dict[il][rechargeScenarioNames[0]][caseStudyNames[0]][ir])) 
        x2_a1.append((list_grad_a_dict[il][rechargeScenarioNames[0]][caseStudyNames[1]][ir]))    
    x1_a1_sum_array = x1_a1_sum_array + x1_a1
    x2_a1_sum_array = x2_a1_sum_array + x2_a1   

# Now take the mean of all layers
x1_a1_mean_array = x1_a1_sum_array/nlay
x2_a1_mean_array = x2_a1_sum_array/nlay
   
# ----- Setting up the data for Recharge Scenario 2 (Upgradient recharge)

# Start with an empty array
x1_a2_sum_array = np.zeros(len(list_grad_a_dict[0][rechargeScenarioNames[1]][caseStudyNames[1]]))
x2_a2_sum_array = np.zeros(len(list_grad_a_dict[0][rechargeScenarioNames[1]][caseStudyNames[1]]))

# No norm(alisation) are for the discussion section on velocities (below)
x1_a2_sum_array_nonorm = np.zeros(len(list_grad_a_dict[0][rechargeScenarioNames[1]][caseStudyNames[1]]))
x2_a2_sum_array_nonorm = np.zeros(len(list_grad_a_dict[0][rechargeScenarioNames[1]][caseStudyNames[1]]))

# Now add the array for each layer so that you get a sum of all of the layers
for il in range(nlay):
    x1_a2 = []
    x2_a2 = [] 
    x1_a2_nonorm = []
    x2_a2_nonorm = [] 
    for ir in range(len(list_grad_a_dict[il][rechargeScenarioNames[1]][caseStudyNames[1]])):
        x1_a2.append((list_grad_a_dict[il][rechargeScenarioNames[1]][caseStudyNames[0]][ir]))  
        x2_a2.append((list_grad_a_dict[il][rechargeScenarioNames[1]][caseStudyNames[1]][ir]))  
        x1_a2_nonorm.append((list_grad_a_dict_nonorm[il][rechargeScenarioNames[1]][caseStudyNames[0]][ir]))  
        x2_a2_nonorm.append((list_grad_a_dict_nonorm[il][rechargeScenarioNames[1]][caseStudyNames[1]][ir]))  
    x1_a2_sum_array = x1_a2_sum_array + x1_a2
    x2_a2_sum_array = x2_a2_sum_array + x2_a2 
    x1_a2_sum_array_nonorm = x1_a2_sum_array_nonorm + x1_a2_nonorm
    x2_a2_sum_array_nonorm = x2_a2_sum_array_nonorm + x2_a2_nonorm   
    
# Now take the mean of all layers
x1_a2_mean_array = x1_a2_sum_array/nlay
x2_a2_mean_array = x2_a2_sum_array/nlay
x1_a2_mean_array_nonorm = x1_a2_sum_array_nonorm/nlay
x2_a2_mean_array_nonorm = x2_a2_sum_array_nonorm/nlay    
    
# Setting up HEAD data for BURIED BARRIER (c)

x_axes_r1_cs3_h = []
x_axes_r2_cs3_h = []
x_axes_r1_cs1_h = []
x_axes_r2_cs1_h = []

y_axes_depth = np.linspace(25/2, 300-25/2, 12)
#y_axes_depth = 300 - y_axes_depth


x_solid_z = [0]*len(y_axes_depth)

# These are a mean across the rows
for i in range(nlay):
    x_axes_r1_cs3_h.append((np.mean(list_grad_h_dict[i][rechargeScenarioNames[0]][caseStudyNames[2]]))) 
    x_axes_r2_cs3_h.append((np.mean(list_grad_h_dict[i][rechargeScenarioNames[1]][caseStudyNames[2]]))) 
    x_axes_r1_cs1_h.append((np.mean(list_grad_h_dict[i][rechargeScenarioNames[0]][caseStudyNames[0]]))) 
    x_axes_r2_cs1_h.append((np.mean(list_grad_h_dict[i][rechargeScenarioNames[1]][caseStudyNames[0]]))) 

# Setting up the AGE data for BURIED BARRIER (d)   

x_axes_r1_cs3_a = []
x_axes_r1_cs3_a_nonorm = []
x_axes_r2_cs3_a = []
x_axes_r2_cs3_a_nonorm = []
x_axes_r1_cs1_a = []
x_axes_r1_cs1_a_nonorm = []
x_axes_r2_cs1_a = []
x_axes_r2_cs1_a_nonorm = []

# These are a mean across the rows
for i in range(nlay):
    x_axes_r1_cs3_a.append(np.mean(list_grad_a_dict[i][rechargeScenarioNames[0]][caseStudyNames[2]])) 
    x_axes_r1_cs3_a_nonorm.append(np.mean(list_grad_a_dict_nonorm[i][rechargeScenarioNames[0]][caseStudyNames[2]])) 

    x_axes_r2_cs3_a.append(np.mean(list_grad_a_dict[i][rechargeScenarioNames[1]][caseStudyNames[2]])) 
    x_axes_r2_cs3_a_nonorm.append(np.mean(list_grad_a_dict_nonorm[i][rechargeScenarioNames[1]][caseStudyNames[2]])) 

    x_axes_r1_cs1_a.append(np.mean(list_grad_a_dict[i][rechargeScenarioNames[0]][caseStudyNames[0]])) 
    x_axes_r1_cs1_a_nonorm.append(np.mean(list_grad_a_dict_nonorm[i][rechargeScenarioNames[0]][caseStudyNames[0]])) 

    x_axes_r2_cs1_a.append(np.mean(list_grad_a_dict[i][rechargeScenarioNames[1]][caseStudyNames[0]])) 
    x_axes_r2_cs1_a_nonorm.append(np.mean(list_grad_a_dict_nonorm[i][rechargeScenarioNames[1]][caseStudyNames[0]])) 
    
# PLOTTING ####################################################################

#=========================
plt.figure(figsize=(9,10))
gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1]) 
#=========================
shade_colour = "yellow"
#---------------------------#--------------------------------------------------
# A
plt1 = plt.subplot(gs[0]) 
#---------------------------#--------------------------------------------------

# Show barrier location
plt.axvspan(-3., 0, ymin=(500/5000), ymax=(4500/5000), color =shade_colour, alpha = 0.2) 

# Plot Uniform recharge
l1 = plt.plot(x2_h1_mean_array, y_y, c=c2, lw=lw_value , alpha=alpha_value, ls="--", label = '_nolegend_')
l1_baseline = plt.plot(x1_h1_mean_array, y_y, c=c1, lw=lw_value , alpha=alpha_value, ls="--", label = '_nolegend_')
  
# Plot Upgradient recharge 

l2 = plt.plot(x2_h2_mean_array, y_y, c=c5, lw=lw_value ,ls="-", alpha=alpha_value, label = '_nolegend_')
l2_baseline = plt.plot(x1_h2_mean_array, y_y, c=c4, lw=lw_value ,ls="-", alpha=alpha_value, label = '_nolegend_')

axes = plt.gca()
axes.set_xlim([-3., 0])
plt.ylabel("$y$ [m]")

#axes.set_xticklabels([])
#plt.xlabel("Head gradient (m/m)")
axes.set_xticklabels([])

axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

axes.text(-2.8, 4700, (r"(a) Fully-penetrating, $\delta^{\mathrm{h}}_{y,\bar{z}}$"), zorder=15) 
axes.text(-2.8, 5050, ("Plan view"), zorder=15)

#---------------------------#--------------------------------------------------
# B
plt2 = plt.subplot(gs[1]) 
#---------------------------#--------------------------------------------------

# Show barrier location
plt.axvspan(-.5, .2, ymin=(500/5000), ymax=(4500/5000), color =shade_colour, alpha = 0.2) # PUT BACK

# Plot Uniform recharge
plt.plot(x2_a1_mean_array, y_y, c=c2, lw=lw_value, alpha=alpha_value, 
         ls="--", label = "Fully-penetrating, uniform recharge")
plt.plot(x1_a1_mean_array, y_y, c=c1, lw=lw_value, alpha=alpha_value, 
         ls="--", label = "Baseline, uniform recharge")
   
# Plot Upgradient recharge
plt.plot(x2_a2_mean_array, y_y, c=c5, lw=lw_value ,ls="-", alpha=alpha_value, 
         label = 'Fully-penetrating, Upgradient recharge')

plt.plot(x1_a2_mean_array, y_y, c=c4, lw=lw_value ,ls="-", alpha=alpha_value, 
         label = 'Baseline, Upgradient recharge')

axes = plt.gca()
axes.set_xlim([-.5, .2])
axes.set_ylim([0,5000])

axes.set_yticklabels([])
axes.set_xticklabels([])

axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

axes.text(-.47, 4700, (r"(b) Fully-penetrating,"), zorder=15)
axes.text(-.41, 4500, (r"$\delta^{\mathrm{a}}_{y,\bar{z}}$"), zorder=15)
#---------------------------#--------------------------------------------------
# C
plt3 = plt.subplot(gs[2]) 
#---------------------------#--------------------------------------------------

# Show barrier location
plt.axvspan(-3.0, 0, ymin=(0), ymax=(150/300), color =shade_colour, alpha = 0.2)

# Plot Uniform recharge
plt.plot(x_axes_r1_cs1_h, (300-y_axes_depth),lw=lw_value, 
         alpha=alpha_value,ls = "--", color = c1,label = '_nolegend_') 
plt.plot(x_axes_r2_cs1_h, (300-y_axes_depth), lw=lw_value, 
         alpha=alpha_value, ls = "-", color = c4, label = '_nolegend_') 

# Plot Upgradient recharge  
plt.plot(x_axes_r1_cs3_h, (300-y_axes_depth),lw=lw_value, 
         alpha=alpha_value,ls = "--", color = buried_colour,label = '_nolegend_') 
plt.plot(x_axes_r2_cs3_h, (300-y_axes_depth), lw=lw_value, 
         alpha=alpha_value, ls = "-", color = buried_colour, label = '_nolegend_')   
  
axes = plt.gca()
axes.invert_yaxis()
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

axes.set_xlim([-3., 0.0])
plt.xlabel("Head change across barrier [m]")
#every_nth = 3
#for n, label in enumerate(axes.xaxis.get_ticklabels()):
#    if n % every_nth != 0:
#        label.set_visible(False)

axes.set_ylim([25, 275])
plt.ylabel("$z$ [m]")

axes.text(-2.8, (300-75), (r"(c) Buried, $\delta^{\mathrm{h}}_{\bar{y},z}$"), zorder=15)
axes.text(-2.8, (300-19), ("Vertical section"), zorder=15)

#---------------------------#---------------------------
# D  
plt4 = plt.subplot(gs[3]) 
#---------------------------#---------------------------

# Show barrier location
plt.axvspan(-.5, .2, ymin=(0), ymax=(150/300), color =shade_colour, alpha = 0.2)

# Plot Uniform recharge
plt.plot(x_axes_r1_cs1_a, (300-y_axes_depth), 
         alpha=alpha_value,ls = "--", color = c1, lw=lw_value, label = "Baseline, Uniform") 
plt.plot(x_axes_r2_cs1_a, (300-y_axes_depth),
         alpha=alpha_value, ls = "-", color = c4, lw=lw_value, label = "Baseline, Upgradient")   

# Plot Upgradient recharge   
plt.plot(x_axes_r1_cs3_a, (300-y_axes_depth), 
         alpha=alpha_value,ls = "--", color = buried_colour, lw=lw_value, label = "Buried, Uniform") 
plt.plot(x_axes_r2_cs3_a, (300-y_axes_depth),
         alpha=alpha_value, ls = "-", color = buried_colour, lw=lw_value, label = "Buried, Upgradient")   

# For legend only #############################################
plt.plot([100, 101], [-10, -11], lw = lw_value, color = horizontal_colour, ls="--", alpha=alpha_value,
         label = "Fully-penetrating, Uniform")

plt.plot([100, 101], [-10, -11], lw = lw_value, color = horizontal_colour, ls="-", alpha=alpha_value,
         label = "Fully-penetrating, Upgradient")
#############################################

axes = plt.gca()
axes.invert_yaxis()
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

plt.xlabel("Normalised age change across barrier [-]")
axes.set_xlim([-.5, 0.2])
axes.set_yticklabels([])
axes.set_ylim([25, 275])

axes.text(-.47, (300-75), (r"(d) Buried, $\delta^{\mathrm{a}}_{\bar{y},z}$"), zorder=15)

# LEGEND #############################################
mpl.rc('legend', labelspacing=0.1,  frameon=False, facecolor="white", edgecolor = "white")  # , framealpha=0.5
plt.subplots_adjust(bottom=0.2, wspace=0.15, hspace=0.1)
legend = plt.legend(frameon=True, ncol=2, bbox_to_anchor=(1.03, -0.35))


# SAVING THE FIGURE # 
name = "_fig_9"
name2 = os.path.join(r"C:\SMarshall_PhD\Papers\2_Paper_Two\Figures", name)
plt.savefig(name2, format="pdf", dpi = dpi_value)

            
#------------------------------------------------------------------------------
### VALUES FOR DISCUSSION - on groundwater velocity
#------------------------------------------------------------------------------

### For Fully-penetrating barrier ---------------------------------------------

# Measurement at y = 2500 m  
y_measurement = 2500
cell_measurement_y = int(y_measurement/delc)

# Upgradient recharge --> head

fp_baseline_head_change_r2 = x1_h2_mean_array[cell_measurement_y]
fp_barrier_head_change_r2 = x2_h2_mean_array[cell_measurement_y]

fp_baseline_h_velocity_est_r2 = (-hk_aquifer*(fp_baseline_head_change_r2/length_gradient))/prsity # m/d
fp_barrier_h_velocity_est_r2 = (-hk_aquifer*(fp_barrier_head_change_r2/length_gradient))/prsity # m/d

# Upgradient recharge --> age

fp_baseline_age_change_r2 = x1_a2_mean_array_nonorm[cell_measurement_y]
fp_barrier_age_change_r2 = x2_a2_mean_array_nonorm[cell_measurement_y]

fp_baseline_a_velocity_est_r2 = length_gradient/(fp_baseline_age_change_r2*365) # m/d
fp_barrier_a_velocity_est_r2 = length_gradient/(fp_barrier_age_change_r2*365) # m/d

### For Buried barrier ---------------------------------------------

# Measurement at z = 212.5 m  
y_measurement = 212.5 # m
cell_measurement_z = 8

# Upgradient recharge --> head
b_baseline_head_change_r2 = x_axes_r2_cs1_h[cell_measurement_z]
b_barrier_head_change_r2 = x_axes_r2_cs3_h[cell_measurement_z]

b_baseline_h_velocity_est_r2 = (-hk_aquifer*(b_baseline_head_change_r2/length_gradient))/prsity # m/d
b_barrier_h_velocity_est_r2 = (-hk_aquifer*(b_barrier_head_change_r2/length_gradient))/prsity # m/d

# Upgradient recharge --> age

b_baseline_age_change_r2 = x_axes_r2_cs1_a_nonorm[cell_measurement_z]
b_barrier_age_change_r2 = x_axes_r2_cs3_a_nonorm[cell_measurement_z]

b_baseline_a_velocity_est_r2 = length_gradient/(b_baseline_age_change_r2*365) # m/d
b_barrier_a_velocity_est_r2 = length_gradient/(b_barrier_age_change_r2*365) # m/d

# Now what happens if we took the mean with depth, very much an approximation to 
# an open well (without considering effects of flux differences/mixing)

b_baseline_age_change_r2_m = np.mean(x_axes_r2_cs1_a_nonorm)
b_barrier_age_change_r2_m = np.mean(x_axes_r2_cs3_a_nonorm)

b_baseline_a_velocity_est_r2_m = length_gradient/(b_baseline_age_change_r2_m*365) # m/d
b_barrier_a_velocity_est_r2_m = length_gradient/(b_barrier_age_change_r2_m*365) # m/d

# Now other recharge scen to check # - - - - - - - - - - - - - - - - - - - -

# Uniform recharge --> head
b_baseline_head_change_r1 = x_axes_r1_cs1_h[cell_measurement_z]
b_barrier_head_change_r1 = x_axes_r1_cs3_h[cell_measurement_z]

b_baseline_h_velocity_est_r1 = (-hk_aquifer*(b_baseline_head_change_r1/length_gradient))/prsity # m/d
b_barrier_h_velocity_est_r1 = (-hk_aquifer*(b_barrier_head_change_r1/length_gradient))/prsity # m/d

# Uniform recharge --> age

b_baseline_age_change_r1 = x_axes_r1_cs1_a_nonorm[cell_measurement_z]
b_barrier_age_change_r1 = x_axes_r1_cs3_a_nonorm[cell_measurement_z]

b_baseline_a_velocity_est_r1 = length_gradient/(b_baseline_age_change_r1*365) # m/d
b_barrier_a_velocity_est_r1 = length_gradient/(b_barrier_age_change_r1*365) # m/d

# Now what happens if we took the mean with depth, very much an approximation to 
# an open well (without considering effects of flux differences/mixing)

b_baseline_age_change_r1_m = np.mean(x_axes_r1_cs1_a_nonorm)
b_barrier_age_change_r1_m = np.mean(x_axes_r1_cs3_a_nonorm)

b_baseline_a_velocity_est_r1_m = length_gradient/(b_baseline_age_change_r1_m*365) # m/d
b_barrier_a_velocity_est_r1_m = length_gradient/(b_barrier_age_change_r1_m*365) # m/d

#------------------------------------------------------------------------------
### FIGURE 10
#------------------------------------------------------------------------------

# Run 114_PP Then...

param_list = [0, 11.25, 22.5, 33.75, 45, 56.25, 67.5, 78.75, 90.]
colour_scale = np.linspace(0.0, 1.0, len(param_list), endpoint=False)
colour_scale=cm.inferno(np.linspace(0,1,len(param_list))) # cool_r

display_value = "ratio_xsectn" # "ratio_length"
xlim_values = [0.1, 15]
shade_colour = "0.5" 

############# HEAD: Plot R 1 and R 2 together on the same plot - but in four sub-plots
lw_spatial = 2
# ========================== 
plt.figure(figsize=(9,9))
# ==========================

ax = plt.subplot(2, 2, 1)


plt.plot(list_spatial_head[7][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_head[7][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[8]), ls="--",lw=lw_spatial, label = str(90-param_list[8])+ u'\N{DEGREE SIGN}')   
plt.plot(list_spatial_head[6][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_head[6][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[7]), ls="--",lw=lw_spatial, label = str(90-param_list[7])+ u'\N{DEGREE SIGN}')   
plt.plot(list_spatial_head[5][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_head[5][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[6]), ls="--",lw=lw_spatial, label = str(90-param_list[6])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_head[4][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_head[4][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[5]), ls="--",lw=lw_spatial, label = str(90-param_list[5])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_head[3][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_head[3][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[4]), ls="--",lw=lw_spatial, label = str(90-param_list[4])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_head[2][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_head[2][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[3]), ls="--",lw=lw_spatial, label = str(90-param_list[3])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_head[1][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_head[1][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[2]), ls="--",lw=lw_spatial, label = str(90-param_list[2])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_head[0][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
           list_spatial_head[0][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[1]), ls="--",lw=lw_spatial, label = str(90-param_list[1])+ u'\N{DEGREE SIGN}')
plt.plot(spatAbsHeadDict_Standard[rechargeScenarioNames[0]][caseStudyNames[1]][display_value],
         spatAbsHeadDict_Standard[rechargeScenarioNames[0]][caseStudyNames[1]]["brackets_change_in_age"],
         ls="--", color = (colour_scale[0]), lw=lw_spatial, label = str(90-param_list[0])+ u'\N{DEGREE SIGN}')

plt.ylabel(r"$|\Delta h|$ [m]")
#plt.legend(loc=1)
axes = plt.gca()
#axes.set_yscale("log")
axes.set_ylim([0.1, 2.5])
#    axes.set_xscale("log")
#    axes.set_ylim([0.005, 2])
axes.set_xticklabels([])
axes.set_xlim(xlim_values)
plt.axvspan(0, xlim_values[1], ymin=0, ymax=((0.5-0.1)/(2.5-0.1)), color=shade_colour, alpha=0.4)
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

axes.text(0.25, 2.25, ("(a) Hydraulic head difference"), zorder=15)

pos1 = axes.get_position()
axes.text(0.125, 2.55, "Uniform recharge")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    
ax = plt.subplot(2, 2, 2)

plt.plot(list_spatial_head[7][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_head[7][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[8]), lw=lw_spatial,ls = "-",  label = "_nolegend_")
plt.plot(list_spatial_head[6][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_head[6][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[7]), lw=lw_spatial, ls = "-", label = "_nolegend_")
plt.plot(list_spatial_head[5][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_head[5][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[6]), ls = "-", lw=lw_spatial, label = "_nolegend_")
plt.plot(list_spatial_head[4][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_head[4][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[5]), lw=lw_spatial, ls = "-",  label = "_nolegend_")
plt.plot(list_spatial_head[3][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_head[3][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[4]), lw=lw_spatial, ls = "-", label = "_nolegend_")
plt.plot(list_spatial_head[2][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_head[2][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[3]), ls = "-",  lw=lw_spatial, label = "_nolegend_")
plt.plot(list_spatial_head[1][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_head[1][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[2]), ls = "-", lw=lw_spatial, label = "_nolegend_")
plt.plot(list_spatial_head[0][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
           list_spatial_head[0][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[1]), ls = "-", lw=lw_spatial, label = "_nolegend_")
plt.plot(spatAbsHeadDict_Standard[rechargeScenarioNames[1]][caseStudyNames[1]][display_value],
         spatAbsHeadDict_Standard[rechargeScenarioNames[1]][caseStudyNames[1]]["brackets_change_in_age"],
         color = (colour_scale[0]), lw=lw_spatial, label = "_nolegend_")

axes = plt.gca()
#axes.set_yscale("log")
axes.set_ylim([0.1, 2.5])
axes.set_xlim(xlim_values)
plt.axvspan(0, xlim_values[1], ymin=0, ymax=((0.5-0.1)/(2.5-0.1)), color=shade_colour, alpha=0.4)
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

axes.set_xticklabels([])
axes.set_yticklabels([])
axes.text(.25, 2.25, ("(b) Hydraulic head difference"), zorder=15)

pos1 = axes.get_position()
axes.text(0.125, 2.55, "Upgradient recharge")

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    

ax = plt.subplot(2, 2, 3)

#plt.title("Age, R1 and R2")
plt.plot(list_spatial_age[7][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_age[7][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[8]), ls="--",lw=lw_spatial, label = str(param_list[8])) 
plt.plot(list_spatial_age[6][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_age[6][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[7]), ls="--",lw=lw_spatial, label = str(param_list[7]))
plt.plot(list_spatial_age[5][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_age[5][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[6]), ls="--",lw=lw_spatial, label = str(param_list[6]))
plt.plot(list_spatial_age[4][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_age[4][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[5]), ls="--",lw=lw_spatial, label = str(param_list[5]))
plt.plot(list_spatial_age[3][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_age[3][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[4]), ls="--",lw=lw_spatial, label = str(param_list[4]))
plt.plot(list_spatial_age[2][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_age[2][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[3]), ls="--",lw=lw_spatial, label = str(param_list[3]))
plt.plot(list_spatial_age[1][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
         list_spatial_age[1][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[2]), ls="--",lw=lw_spatial, label = str(param_list[2]))
plt.plot(list_spatial_age[0][rechargeScenarioNames[0]][caseStudyNames[3]][display_value],
           list_spatial_age[0][rechargeScenarioNames[0]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[1]), ls="--", lw=lw_spatial, label = str(param_list[1])) 
plt.plot(spatAbsAgeDict_Standard[rechargeScenarioNames[0]][caseStudyNames[1]][display_value],
         spatAbsAgeDict_Standard[rechargeScenarioNames[0]][caseStudyNames[1]]["brackets_change_in_age"],
         ls="--", color = (colour_scale[0]), lw=lw_spatial, label = str(param_list[0]))
   
plt.xlabel("Effective distance from the barrier [km]")
plt.ylabel("$|\Gamma^\mathrm{a}|$ [-]")

axes = plt.gca()
axes.set_yscale("log")
axes.set_xlim(xlim_values)
axes.set_ylim([.01, 20])
plt.axvspan(0, xlim_values[1], ymin=0, ymax=((np.log(0.2/0.01))/np.log(20/0.01)), color=shade_colour, alpha=0.4)
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

axes.text(.25, 9, ("(c) Normalised age difference"), zorder=15)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    

plt.subplot(2, 2, 4) 

#plt.figure()
plt.plot(list_spatial_age[7][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_age[7][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[8]), lw=lw_spatial,ls = "-",  label = str(90-param_list[8])+ u'\N{DEGREE SIGN}')  
plt.plot(list_spatial_age[6][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_age[6][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[7]), lw=lw_spatial, ls = "-", label = str(90-param_list[7])+ u'\N{DEGREE SIGN}') 
plt.plot(list_spatial_age[5][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_age[5][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[6]), ls = "-", lw=lw_spatial, label = str(90-param_list[6])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_age[4][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_age[4][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[5]), lw=lw_spatial, ls = "-",  label = str(90-param_list[5])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_age[3][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_age[3][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
        color = (colour_scale[4]), lw=lw_spatial, ls = "-", label = str(90-param_list[4])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_age[2][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_age[2][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[3]), ls = "-",  lw=lw_spatial, label = str(90-param_list[3])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_age[1][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
         list_spatial_age[1][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[2]), ls = "-", lw=lw_spatial, label = str(90-param_list[2])+ u'\N{DEGREE SIGN}')
plt.plot(list_spatial_age[0][rechargeScenarioNames[1]][caseStudyNames[3]][display_value],
           list_spatial_age[0][rechargeScenarioNames[1]][caseStudyNames[3]]["brackets_change_in_age"],
            color = (colour_scale[1]), ls = "-", lw=lw_spatial, label = str(90-param_list[1])+ u'\N{DEGREE SIGN}')
plt.plot(spatAbsAgeDict_Standard[rechargeScenarioNames[1]][caseStudyNames[1]][display_value],
         spatAbsAgeDict_Standard[rechargeScenarioNames[1]][caseStudyNames[1]]["brackets_change_in_age"],
         color = (colour_scale[0]), lw=lw_spatial, label = str(90.0-param_list[0])+ u'\N{DEGREE SIGN}')


plt.xlabel("Effective distance from the barrier [km]")
axes = plt.gca()
axes.set_yscale("log")
axes.set_xlim(xlim_values)
axes.set_ylim([.01, 20])
plt.axvspan(0, xlim_values[1], ymin=0, ymax=((np.log(0.2/0.01))/np.log(20/0.01)), color=shade_colour, alpha=0.4)
axes.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)

axes.set_yticklabels([])
axes.text(.25, 9, ("(d) Normalised age difference"), zorder=15)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =    

# LEGEND
mpl.rc('legend', labelspacing=0.1,  frameon=False, facecolor="white", edgecolor = "white")  # , framealpha=0.5
plt.subplots_adjust(bottom=0.2)
legend = plt.legend(frameon=True, bbox_to_anchor=(1.02, -0.2), ncol= 3)# , bbox_to_anchor=(0.5, 0., 0.5, 0.5))

plt.subplots_adjust(wspace=0.1, hspace = 0.1)

name = "fig_10"
plt.savefig((figureDirectory + "\_" + str(name)), format = "pdf", dpi = dpi_value)     

#------------------------------------------------------------------------------
### FIGURE 11
#------------------------------------------------------------------------------

# What is the max value if we don't consider the barrier cells?

print(need to run other model parts to get arrays)
## This is for filtering the data, below to calculate min and max, but I might not need to do this.

#Run 14_b --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14b = hk_array_4

#Run 14_c --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14c = hk_array_4

#Run 14_d --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14d = hk_array_4

#Run 14_e --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14e = hk_array_4

#Run 14_f --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14f = hk_array_4

#Run 14_g --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14g = hk_array_4

#Run 14_h --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14h = hk_array_4

#Run 14_i --> "BARRIER CHARACTERISTICS" & "HYDRAULIC CONDUCTIVITY ARRAYS"
hk_array_14i = hk_array_4

list_orientation_hks = [hk_array_14b, hk_array_14c, hk_array_14d, hk_array_14e, 
                        hk_array_14f, hk_array_14g, hk_array_14h, hk_array_14i]

list_max_values_orient = []
list_min_values_orient = []   
maxs_no_filter = []
mins_no_filter = []                     
                        
i = 0
for i in range(len(list_orientation_hks)):                   
    age_array_orient = all_models_gamma_age[i][rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]]
    hk_array_orient = list_orientation_hks[i]               
    list_max_values_orient.append(np.amax(age_array_orient[hk_array_orient>0.09]))
    list_min_values_orient.append(np.amin(age_array_orient[hk_array_orient>0.09]))
    
    maxs_no_filter.append(np.amax(age_array_orient))
    mins_no_filter.append(np.amin(age_array_orient))   

"Max overall with filter is: " + str(np.amax(list_max_values_orient))
# 'Max overall is: 1.5957
"Min overall with filter is: " + str(np.amin(list_min_values_orient))
# 'Min overall is: -0.676625'

"Max overall without filter is: " + str(np.amax(maxs_no_filter))
# 'Max overall is: 75.3168'
"Min overall without filter is: " + str(np.amin(mins_no_filter))
# 'Min overall is: -0.676625'

gamma_age_min = -1.6
gamma_age_max = 1.6

# -------------------- Plotting
cs_labels = [0]
levels = np.arange(-4, 4, 1) # Plotting contour set

colormap_options = ["gist_ncar", "RdGy_r", "gist_stern", "hsv", 
                    "nipy_spectral", 'gist_rainbow', "seismic",
                    "PiYG_r", "PuOr", "Spectral"]

contour_colors = "k"

gamma_cmap = colormap_options[6] 

bbox_props = {'facecolor':'white', 'alpha':0.0, 'lw':0}
subplot_label = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", "(h)"]
                 

## The h and i plot contours are SO messy, I need to filter the data a little bit to make them look better.
data_h = all_models_gamma_age[6][rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]
data_i = all_models_gamma_age[7][rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]

sigma = 4.0 # this depends on how noisy your data is, play with it!

filter_h = gaussian_filter(data_h, sigma)
filter_i = gaussian_filter(data_i, sigma)
##

# ==========================
fig = plt.figure(figsize=(10,12)) # figsize=(10,13)
# ==========================
for i in range(len(all_models_gamma_age)):    
    print(i)
    mean_value = list_statsDict[i][rechargeScenarioNames[1]].loc["CS_4", "mean_diff_age_gamma"]
    max_value = list_statsDict[i][rechargeScenarioNames[1]].loc["CS_4", "max_diff_age_gamma"]
    std_value = list_statsDict[i][rechargeScenarioNames[1]].loc["CS_4", "std_diff_age_gamma"]
    min_value = list_statsDict[i][rechargeScenarioNames[1]].loc["CS_4", "min_diff_age_gamma"]
    plt.subplot(4,2,(i+1))
    for rechargeScenario in range(1,2): # number_of_recharge_scenarios
        for caseStudy in range(3,4):
            for nLayer in range(whichLayerSave, (whichLayerSave+1)): # nlay

                if i == 6:
                    raster = filter_h                
#                elif i == 7:
#                    raster = filter_i
#
                else:
                    raster = all_models_gamma_age[i][rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]

    

                extent = (delr/2., Lx - delr/2., Ly - delc/2., delc/2.)
#                extent = (2000.0, 9995.0, 4995.0, 5.0)

                img = plt.imshow(np.flipud(raster), cmap=gamma_cmap, extent=extent, 
                                 vmin = gamma_age_min, vmax = gamma_age_max) #  cmap='RdGy_r'
#                cbr = plt.colorbar(fraction=0.016, pad=0.02)
#                cbr.set_label("Rel diff age (%)")
                plt.ylabel('$y$ [m]')
                plt.xlabel('$x$ [m]')
                axes = plt.gca()
                               
                axes.text(500, 1200, (r"$\mu$ = %.2f" % (mean_value)),
                bbox=bbox_props, zorder=15)
                axes.text(500, 1800, (r"$\sigma$ = %.2f" % (std_value)),
                bbox=bbox_props, zorder=15)
#                axes.text(500, 1500, (r"max = %.2f, min = %.2f" % (max_value, min_value)))
                axes.text(500, 600, (str(90-param_list[i+1])+ u'\N{DEGREE SIGN}'),
                bbox=bbox_props, zorder=15)
                if (i != 6 and i != 7):
                    axes.axes.get_xaxis().set_visible(False)
                else:
                    pass
                if (i == 1 or i == 3 or i == 5 or i == 7):
                    axes.axes.get_yaxis().set_visible(False)
                else:
                    pass
                cs = plt.contour(raster, antialiased = True, levels=[0], extent=extent, zorder=10, colors=contour_colors) 
                
                if i == 6:
                    plt.clabel(cs, cs_labels, manual = False, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)
                else:
                    plt.clabel(cs, cs_labels, manual = False, inline=1, fontsize=10, fmt='%1i', zorder=11, colors=contour_colors)
                axes.text(8900, 600, subplot_label[i])

plt.subplots_adjust(bottom=0.16) # plt.subplots_adjust(bottom=0.2)

cbaxes = fig.add_axes([0.2, 0.1, 0.6, 0.02])  # [left, bottom, width, height] ([0.2, 0.1, 0.6, 0.02])
cb = plt.colorbar(img, cax = cbaxes, orientation="horizontal")  
cb.set_label(r"Normalised age difference, $\Gamma^\mathrm{a}$ [-]")
#tick_locator = ticker.MaxNLocator(nbins = 5)
#cb.locator = tick_locator
#cb.update_ticks()

plt.subplots_adjust(wspace=0.05, hspace = 0.01)

figureDirectory = r"C:\SMarshall_PhD\Papers\2_Paper_Two\Figures" 
name = "fig_11"
plt.savefig((figureDirectory + "\_" + str(name)), format = "pdf",  dpi = dpi_value) # format="pdf",

i = 8
print(np.max(all_models_gamma_age[i][rechargeScenarioNames[rechargeScenario]][caseStudyNames[caseStudy]][nLayer, :, :]))
#------------------------------------------------------------------------------
### FIGURE 12 - head and age sensitivity analysis
#------------------------------------------------------------------------------
hk_1_head_x 
hk_1_head_y
hk_1_age_x
hk_1_age_y 
#hk_param_list = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 
#                 9.0, 10.0, 11.0, 12.0]
                 
hk_param_list = [2.0, 5.0, 10.0]
                                  
hk_param_list_string = [r"$1.0$", r"$2.0$", r"$3.0$", r"$4.0$", r"$5.0$", 
                        r"$6.0$", r"$7.0$", r"$8.0$", 
                        r"$9.0$", r"$10.0$", r"$11.0$", r"$12.0$"]
                        
hk_param_list_string = [r"$2.0$", r"$5.0$", r"$10.0$"]

              
ratio_param_list = [10, 100, 1000, 10000, 100000, 1000000]
              
ratio_param_list_string = [r"$10^1$", r"$10^2$", r"$10^3$", r"$10^4$", 
                           r"$10^5$", r"$10^6$"] 

#width_param_list = [50, 150, 250, 350, 450, 550, 650, 750, 850, 950, 1050,
#              1250, 1550, 1750, 2050, 3050, 4050, 5050]
              
#width_param_list = [50, 150, 250, 350, 450, 550, 650, 750, 850, 950, 1050]
width_param_list = [50, 250, 500, 750, 1000]

#width_param_list_string = [r"$50$", r"$150$", r"$250$", r"$350$", r"$450$", r"$550$", r"$650$", 
#                          r"$750$", r"$850$", r"$950$", r"$1050$"]
                          
width_param_list_string = [r"$50$", r"$250$", r"$500$", r"$750$", r"$1000$"]

porosity_param_list = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35]

aqthickness_param_list = [100, 300, 500, 700, 900, 1100]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

hk_spat_list = []
hk_spat_list_h = []

#dataDirectory = r'C:\workspace\Proj2_TracersBarriers\TB_18_a_3\Data'
dataDirectory = r'C:\workspace\Proj2_TracersBarriers\TB_118_a\Data'
model_run = 0
for model_run in range(len(hk_param_list)): 
    name = "spat_age_" + str(model_run)
    fileName = os.path.join(dataDirectory, (str(name) + '.csv'))
    spat_data = pd.read_csv(fileName)
    hk_spat_list.append(spat_data)

    name = "spat_head_" + str(model_run)
    fileName = os.path.join(dataDirectory, (str(name) + '.csv'))    
    spat_data = pd.read_csv(fileName)
    hk_spat_list_h.append(spat_data)

# - - - - - - - - - - - - - - - - - - - - 
ratio_spat_list = []
ratio_spat_list_h = []

#dataDirectory = r'C:\workspace\Proj2_TracersBarriers\TB_18_b_3\Data'
dataDirectory = r'C:\workspace\Proj2_TracersBarriers\TB_18_b_4\Data'

for model_run in range(1, (len(ratio_param_list)+1)): 
    name = "spat_age_" + str(model_run)
    fileName = os.path.join(dataDirectory, (str(name) + '.csv'))    
    spat_data = pd.read_csv(fileName)
    ratio_spat_list.append(spat_data)
    
    name = "spat_head_" + str(model_run)
    fileName = os.path.join(dataDirectory, (str(name) + '.csv'))    
    spat_data = pd.read_csv(fileName)
    ratio_spat_list_h.append(spat_data)

# - - - - - - - - - - - - - - - - - - - -     
width_spat_list = []
width_spat_list_h = []

#dataDirectory = r'C:\workspace\Proj2_TracersBarriers\TB_18_c_4\Data'
dataDirectory = r'C:\workspace\Proj2_TracersBarriers\TB_18_c_5\Data'

for model_run in range(len(width_param_list)): 
    name = "spat_age_" + str(model_run)
    fileName = os.path.join(dataDirectory, (str(name) + '.csv'))    
    spat_data = pd.read_csv(fileName)
    width_spat_list.append(spat_data)

    name = "spat_head_" + str(model_run)
    fileName = os.path.join(dataDirectory, (str(name) + '.csv'))    
    spat_data = pd.read_csv(fileName)
    width_spat_list_h.append(spat_data)

       
###########################
plt.figure(figsize=(15, 15))
###########################
legloc = "upper right"
x_limit = [0, 15]
value_used = "ratio_xsectn" # "ratio_length"
shade_colour = "0.5"
head_abc_loc = [.5, 2.55]
age_abc_loc= [.5, 57]
font_size_leg = 17
# -----------------------------------------------------------------------------
# # # HYDRAULIC CONDUCTIVITY OF AQUIFER 
# Head
plt.subplot(3, 2, 1)

color=cm.rainbow(np.linspace(0,1,len(hk_param_list)))

i = 0

plt.plot(hk_1_head_x, hk_1_head_y, lw=5, color="k", label=r"$1.0$")

for i in range(len(hk_param_list)):
    label = str(hk_param_list_string[i])

    plt.plot(hk_spat_list_h[i][value_used], # hk_spat_list_h[i]["ratio_length"]
               hk_spat_list_h[i]["brackets_change_in_age"],
                lw=5, color=color[i], label = label)
     
plt.ylabel("")
axes1 = plt.gca()
#axes1.set_yscale("log")
axes1.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
axes1.set_xticklabels([])
axes1.set_ylim([0.1, 2.5])
axes1.set_xlim(x_limit)

plt.axvspan(x_limit[0], x_limit[1], ymin=0, ymax=((0.5-0.1)/(2.5-0.1)), color=shade_colour, alpha=0.4)
axes1.text(head_abc_loc[0], head_abc_loc[1], ("(a)"), zorder=15)

#leg1 = axes1.legend(loc=legloc, ncol=2, fontsize = 8, framealpha = 0.7, 
#             title = "(a) $K_{aquifer}$ [m/d]", frameon=True)
#leg1.get_frame().set_linewidth(0.0)

###########################################################
# Age
plt.subplot(3, 2, 2)

color=cm.rainbow(np.linspace(0,1,len(hk_param_list)))

plt.plot(hk_1_age_x, hk_1_age_y, lw=5, color="k", label=r"$1.0$")

i = 2
for i in range(len(hk_param_list)):
    label = str(hk_param_list_string[i])

    plt.plot(hk_spat_list[i][value_used],
               hk_spat_list[i]["brackets_change_in_age"],
                lw=5, color=color[i], label = label)
    
plt.plot(hk_1_age_x, hk_1_age_y, lw=5, color="k", label='_nolegend_')

plt.ylabel("")
axes1 = plt.gca()
axes1.set_yscale("log")
axes1.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
axes1.set_xticklabels([])
axes1.set_ylim([0.01, 50])
axes1.set_xlim(x_limit)
axes1.yaxis.tick_right()

plt.axvspan(x_limit[0], x_limit[1], ymin=0, ymax=((np.log(0.2/0.01))/np.log(50/0.01)), color=shade_colour, alpha=0.4)
axes1.text(age_abc_loc[0], age_abc_loc[1], ("(b)"), zorder=15)

leg1 = axes1.legend(loc=legloc, ncol=2, fontsize = font_size_leg, framealpha = 0.7, 
             title = "$K_{\mathrm{a}}$ [m/d]", frameon=True)
leg1.get_frame().set_linewidth(0.0)
leg1.get_title().set_fontsize(font_size_leg)

# -----------------------------------------------------------------------------
# # # RATIO
# Head
plt.subplot(3, 2, 3)

color=cm.rainbow(np.linspace(0,1,len(ratio_param_list)))

plt.plot(hk_1_head_x, hk_1_head_y, lw=5, color="k", label=r"$10^3$")

for i in [0, 1, 3, 4, 5]:#range(len(ratio_param_list)):
    label = str(ratio_param_list_string[i])

    plt.plot(ratio_spat_list_h[i][value_used],
               ratio_spat_list_h[i]["brackets_change_in_age"],
                lw=5, color=color[i], label = label)
    
plt.plot(hk_1_head_x, hk_1_head_y, lw=5, color="k", label='_nolegend_')

     
plt.ylabel("Hydraulic head difference, $|\Delta^h|$ [m]", fontsize = font_size_leg)
axes2 = plt.gca()
#axes2.set_yscale("log")
axes2.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
axes2.set_xticklabels([])
axes2.set_ylim([0.1, 2.5])
axes2.set_xlim(x_limit)


plt.axvspan(x_limit[0], x_limit[1], ymin=0, ymax=((0.5-0.1)/(2.5-0.1)), color=shade_colour, alpha=0.4)
axes2.text(head_abc_loc[0], head_abc_loc[1], ("(c)"), zorder=15)


##############################################################
# Age
plt.subplot(3, 2, 4)

color=cm.rainbow(np.linspace(0,1,len(ratio_param_list)))


for i in [0, 1]:#range(len(ratio_param_list)):
    label = str(ratio_param_list_string[i])

    plt.plot(ratio_spat_list[i][value_used],
               ratio_spat_list[i]["brackets_change_in_age"],
                lw=5, color=color[i], label = label)
    
plt.plot(hk_1_age_x, hk_1_age_y, lw=5, color="k", label=r"$10^3$")

for i in [3, 4, 5]:#range(len(ratio_param_list)):
    label = str(ratio_param_list_string[i])

    plt.plot(ratio_spat_list[i][value_used],
               ratio_spat_list[i]["brackets_change_in_age"],
                lw=5, color=color[i], label = label)

plt.plot(hk_1_age_x, hk_1_age_y, lw=5, color="k", label="_nolegend_")
     
axes2 = plt.gca()
axes2.set_yscale("log")
axes2.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
axes2.set_xticklabels([])
axes2.set_ylim([0.01, 50])
axes2.set_xlim([0, 15])
#axes2.text(.25, 20, ("(b)"), zorder=15)
axes2.yaxis.tick_right()
plt.ylabel("Normalised age difference, $|\Gamma^a|$ [-]", fontsize = font_size_leg)
axes2.yaxis.set_label_position("right")


plt.axvspan(x_limit[0], x_limit[1], ymin=0, ymax=((np.log(0.2/0.01))/np.log(50/0.01)), color=shade_colour, alpha=0.4)
axes2.text(age_abc_loc[0], age_abc_loc[1], ("(d)"), zorder=15)

leg2 = axes2.legend(loc=legloc, ncol=2, fontsize = font_size_leg, framealpha = 0.7, 
             title = "$K_{\mathrm{a}}$ : $K_{\mathrm{b}}$ [-]", frameon=True)
leg2.get_frame().set_linewidth(0.0)
leg2.get_title().set_fontsize(font_size_leg)

# -----------------------------------------------------------------------------
# # # WIDTH
# Head
plt.subplot(3, 2, 5)

color=cm.rainbow(np.linspace(0,1,len(width_param_list)))

plt.plot(hk_1_head_x, hk_1_head_y, lw=5, color="k", label=r"$50$")

for i in range(1, len(width_param_list)):
    label = str(width_param_list_string[i])

    plt.plot(width_spat_list_h[i][value_used],
               width_spat_list_h[i]["brackets_change_in_age"],
                lw=5, color=color[i], label = label)

plt.plot(hk_1_head_x, hk_1_head_y, lw=5, color="k", label="_nolegend_")
     
plt.ylabel("")
plt.xlabel("Effective distance from barrier [km]", fontsize = font_size_leg)
axes3 = plt.gca()
#axes3.set_yscale("log")
axes3.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
axes3.set_ylim([0.1, 4.5])
axes3.set_xlim([0, 15])

plt.axvspan(x_limit[0], x_limit[1], ymin=0, ymax=((0.5-0.1)/(4.5-0.1)), color=shade_colour, alpha=0.4)
axes3.text(head_abc_loc[0], (head_abc_loc[1]+2.05), ("(e)"), zorder=15)

###############################################################
# Age
plt.subplot(3, 2, 6)

color=cm.rainbow(np.linspace(0,1,len(width_param_list)))

plt.plot(hk_1_age_x, hk_1_age_y, lw=5, color="k", label=r"$50$")

for i in range(1, len(width_param_list)):
    label = str(width_param_list_string[i])

    plt.plot(width_spat_list[i][value_used],
               width_spat_list[i]["brackets_change_in_age"],
                lw=5, color=color[i], label = label)

plt.plot(hk_1_age_x, hk_1_age_y, lw=5, color="k", label="_nolegend_")
     
plt.ylabel("")
plt.xlabel("Effective distance from barrier [km]", fontsize = font_size_leg)
axes3 = plt.gca()
axes3.set_yscale("log")
axes3.grid(color='grey', which ="both", linestyle='-', linewidth=0.5, alpha=0.5)
axes3.set_ylim([0.01, 50])
axes3.set_xlim(x_limit)
axes3.yaxis.tick_right()

plt.axvspan(x_limit[0], x_limit[1], ymin=0, ymax=((np.log(0.2/0.01))/np.log(50/0.01)), color=shade_colour, alpha=0.4)
axes3.text(age_abc_loc[0], age_abc_loc[1], ("(f)"), zorder=15)


leg3 = axes3.legend(loc=legloc, ncol=2, fontsize = font_size_leg, framealpha = 0.7, 
             title = "$W_\mathrm{b}$ [m]", frameon=True)
leg3.get_frame().set_linewidth(0.0)
leg3.get_title().set_fontsize(font_size_leg)

plt.subplots_adjust(hspace=0.1, wspace = 0.1)

# SAVING THE FIGURE # ---------------------------------------------------------
name = "_fig_12"
name2 = os.path.join(r"C:\SMarshall_PhD\Papers\2_Paper_Two\Figures", name)
plt.savefig(name2, format='pdf', dpi = dpi_value)






























