from SimulationModel.DataImportFunctions import *
from SimulationModel.HelperFunctions import *
import numpy as np
import pandas as pd

# Path Variables
PATH    = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\Leading_EdgeGrid.gdb\Leading_EdgeGrids.xlsx"
PATH1   = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_16green.xls"
PATH2   = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_17green.xls"
PATH3   = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_17greenPoissonCount.xls"
PATH4   = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\Hinton_Climate_Data\GrowthRates_All.xlsx"
PATH5   = r'C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\Hinton_Climate_Data\StochasticGrowthRates\Iters3_Test.txt'
# Path variables for the juvenile stem data
PATH6   = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\LeadingEdge_Juveniles.gdb\OAge_grid.xlsx"
PATH7   = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\LeadingEdge_Juveniles.gdb\UAge_grid.xlsx"
PATH8   = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\LeadingEdge_Juveniles.gdb\OStem_grid.xlsx"
PATH9   = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\LeadingEdge_Juveniles.gdb\UStem_grid.xlsx"
PATH10  = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\LeadingEdge_Juveniles.gdb\OJuvi_grid.xlsx"
PATH11  = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\LeadingEdge_Juveniles.gdb\UJuvi_grid.xlsx"
PATH12  = r'C:\Users\djgra\Jupyter_Notebooks\HintonRelated\GDD5_data.txt'
######################################################################################################################################################
### Mutable Parameters ###############################################################################################################################
######################################################################################################################################################

# Simulation Parameters
Iterations = 200

# Model Parameters
c0, c1 = 0.0002, 3
beta = np.array(10.8 * 10**-6, dtype=np.float)
rho = 5000
price, tau = 177.54 , 0.05
f_healthy, f_dead = 6.58, 0.95
f = 5.63
c = 101.00
discount = 0.015
SI_26 = 20
SI_45 = 22
SI_85 = 26
LC = np.arange(0,Iterations)

# Important For Convolution2d - Keep as is.
mask = np.ones((3,3))
mask[1,1] = 0

######################################################################################################################################################
### Hard Coded Data  #################################################################################################################################
######################################################################################################################################################

# Set Definitions
Imax, Jmax, Tmax =  115, 115, 82 #(115,115,81)
Iset, Jset, Tset = range(0,Imax), range(0,Jmax), range(0,Tmax )
dimensions = (Tmax, Imax, Jmax)

# Variable Definitions - Preallocation Stage

obj_expr, obj_expr_summ, obj_func_value = np.zeros(dimensions), np.zeros(Tmax), np.zeros(Tmax)

J = np.zeros(dimensions)
Y = np.zeros(dimensions)
M = np.zeros(dimensions)
I = np.zeros(dimensions)
B = np.zeros(dimensions)
B_SDD = np.zeros(dimensions)
damages = np.zeros(dimensions)
costs = np.zeros(dimensions)

incell_weight = np.zeros(dimensions)
Redsnag, Greysnag = np.zeros(dimensions), np.zeros(dimensions)
dfactor = np.zeros(Tmax)


OYoung, UYoung = np.zeros(dimensions), np.zeros(dimensions)
OMature, UMature = np.zeros(dimensions), np.zeros(dimensions)
OAge, UAge = np.zeros(dimensions), np.zeros(dimensions)
OBhage, UBhage = np.zeros(dimensions), np.zeros(dimensions)

O_height, U_height = np.zeros(dimensions), np.zeros(dimensions)
O_diameter, U_diameter = np.zeros(dimensions), np.zeros(dimensions)
O_volumes, U_volumes = np.zeros(dimensions), np.zeros(dimensions)
young_alpha = np.zeros(dimensions)

# Variable Definitions - Assignments
PostTreatment_Initial = np.array(mpbdata(1, Imax, 1, Jmax, 3, PATH))
OAge_Initial = np.array(mpbdata(1, Imax, 1, Jmax, 0, PATH6))
UAge_Initial = np.array(mpbdata(1, Imax, 1, Jmax, 0, PATH7))
OStem_Initial = np.array(mpbdata(1, Imax, 1, Jmax, 0, PATH8))
UStem_Initial = np.array(mpbdata(1, Imax, 1, Jmax, 0, PATH9))
OJuvi_Initial = np.array(mpbdata(1, Imax, 1, Jmax, 0, PATH10))
UJuvi_Initial = np.array(mpbdata(1, Imax, 1, Jmax, 0, PATH11))

OBhage_Initial = Bhage(OAge_Initial, Iterations, Imax, Jmax)
UBhage_Initial = Bhage(UAge_Initial, Iterations, Imax, Jmax)

ones = np.ones((Imax,Imax))

# Dividing the OStem and UStem into young ( 60 <= age <= 120) and mature (age>120)
    # Adding the young and mature stems should be True for all entries
OYoung0 = np.where( OAge_Initial >= 60, OStem_Initial, 0)
OMature_Initial = np.where( OAge_Initial > 120, OStem_Initial, 0)
OYoung_Initial = OYoung0 - OMature_Initial


UYoung0 = np.where( UAge_Initial >= 60, UStem_Initial, 0)
UMature_Initial = np.where( UAge_Initial > 120, UStem_Initial, 0)
UYoung_Initial = UYoung0 - UMature_Initial

J_Initial = OJuvi_Initial + UJuvi_Initial
Y_Initial = OYoung_Initial + UYoung_Initial
M_Initial = OMature_Initial + UMature_Initial
I_Initial = PostTreatment_Initial




######################################################################################################################################################
### Growth Rate Data #################################################################################################################################
######################################################################################################################################################

# Non-Stochastic RCP Data (Order ~ 2.6, 4.5, 8.5) (From 2020 to 2100)
R_26 = mpb_growthrates(rowstart = 21, rowend = 101, column = 1, file_location=PATH4)
R_45 = mpb_growthrates(rowstart = 21, rowend = 101, column = 2, file_location=PATH4)
R_85 = mpb_growthrates(rowstart = 21, rowend = 101, column = 3, file_location=PATH4)

# Stochastic RCP Data (Order ~ 2.6, 4.5, 8.5) (From 2020 to 2100)...Note: Note that PATH5 will need to change depending on what type of std. dev is being used.
SR_data = mpb_stochasticgrowthrates(PATH5)
SR_26 = SR_data['rcp26']
SR_45 = SR_data['rcp45']
SR_85 = SR_data['rcp85']
######################################################################################################################################################
### Probability Distributions + Grids ################################################################################################################
######################################################################################################################################################
df_grid1,df_grid2,df_pgrid = pd.read_excel(PATH1), pd.read_excel(PATH2), pd.read_excel(PATH3)
I_x_J = Imax * Jmax
grid1 = df_grid1['Infest'].iloc[0:I_x_J].to_numpy().reshape((Imax,Jmax))
grid2 = df_grid2['Infest2'].iloc[0:I_x_J].to_numpy().reshape((Imax,Jmax))
grid3 = grid1 + grid2
pgrid = df_pgrid['inf_count'].iloc[0:I_x_J].to_numpy().reshape((Imax,Jmax))

ldd_prob, ldd_grid, counter = ldd_probability(grid3,Imax)
p_lambda = poisson_lambda(ldd_grid,pgrid,Imax)

size_of_sample = 1
infestation_probability = ldd_prob
number_of_samples = Tmax*Imax*Jmax


LDD_assignment = np.random.binomial(size_of_sample, infestation_probability, number_of_samples).reshape((Tmax,Imax,Imax))
LDD_deposit = np.random.poisson(p_lambda, (Tmax,Imax, Imax))
I_LDD = np.multiply(LDD_assignment, LDD_deposit)

print(I_LDD[0].sum())
######################################################################################################################################################
### Preallocating Solution Arrays ####################################################################################################################
######################################################################################################################################################

Y_solutions     = np.zeros((Iterations,Tmax))
J_solutions = np.zeros((Iterations,Tmax))
M_solutions     = np.zeros((Iterations,Tmax))
B_solutions     = np.zeros((Iterations,Tmax))
B_SDDsolutions = np.zeros((Iterations,Tmax))
I_solutions     = np.zeros((Iterations,Tmax))
I_LDDsolutions  = np.zeros((Iterations, Tmax))
Redsnag_solutions = np.zeros((Iterations,Tmax))
Greysnag_solutions     = np.zeros((Iterations,Tmax))
Obj_solutions   = np.zeros((Iterations,Tmax))
cost_solutions = np.zeros((Iterations,Tmax))
damage_solutions = np.zeros((Iterations,Tmax))

dme_solutions = np.zeros((Iterations,Tmax))