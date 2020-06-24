from SimulationModel.DataImportFunctions import *
from SimulationModel.HelperFunctions import *
import numpy as np

######################################################################################################################################################
### Mutable Parameters ###############################################################################################################################
######################################################################################################################################################

# Simulation Parameters
Iterations = 100000

# Model Parameters
c0, c1 = 0.000032, 3
beta = np.array(10.8 * 10**-6, dtype=np.float)
rho = 5000
price, tau = 177.54 , 0.05
f_healthy, f_dead = 6.58, 0.95
f = 5.63
c = 101.00
d = 0.02
discount = 0.05

LC = np.arange(0,10)

# Important For Convolution2d - Keep as is.
mask = np.ones((3,3))
mask[1,1] = 0

######################################################################################################################################################
### Hard Coded Data  #################################################################################################################################
######################################################################################################################################################

# Path Variables
PATH = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\Leading_EdgeGrid.gdb\Leading_EdgeGrids.xlsx"
PATH1 = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_16green.xls"
PATH2 = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_17green.xls"
PATH3 = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_17greenPoissonCount.xls"
PATH4 = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\Hinton_Climate_Data\GrowthRates_All.xlsx"

# Set Definitions
Imax, Jmax, Tmax = 115, 115, 11
Iset, Jset, Tset = range(0,Imax), range(0,Jmax), range(0,Tmax)
dimensions = (Tmax, Imax, Jmax)

# Spatial Data
TotalVol_data = mpbdata(1,Imax,1,Jmax,0,PATH)
TotalVol = np.array(TotalVol_data).sum()
alpha_data = np.array(mpbdata(1,Imax,1,Jmax,2,PATH))
alpha = np.zeros(dimensions)
alpha[0:Tmax] = alpha_data

# RCP Data (Order ~ 2.6, 4.5, 8.5)
R_data1_2030 = mpb_growthrates(rowstart=21,rowend=31,column=1, file_location=PATH4) # Row Start offset bythe column titles in Excel
R_data1_5060 = mpb_growthrates(rowstart=51,rowend=61,column=1, file_location=PATH4)
R_data1_8090 = mpb_growthrates(rowstart=81,rowend=91,column=1, file_location=PATH4)

R_data2_2030 = mpb_growthrates(rowstart=21,rowend=31,column=2, file_location=PATH4) # Row Start offset bythe column titles in Excel
R_data2_5060 = mpb_growthrates(rowstart=51,rowend=61,column=2, file_location=PATH4)
R_data2_8090 = mpb_growthrates(rowstart=81,rowend=91,column=2, file_location=PATH4)

R_data3_2030 = mpb_growthrates(rowstart=21,rowend=31,column=3, file_location=PATH4) # Row Start offset bythe column titles in Excel
R_data3_5060 = mpb_growthrates(rowstart=51,rowend=61,column=3, file_location=PATH4)
R_data3_8090 = mpb_growthrates(rowstart=81,rowend=91,column=3, file_location=PATH4)

# Variable Definitions - Initial Value Assignment
Susceptible_Initial = mpbdata(1, Imax, 1, Jmax, 1, PATH)
PostTreatment_Initial = mpbdata(1, Imax, 1, Jmax, 3, PATH)

# Variable Definitions - Empty Arrays
S, G_NewGrowth, incell_weight = np.zeros(dimensions), np.zeros(dimensions), np.zeros(dimensions)
G_PreDispersal, G_PostDispersal, G_PostTreatment = np.zeros(dimensions), np.zeros(dimensions), np.zeros(dimensions)
Redsnag, Greysnag = np.zeros(dimensions, dtype=np.float64), np.zeros(dimensions, dtype = np.float64)
TotalValue = np.zeros(Tmax)
dfactor = np.zeros(Tmax)
obj_expr = np.zeros(dimensions)
obj_expr_summ = np.zeros(Tmax)
obj_func_value = np.zeros(Tmax)
ones = np.ones((Imax,Imax))

######################################################################################################################################################
### Probability Distributions + Grids ################################################################################################################
######################################################################################################################################################
df_grid1,df_grid2,df_pgrid = pd.read_excel(PATH1), pd.read_excel(PATH2), pd.read_excel(PATH3)
I_x_J = Imax * Jmax
grid1 = df_grid1['Infest'].iloc[0:I_x_J].to_numpy().reshape((Imax,Jmax))
grid2 = df_grid2['Infest2'].iloc[0:I_x_J].to_numpy().reshape((Imax,Jmax))
grid3 = grid1 + grid2
pgrid = df_pgrid['inf_count'].iloc[0:I_x_J].to_numpy().reshape((Imax,Jmax))


ldd_prob, ldd_grid = ldd_probability(grid3,Imax)
p_lambda = poisson_lambda(ldd_grid,pgrid,Imax)

size_of_sample = 1
detection_probability = 0.90
infestation_probability = ldd_prob
number_of_samples = Tmax*Imax*Jmax

LDD_assignment = np.random.binomial(size_of_sample, infestation_probability, number_of_samples).reshape((Tmax,Imax, Imax))
LDD_deposit = np.random.poisson(p_lambda, (Tmax,Imax, Imax))
Z = np.multiply(LDD_assignment, LDD_deposit)
Detection = np.random.binomial(size_of_sample, detection_probability, number_of_samples).reshape((Tmax,Imax,Imax))

######################################################################################################################################################
### Probability Distributions + Grids ################################################################################################################
######################################################################################################################################################

