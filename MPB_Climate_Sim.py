import time
import concurrent.futures
import numpy as np
import pandas as pd
from math import exp
from SimulationModel.DataImportFunctions import *
from SimulationModel.HelperFunctions import *
import SimulationModel.Simulation_Equations as eq
# Path Variables
PATH = r"C:\Users\djgra\Desktop\Thesis_MPB\MPB_Data\GridExperiments\Leading_EdgeGrid.gdb\Leading_EdgeGrids.xlsx"

# Set Path To Probability Distribution Data
PATH1 = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_16green.xls"
PATH2 = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_17green.xls"
PATH3 = r"C:\Users\djgra\Jupyter_Notebooks\HintonRelated\LongDistanceDispersal\grid_17greenPoissonCount.xls"



######################################################################################################################################################
### Sets, Parameters, Variables ######################################################################################################################
######################################################################################################################################################

# Set Definitions
Imax, Jmax, Tmax = 115, 115, 15
Iset, Jset, Tset = range(0,Imax), range(0,Jmax), range(0,Tmax)
dimensions = (Tmax, Imax, Jmax)

# Non-Spatial Parameters
c0, c1 = 0.000032, 3
beta = 10.8 * 10**-6
rho = 5000
price, tau = 177.54 , 0.05
f_healthy, f_dead = 6.58, 0.95
f = 5.63
c = 97.00
d = 0.02
discount = 0.05

# Spatial Parameters
TotalVol_data = mpbdata(1,Imax,1,Jmax,0,PATH)
TotalVol = np.array(TotalVol_data).sum()
alpha_data = mpbdata(1,Imax,1,Jmax,2,PATH)
alpha = np.array(alpha_data)

# RCP Data (Order ~ 2.6, 4.5, 8.5)
R_data1_2030 = mpb_growthrates(Imax,Jmax,rowstart=21,rowend=31,column=1) # Row Start offset bythe column titles in Excel
R_data1_5060 = mpb_growthrates(Imax,Jmax,rowstart=51,rowend=61,column=1)
R_data1_8090 = mpb_growthrates(Imax,Jmax,rowstart=81,rowend=91,column=1)

R_data2_2030 = mpb_growthrates(Imax,Jmax,rowstart=21,rowend=31,column=2) # Row Start offset bythe column titles in Excel
R_data2_5060 = mpb_growthrates(Imax,Jmax,rowstart=51,rowend=61,column=2)
R_data2_8090 = mpb_growthrates(Imax,Jmax,rowstart=81,rowend=91,column=2)

R_data3_2030 = mpb_growthrates(Imax,Jmax,rowstart=21,rowend=31,column=3) # Row Start offset bythe column titles in Excel
R_data3_5060 = mpb_growthrates(Imax,Jmax,rowstart=51,rowend=61,column=3)
R_data3_8090 = mpb_growthrates(Imax,Jmax,rowstart=81,rowend=91,column=3)

# Variable Definitions - Empty Arrays

S, G_NewGrowth, incell_weight = np.zeros(dimensions), np.zeros(dimensions), np.zeros(dimensions)
G_PreDispersal, G_PostDispersal, G_PostTreatment = np.zeros(dimensions), np.zeros(dimensions), np.zeros(dimensions)
Redsnag, Greysnag = np.zeros(dimensions), np.zeros(dimensions)
TotalValue = np.zeros((Tmax))
dfactor = np.zeros((Tmax))
obj_expr = np.zeros(dimensions)
obj_expr_summ = np.zeros((Tmax))
obj_func_value = np.zeros((Tmax))

# Variable Definitions - Initial Value Assignment
Susceptible_Initial = mpbdata(1, Imax, 1, Jmax, 1, PATH)
PostTreatment_Initial = mpbdata(1, Imax, 1, Jmax, 3, PATH)

# Set t=0 Values:
S[0] = np.array(Susceptible_Initial)
G_NewGrowth[0] = np.array(PostTreatment_Initial)
G_PostTreatment[0] = G_NewGrowth[0].copy()
obj_func_value[0] = (price*tau + f_healthy)*TotalVol
L = 0

######################################################################################################################################################
### Dictionary Assignments For Indexing ##############################################################################################################
######################################################################################################################################################
neighbour_dict = {(a,b,c):neighbours(a,b,c,Imax) for a in Tset for b in Iset for c in Jset}

n_of_n_for_ijt_dict = {}
for t in Tset:
    if t == 10:
        break
    for i in Iset:
        for j in Jset:
            neighbour = neighbour_dict.get((t + 1, i, j))
            test2 = [neighbour_dict.get(entry) for entry in neighbour]  # Get the Neighbours of neighbours for (t+1,i,j)
            n_of_n_for_ijt_dict[(t + 1, i, j)] = test2

target_index_dict = {}
for t in Tset:
    if t == 10:
        break
    for i in Iset:
        for j in Jset:
            target_index_dict[(t + 1, i, j)] = [entry.index((t + 1, i, j)) for entry in n_of_n_for_ijt_dict.get((t + 1, i, j))]

######################################################################################################################################################
### Grids & Probability Distributions ################################################################################################################
######################################################################################################################################################

# grid1 == Green Attack 16 MPB data, grid 2 == Green Attack 17 MPB data
df_grid1,df_grid2,df_pgrid = pd.read_excel(PATH1), pd.read_excel(PATH2), pd.read_excel(PATH3)

grid1 = df_grid1['Infest'].to_numpy().reshape((Imax,Jmax))
grid2 = df_grid2['Infest2'].to_numpy().reshape((Imax,Jmax))
grid3 = grid1 + grid2
pgrid = df_pgrid['inf_count'].to_numpy().reshape((Imax,Jmax))

ldd_prob, ldd_grid = ldd_probability(grid3,115)
p_lambda = poisson_lambda(ldd_grid,pgrid,115)

size_of_sample = 1
detection_probability = 0.90
infestation_probability = ldd_prob
number_of_samples = 13225
######################################################################################################################################################
### Simulation Equations  ############################################################################################################################
######################################################################################################################################################

def susceptible_rule(t, i, j):
    S[t + 1, i, j] = S[t, i, j] - G_NewGrowth[t, i, j]
    # Enforce Non-Negativity
    if S[t + 1, i, j] < 0:
        S[t + 1, i, j] = 0
    else:
        pass

def incell_weight_rule(t, i, j):
    if S[t + 1, i, j] > 60000:
        incell_weight[t + 1, i, j] = 1
    else:
        incell_weight[t + 1, i, j] = 1 - 1 / exp((c0 * S[t + 1, i, j]) ** c1)

def post_dispersal_rule(t, i, j, local_dict):
    # Compute the first half of Equation 4: G_PostDispersal[t+1,i,j] = G_PostTreatment[t,i,j] - G_PostTreatment[t,i,j]*emmigration_weight_sum

    # local_dict = {(t+1,b,c):neighbour_summation(t+1,b,c,dim,neighbour_dict) for b in Iset for c in Jset }
    emmigration_values = local_dict.get((t + 1, i, j))

    # Eq. 3 == emmigration_weights
    emmigration_weights = [element * (1 - incell_weight[t + 1, i, j]) for element in emmigration_values]

    # Emmigration summation of Eq.4
    emmigration_weight_sum = sum(emmigration_weights)

    # Compute last half of Eq. 4
    immigration_list = [G_PostTreatment[a, b, c] for (a, b, c) in neighbour_dict.get((t, i, j))]

    test_list = [1 - incell_weight[a, b, c] for (a, b, c) in neighbour_dict.get((t + 1, i, j))]

    # x = func3(t+1,i,j,local_dict)
    coord = (t + 1, i, j)
    neighbours = neighbour_dict.get(coord)
    target_indicies = target_index_dict.get(coord)

    x = [local_dict[neighbours[i]][target_indicies[i]] for i, _ in enumerate(neighbours)]
    # x = [local_dict[neighbour_dict.get((t+1,i,j))[i]][target_index_dict.get((t+1,i,j))[i]] for i,_ in enumerate(neighbour_dict.get((t+1,i,j)))

    b1 = [test_list[i] * x[i] * immigration_list[i] for i in range(len(x))]

    G_PostDispersal[t + 1, i, j] = G_PostTreatment[t, i, j] - G_PostTreatment[t, i, j] * emmigration_weight_sum + sum(b1)

def g_newgrowth_rule(t, i, j, R):
    G_NewGrowth[t + 1, i, j] = R[t + 1] * G_PostDispersal[t + 1, i, j] * exp(-beta * (Redsnag[t + 1, i, j] + Greysnag[t + 1, i, j]))

def g_posttreatment_rule(t, i, j):
    G_PostTreatment[t + 1, i, j] = G_NewGrowth[t + 1, i, j] * (1 - L)

def redsnag_rule(t, i, j):
    Redsnag[t + 1, i, j] = G_PostTreatment[t, i, j]

def greysnag_rule(t, i, j):
    Greysnag[t + 1, i, j] = (1 - d) * Greysnag[t, i, j] + Redsnag[t, i, j]

def objective_function(t, i, j):
    obj_expr[t + 1, i, j] = (
                (c + price * tau * alpha[i, j] + f_healthy) * L * G_NewGrowth[t + 1, i, j] + alpha[i, j] * (f + price * tau) * G_PostTreatment[
            t + 1, i, j])
    obj_expr_summ[t + 1] = (1 / 1.05 ** t) * obj_expr[t + 1].sum()

######################################################################################################################################################
### Simulation  ######################################################################################################################################
######################################################################################################################################################

def simulation(RCP_data):
    Time_Out = 10
    simulation_count = 0
    simulation_cieling = 5
    simulation_solutions = []
    RCP_array = np.array(RCP_data)

    while simulation_count < simulation_cieling:
        # Set t=0 Values:
        S[0], G_NewGrowth[0], G_PostTreatment[0] = np.array(Susceptible_Initial), np.array(PostTreatment_Initial), G_NewGrowth[0].copy()

        for t in Tset:
            if t == Time_Out:
                break
            else:
                # print(f"Period:{t+1}")
                for i in Iset:
                    for j in Jset:
                        susceptible_rule(t, i, j)
                        redsnag_rule(t, i, j)
                        greysnag_rule(t, i, j)
                        incell_weight_rule(t, i, j)

                # Need New Loop - G_PostDispersal is dependent on all of the incells being computed first for the period.
                # Generate Probability Grids To Be Utilized For Next Looping Sequence
                LDD_assignment = np.random.binomial(size_of_sample,
                                                    infestation_probability,
                                                    number_of_samples).reshape((115, 115))

                LDD_deposit = np.random.poisson(p_lambda, (115, 115))

                Z = np.multiply(LDD_assignment, LDD_deposit)

                local_dict = {(t + 1, b, c): neighbour_summation(t + 1, b, c, S, neighbour_dict) for b in Iset for c in Jset}

                for i in Iset:
                    for j in Jset:
                        # Implement LDD by Adding to the G_t,PostDispersal variable
                        ### Compute PostDispersal & then add to it.
                        post_dispersal_rule(t, i, j, local_dict)

                G_PostDispersal[t + 1] = G_PostDispersal[t + 1] + Z

                for i in Iset:
                    for j in Jset:
                        g_newgrowth_rule(t, i, j, RCP_array)

                        # No Need to implement detection probability for L = 0
                        # if Detection[i,j] == 1:
                        g_posttreatment_rule(t, i, j)
                        objective_function(t, i, j)

        simulation_solutions.append(obj_expr_summ.sum())
        simulation_count += 1

    average_obj_value = np.array(simulation_solutions).sum() / len(simulation_solutions)
    return average_obj_value


solutions1 = {'2020-2030':0, '2050-2060':0, '2080-2090':0}
solutions2 = {'2020-2030':0, '2050-2060':0, '2080-2090':0}
solutions3 = {'2020-2030':0, '2050-2060':0, '2080-2090':0}



if __name__ == '__main__':
    t0 = time.time()
    with concurrent.futures.ProcessPoolExecutor() as executor:
        x1 = executor.submit(simulation, R_data1_2030)
        x2 = executor.submit(simulation, R_data1_5060)
        x3 = executor.submit(simulation, R_data1_8090)
        y1 = executor.submit(simulation, R_data2_2030)
        y2 = executor.submit(simulation, R_data2_5060)
        y3 = executor.submit(simulation, R_data2_8090)
        z1 = executor.submit(simulation, R_data3_2030)
        z2 = executor.submit(simulation, R_data3_5060)
        z3 = executor.submit(simulation, R_data3_8090)
        solutions1['2020-2030'] = x1.result()
        solutions1['2050-2060'] = x2.result()
        solutions1['2080-2090'] = x3.result()
        solutions2['2020-2030'] = y1.result()
        solutions2['2050-2060'] = y2.result()
        solutions2['2080-2090'] = y3.result()
        solutions3['2020-2030'] = z1.result()
        solutions3['2050-2060'] = z2.result()
        solutions3['2080-2090'] = z3.result()

        print(solutions1)
        print(solutions2)
        print(solutions3)

        t1 = time.time()
        print(t1 - t0)




