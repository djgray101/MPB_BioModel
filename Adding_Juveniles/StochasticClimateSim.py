import time
import concurrent.futures
from SimulationModel.Adding_Juveniles.Simulation_Equations import *


def simulation(RCP_data,Stochastic_data, SiteIndex, Iterations, policy_number):
    simulation_count = 0
    Time_Out = Tmax - 1
    #Time_Out = Tmax
    #RCP_array = np.array(RCP_data)
    L = climate_policy(RCP_data, strategy_number=policy_number)
    while simulation_count < Iterations:

        SRCP_array = np.array(Stochastic_data)[simulation_count*81:(simulation_count+1)*81]

        # Set t=0 Values:
        J[0],Y[0],M[0],I[0] = J_Initial, Y_Initial, M_Initial, PostTreatment_Initial.copy()
        B[0] = rho * PostTreatment_Initial
        Redsnag[0],Greysnag[0] = 0, 0
        OAge[0], UAge[0] = OAge_Initial, UAge_Initial
        OBhage[0], UBhage[0] = OBhage_Initial[simulation_count], UBhage_Initial[simulation_count]
        alpha_volumes(SiteIndex,OBhage,UBhage,0)   # Initialize the first period alpha's.

        for t in Tset:
            if t == Time_Out:
                break
            else:
                #print(f'Simulation Run:{simulation_count}')
                # Demographic Updates Before Model Dynamics
                Age_Updates(t)
                alpha_volumes(SiteIndex,OBhage,UBhage,t+1)
                juvenile_rule(t)
                young_rule(t)
                mature_rule(t)
                redattack_rule(t,L[t])
                greyattack_rule(t)
                B_SDD_rule(t)
                infection_rule(t,simulation_count)
                B_rule(t, SRCP_array, L[t], simulation_count)
                objective_function_noclimate(t, L[t])

                # Solution Vectors
                Obj_solutions[simulation_count][t] = obj_expr_summ[t].sum()
                cost_solutions[simulation_count][t] = costs[t].sum()
                damage_solutions[simulation_count][t] = damages[t].sum()
                Y_solutions[simulation_count][t] = Y[t].sum()
                J_solutions[simulation_count][t] = J[t].sum()
                M_solutions[simulation_count][t] = M[t].sum()
                I_solutions[simulation_count][t] = I[t].sum()
                Redsnag_solutions[simulation_count][t] = Redsnag[t].sum()
                Greysnag_solutions[simulation_count][t] = Greysnag[t].sum()

        simulation_count += 1


    soln_list = np.array([Obj_solutions, cost_solutions, damage_solutions,J_solutions, Y_solutions, M_solutions,I_solutions, Redsnag_solutions,Greysnag_solutions])
    col_names = ['Objective','Costs', 'Damages','J','Y','M','I', 'Redsnag', 'Greysnag'];
    df = pd.DataFrame()
    years = np.arange(2020,2102,1)
    df['years'] = years

    for name, soln in zip(col_names,soln_list):
        avg_solutions = np.zeros(Tmax)
        for t in Tset:
            avg_solutions[t] = np.average(soln[:,t])

        df[name] = avg_solutions.round()

    df['control'] = L
    return df


if __name__ == '__main__':
    t0 = time.time()
    R_list = [R_26, R_45, R_85]
    SR_list = [SR_26, SR_45, SR_85]
    SI_list = [SI_26, SI_45, SI_85]
    print("Starting Simulation...")
    print("")

    # Cycle through each control policy
    for num in range(1,4):
        solutions_df = pd.DataFrame()
        # Cycle through each scenario
        for x,y,z in zip(R_list, SR_list, SI_list):

            solutions = simulation(x, y, z, Iterations, policy_number=num)
            solutions_df = pd.concat([solutions_df, solutions], ignore_index=True)

            solutions_df.to_csv(f'Results\StochasticControl\S_climatetestsolns{num}.txt',index=False)

    t1 = time.time()
    print(f'Time To Complete {Iterations} Simulations:', (t1 - t0)/3600)