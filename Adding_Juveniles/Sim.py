import time
from SimulationModel.Adding_Juveniles.Simulation_Equations import *

def simulation(RCP_data, SiteIndex, Iterations, L):
    simulation_count = 0
    Time_Out = Tmax - 1
    RCP_array = np.array(RCP_data)

    while simulation_count < Iterations:
        # Set t=0 Values:
        J[0],Y[0],M[0],I[0] = J_Initial, Y_Initial, M_Initial, PostTreatment_Initial.copy()
        B[0] = rho * PostTreatment_Initial
        B_SDD[0] = B[0].copy()
        Redsnag[0],Greysnag[0] = 0, 0
        OAge[0], UAge[0] = OAge_Initial, UAge_Initial
        OBhage[0], UBhage[0] = OBhage_Initial[simulation_count], UBhage_Initial[simulation_count]
        alpha_volumes(SiteIndex,OBhage,UBhage,0)   # Initialize the first period alpha's.
        for t in Tset:
            if t == Time_Out:
                break
            else:
                # Demographic Updates Before Model Dynamics
                Age_Updates(t)
                alpha_volumes(SiteIndex,OBhage,UBhage,t+1)
                juvenile_rule(t)
                young_rule(t)
                mature_rule(t)
                redattack_rule(t,L)
                greyattack_rule(t)
                B_rule(t, RCP_array, L, simulation_count)
                B_SDD_rule(t)
                infection_rule(t,simulation_count)
                objective_function_noclimate(t, L)

                # Solution Vectors
                Obj_solutions[simulation_count][t] = obj_expr_summ[t].sum()
                cost_solutions[simulation_count][t] = costs[t].sum()
                damage_solutions[simulation_count][t] = damages[t].sum()
                Y_solutions[simulation_count][t] = Y[t].sum()
                J_solutions[simulation_count][t] = J[t].sum()
                M_solutions[simulation_count][t] = M[t].sum()
                I_solutions[simulation_count][t] = I[t].sum()
                B_solutions[simulation_count][t] = B[t].sum()
                B_SDDsolutions[simulation_count][t] = B_SDD[t].sum()
                I_LDDsolutions[simulation_count][t] = I_LDD[t].sum()
                Redsnag_solutions[simulation_count][t] = Redsnag[t].sum()
                Greysnag_solutions[simulation_count][t] = Greysnag[t].sum()

        simulation_count += 1

        #simulation_solutions.append(obj_expr_summ.sum())
    #average_obj_value = np.array(simulation_solutions).sum() / len(simulation_solutions)

    soln_list = np.array([Obj_solutions, cost_solutions, damage_solutions,J_solutions, Y_solutions, M_solutions,I_solutions, Redsnag_solutions,Greysnag_solutions, B_solutions, B_SDDsolutions, I_LDDsolutions])
    col_names = ['Objective','Costs', 'Damages','J','Y','M','I', 'Redsnag', 'Greysnag', 'B', 'B_SDD','I_LDD']
    df = pd.DataFrame()
    years = np.arange(2020,2102,1)
    df['years'] = years


    for name, soln in zip(col_names, soln_list):
        df[name] = np.average(soln, axis = 0)


    return df


if __name__ == '__main__':

    t0 = time.time()
    R_list = [R_26, R_45, R_85]
    SI_list = [SI_26, SI_45, SI_85]
    txt_labels = ['26','45','85']
    for x,y,z in zip(R_list, SI_list, txt_labels):
        solutions_df = pd.DataFrame()
        # L = 1
        # solutions = simulation(x, y, Iterations, L)
        # solutions_df = pd.concat([solutions_df, solutions], ignore_index=True)
        for L in np.arange(0,1.1,0.1):
            if L == 1.1:
                break
            else:
                print(L)
                print("Starting Simulation...")
                print("")
                #solutions['RCP85'] = simulation(R_26, SI_26, Iterations, L)
                solutions = simulation(x, y, Iterations, L)
                solutions_df = pd.concat([solutions_df, solutions], ignore_index=True)

        solutions_df.to_csv(f'Results\ConstantControl\const{z}_1.txt',index=False)

    t1 = time.time()
    print(f'Time To Complete {Iterations} Simulations:', (t1 - t0)/3600)