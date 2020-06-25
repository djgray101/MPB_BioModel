import time
import concurrent.futures
from SimulationModel.Simulation_Equations import *


def simulation(RCP_data, Iterations, L):
    print('Starting Simulation')
    simulation_count = 0
    Time_Out = 10
    simulation_solutions = []
    RCP_array = np.array(RCP_data, dtype = np.float64)

    # Constant Control L:

    while simulation_count < Iterations:
        # Set t=0 Values:
        S[0], G_NewGrowth[0], G_PostTreatment[0] = np.array(Susceptible_Initial), np.array(PostTreatment_Initial), G_NewGrowth[0].copy()

        for t in Tset:
            if t == Time_Out:
                break
            else:
                susceptible_rule(t)
                redsnag_rule(t)
                greysnag_rule(t)
                incell_weight_rule2(t)
                post_dispersal_rule2(t)
                g_newgrowth_rule(t,RCP_array)
                g_posttreatment_rule_noclimate(t, L)
                objective_function_noclimate(t,L)

        simulation_solutions.append(obj_expr_summ.sum())
        simulation_count += 1

    average_obj_value = np.array(simulation_solutions).sum() / len(simulation_solutions)
    return average_obj_value


solutions1 = {'2020-2030':0, '2050-2060':0, '2080-2090':0}
solutions2 = {'2020-2030':0, '2050-2060':0, '2080-2090':0}
solutions3 = {'2020-2030':0, '2050-2060':0, '2080-2090':0}

if __name__ == '__main__':
    soln_file = open("ResultFile.txt", "w+")
    for L in np.arange(0,0.55, 0.05):
        t0 = time.time()
        with concurrent.futures.ProcessPoolExecutor() as executor:

            x1 = executor.submit(simulation, R_data1_2030, Iterations, L)
            x2 = executor.submit(simulation, R_data1_5060, Iterations, L)
            x3 = executor.submit(simulation, R_data1_8090, Iterations, L)
            y1 = executor.submit(simulation, R_data2_2030, Iterations, L)
            y2 = executor.submit(simulation, R_data2_5060, Iterations, L)
            y3 = executor.submit(simulation, R_data2_8090, Iterations, L)
            z1 = executor.submit(simulation, R_data3_2030, Iterations, L)
            z2 = executor.submit(simulation, R_data3_5060, Iterations, L)
            z3 = executor.submit(simulation, R_data3_8090, Iterations, L)
            solutions1['2020-2030'] = x1.result()
            solutions1['2050-2060'] = x2.result()
            solutions1['2080-2090'] = x3.result()
            solutions2['2020-2030'] = y1.result()
            solutions2['2050-2060'] = y2.result()
            solutions2['2080-2090'] = y3.result()
            solutions3['2020-2030'] = z1.result()
            solutions3['2050-2060'] = z2.result()
            solutions3['2080-2090'] = z3.result()

        soln_file.write(f"Control Effort: {L} \r\n")
        soln_file.write(f"RCP 2.6: {solutions1} \r")
        soln_file.write(f"RCP 4.5: {solutions2} \r")
        soln_file.write(f"RCP 8.5: {solutions3} \r\n")

    soln_file.close()
