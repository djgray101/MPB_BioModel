from SimulationModel.Adding_Juveniles.JuvenileData import *
from scipy.signal import convolve2d
from math import exp

def juvenile_rule(t):
    J[t + 1] = np.where( (OAge[t] == 59) | (UAge[t] == 59), 0, J[t])

def young_rule(t):
    Y[t + 1] = Y[t] - I[t] \
               + np.where( (OAge[t] == 59) | (UAge[t] == 59), J[t], 0) \
               - np.where( (OAge[t] == 119) | (UAge[t] == 119), Y[t] - I[t], 0) \
    # Enforce Non-Negativity
    Y[t + 1] = np.where( Y[t+1] <= 0, 0, Y[t + 1])

def mature_rule(t):
    M[t + 1] = M[t] + np.where( (OAge[t] == 119) | (UAge[t] == 119), Y[t] - I[t], 0)

def redattack_rule(t,L):
    Redsnag[t+1] = (1 - L) * I[t]
    Harvest[t+1] = L*I[t]

def greyattack_rule(t):
    Greysnag[t + 1] = Redsnag[t] + Greysnag[t]

def B_rule(t,R,L,Iteration):
    if t == 0:
        B[t + 1] = (1 - L) * R[t] * (B[0]) * (1 / np.exp(beta * (Redsnag[t] + Greysnag[t] + J[t] + M[t]), dtype=np.float64))
    else:
        B[t + 1] =  (1 - L) * R[t] * (B_SDD[t] + rho * I_LDD[t]) * (1 / np.exp(beta * (Redsnag[t] + Greysnag[t] + J[t] + M[t]), dtype=np.float64))


    # Enforce Non-Negativity
    #B[t+1] = np.where(B[t+1] <= 0, 0, B[t+1])

def B_SDD_rule(t):
    # Computing Incell Dispersal
    foo = np.where(Y[t + 1] > 20000, ones, Y[t + 1])
    incell_weight[t + 1] = np.where(foo == 1, 1, 1 - 1 / np.exp((c0 * foo) ** c1))
    # Computing Outcell Dispersal
    S_C = convolve2d(Y[t + 1], mask, mode='same')
    S_CR = np.divide(ones, S_C, out=np.zeros_like(ones), where = S_C != 0)

    new_matrix = B[t + 1] * (1 - incell_weight[t + 1]) * S_CR
    new_matrix_C = convolve2d(new_matrix, mask, mode='same')

    B_SDD[t + 1] = B[t+1] * incell_weight[t+1] + Y[t + 1] * new_matrix_C

def infection_rule(t,simulation_count):
    I[t+1] = (B_SDD[t+1]/rho) + I_LDD[t+1]

def objective_function_noclimate(t, L):
    obj_expr[t + 1] = 89.7*(L*I[t+1])**1.0889 + (young_alpha[t+1]*(price*tau + f_healthy) ) * L * I[t+1] + young_alpha[t]*(price*tau + f)*Redsnag[t+1]
    obj_expr_summ[t + 1] = (1 / (1 + discount) ** (t + 1)) * obj_expr[t + 1].sum()

    costs[t + 1] = (1 / (1 + discount) ** (t+1)) *  ( 89.7*(L*I[t+1])**1.0889 + ( young_alpha[t + 1] * (price * tau + f_healthy)) * L * I[t + 1] )
    damages[t + 1] = (1 / (1 + discount) ** (t+1)) *  young_alpha[t] * (price * tau + f) * Redsnag[t + 1]

def Age_Updates(t):
    OAge[t+1] = np.where(OAge[t] == 0, 0, OAge[t] + 1)
    UAge[t+1] = np.where(UAge[t] == 0, 0, UAge[t] + 1)
    OBhage[t+1] = np.where(OAge[t] == 0, 0, OBhage[t] + 1)
    UBhage[t+1] = np.where(UAge[t] == 0, 0, UBhage[t] + 1)

def stem_height(SI,Bhage,t):
    b1,b2,b3 = 8.1656,-1.3746,1.3922
    numerator = 1 + np.exp( b1 +  b2*np.log(50 + b3) - np.log(SI - 1.3) )
    denominator = 1 + np.exp( b1 + b2*np.log(Bhage[t] + b3) - np.log(SI - 1.3) )

    height = 1.3 + (SI - 1.3) * np.divide(numerator,denominator)
    height = np.where(height > 37, 37, height)
    return height

def stem_diameter(Height,t):
    a,b,c = 4.2512,-5.7514,-0.4588
    D = np.where(Height[t] < 1.3, 0, ( (np.log(Height[t] - 1.3) - a)/ b)**(1/c))
    D = np.where(D > 100, 100, D)
    return D

def alpha_volumes(SI,OBhage,UBhage,t):
    O_height[t],U_height[t] = stem_height(SI,OBhage,t), stem_height(SI,UBhage,t)
    O_diameter[t], U_diameter[t] = stem_diameter(O_height,t), stem_diameter(U_height,t)

    O_volumes[t] = np.pi * ((O_diameter[t]/(2*100))**2) * (O_height[t] - 0.3)
    U_volumes[t] = np.pi * ((U_diameter[t]/(2*100))**2) * (U_height[t] - 0.3)

    young_alpha[t] = np.where((OAge[t] >= 60) & (OAge[t] <= 120), O_volumes[t], 0) \
                        + np.where((UAge[t] >= 60) & (UAge[t] <= 120), U_volumes[t], 0)


######################################################################################################################################################
### Climate Based Control Policy #####################################################################################################################
######################################################################################################################################################
def climate_policy(R_scenario,strategy_number):
    C = np.zeros(Tmax)
    RCP_array = np.array(R_scenario)
    if strategy_number == 1:
        for t in range(len(RCP_array)):
            if 0 <= RCP_array[t] < 1:
                C[t] = 0.20
            elif 1 <= RCP_array[t] < 2:
                C[t] = 0.35
            elif 2 <= RCP_array[t] < 3:
                C[t] = 0.50
            elif 3 <= RCP_array[t] < 4:
                C[t] = 0.65
            elif 4 <= RCP_array[t] < 5:
                C[t] = 0.80
            elif 5 <= RCP_array[t]:
                C[t] = 0.95
        return C

    elif strategy_number == 2:
        for t in range(len(RCP_array)):
            if 0 <= RCP_array[t] < 1:
                C[t] = 0.45
            elif 1 <= RCP_array[t] < 2:
                C[t] = 0.55
            elif 2 <= RCP_array[t] < 3:
                C[t] = 0.65
            elif 3 <= RCP_array[t] < 4:
                C[t] = 0.75
            elif 4 <= RCP_array[t] < 5:
                C[t] = 0.85
            elif 5 <= RCP_array[t]:
                C[t] = 0.95
        return C

    else:
        for t in range(len(RCP_array)):
            if 0 <= RCP_array[t] < 1:
                C[t] = 0.70
            elif 1 <= RCP_array[t] < 2:
                C[t] = 0.75
            elif 2 <= RCP_array[t] < 3:
                C[t] = 0.80
            elif 3 <= RCP_array[t] < 4:
                C[t] = 0.85
            elif 4 <= RCP_array[t] < 5:
                C[t] = 0.90
            elif 5 <= RCP_array[t]:
                C[t] = 0.95
        return C
#
# def climate_policy2(R_scenario,strategy_number):
#     C = np.zeros(Tmax)
#     RCP_array = np.array(R_scenario)
#     if strategy_number == 1:
#         for t in range(len(RCP_array)):
#             if 0 <= RCP_array[t] < 1:
#                 C[t] = 0.20
#             elif 1 <= RCP_array[t] < 2:
#                 C[t] = 0.35
#             elif 2 <= RCP_array[t] < 3:
#                 C[t] = 0.50
#             elif 3 <= RCP_array[t] < 4:
#                 C[t] = 0.65
#             elif 4 <= RCP_array[t] < 5:
#                 C[t] = 0.80
#             elif 5 <= RCP_array[t]:
#                 C[t] = 0.95
#         return C**2
#
#     elif strategy_number == 2:
#         for t in range(len(RCP_array)):
#             if 0 <= RCP_array[t] < 1:
#                 C[t] = 0.45
#             elif 1 <= RCP_array[t] < 2:
#                 C[t] = 0.55
#             elif 2 <= RCP_array[t] < 3:
#                 C[t] = 0.65
#             elif 3 <= RCP_array[t] < 4:
#                 C[t] = 0.75
#             elif 4 <= RCP_array[t] < 5:
#                 C[t] = 0.85
#             elif 5 <= RCP_array[t]:
#                 C[t] = 0.95
#         return C**2
#
#     else:
#         for t in range(len(RCP_array)):
#             if 0 <= RCP_array[t] < 1:
#                 C[t] = 0.70
#             elif 1 <= RCP_array[t] < 2:
#                 C[t] = 0.75
#             elif 2 <= RCP_array[t] < 3:
#                 C[t] = 0.80
#             elif 3 <= RCP_array[t] < 4:
#                 C[t] = 0.85
#             elif 4 <= RCP_array[t] < 5:
#                 C[t] = 0.90
#             elif 5 <= RCP_array[t]:
#                 C[t] = 0.95
#         return C**2

