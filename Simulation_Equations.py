from SimulationModel.Data import *
from scipy.signal import convolve2d

def susceptible_rule(t):
    S[t + 1] = S[t] - G_NewGrowth[t]
    # Enforce Non-Negativity
    S[t + 1] = np.where(S[t+1] <= 0, 0, S[t+1])


def incell_weight_rule2(t):
    foo = np.where(S[t+1] > 60000, ones, S[t+1])
    incell_weight[t+1] = np.where(foo == 1, 1, 1 - 1/np.exp((c0 * foo) ** c1))


def post_dispersal_rule2(t):
    S_C = convolve2d(S[t+1],mask,mode='same')
    S_CR = np.divide(ones, S_C, out=np.zeros_like(ones), where = S_C != 0)
    new_matrix = G_PostTreatment[t] * (1-incell_weight[t+1]) * S_CR
    new_matrix_C = convolve2d(new_matrix, mask, mode='same')
    G_PostDispersal[t+1] = G_PostTreatment[t]*incell_weight[t+1] + S[t+1]*new_matrix_C + Z[t+1]

def g_newgrowth_rule(t,R):
    G_NewGrowth[t+1] = R[t+1]*G_PostDispersal[t+1] * (1 / np.exp(beta * (Redsnag[t+1] + Greysnag[t+1]),dtype=np.float64))

def g_posttreatment_rule_noclimate(t,L):
    G_PostTreatment[t + 1] = G_NewGrowth[t + 1] * (1 - L)

def g_posttreatment_rule_climate1(t,R,L):
    if 0 <= R[t] < 1:
        L[t] = 0.20
    elif 1 <= R[t] < 2:
        L[t] = 0.35
    elif 2 <= R[t] < 3:
        L[t] = 0.50
    elif 3 <= R[t] < 4:
        L[t] = 0.65
    elif 4 <= R[t] < 5:
        L[t] = 0.80
    elif 5 <= R[t]:
        L[t] = 0.95
    G_PostTreatment[t + 1] = G_NewGrowth[t + 1] * (1 - L[t])

def g_posttreatment_rule_climate2(t,R,L):
    if 0 <= R[t] < 1:
        L[t] = 0.45
    elif 1 <= R[t] < 2:
        L[t] = 0.55
    elif 2 <= R[t] < 3:
        L[t] = 0.65
    elif 3 <= R[t] < 4:
        L[t] = 0.75
    elif 4 <= R[t] < 5:
        L[t] = 0.85
    elif 5 <= R[t]:
        L[t] = 0.95
    G_PostTreatment[t + 1] = G_NewGrowth[t + 1] * (1 - L[t])

def g_posttreatment_rule_climate3(t,R, L):
    if 0 <= R[t] < 1:
        L[t] = 0.70
    elif 1 <= R[t] < 2:
        L[t] = 0.75
    elif 2 <= R[t] < 3:
        L[t] = 0.80
    elif 3 <= R[t] < 4:
        L[t] = 0.85
    elif 4 <= R[t] < 5:
        L[t] = 0.90
    elif 5 <= R[t]:
        L[t] = 0.95
    print(L)
    G_PostTreatment[t + 1] = G_NewGrowth[t + 1] * (1 - L[t])

def redsnag_rule(t):
    Redsnag[t + 1] = G_PostTreatment[t]

def greysnag_rule(t):
    Greysnag[t + 1] = (1 - d) * Greysnag[t] + Redsnag[t]

def objective_function_noclimate(t, L):
    obj_expr[t + 1] = (
                (c + price * tau * alpha[t] + f_healthy) * L * G_NewGrowth[t + 1] + alpha[t] * (f + price * tau) * G_PostTreatment[
            t + 1])
    obj_expr_summ[t + 1] = (1 / 1.05 ** t) * obj_expr[t + 1].sum()

def objective_function_climate(t, L):
    obj_expr[t + 1] = (
            (c + price * tau * alpha[t] + f_healthy) * L[t] * G_NewGrowth[t + 1] + alpha[t] * (f + price * tau) * G_PostTreatment[
        t + 1])
    obj_expr_summ[t + 1] = (1 / 1.05 ** t) * obj_expr[t + 1].sum()