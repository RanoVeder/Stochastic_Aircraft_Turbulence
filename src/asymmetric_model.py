import numpy as np
import control
#Conditions
V = 121.3
m = 5445
W = 53361
S = 24.2
c = 2.022
b = 13.36
h = 9144
rho = 0.4587
mu_c = 209
mu_b = 32
K_X2 = 0.013
K_Z2 = 0.037
KXZ = 0.002
KY2 = 0.950

CL = 0

CY_beta  = -1.3250
CYp  = -0.1320
CYr  = 0.4300
CY_delta_alpha = 0.0000
CY_delta_r = 0.3037 

Cl_beta  =-0.1070
Clp  =-0.3684
Clr  = 0.1750
Cl_delta_alpha =-0.2292
Cl_delta_r = 0.0446

Cn_beta  = 0.1835
Cnp  =-0.0035
Cnr  =-0.1930
Cn_delta_alpha = 0.0071
Cn_delta_r =-0.1261

Clpw = 0.8*Clp
Cnpw = 0.9*Cnp
Clrw = 0.7*Clr
Cnrw = 0.2*Cnr


fact_y = V/(b*2*mu_b)
y_beta = fact_y*CY_beta
y_phi = 0
y_p = fact_y*CYp
y_r = fact_y*(CYr-4*mu_b)
y_delta_r = fact_y*CY_delta_r

fact_l = b*4*mu_b*(K_X2*K_Z2 - KXZ**2)/V
l_beta = (Cl_beta*K_Z2+ Cn_beta * KXZ) / fact_l
l_p = (Clp*K_Z2 + Cnp * KXZ)/ fact_l
l_r = (Clr * K_Z2 + Cnr * KXZ) / fact_l
l_delta_alpha = (Cl_delta_alpha * K_Z2 + Cn_delta_alpha*KXZ) / fact_l
l_delta_r = (Cl_delta_r*K_Z2 + Cn_delta_r*KXZ) / fact_l

n_beta = (Cl_beta * KXZ  + Cn_beta  * K_X2) / fact_l
n_p = (Clp * KXZ + Cnp *K_X2) / fact_l
n_r = (Clr * KXZ + Cnr * K_X2) / fact_l
n_delta_alpha = (Cl_delta_alpha * KXZ + Cn_delta_alpha * K_X2) / fact_l
n_delta_r = (Cl_delta_r * KXZ +  Cn_delta_r * K_X2) / fact_l


A = np.array([
    [y_beta, y_phi, y_p,    y_r],
    [0,      0,     2*V/b,  0],
    [l_beta, 0,     l_p,    l_r],
    [n_beta, 0,     n_p,    n_r]
])

B = np.array([
    [0,             y_delta_r],
    [0,             0        ],
    [l_delta_alpha, l_delta_r],
    [n_delta_alpha, n_delta_r]
])

C = np.eye(4)
D = np.zeros((4,2))
sys = control.StateSpace(A, B, C, D)

t, y = control.step(sys)
print(np.linalg.eigvals(A))