import sympy
import math

Vout = sympy.Symbol("Vout")
Vout_aux = 0
hl1 = 0
hl2 = 0
rho = 1000
mi = 0.001
L1 = 10
L2 = 30
g = 9.81

ex = 1000*Vout*((math.pi/4)*(0.04)**2)*(-(Vout**2/2)+20-hl1-hl2) - 500
def Rey(Vel, diam):
    global mi
    global rho
    Re = (rho*Vel*diam)/mi
    return Re
def coeficiente(Reynolds):
    alpha = 6.9/Reynolds
    delta = (math.log10(alpha**(-1.8)))
    coe = (1/(delta**2))
    return coe
def tau(Vel, coef_f):
    global rho
    ans = rho*(Vel**2)*(coef_f/8)
    return ans
def HL(T,diam,comp):
    hl = (4*comp*T)/(rho*g*diam)
    return hl

Vin = (4/9)*Vout

while True:
    Vout_aux1 = Vout_aux
    solution = sympy.solveset(1000 * Vout * ((math.pi / 4) * (0.04) ** 2) * (-(Vout ** 2 / 2) + 20 - hl1 - hl2) + 500, Vout)
    Vout_aux, dontmatter1, dontmatter2 = solution
    Vin = (4 / 9) * Vout_aux
    Rout =Rey(Vout_aux,0.04)
    Rin = Rey(Vin,0.06)
    Cout = coeficiente(Rout)
    Cin = coeficiente(Rin)
    tauout = tau(Vout_aux,Cout)
    tauin = tau(Vin,Cin)
    hl2 = HL(tauout, 0.04, 30)
    hl1 = HL(tauin,0.06,10)
    Vout_aux2 = Vout_aux
    if round(Vout_aux1,2) == round(Vout_aux2,2):
        break
print("Velocidade de saida:",Vout_aux)