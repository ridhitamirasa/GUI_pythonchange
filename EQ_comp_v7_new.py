import numpy as np
import itertools
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

R = 8.314

#Written by Pilsun Yoo
def G_A_a_(dH_A_a, Tm_A, T):
    G_A_a = dH_A_a * (1. - (T / Tm_A))
    return G_A_a
def G_A_b_(dH_A_b, Tm_A, T):
    G_A_b = dH_A_b * (1. - (T / Tm_A))
    return G_A_b
def G_B_a_(dH_B_a, Tm_B, T):
    G_B_a = dH_B_a * (1. - (T / Tm_B))
    return G_B_a
def G_B_b_(dH_B_b, Tm_B, T):
    G_B_b = dH_B_b * (1. - (T / Tm_B))
    return G_B_b

#Written by Pilsun Yoo
def Gibbs_phase(X, T, G_A, G_B, E_phase):
    #G_A and G_B should be J/mol unit
    #E_phase is interaction parameter in J/mol
    G_phase = (1.-X) * G_A + X * G_B
    G_phase += R * T * (X * np.log(X) + (1.-X) * np.log(1.-X))
    G_phase += E_phase * X * (1.-X)
    return G_phase

def dG_dx(X, T, G_A, G_B, E_phase):
    dG_dx = -G_A + G_B + R * T * np.log(X/(1.-X)) + E_phase * (1. - 2. * X)
    return dG_dx

# Written by Sritapaswi Nori
def Linear(X, T, G_A, G_B, E_phase):
    x_list = []
    y_list = []
    
    slope = dG_dx(X, T, G_A, G_B, E_phase)
    x0 = X
    y0 = Gibbs_phase(X, T, G_A, G_B, E_phase)
    for i in range(1,100):
        x = i * 0.01
        y = slope * (x - x0) + y0
        x_list.append(x)
        y_list.append(y)
    return x_list, y_list, slope

# Written by Sritapaswi Nori
def Ga_Gb(X, T, G_A_a, G_B_a, G_A_b, G_B_b, Ea, Eb):
    G_alpha = Gibbs_phase(X=X, 
                          T=T, 
                          G_A=G_A_a, 
                          G_B=G_B_a, 
                          E_phase=Ea)
    G_beta = Gibbs_phase(X=X, 
                         T=T, 
                         G_A=G_A_b, 
                         G_B=G_B_b, 
                         E_phase=Eb)
    return G_alpha - G_beta

# Written by Jiaqi Yang
def dGa_dGb(X, T, G_A_a, G_B_a, G_A_b, G_B_b, Ea, Eb):
    dGa = dG_dx(X, T, G_A_a, G_B_a, Ea)
    dGb = dG_dx(X, T, G_A_b, G_B_b, Eb)
    return dGa - dGb

# Written by Sritapaswi Nori and Pilsun Yoo
# Fisrt database for binary system
def G_min_AgCu(T):
    x_list = []
    G_phase1_list = []
    G_liquid_list = []
    G_min_list = []
    for x_ in range(1, 9999):
        X = x_*0.0001
        G_phase1 = (1.-X) * (-11945 + 9.67*T) + X * (-13054 + 9.62*T)
        G_phase1 += R * T * (X * np.log(X) + (1.-X) * np.log(1.-X))
        G_phase1 += 25066*X*(1.-X)
        G_liquid = R * T * (X * np.log(X) + (1.-X) * np.log(1.-X))
        G_liquid += 14903*X*(1.-X)
        x_list.append(X)
        G_phase1_list.append(G_phase1)
        G_liquid_list.append(G_liquid)
        G_min_list.append(min([G_phase1, G_liquid]))
    
    return x_list, G_min_list

# Written by Jiaqi Yang
def Twopoint_Linear(x, y, x0, y0, x1, y1):
    slope = (y1-y0) / (x1-x0)
    t = y-y0 - slope * (x-x0)
    return t

# Written by Sritapaswi Nori and Pilsun Yoo
def Eq_comp_ConvexHull(x_list, G_sum_list):
    Eq_comp = []
    points = np.column_stack([x_list, G_sum_list])
    hull = ConvexHull(points)
    xhullvertices = points[hull.vertices, 0]
    sx = np.sort(xhullvertices, axis=None)
    
    #plt.plot(x_list, G_sum_list, 'b--')
    #plt.xlim(0, 1.)
    #plt.scatter(points[hull.vertices,0], points[hull.vertices,1], s=10, color='red')
    #plt.ylabel('Gibbs Free energy (J/mol)')
    #plt.xlabel('Fraction of B')
    #plt.title("Temperature: %d Kelvin" % T)
    
    sorted_xhullvertices = []
    for i in range(len(sx)-1):
        diff = sx[i+1] - sx[i]
        if diff > 0.00012:
            sorted_xhullvertices.append(sx[i])
            sorted_xhullvertices.append(sx[i+1])

    G = []
    for index_x in range(len(sorted_xhullvertices)):
        for index_x0 in range(len(x_list)):
            if x_list[index_x0] == sorted_xhullvertices[index_x]:
                G.append(G_sum_list[index_x0])
    for a in range(len(sorted_xhullvertices)):
        for b in range(a+1, len(sorted_xhullvertices[(a+1):])+(a+1)):
            x0 = sorted_xhullvertices[a]
            y0 = G[a]
            x1 = sorted_xhullvertices[b]
            y1 = G[b]

            test = [round(Twopoint_Linear(sorted_xhullvertices[index], G[index], x0, y0, x1, y1), 3)
                    for index in range(len(sorted_xhullvertices))]
            test2 = list(map(lambda x: x if x < 0 else 0, test))
            if len(list(set(test2))) <= 1:
                if abs(x0 - x1) > 0.02:
                    Eq_comp.append([x0, x1])
                    #plt.scatter(x0, y0, s=20, color='black')
                    #plt.scatter(x1, y1, s=20, color='black')
    #plt.show()
    #print(Eq_comp)
    return Eq_comp


# Written by Pilsun Yoo and Jiaqi Yang
def Binary_Phase(T1, T2, system_string, step=10):
    #T1 lowest temperature
    #T2 highest temperature
    #system_string : the number of the system ex) 'AgCu'
    Tem = []
    Compositions = []
    for T in range(T1, T2, step):
        x_list, G_sum_list = eval("G_min_%s(T)" % system_string)
        Eq = Eq_comp_ConvexHull(x_list, G_sum_list)
        for j in range(len(Eq)):
            for k in range(len(Eq[j])):
                Tem.append(round(T,2))
                Compositions.append(Eq[j][k])
    return Tem, Compositions

# Written by Jiaqi Yang
if __name__=='__main__':
    
    #Temperature in Degree C
    T1=200
    T2=1300
    step=5
    system_string = 'AgCu'
    
    #Convert from Degree C to K
    T1 = T1 + 273.15
    T2 = T2 + 273.15
    Tem, Compositions = Binary_Phase(int(T1), int(T2), system_string, step=step)
    Tem = list(np.array(Tem) - 273.15)
    
    """
    Tem_count = {round(i,2):Tem.count(i) for i in Tem}
    parallel = []
    k = list(Tem_count.keys())
    v = list(Tem_count.values())
    for i in range(len(Tem_count.keys())-1):
    pn = v[i+1] - v[i]
    if pn > 0:
    parallel.append(k[i+1])
    D = {}
    for i in range(len(parallel)):
    tl = []
    t = parallel[i]
    for j in range(len(Tem)):
    if t == Tem[j]:
    tl.append(Compositions[j])
    D[t] = min(tl), max(t1)
    """
    
    plt.scatter(Compositions, Tem, s=3, color=(0,0,0), alpha=1)
    plt.xlabel('Molar Fraction of B')
    plt.xlim(0, 1)
    plt.ylabel('Temperature (°C)')
    figname = E1 + "_" + E2 + "_" + str(T1) + "_" + str(T2) + ".png"
    plt.savefig(figname)
    plt.show()

