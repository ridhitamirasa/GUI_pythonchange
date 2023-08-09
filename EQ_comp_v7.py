import numpy as np
import itertools
from scipy.signal import argrelextrema
# import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt

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
    for x_ in range(1, 100):
        X = x_*0.01
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
    sorted_xhullvertices = np.sort(xhullvertices, axis=None)

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
                Tem.append(T-273.15)
                Compositions.append(Eq[j][k])
    return Tem, Compositions

# Written by Jiaqi Yang
if __name__=='__main__':
    T1=400
    T2=900
    step=10
    system_string = 'AgCu'

    Tem, Compositions = Binary_Phase(T1, T2, system_string, step=step)
    plt.plot(Compositions, Tem, s=1, c=(0,0,0), alpha=1)
    plt.xlabel('Molar Fraction of B')
    plt.xlim(0, 1)
    plt.ylabel('Temperature (°C)')
    figname = E1 + "_" + E2 + "_" + str(T1) + "_" + str(T2) + ".png"
    plt.savefig(figname)
    plt.show()

################### Previous Versions

# Written by Pilsun Yoo
def G_Intersection(T, G_info):

    #Find intersections
    #return list of intersections
    
    G_comb = list(itertools.combinations(G_info, 2))
    N_comb = len(G_comb)
    intersections = []
    n_inter = []
    
    for i in range(N_comb):
        a, b = G_comb[i]
        G_A_a, G_B_a, Ea = a[0], a[1], a[2]
        G_A_b, G_B_b, Eb = b[0], b[1], b[2]
        Xs = []
        Xs2 = []
        for i in range(1,10000):
            X = i * 0.0001
            Ga_Gb1 = Ga_Gb(X, T, G_A_a, G_B_a, G_A_b, G_B_b, Ea, Eb)
            if abs(Ga_Gb1) < 1.0:
                Xs.append(X)
        #print(Xs)
        for X2 in Xs:
            intersections.append(round(X2,4))

    intersections = sorted(intersections)
    remove_list = []
    for c1 in range(len(intersections)):
        for c2 in range(c1+1, len(intersections)):
            if abs(intersections[c2] - intersections[c1]) < 0.001:
                remove_list.append(intersections[c2])
                
    remove_list = list(set(remove_list))
    for r in remove_list:
        intersections.remove(r)
            
    return intersections

# Written by Pilsun Yoo
#Find local minimas of Gibb free energy curves
def Local_min_Gibbs(T, G_info):
    
    G_comb = list(itertools.combinations(G_info, 2))
    N_comb = len(G_comb)
    
    local_min = []
    
    for i in range(N_comb):
        a, b = G_comb[i]
        G_A_a, G_B_a, Ea = a[0], a[1], a[2]
        G_A_b, G_B_b, Eb = b[0], b[1], b[2]    
    
        x_list = []
        y1_list = []
        y2_list = []
        for i in range(1,1000):
            xt = i * 0.001
            x_list.append(xt)
            y1_list.append(Gibbs_phase(xt, T, G_A_a, G_B_a, Ea))
            y2_list.append(Gibbs_phase(xt, T, G_A_b, G_B_b, Eb))

        x_list = np.array(x_list)
        y1_list = np.array(y1_list)
        y2_list = np.array(y2_list)
        minm1 = argrelextrema(y1_list, np.less)
        minm2 = argrelextrema(y2_list, np.less)
        Ga_min = x_list[minm1]
        Gb_min = x_list[minm2]
        
        for Ga in Ga_min:
            local_min.append(Ga)
        for Gb in Gb_min:
            local_min.append(Gb)
    local_min = list(set(sorted(local_min)))
    return local_min

#2nd version
#Written by Pilsun Yoo and Sritapaswi Nori
def Linear2(m, x, T, G_A, G_B, E_phase):
    
    slope = dG_dx(m, T, G_A, G_B, E_phase)
    x0 = m
    y0 = Gibbs_phase(m, T, G_A, G_B, E_phase)        
    y = slope * (x - x0) + y0
    return y, slope

#Written by Pilsun Yoo and Sritapaswi Nori
def eqcomp_convexhull(T, G_A_a, G_B_a, G_A_b, G_B_b, Ea, Eb):
    
    G_info = [[G_A_a, G_B_a, Ea], [G_A_b, G_B_b, Eb]]
    intersections = G_Intersection(T, G_info)
    minima = Local_min_Gibbs(T, G_info)
    
    intersections.insert(0, 0.0000001)
    intersections.append(0.9999999)
    a = []
    b = []
    xt = []
    for r in range(len(intersections)-1):
        
        lower = intersections[r]
        upper = intersections[r+1]
        mid = (upper - lower) / 2
        
        G_alpha = Gibbs_phase(X=mid, 
                              T=T, 
                              G_A=G_A_a, 
                              G_B=G_B_a, 
                              E_phase=Ea)
        G_beta = Gibbs_phase(X=mid, 
                             T=T, 
                             G_A=G_A_b, 
                             G_B=G_B_b, 
                             E_phase=Eb)
        Gd = G_alpha - G_beta
        step = 0.0005
        #step = (upper - lower) / 1000
        x = lower
        if Gd > 0:
            mu_a = []
            mu_b = []
            xt0 = []
            while x <= upper:
                mu_a1, sa1 = Linear2(x, 0.0000001, T, G_A_b, G_B_b, Eb)
                mu_b1, sb1 = Linear2(x, 0.9999999, T, G_A_b, G_B_b, Eb)
                
                mu_a.append(mu_a1)
                mu_b.append(mu_b1)
                xt0.append(x)
                x += step
            xt.append(xt0)
            a.append(mu_a)
            b.append(mu_b)
            
        elif Gd < 0:
            x = lower
            mu_a = []
            mu_b = []
            xt0 = []
            while x <= upper:
                mu_a1, sa1 = Linear2(x, 0.0000001, T, G_A_a, G_B_a, Ea)
                mu_b1, sb1 = Linear2(x, 0.9999999, T, G_A_a, G_B_a, Ea)
                
                mu_a.append(mu_a1)
                mu_b.append(mu_b1)
                xt0.append(x)
                x += step
            xt.append(xt0)
            a.append(mu_a)
            b.append(mu_b)
    tmp = []
    for k in range(len(a[0])):
        a1 = a[0][k]
        for l in range(len(a[1])):
            a2 = a[1][l]
            if abs(a1-a2) < 10:
                #print(xt[0][k])
                tmp.append([xt[0][k], xt[1][l]])
    
    com = []
    slope = []
    for t in range(len(tmp)):
        mu_A, sA = Linear2(tmp[t][0], 0.0000001, T, G_A_a, G_B_a, Ea)
        mu_A, sB = Linear2(tmp[t][1], 0.0000001, T, G_A_b, G_B_b, Eb)
        com.append(tmp[t])
        slope.append(abs(sA-sB))
    index = slope.index(min(slope))
    Eq = [com[index][0], com[index][1]]
    return Eq
