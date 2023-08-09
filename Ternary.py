import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
# from mpl_toolkits.mplot3d import Axes3D
# from shapely.geometry import Point
# from shapely.geometry.polygon import Polygon
# import random
# import itertools

R = 8.3145


# Real example 1
#Developed by Sritapaswi Nori & Pilsun Yoo
def Gibbs_phase_1(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC):
    # G_A,G_B, and G_C should be J/mol unit
    G_phase = XA * G_A + XB * G_B + (1 - XA - XB) * G_C
    G_phase += R * T * (XA * np.log(XA) + XB * np.log(XB) + (1 - XA - XB) * np.log(1 - XA - XB))
    G_phase += LAB * XA * XB + LBC * XB * (1 - XA - XB) + LAC * XA * (1 - XA - XB) + XA * XB * (1 - XA - XB) * (LABC)
    return G_phase


def Gibbs_phase_2(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC):
    T = 50.33
    a = 2510
    G_phase = a * (XA * XB + XB * (1 - XA - XB) + XA * (1 - XA - XB))
    G_phase += R * T * (XA * np.log(XA) + XB * np.log(XB) + (1 - XA - XB) * np.log(1 - XA - XB))
    return G_phase


def Gibbs_phase_3(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC):
    T = 1900
    b = 2 * 10 ** 5
    G_phase = b * (XA * XB * (1 - XA - XB))
    G_phase += R * T * (XA * np.log(XA) + XB * np.log(XB) + (1 - XA - XB) * np.log(1 - XA - XB))
    return G_phase


def Gibbs_phase(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC, tag):
    if tag == 1:
        return Gibbs_phase_1(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC)
    elif tag == 2:
        return Gibbs_phase_2(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC)
    elif tag == 3:
        return Gibbs_phase_3(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC)


def Gliquid(XA, XB, T, tag):
    LAB = 6200 - 0.418 * T + (XA - XB) * (790 - 1.914 * T)
    LAC = 110 - (2.5 * T) + (XA - (1 - XA - XB)) * (420 - 1.05 * T) + ((XA - (1 - XA - XB)) ** 2) * (0.36 * T)
    LBC = - 5695 - (1.71 * T) + (1 - XA - 2 * XB) * 780 + ((1 - XA - 2 * XB) ** 2) * 1840
    LABC = (4910 - 34 * T) * XA - 9240 * XB - 13120 * (1 - XA - XB)
    if 298 < T <= 600:
        G_A = 4672.157 - (7.750257 * T) - (6.0144 * 10 ** -19) * (T ** 7)
    elif 600 < T < 5000:
        G_A = 4853.112 - (8.066587 * T) - (8.05644 * 10 ** 25) * (T ** -9)
    if 298 < T <= 505:
        G_B = 7104.38 - (14.0895659 * T) + (1.49503 * 10 ** -18) * (T ** 7)
    elif 505 < T < 3000:
        G_B = 6970.584 - (13.811447 * T) + (1.25305 * 10 ** 25) * (T ** -9)
    if 298 < T < 904:
        G_C = 19822.595 - (21.920597 * T) - (1.73785 * 10 ** -20) * (T ** 7)
    elif 904 <= T < 2000:
        G_C = 19913.982 - (22.026755 * T) - (1.610442 * 10 ** 27) * (T ** -9)
    return Gibbs_phase(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC, tag)


def Galpha(XA, XB, T, tag):
    LAB = 10000
    LBC = -5760 + 11.834 * T
    LAC = 21360 - 5.66 * T
    LABC = 0
    G_A = 9900 - 12.5 * T
    G_B = 2035
    G_C = 0
    return Gibbs_phase(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC, tag)


def Ggamma(XA, XB, T, tag):
    LAB = 19700 - 15.89 * T
    LBC = 25400 - 36.23 * T + (1 - XA - 2 * XB) * (-10400 + 36.23 * T)
    LAC = 10000
    LABC = 0
    G_A = 489 + 3.52 * T
    G_B = 0
    G_C = 1000
    return Gibbs_phase(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC, tag)


def Gepsilon(XA, XB, T, tag):
    LAB = 16510 - 19.9 * T + (XA - XB) * 1850
    LBC = 0
    LAC = 0
    LABC = 0
    G_A = 300 + T
    G_B = 2400 - 3.1 * T
    G_C = 0
    return Gibbs_phase(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC, tag)


def Gdelta(XA, XB, T, tag):
    LAB = 7860 - 4.94 * T
    LBC = 0
    LAC = 11400 - 22.66 * T
    LABC = 0
    G_A = 0
    G_B = 4150 - 5.2 * T
    G_C = 19874 - 13.7 * T
    return Gibbs_phase(XA, XB, T, G_A, G_B, G_C, LAB, LBC, LAC, LABC, tag)


# T = 350

#Convex hull method and boundary screening
#Developed by Sritapaswi Nori & Pilsun Yoo & Jiaqi Yang
def cal_ternary(T, tag):
    XA_list = []
    XB_list = []
    XC_list = []
    G_liquid_list = []
    G_alpha_list = []
    G_gamma_list = []
    G_epsilon_list = []
    G_delta_list = []
    G_min_list = []
    for xA_ in range(1, 50):
        XA = xA_ * 0.02
        for xB_ in range(1, 50):
            XB = xB_ * 0.02
            if (XA + XB) < 1:
                G_liquid = Gliquid(XA=XA, XB=XB, T=T, tag=tag)
                G_alpha = Galpha(XA=XA, XB=XB, T=T, tag=tag)
                G_gamma = Ggamma(XA=XA, XB=XB, T=T, tag=tag)
                G_epsilon = Gepsilon(XA=XA, XB=XB, T=T, tag=tag)
                G_delta = Gdelta(XA=XA, XB=XB, T=T, tag=tag)
                XA_list.append(XA)
                XB_list.append(XB)
                XC_list.append(1 - XA - XB)
                G_liquid_list.append(G_liquid)
                G_alpha_list.append(G_alpha)
                G_gamma_list.append(G_gamma)
                G_epsilon_list.append(G_epsilon)
                G_delta_list.append(G_delta)
                G_min_list.append(min([G_liquid, G_alpha, G_gamma, G_epsilon, G_delta]))
    points = np.column_stack([XA_list, XB_list, G_min_list])
    points0 = np.column_stack([XA_list, XB_list])
    points0_t = tuple(tuple(x) for x in points0)

    hull = ConvexHull(points)
    xhullvertices = [points[hull.vertices, 0], points[hull.vertices, 1], points[hull.vertices, 2]]
    sorted_xhullvertices = np.sort(xhullvertices, axis=None)
    data = np.column_stack([points[hull.vertices, 0], points[hull.vertices, 1], points[hull.vertices, 2]])

    fig = plt.figure()
    l = plt.plot(data[:, 0], data[:, 1], 'bo')
    plt.setp(l, markersize=3)

    xa = points[hull.vertices, 0]
    xb = points[hull.vertices, 1]
    points2 = np.column_stack([xa, xb])

    points2_t = tuple(tuple(x) for x in points2)

    com = set(points0_t) - set(points2_t)
    points3 = np.array(list(list(x) for x in com))

    l = plt.plot(points3[:, 0], points3[:, 1], 'yo')
    plt.setp(l, markersize=3)

    # Only work with points2 and points3 to find boundary

    t3 = []

    for i in range(len(points2) - 1):
        if points2[i + 1][1] - points2[i][1] > 0.03:
            t3.append(list(points2[i]))
            t3.append(list(points2[i + 1]))
        if points2[i + 1][1] - points2[i][1] < 0.0 and \
                points2[i][0] + points2[i][1] < 0.98:
            t3.append(list(points2[i]))

    for i in range(1, len(points2)):
        if points2[i][0] + points2[i][1] >= 0.98:
            if points2[i][0] != points2[i - 1][0]:
                t3.append(points2[i])
            else:
                if points2[i][1] - points2[i - 1][1] > 0.03:
                    t3.append(points2[i])

    for i in range(1, len(points2)):
        if points2[i - 1][0] != points2[i][0] and points2[i][1] > 0.02:
            t3.append(points2[i])
    points2_l = list(points2)

    points2_l.sort(key=lambda x: x[1])
    points2_l = np.array(points2_l)
    for i in range(len(points2_l) - 1):
        if points2_l[i + 1][0] - points2_l[i][0] > 0.03:
            t3.append(list(points2_l[i]))
            t3.append(list(points2_l[i + 1]))

    t3 = np.array(t3)

    l = plt.plot(t3[:, 0], t3[:, 1], 'ro')
    plt.setp(l, markersize=3)
    # plt.savefig("test.png")
    t4 = []
    for x, y in t3:
        z = 1 - x - y
        t4.append([x, y, z])
    # t4 = np.array(t4)

    return t4

#Ternary system plotting funcion
#Developed by Pilsun Yoo & Jiaqi Yang
def ternary_main(tag,T):
    t4 = []
    t4 = cal_ternary(T, tag)

    # t4=np.array(t4)
    # %matplotlib inline
    import ternary
    figname = "ternary_model_" + str(tag) +'_'+str(T)+ ".png"
    # print(figname)
    ## Boundary and Gridlines
    scale = 1.0
    figure, tax = ternary.figure(scale=scale)
    figure.set_size_inches(2, 2)
    figure.set_dpi(300)

    # Draw Boundary and Gridlines
    tax.boundary(linewidth=0.1)
    tax.gridlines(color="black", multiple=2)
    tax.gridlines(color="blue", multiple=2, linewidth=0.01)

    # Set Axis labels and Title
    fontsize = 5
    tax.set_title("Simplex Boundary and Gridlines", fontsize=fontsize)
    tax.left_axis_label("Left label $\\alpha^2$", fontsize=fontsize)
    tax.right_axis_label("Right label $\\beta^2$", fontsize=fontsize)
    tax.bottom_axis_label("Bottom label $\\Gamma - \\Omega$", fontsize=fontsize)

    # Set ticks
    tax.ticks(axis='lbr', linewidth=1, fontsize=1)

    # Remove default Matplotlib Axes
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    tax.scatter(t4, marker='o', color='black', label="Test", s=0.1)

    ternary.plt.savefig(figname)
    # ternary.plt.show()


if __name__ == '__main__':
    ternary_main(1,400)
