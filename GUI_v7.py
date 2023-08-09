# GUI design
# Written by Jiaqi Yang

import tkinter as tk
from PIL import Image, ImageTk
from EQ_comp_v7_new import *
import os
from tkinter import ttk
from Ternary import ternary_main

#############constant input###############
R = 8.314


# Picture resize and insert
# This part can be funcionalized later to fits the requirement of whole project
def resize(w_box, h_box, pil_image):
    w, h = pil_image.size
    f1 = 1.0 * w_box / w
    f2 = 1.0 * h_box / h
    factor = min([f1, f2])
    width = int(w * factor)
    height = int(h * factor)
    return pil_image.resize((width, height), Image.ANTIALIAS)


# binary system show function
def show_binary():
    # pic_name = e_1.get() + "_" + e_2.get() + ".png"
    # print(pic_name)
    # Show picture
    w_box = 700
    h_box = 600
    if len(t_1.get()) != 0:
        T1 = int(t_1.get())
    else:
        T1 = 400
    if len(t_2.get()) != 0:
        T2 = int(t_2.get())
    else:
        T2 = 800
    if len(name.get()) != 0:
        e = str(name.get())
    else:
        e = "Ag_Cu"
    figname = e + "_" + str(T1) + "_" + str(T2) + ".png"
    print(figname)
    if not os.path.exists(figname):
        # Show picture
        step = 10
        E_symbol = click_get_element()
        E1 = E_symbol[0]
        E2 = E_symbol[1]
        system = sorted([E1, E2])
        system_string = ""
        for ele in system:
            system_string += ele

        print(system_string)

        Tem, Compositions = Binary_Phase(T1, T2, system_string, step=step)
        # print(Tem)
        # print(Compositions)
        plt.scatter(Compositions, Tem, s=1, c=(0, 0, 0), alpha=1)
        plt.xlabel('Molar Fraction of B')
        plt.xlim(0, 1)
        plt.ylabel('Temperature (°C)')
        figname = E1 + "_" + E2 + "_" + str(T1) + "_" + str(T2) + ".png"
        plt.savefig(figname)

    img_open = Image.open(figname)

    img_resize = resize(w_box, h_box, img_open)
    tk_image = ImageTk.PhotoImage(img_resize)

    L_p = tk.Label(root, image=tk_image, width=w_box, height=h_box)
    L_p.place(x=80, y=150)
    root.mainloop()


# ternary system show function
def show_ternary():
    T=int(t_1.get())
    E_symbol = click_get_element()
    tag = int(E_symbol[2])
    print(tag)
    print(type(tag))
    ternary_main(tag,T)
    figname = "ternary_model_" + str(E_symbol[2]) +'_'+str(T)+ '.png'
    # print(figname)
    w_box = 700
    h_box = 600
    img_open = Image.open(figname)

    img_resize = resize(w_box, h_box, img_open)
    tk_image = ImageTk.PhotoImage(img_resize)

    L_p = tk.Label(root, image=tk_image, width=w_box, height=h_box)
    L_p.place(x=80, y=150)
    root.mainloop()


# Function to get element information

def click_get_element():
    text = "Alloy System is : " + str(name.get())
    sys = name.get()
    var.set(text)
    if len(sys.split("_")) == 2:
        E1 = sys.split("_")[0]
        E2 = sys.split("_")[1]
        E_symbol = [E1, E2]
    elif len(sys.split("_")) == 3:
        E1 = sys.split("_")[0]
        E2 = sys.split("_")[1]
        E3 = sys.split("_")[2]
        E_symbol = [E1, E2, E3]
    return E_symbol


# Show function combined binary and ternary system
def sys_show():
    E_symbol = click_get_element()
    if len(E_symbol) == 2:
        show_binary()
    elif len(E_symbol) == 3:
        show_ternary()


# create the GUI
root = tk.Tk()

root.geometry("1000x800+10+10")
root.title('ProSE')
# Label
L_0 = tk.Label(root, text="ProSE Group2")
L_0.config(font='Helvetica -15 bold', fg='black')
L_0.place(x=700, y=780, anchor='center')

L_a = tk.Label(root, text="Author:\tPilsun,Nori,Jiaqi")
L_a.config(font='Helvetica -12', fg='black')
L_a.place(x=800, y=770)

# Temperature entry

l_t_1 = tk.Label(root, text="Input Starting Temperature").grid(row=3)
l_t_2 = tk.Label(root, text="Input Ending Temperature").grid(row=4)
t_1 = tk.Entry(root)
t_2 = tk.Entry(root)

t_1.grid(row=3, column=1)
t_2.grid(row=4, column=1)

# Element entry by combobox
la_e = tk.Label(root, text='Choose Element').grid(row=1)
var = tk.StringVar()
la = tk.Label(root, textvariable=var)
la.place(x=500, y=60)

name = tk.StringVar()
nameChosen = ttk.Combobox(root, width=12, textvariable=name)
nameChosen['values'] = ("Ag_Cu", "Pd_Ni", "ternary_model_1", "ternary_model_2", "ternary_model_3")
nameChosen.grid(row=1, column=1)
nameChosen.current = 0
# Button to get element information
B_1 = tk.Button(root, text='ELEMENT', command=click_get_element)
B_1.place(x=500, y=20)

# Button
B_text = tk.Label(root, text="Push this button to calculate phase diagram")
B_text.place(x=500, y=0)

B_0 = tk.Button(root, text="TEST", command=sys_show)
B_0.place(x=600, y=20)

# Menu construction
menuBar = tk.Menu(root)
root.config(menu=menuBar)


# Menu function construction
def _quit():
    root.quit()
    root.destroy()
    exit()


# Menu option setting
op_file = tk.Menu(menuBar, tearoff=0)
menuBar.add_cascade(label='File', menu=op_file)
op_file.add_command(label='Exit', command=_quit)

# Enter the main loop of tk
root.mainloop()
