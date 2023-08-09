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


# Function for Tab 1
def tab1_function():
    print("Tab 1 Function")
root = tk.Tk()

# Create a Notebook widget
notebook = ttk.Notebook(root)
notebook.grid(row=0, column=0, sticky="nsew")
# Create Tab 1
tab1 = ttk.Frame(notebook)
notebook.add(tab1, text="Tab 1")

# Add content to Tab 1
tab1_label = tk.Label(tab1, text="Tab 1 Content")
tab1_label.grid(row=0, column=0)
tab1_button = tk.Button(tab1, text="Tab 1 Button", command=tab1_function)
tab1_button.grid(row=1, column=0)

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




# Define the function to handle element click event
def element_clicked(symbol):
    print("Element clicked:", symbol)


# Load the periodic table image
image = Image.open("C:/Users/Ridhi Tamirasa/Downloads/GUI_final_1/GUI_final_1/periodic_table.png")
# Resize the image to a smaller size
new_width = 200
new_height = 200
image = image.resize((new_width, new_height))

photo = ImageTk.PhotoImage(image)

# Create a Canvas widget
canvas = tk.Canvas(root, width=image.width, height=image.height)
canvas.grid()

# Display the periodic table image
canvas.create_image(0, 0, image=photo, anchor=tk.NW)
image_label = tk.Label(root, image=photo)
image_label.grid()

# Define the positions and colors of element tiles
element_positions = {
    "H": {"position": (90, 82), "color": "lightblue"},
    "He": {"position": (78, 82), "color": "lightgreen"},
    "Li": {"position": (46, 114), "color": "pink"},
    "Be": {"position": (78, 114), "color": "lightyellow"},
    # Add more elements and their positions...
    "Ar": {"position": (204, 210), "color": "lightcoral"}
}

# Create buttons for each element at their respective positions with colors
for symbol, info in element_positions.items():
    position = info["position"]
    color = info["color"]
    x, y = position
    element_label = tk.Label(root, text=symbol, bg="lightblue")
    element_label.place(x=x, y=y, width=32, height=32)
    element_label.bind("<Button-1>", lambda event, s=symbol: element_clicked(s))
#Start the Tkinter event loop

# Enter the main loop of tk

# Function for Tab 2
def tab2_function():
    print("Tab 2 Function")
root.title("Multiple Tabs Example")

# Create Tab 2
tab2 = ttk.Frame(notebook)
notebook.add(tab2, text="Tab 2")

# Add content to Tab 2
tab2_label = tk.Label(tab2, text="Tab 2 Content")
tab2_label.grid(row=0, column=0)
tab2_button = tk.Button(tab2, text="Tab 2 Button", command=tab2_function)
tab2_button.grid(row=1, column=0)

root.mainloop() 


