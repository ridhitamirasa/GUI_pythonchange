import tkinter as tk
from tkinter import ttk

root = tk.Tk()

# Create two frames
frame1 = ttk.Frame(root)
frame1.pack()

frame2 = ttk.Frame(root)
frame2.pack()

# Add content to frame 1
label1 = ttk.Label(frame1, text="Content in Frame 1")
label1.pack()

# Add content to frame 2
label2 = ttk.Label(frame2, text="Content in Frame 2")
label2.pack()

root.mainloop()
