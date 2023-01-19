from tkinter import *
from tkinter import *
from tkinter.filedialog import *
from tkinter.messagebox import *
from tkinter.font import Font
from tkinter.scrolledtext import *
import file_menu
import edit_menu
import format_menu
import help_menu
from PIL import ImageTk,Image
import sqlite3
import os

root = Tk()
root.title("SW Launcher")
root.geometry("300x250")
root.minsize(width=400, height=400)

#pad inside the frame is the size of the frame
#run_bar = LabelFrame(root,padx=20,pady=2)
#run_bar.grid(row=0,column=0,columnspan=2)

text = ScrolledText(root, state='normal', height=100, width=250, wrap='word', pady=2, padx=3, undo=True)
text.grid(row=1, column=0,sticky=W+E)
text.focus_set()

# creating a menubar
menubar = Menu(root)

# adding our menus to the menubar
file_menu.main(root, text, menubar)
edit_menu.main(root, text, menubar)
format_menu.main(root, text, menubar)
help_menu.main(root, text, menubar)


'''
def show_Out():
    btn = Button(run_bar,text="End",command=root.destroy)
    btn.grid(row=0,column=1)

def run_SWM():
    print("dentro")
    global prueba
    
    path_run = "/home/jorge/M2/Project/Molecular-Liquids/program.out"
    #path_run = "ls"
    #lab = Label(run_bar,text=path_run)
    #lab.grid(row=1,column=0)
    
    os.system(path_run)

    prueba=1

    if prueba==1:
        show_Out()

run_btn = Button(run_bar,text="Press for SWM simulation",command=run_SWM)
#run_btn.grid_columnconfigure(4, weight = 1)
run_btn.grid(row=0,column=0)
'''

root.mainloop()