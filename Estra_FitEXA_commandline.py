from Tkinter import *
from tkFileDialog import *
import subprocess, os
my_env = os.environ.copy()
my_env["PATH"] = os.getcwd()+';' + my_env["PATH"]
root = Tk()
root.withdraw() # hide root
path = askdirectory(initialdir='%HOMEPATH%')
os.chdir(path)
subprocess.call(['cmd'], env = my_env)



