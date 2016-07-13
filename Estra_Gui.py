# A cx_freeze setup script to create PyMca executables
#
# Use "python cx_setup.py install"
#
# It expects a properly configured compiler.
#
# Under windows you may need to set MINGW = True (untested) if you are
# not using VS2003 (python 2.5) or VS2008 (python 2.6)
#
# If everything works well one should find a directory in the build
# directory that contains the files needed to run the PyMca without Python

import text
import sys
import os
import ConfigParser
from Tkinter import *
from ttk import *
import tkFileDialog

from numpy import savetxt      as np_savetxt
from numpy import loadtxt      as np_loadtxt
from numpy import array        as np_array
from numpy import column_stack as np_column_stack 
from numpy import arange as np_arange 
from numpy import zeros        as np_zeros
from numpy import log        as np_log
from numpy import sqrt        as np_sqrt
from numpy import sign        as np_sign
from numpy import hsplit      as np_hsplit
from numpy import argmax      as np_argmax
from numpy import gradient     as np_gradient
from numpy import interp     as np_interp
from numpy import linspace    as np_linspace
from numpy import genfromtxt    as np_genfromtxt

global __verbose__                                                                    
__verbose__=False#True
global __version__
__version__= "a0.0.0"
global inivar
inivar=ConfigParser.ConfigParser()

import warnings

global estra_output
global paste




class input_exafs():
    """
    """
    def __init__(self, key_dict={}):
        self.key_dict=key_dict
    def read_inp(self,filename):
        self.key_dict["filename"]=filename.next().rstrip()
        self.key_dict["colums"]=filename.next().split("#")[0].split()
        self.key_dict["pre_edge"]=map(float,filename.next().split("#")[0].split())
        self.key_dict["edge"]=map(float,filename.next().split("#")[0].split())
        self.key_dict["glitches_n"]=int(filename.next().split("#")[0])
        #print "number of glitches" ,self.key_dict["glitches_n"]
        if self.key_dict["glitches_n"]!=0:
           print "leggo i glich" 
           self.key_dict["glitches_lim"] =[]
           for i in range(self.key_dict["glitches_n"]):
               stringa=filename.next().replace(","," ").split("#")[0]
               self.key_dict["glitches_lim"].append(map(
                                float,stringa.split()))
        #print "\n",  self.key_dict["glitches_lim"]     
        self.key_dict["splines_n"]=int(filename.next().split("#")[0]) 
        #print "\n\n\n spline_n",  self.key_dict["splines_n"]
        self.key_dict["splines_int"]=[]
        buffero=map(float,filename.next().split("#")[0].split())
        if len(buffero)==self.key_dict["splines_n"]+1:
            self.key_dict["splines_int"]=buffero
        else:
            self.key_dict["splines_int"].extend(buffero)
            for item in range(self.key_dict["splines_n"]):
                self.key_dict["splines_int"].append(float(filename.next().split("#")[0]))
        #print self.key_dict["splines_int"]        
        self.key_dict["splines_degree"]=map(int, filename.next().split("#")[0].split())   
        self.key_dict["sec_edge_n"]=int(filename.next().split("#")[0])
        if self.key_dict["sec_edge_n"]!=0:
            self.key_dict["sec_edge_lim"]=[] 
            for i in range(self.key_dict["sec_edge_n"]):   
                self.key_dict["sec_edge_lim"].append(map(
                                float,filename.next().split("#")[0].split())) 
        self.key_dict["lim_Jump_eval"]=map(float, filename.next().split("#")[0].split())   
        self.key_dict["Const_Edge_eval"]= filename.next().split("#")[0].split()
        self.key_dict["Const_Edge_eval"][0]=str(self.key_dict["Const_Edge_eval"][0])
        try:
            self.key_dict["Const_Edge_eval"][1]=str(self.key_dict["Const_Edge_eval"][1])   
        except IndexError:    
            self.key_dict["Const_Edge_eval"].append("E")
        try:    
            self.key_dict["Const_Edge_eval"][2]=int(float(self.key_dict["Const_Edge_eval"][2]))
        except IndexError:
            self.key_dict["Const_Edge_eval"].append(0)
        #print "\n",  self.key_dict["Const_Edge_eval"]
        if self.key_dict["Const_Edge_eval"][0]!="n":
            self.key_dict["Edge_jump"]=float(filename.next().split("#")[0].rstrip())
  
            
    def write_inp(self,filename) :
        if self.key_dict["pre_edge"][0]>=self.key_dict["pre_edge"][1]:
            print "\n"*5+"ERROR-----"*10+"\n\nERROR  pre_edge lim1 >= pre_edge lim2\n\n"+"ERROR-----"*10+"\n"*5
            raise ValueError() 
            return
        if self.key_dict["pre_edge"][1]>=abs(self.key_dict["edge"][0]):
            print "\n"*5+"ERROR-----"*10+"\n\nERROR  pre_edge lim2 >= edge_max\n\n"+"ERROR-----"*10+"\n"*5
            raise ValueError() 
            return
        filename.write(self.key_dict["filename"]+ "\n")
        filename.write("{0:s}   {1:s}".format(*self.key_dict["colums"]).ljust(30)+"# columns for E (eV/keV) and alpha?\n")
        filename.write("{0:5.1f}    {1:5.1f}".format(*self.key_dict["pre_edge"]).ljust(30)+"# limits (eV) for PRE EDGE (linear)fit\n")
        filename.write("{0:5.1f}    {1:1.2f}".format(*self.key_dict["edge"]).ljust(30)+"# Max. energy and tolerance for edge search\n")
        filename.write("{0:d}".format(self.key_dict["glitches_n"]).ljust(30)+"# Number of glitches\n")   
        if self.key_dict["glitches_n"]!=0:
           for i in range(self.key_dict["glitches_n"]):  
               filename.write("{0:5.1f}   {1:5.1f}".format(*self.key_dict["glitches_lim"][i][:2]).ljust(30)+"#glitch {0:d}  limits \n".format(i))  
               
        filename.write("{0:d}".format(self.key_dict["splines_n"]).ljust(30)+"# Number of spline intervals\n")   
        for item in self.key_dict["splines_int"]: filename.write("{0:5.2f}\n".format(item))
        #filename.write("      #spline intervals\n")
        buffero=""
        for item in self.key_dict["splines_degree"]:
            buffero+="{0:d}  ".format(item)    
        filename.write(buffero.ljust(30)+"# degrees for splines\n")
        filename.write("{0:d}".format(self.key_dict["sec_edge_n"]).ljust(30)+"# Number of secondary edges Ede, Ade, Wde\n")
        
        if self.key_dict["sec_edge_n"]!=0:
           for i in range(self.key_dict["sec_edge_n"]):  
               filename.write("{1:5.2f}   {2:5.3f}   {3:5.3f}".format(*self.key_dict["sec_edge_lim"][i][:3]).ljust(30)+
                              "#sec edge {0:d}  Ede, Ade, Wde \n".format(i))            
        filename.write("{0:5.1f}    {1:1.2f}".format(*self.key_dict["lim_Jump_eval"]).ljust(30)+"# Limits for Jump evaluation (eV/k)\n") 
        filename.write("{0:s}  {1:s}  {2:d}".format(*self.key_dict["Const_Edge_eval"]).ljust(30)+"# constraints at the edge (n/J  E/k  wgt)\n")      
        if self.key_dict["Const_Edge_eval"][0]!='n':
            filename.write("{0:1.5f}".format(self.key_dict["Edge_jump"]).ljust(30)+"# Edge Jump\n")
        filename.write("# ******  EXAFS done **********\n")    
            
         
        
class input_fourier():
    """
    """
    def __init__(self, key_dict={}):
        self.key_dict=key_dict
    def read_inp(self,filename):
        self.key_dict["k_lim"]=map(float,filename.next().split("#")[0].split())
        self.key_dict["k_w"]=int(filename.next().split("#")[0])
        self.key_dict["window"]=filename.next().split("#")[0].strip()
        self.key_dict["window_apodization"]=float(filename.next().split("#")[0].strip())
        self.key_dict["R_lim"]=map(float,filename.next().split("#")[0].strip().split())
        
    def write_inp(self,filename) :
        filename.write("{0:5.1f}  {1:5.1f}".format(*self.key_dict["k_lim"]).ljust(30)+"# K_min, K_max for the FT\n")  
        filename.write("{0:d}".format(self.key_dict["k_w"]).ljust(30)+"# K-Weight for FT\n")          
        filename.write("{0:s}".format(self.key_dict["window"]).ljust(30)+"# window: [G]auss or [H]anning or [N]one\n") 
        filename.write("{0:5.1f}".format(self.key_dict["window_apodization"]).ljust(30)+"# Gaussian width\n")    
        filename.write("{0:5.1f}  {1:5.5f}".format(*self.key_dict["R_lim"]).ljust(30)+"# R max and DR for FT\n")
        filename.write("# ******  Fourier done **********\n")            
  
       
        
class input_back_fourier():
    """
    """
    def __init__(self, key_dict={}):
        self.key_dict=key_dict
    def read_inp(self,filename):
        self.key_dict["R_lim"]=map(float,filename.next().split("#")[0].split())  
        self.key_dict["R_exp"]=float(filename.next().split("#")[0])
    def write_inp(self,filename) :
        filename.write("{0:5.2f}  {1:5.2f}".format(*self.key_dict["R_lim"]).ljust(30)+"# BF limits\n")  
        filename.write("{0:5.2f}".format(self.key_dict["R_exp"]).ljust(30)+"# expected <R> for the BF\n")          
        filename.write("# ******  Back Fourier done **********\n")  



class class_estra_input():
    """
    classe per l'imput file di Carlo
    """
    def __init__(self,  input_filename, label="A10"):
        self.input_filename =input_filename
        self.label =label 
        pass
    
    def read_input(self):
        with open(self.input_filename,'r') as realfile:
            filename=realfile.readlines()
            realfile.close()
        for i,item in enumerate(filename): filename[i].replace(","," ")
        
        filename=iter(filename)

                
        self.label=filename.next().split("#")[0].strip()
        for i in range(50):
            o=filename.next().split("#")[0].strip()
            if o.strip()=="e":
                #print "###############  exafs input"
                self.inp_exa=input_exafs()
                self.inp_exa.read_inp(filename)
            elif  o.strip()=="f" : 
                # "###############  fourier input"
                self.inp_fou=input_fourier()
                self.inp_fou.read_inp(filename)  
            elif  o.strip()=="b":  
                #print "###############  back fourier input"
                self.inp_bfou=input_back_fourier()
                self.inp_bfou.read_inp(filename)
            elif  o.strip()=="q":
                break
        if os.access(self.inp_exa.key_dict["filename"], os.F_OK):
            self.filedata=self.inp_exa.key_dict["filename"]
        else :self.filedata=""
          

    def write_input(self,filename="estra.inp"):
        with open(filename,'w') as filename:
            filename.write("{0:s}".format(self.label).ljust(30)+"# SAMPLE \n")          
            filename.write("e".ljust(30)+"# EXAFS option\n")         
            self.inp_exa.write_inp(filename)
            filename.write("f".format(self.label).ljust(30)+"# fourier option\n")                         
            self.inp_fou.write_inp(filename)
            filename.write("b".format(self.label).ljust(30)+"# back fourier option\n")                         
            self.inp_bfou.write_inp(filename)
            filename.write("p".format(self.label).ljust(30)+"# gnuplot macro\n")        
            filename.write("q".format(self.label).ljust(30)+"# quit\n")
            filename.close()
            
    def default_input(self, Datafile,E,Mu,label=None):
        e2k=lambda x: 0.512317 * np_sqrt(x - E_edge)
        E = E*1000 if E[0]<90 else E 
        E_edge=E[np_argmax(np_gradient(Mu))] 
        if label:
            self.label=label
        ##### EXAFS DICT ########################################
        key_dict={}
        key_dict["filename"]= Datafile
        key_dict["colums"]=["1" ,"2"]
        key_dict["pre_edge"]=[E[0], E_edge-30]
        key_dict["edge"]=[E_edge, 0.7]
        key_dict["glitches_n"]=0
        key_dict["splines_n"]= int((e2k(E[-1])-3)//3)
        roundmore=lambda x: round(x,2)
        key_dict["splines_int"]=map(roundmore,np_linspace(3,e2k(E[-1]),key_dict["splines_n"]+1))
        key_dict["splines_degree"]=[3]*key_dict["splines_n"]
        key_dict["sec_edge_n"]=0
        key_dict["lim_Jump_eval"]=[3, 6]
        key_dict["Const_Edge_eval"]=["n", "E",  0]
        self.inp_exa=input_exafs(key_dict)
        ##### Fourier DICT ########################################
        key_dict={}
        fin =self.inp_exa.key_dict["splines_int"][-1] if self.inp_exa.key_dict[
                                   "splines_int"][-1]<90 else e2k(self.inp_exa.key_dict["splines_int"][-1])
        key_dict["k_lim"]=[3.0, fin]
        key_dict["k_w"]=3
        key_dict["window"]="g"
        key_dict["window_apodization"]=0.2
        key_dict["R_lim"]=[6.0, .01]
        self.inp_fou=input_fourier(key_dict)
        ##### back Fourier DICT ####################################
        key_dict={}
        key_dict["R_lim"]=[1.0, 2.0]
        key_dict["R_exp"]=1.80
        self.inp_bfou=input_back_fourier(key_dict)
        
    def default_spline_intervals(self,n_splines):                      
        """serve a ricalcolare gli intervalli equispaziati quando cambi il numero di spline
        """
        roundmore=lambda x: round(x,2)
        self.inp_exa.key_dict["splines_n"]=n_splines
        self.inp_exa.key_dict["splines_int"]=map(roundmore,np_linspace(self.inp_exa.key_dict["splines_int"][0],
                                                                       self.inp_exa.key_dict["splines_int"][-1],
                                                                       self.inp_exa.key_dict["splines_n"]+1))
        self.inp_exa.key_dict["splines_degree"]=[3]*self.inp_exa.key_dict["splines_n"]

        
        
        
        
        

################################################################################

class class_estra_output():
    def __init__(self,  radix ):
        self.radix =radix 
        pass    
    def read_output(self):
        vaff=self.radix+".pre"
        with open(vaff,'r') as filename:
            while True:
                line=filename.readline()
                if "Edge Energy" in line :
                    self.edge= float(line.split(":")[1])
                elif "White line energy" in line:    
                    self.wline= float(line.split(":")[1])
                elif "#   1_E" in line:
                    break
            data= np_loadtxt(filename)
            filename.close
        self.Ep, self.ABSpre, self.der1, self.der2, self.noiseE=np_hsplit(data,5)
        #------------------------------------------------------------------
        vaff=self.radix+".jmp"
        #os.system("sed \'s/\-0\./ \-0\./g\' " +vaff+" >"+vaff+"o")
        with open(vaff,'r') as filename:
            data= np_loadtxt(filename, usecols=(0,1,2))
            self.Ej, self.A_nor, self.A_0,  =np_hsplit(data,3)
            filename.seek(0)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                data=np_genfromtxt(filename, invalid_raise=False)
            self.Ejsp,a,b, self.SPL_nor, self.SPL_ME_nor = np_hsplit(data,5)
            filename.close
        # perche Ep e' differente da Ej?????
        #------------------------------------------------------------------
        vaff=self.radix+".sta"
        with open(vaff,'r') as filename:
            while True:
                line=filename.readline()
                if "sigma^2 Noise" in line :
                    self.s2= float(line.split(":")[1])
                elif "sigma^2 k*Noise" in line:    
                    self.s2k= float(line.split(":")[1])
                elif "K              data" in line:
                    break
            data= np_loadtxt(filename)    
            filename.close
        self.Ks, self.der, self.Discont, self.Spectrum_noise, self.Average_noise = np_hsplit(data,5)
        #------------------------------------------------------------------
        vaff=self.radix+".exa"
        with open(vaff,'r') as filename:
            while True:
                line=filename.readline()
                if "# Jump" in line :
                    self.Jump= float(line.split(":")[1])
                elif "EWL" in line:    
                    self.EWL= float(line.split(":")[1])
                elif "hwl" in line:    
                    self.hwl= float(line.split(")")[1])                   
                elif "#     k " in line:
                    break
            data= np_loadtxt(filename, usecols=(0,1))
            filename.close
        self.k, self.kchi = np_hsplit(data,2)
        #------------------------------------------------------------------  
        vaff=self.radix+".fou"
        with open(vaff,'r') as filename:
            data= np_loadtxt(filename)
            filename.close
        self.R, self.Amp, self.Re, self.Im, self.kf, self.knchi, self.wind = np_hsplit(data,7)
        #------------------------------------------------------------------    
        vaff=self.radix+".bf"
        with open(vaff,'r') as filename:
            data= np_loadtxt(filename)
            filename.close
        self.q, self.qchi, self.qnchi, self.qAmplitude, self.qAmpR2, self.Phase, self.Phase2kr=np_hsplit(data,7)
        #------------------------------------------------------------------ 
        
        
  




################################################################################
################################################################################
################################################################################
#######                 GRAPHYCAL INTERFACE                    #################
################################################################################
################################################################################

####################        Utility          ###################################
def string_range(string):
    """define a listi fo selected number starting from 1
       the character allowed are , to define different ranges 
       and - to define a ranges
    """
    select=[]
    ranges=string.split()
    for item in ranges:
        if item.find('-')<0:
            select.append(int(item)-1)
        else:
            l1,l2=map(int,item.split('-'))
            select.extend(range(l1-1,l2))
    return select 
def datalize(*arg):
    return np_column_stack((arg))
    #return transpose((vstack((arg))))

def filewrite(filename,  newdata, comment=None, footers=None,fmt='%1.10f '):
    	"""function to write a bm29 file, define a series of array with the same name
    	of column defined in the file
        """
        if os.path.exists(filename):
        	filename += ".1"
        outFile = open(filename, 'w')
        if comment is None:
		   np_savetxt(outFile, newdata, fmt=fmt)
        else:
           outFile.writelines(comment)
           np_savetxt(outFile, newdata, fmt=fmt)#%1.10f %1.8f %1.8f %d  %d %d %d %d 
        if 	not(footers is None):
            outFile.writelines(footers) 
        outFile.close
        return


def browse_single():
        filename = tkFileDialog.askopenfilename()   ###-defaultextension extension
        return filename

        
def browse_multiple():
        filenames = tkFileDialog.askopenfilenames()   ###-defaultextension extension
        if "}" in filenames:
            filenames=filenames.split("}")
            filenames=filenames[:-1]
            filenames=[item.replace("{","") for item in filenames]
            filenames=[item.strip() for item in filenames]
        #else :filenames=filenames.split(" ")
        #filenames=filenames.split()
        #filenames = Tk._default_root.tk.splitlist(filenames)
        #filenames = sorted(filenames)
        return filenames        


####################Class used to import data###################################
class bm29file():
    """classe file di BM29 dovrebbe aprire un file leggerlo leggere le
    informazioni dei commenti. ha delle funzioni per il calcolo della derivata prima
    e dell angolo del monocromatore per ogni punto di enegia \n
    datinput: should be or a string or an numpy array with first column an Energy(eV or keV) 
    and as second column a Mu
    datinput =could be a strin with filenames or a numpy array or  a list of two array
    All_Column= True all column, False 2 column, Minimum E Mu I0 Ref 
    """
    def __init__(self,  datinput,  All_Column=True):
          self.All_Column=All_Column
          tipo =str(type(datinput)).split("'")[-2] 
          if tipo == "list":
              if all(datinput[0]<90): 
                  self.E =np_array(datinput[0])*1000 
              else: 
                  self.E = np_array(datinput[0])
              self.Mu=np_array(datinput[1])
              self.All_Column = False

          else:
              print "porca troia"
 




####################Class used to import data###################################
class Col_line_Gen:
    def __init__(self, genitore, label, array, row): 
      #-----------------------------      Declare      --------------------------------------------------
        self._check= IntVar()
        self._label=label
        self._position= IntVar()
        self._compo=IntVar()

      #-----------------------------      Define       --------------------------------------------------  
        self.array=array 
        row= row
        self._check.set(False)
      #-----------------------------      Structure    --------------------------------------------------
        self.check = Checkbutton(genitore, variable=self._check, width=15, text=self._label)#, command=self.underline )
        self.check.grid(row= row, column=0)
        
        
        self.column= Spinbox(genitore, from_ = 1, to = self.array.shape[1], 
                             textvariable= self._position, width = 3)      
        self.column.grid(row= row, column=2)
        
        self.pulsante_Plot = Button(genitore ,
                                      command = self.plot,
                                      text = "Plot",
                                      width = 7)
        self.pulsante_Plot.grid(row= row, column=3)   
        
      #-----------------------------      Functios    --------------------------------------------------
    def plot(self): 
        title= "{0:s}  column {1:d}".format(self._label, int(self._position.get()))
        x_array=[np_arange(self.array.shape[0])]
        y_array=[self.array[:,int(self._position.get())-1]]
        graph = Graph()
        graph.plot(x_array, y_array, title= title)

class Col_line_Gen_2(Col_line_Gen):
    def __init__(self, genitore, label, array, row): 
      #-----------------------------      Declare      --------------------------------------------------
        self._check= IntVar()
        self._label=label
        self._position= StringVar()
        self._compo=IntVar()

      #-----------------------------      Define       --------------------------------------------------  
        self.array=array 
        self._position.set("4,5,6-9")
        row= row
        self._check.set(False)
      #-----------------------------      Structure    --------------------------------------------------
        self.check = Checkbutton(genitore, variable=self._check, width=15, text=self._label)#, command=self.underline )
        self.check.grid(row= row, column=0)
        
        self.column= Entry(genitore, textvariable= self._position, width = 8)      
        self.column.grid(row= row, column=2)
        
        self.pulsante_Plot = Button(genitore ,
                                      command = self.plot,
                                      text = "Plot",
                                      width = 7)
        self.pulsante_Plot.grid(row= row, column=3)   
      #-----------------------------      Functios    --------------------------------------------------
    def plot(self): 
        title= "{0:s}  column {1:d}".format(self._label, int(self._position.get()))
        x_array=[np_arange(self.array.shape[0])]
        index= string_range(self._position.get())
        y_array=np_zeros(self.array.shape[0])
        for item in index:
            y_array+=self.array[:,item-1]
        y_array=[y_array]  
        graph = Graph()
        graph.plot(x_array, y_array, title= title)













class Column_Window:
    def __init__(self,filenames):
        global estra_input
      #--------------------------   Declare-------------------------------------------------
        self.filenames=filenames
        print type(self.filenames)
        self._ChaCom=StringVar()
        self._mode=StringVar()
      #--------------------------   Define--------------------------------------------------  
        self.column_names=["E", "Mu", "Ref", "I0", "I1", "I2"]
        self.column_list=[]
        self._ChaCom.set("#")
        self._mode.set("transmission")
      #--------------------------  Top level + reload--------------------------------------------------        
        #self.top_txt = Toplevel()
        #self.top_txt.title("view  first file")         

        
        s = Style()
        s.configure('Green.TButton',  background="green")

        
        self.top = Toplevel(takefocus=True)
        self.top.title("define column of interest")    

        


        self.win_text = Frame(self.top) 
        self.win_text.pack(side=LEFT,expand=YES, fill=BOTH)
        text.ScrolledText( parent=self.win_text, file=self.filenames[0],hor=True, active=False)

        self.top_con= Frame(self.top) 
        self.top_con.pack(side=LEFT,expand=YES)

        self.win_comment = Frame(self.top_con)
        self.win_comment.pack(side=TOP,expand=YES)
        self.win_comment.focus()

        Label(self.win_comment,text="Comment character").grid(row= 0, column=0) 
        self.Entry_Value =Entry(self.win_comment, textvariable =  self._ChaCom, width=8, justify=CENTER )
        self.Entry_Value.grid(row= 0, column=1) 
        self.Reload_B = Button(self.win_comment, text="Reload array" ,
                                      command = self.load,    
                                      width = 13 ).grid(row= 0, column=2)
      #--------------------------   Params  Entries--------------------------------------------------
        self.win_column = Frame(self.top_con)
        self.win_column.pack(side=TOP,expand=YES)
      #--------------------------   Mode array --------------------------------------------------        

        self.quadro_mode = Frame(self.win_column)
        self.quadro_mode.pack(side = TOP, anchor=W, fill = X, pady= 0, ipadx = 0, ipady = 0, expand = N)
        
        #Label(self.quadro_mode,text="Mode    ").pack(side = LEFT)
        self.combo_type= Combobox(self.quadro_mode , textvariable=self._mode, state ="readonly",
                    values=("transmission", "fluorescence"))
        self.combo_type.pack(side = LEFT) 
 
        
        self.quadro_column = Frame(self.win_column)
        self.quadro_column.pack(side = TOP,  fill = X)
        
      #--------------------------   Header-------------------------------------------------- 
      
        Label(self.quadro_column, text="Use       ").grid(row=0, column=0) 
        Label(self.quadro_column, text="   ").grid(row=0, column=1)        
        Label(self.quadro_column, text="Column   ").grid(row=0, column=2)
        Label(self.quadro_column, text="    ").grid(row=0, column=3)  
        Label(self.quadro_column, text="   ").grid(row=0, column=4)  
      #--------------------------   Set lines --------------------------------------------------   
        try:self.load()
        except AttributeError: pass

      #--------------------------   Button --------------------------------------------------
        self.quadro_buttonp = Frame(self.win_column)      
        self.quadro_buttonp.pack(side = TOP,  fill = X) 
        self.buttonMu = Button(self.quadro_buttonp,
                                      command = self.Muplot,
                                      text = "plot Mu",
                                      width = 13)

        self.buttonMu.pack(side = LEFT, anchor = W,pady = 10, padx=5)
        self.buttonRef = Button(self.quadro_buttonp,
                                      command = self.Refplot,
                                      text = "plot Ref",
                                      width = 13)
        self.buttonRef.pack(side = LEFT, anchor = W,pady = 10, padx=5) 
      #--------------------------   Button --------------------------------------------------        
        self.quadro_button = Frame(self.win_column)      
        self.quadro_button.pack(side = TOP,  fill = X) 
        self.save_cons = Button(self.quadro_button,
                                      command = self.opens,
                                      text = "Define Column",
                                      style='Green.TButton',
                                      width = 13)

        self.save_cons.pack(side = LEFT, anchor = W,pady = 10, padx=5) 
      #--------------------------   Wait --------------------------------------------------    
        self.top_con.wait_window()
      #--------------------------   Function --------------------------------------------------
    
    def read_comment(self,n_file):
        infile=open(n_file)
        comment=[]
        CC="1"
        while CC==self._ChaCom.get():
            comment.append(infile.readline())
            CC=comment[-1][0]
        comment=comment[:-1]
        comment=comment if comment!=[] else ["# Generic\n"]
        comment.append("# spectra  "+n_file+ "\n")
        comment.append("#  ---------------------------------"+ "\n")
        comment.append("#L E  Mu"+ "\n") 
        return comment         
            
    

        
    def load(self):
        if self.column_list==[]:
            try:
                self.array=np_loadtxt(fname=self.filenames[0],comments=self._ChaCom.get())
            except :
                print  "\n"*5+"ERROR-----"*10+"\n\nERROR  change comment character\n\n"+"ERROR-----"*10+"\n"*5 
                return
                
            for i,item in enumerate(self.column_names):
                    self.column_list.append(Col_line_Gen(self.quadro_column, label=item, array=self.array,
                                         row=i+1))     
            self.column_list.append(Col_line_Gen_2(self.quadro_column, label="Detectors Sum",
                                      array=self.array, row=i+2))     
        return

    def Muplot(self):
        E=self.array[:,int(self.column_list[0]._position.get())-1]
        if self.column_list[1]._check.get():
            Mu=self.array[:,int(self.column_list[1]._position.get())-1]
        elif self.column_list[3]._check.get() and self.column_list[4]._check.get():
            I0=self.array[:,int(self.column_list[3]._position.get())-1] 
            I1=self.array[:,int(self.column_list[4]._position.get())-1]           
            if self._mode.get()=="transmission":Mu=np_log(I0/I1)
            elif self._mode.get()=="fluorescence":Mu=I1/I0
        elif self.column_list[5]._check.get() and self._mode.get()=="fluorescence":   
            index= string_range(self._position.get())
            I0=self.array[:,int(self.column_list[3]._position.get())-1] 
            I1=np_zeros(self.array.shape[0])
            for item in index:
                I1+=self.array[:][item]
            Mu=I1/I0
        else: raise ValueError("\n\nneither Mu nor I0-I1 defined\n\n") 
        graph = Graph()
        graph.plot( [E],[Mu], title= "Mu")
        
   
        
    def Refplot(self):    
        E=self.array[:,self.column_list[0]._position.get()-1]
        if self.column_list[2]._check.get():
            Ref=self.array[:,int(self.column_list[2]._position.get())-1] 
        elif self.column_list[5]._check.get() and self.column_list[4]._check.get():
            I1=self.array[:,int(self.column_list[4]._position.get())-1] 
            I2=self.array[:,int(self.column_list[5]._position.get())-1]
            Ref=np_log(I1/I2)
        graph = Graph()
        graph.plot( [E], [Ref], title= "Ref")


    def opens(self):
        global estra_input
        filesel_spectra=[]
        for item in self.filenames:
            self.array=np_loadtxt(fname=item,comments=self._ChaCom.get())
            
            E_col=self.column_list[0]._position.get()-1
            comment_Ecol=str(E_col+1)
            self.array[:,E_col].argsort()
            #self.array=self.array[self.array[:,E_col].argsort()]
            E=self.array[:,E_col]
            if self.column_list[1]._check.get():
                Mu=self.array[:,int(self.column_list[1]._position.get())-1]
                comment_I1= str(self.column_list[1]._position.get())+" for Mu"
            elif self.column_list[3]._check.get() and self.column_list[4]._check.get():
                I0=self.array[:,int(self.column_list[3]._position.get())-1] 
                I1=self.array[:,int(self.column_list[4]._position.get())-1]           
                if self._mode.get=="transmission":Mu=np_log(I0/I1)
                elif self._mode.get=="fluorescence":Mu=I1/I0
                comment_I1= str(int(self.column_list[3]._position.get()))+ " for I0"+str(
                                    int(self.column_list[3]._position.get()))+ " for I1"+self._mode.get()   
            elif self.column_list[5]._check.get() and self._mode.get()=="fluorescence":   
                index= string_range(self._position.get())
                I0=self.array[:,int(self.column_list[3]._position.get())-1] 
                I1=np_zeros(self.array.shape[0])
                for item in index:
                    I1+=self.array[:][item]
                Mu=I1/I 
                comment_I1= str(int(self.column_list[3]._position.get())-1)+ " for I0"+str(
                                int(self._position.get())-1)+ " for I1"+self._mode.get()  
            else: 
                stringa="\n"*5+"ERROR-----"*10+"\n\nERROR  neither Mu nor I0-I1 defined\n\n"+"ERROR-----"*10+"\n"*5
                raise ValueError(stringa) 
            
            filesel_spectra.append(bm29file([E,Mu]))
            filesel_spectra[-1].comments=self.read_comment(item)

            if self.column_list[3]._check.get():
                filesel_spectra[-1].I0=self.array[:,int(self.column_list[3]._position.get())-1]
            if self.column_list[2]._check.get():
                filesel_spectra[-1].ref=self.array[:,int(self.column_list[2]._position.get())-1]   
            elif self.column_list[5]._check.get() and self.column_list[4]._check.get():
                I2=self.array[:,int(self.column_list[5]._position.get())-1] 
                I1=self.array[:,int(self.column_list[4]._position.get())-1]
                filesel_spectra[-1].ref=np_log(I1/I2)
            
        pass
        ##### average the files###########################################################
        try:
           E= np_zeros(filesel_spectra[0].E.shape[0])
           for item in  filesel_spectra:
               E+= item.E
           E/=len(filesel_spectra)
        except(IndexError):
            print "error in averaging data"
            E=np_zeros(filesel_spectra[0].E.shape[0])
        Mu= np_zeros(len(E))  
        for item in filesel_spectra:
            Mu+=np_interp(E, item.E, item.Mu)
        Mu/=len(filesel_spectra)    
        comment="# Mu file produced by Estra_GUI\n"
        comment+="# using columns %s , %s\n"%(comment_Ecol, comment_I1)
        comment+="# for files:\n# "
        comment+="\n# ".join(self.filenames)
        comment+="\n# E Mu\n"
        radix = tkFileDialog.asksaveasfilename(title="Opendatafile" ,filetypes=[("No extension","*")]) 
        radix=radix+".abs"
        data=datalize(E,Mu)
        filewrite(radix,  data, comment)
        os.chdir(os.path.dirname(radix))
        
        
        if globals().has_key("estra_input"):
            if hasattr(estra_input,"label"):            
               radix =os.path.basename(radix)
               estra_input.label=radix[:-4]
               estra_input.inp_exa.key_dict["filename"]=radix
               estra_input.inp_exa.key_dict["colums"]=["1","2"]
            print "none"
        else:
            radix =os.path.basename(radix)
            estra_input=class_estra_input("estra.inp", label=radix[:-4])  
            estra_input.default_input(radix,E,Mu,label=None)
        self.top.destroy()
        del self
      

            
            
                
                
      
      

################################################################################
global paste
import matplotlib
matplotlib.interactive(False)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, cursors
class Graph:
    def __init__(self,title=None):
        #self.sh = Tk.StringVar()
        self.top = Toplevel()
        self.top.title(title)
        #self.top.protocol("WM_DELETE_WINDOW", self.topcallback)
        self.fig = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.top)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,  self.top )
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=LEFT, fill=BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=LEFT, fill=BOTH, expand=1)
        self.figsub = self.fig.add_subplot(111)
        # to enable logaritmic plotinser the next two line
        self.canvas.mpl_connect("key_press_event", self.on_key)
        #self.log=0
        
    def on_key(self,event):
        """to change scale from linear to logaritm
        """
        if str(event.key) =='l':
            if self.figsub.get_yscale()=="log":
                self.figsub.set_yscale("linear")
            else:
                self.figsub.set_yscale("log")
        try:
            self.figsub.set_autoscaley_on(True)
            self.canvas.draw() 
            #self.figsub.set_ylim(ymax = max(self.curves[0]._y),#*+.1*abs(max(self.curves[0]._y)),
                                 #ymin = min(self.curves[0]._y))#*-.1*abs(max(self.curves[0]._y)))  
            #self.canvas.draw() 
        except:
            print "pippo"
        pass    
                

    def errorbar(self, x_array, y_array, z_array, comment= None,title=None, ylabel= "Mu (a.u.)", xlabel="Energy (eV)"):
       if (comment) is None:
            comment = [None for i in x_array]
       self.errcurves = self.figsub.errorbar( x_array[0], y_array[0], z_array[0], label= comment[0] )  #,
       self.figsub.set_ylabel(ylabel, fontsize = 8)
       self.figsub.set_xlabel(xlabel, fontsize = 8)
       if any(comment): 
            self.figsub.legend()
       if (title): self.figsub.set_title(title)
       self.toolbar.update()
       step=x_array[0][1]-x_array[0][1]
       print step
       self.figsub.set_xlim(xmin=x_array[0][0]-step, xmax=x_array[0][-1]+step)
       self.canvas.draw()
       self.figsub.set_autoscale_on(False)

    def plot(self, x_array, y_array, comment= None,title=None, ylabel= "", xlabel=""):
       """ycalcurves = array, calcurves = line!!!!!!!"""
       #print len(x_array)
       #print "*******************"
       #print len(y_array[0])
       if (comment) is None:
            comment = [None for i in x_array]
       if hasattr(self, "curves") and hasattr(self, "calcurves"):
            pass
       elif hasattr(self, "curves"):
            self.ycalcurves  =  y_array
            self.calcurves = self.figsub.plot( x_array[0], y_array[0], label= comment[0] )  #,
            for i in  range(len(x_array)-1):
                self.calcurves += self.figsub.plot(x_array[i+1], y_array[i+1], label = comment[i+1])
            for item in self.curves: item.set_marker('+')
       else:
            self.ycurves  =  y_array
            if (len(y_array)>1)and(not(hasattr(self, "slider"))) :
                    self.slider = Scale(self.top, from_= 0, to=1,       #
                                                     command= self.scale,   #variable= self.sh, 
                                                     orient=VERTICAL,
                                                     label= "Shift"
                                                     )
                    self.slider.pack(side = LEFT,fill = BOTH, anchor = W,pady = 15, ipady = 0)
            self.curves = self.figsub.plot( x_array[0], y_array[0], label= comment[0] )  #,
            for i in  range(len(x_array)-1):
                self.curves += self.figsub.plot(x_array[i+1], y_array[i+1], label = comment[i+1])
       self.figsub.set_ylabel(ylabel, fontsize = 8)
       self.figsub.set_xlabel(xlabel, fontsize = 8)
       if any(comment): 
            self.figsub.legend()
       if (title): self.figsub.set_title(title)
       self.toolbar.update()
       self.step=max(y_array[0])-min(y_array[0])
       self.figsub.set_ylim(ymin=(min(self.curves[0]._y)-self.step/10))
       self.figsub.set_autoscale_on(False)
       if len(self.curves)>1:
           self.slider.configure(to = self.step, resolution =self.step/50)
       self.canvas.draw()        
       
    def scale(self,event):
       if  hasattr(self, "calcurves"):
           for i,item in enumerate(self.calcurves):
               item.set_ydata(np_array(self.ycalcurves[i]) + i*float(event))
       if hasattr(self, "curves"):
           for i,item in enumerate(self.curves):
               item.set_ydata(np_array(self.ycurves[i]) + i*float(event))
               
           if len(self.curves)<10:
               self.figsub.set_ylim(ymax = max(map(max,[item._y for item in self.curves]))
                                +self.step/10)  
           else:    
               self.figsub.set_ylim(ymax = max(self.curves[-1]._y)+self.step/10)    
           self.canvas.draw()
       else: pass

               
    def clear(self):
        self.figsub.clear()
        if hasattr(self, "curves"): del self.curves
        if hasattr(self, "curves"): del self.calcurves
        self.canvas.draw() 
   
class EstraGraph:
    """
    class to have a  graph windows interface with some line to obtain some values
    genitore = quadro in wich insert the graph, (often is a top windows)
    plotting_list = list with all the object containg the data
    xattr= string that defining the attribute in wich is contained the abscissa
    yattr= list of strin defining the attributes for ordinates
    """
    def __init__(self, genitore, plotting_set):    
        """Plotting set is a dictionary with key(yattrib) link to a 
        string xattribute
        boolean   plot active
        stype line style
        """
        
        global estra_output
        self.curves={}
        self.check={}
        self.check_variable={}        

        self._check_deriv=1 #IntVar()
        #self._check_deriv.set(1)
        self.genitore=Frame(genitore)
        self.genitore.pack(side=LEFT, expand =Y , fill= BOTH)
        self.plotting_set = plotting_set
        self.fig = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.genitore)
        toolbar = NavigationToolbar2TkAgg(self.canvas,  self.genitore)
        toolbar.update()
        toolbar.pack(side=TOP, fill=X, expand=0)        
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        self.figsub = self.fig.add_subplot(111) 
        
        self.Quadro_sel=Frame(genitore)
        self.Quadro_sel.pack(side=LEFT, expand =N , fill= BOTH)
        for key in self.plotting_set.iterkeys():
            self.check_variable[key]=IntVar() 
            self.check_variable[key].set(self.plotting_set[key][1])   
            self.check[key] = Checkbutton(self.Quadro_sel, 
                                          variable= self.check_variable[key],
                                          text=key,
                                          width=15,
                                          command=self.check_command )
            self.check[key].pack(side=BOTTOM,expand=N,fill=BOTH)
        #self.link = self.canvas.mpl_connect('key_press_event', self.onspace)
        self.link = self.canvas.mpl_connect('button_press_event', self.onspace)
        pass    
            
            
    def onspace(self,event):
        global paste
        #print 'key=\'%s\', x=%d, y=%d, xdata=%f, ydata=%f'%(
        #    event.key, event.x, event.y, event.xdata, event.ydata)
        #if event.key==" ":
        if event.dblclick:
            paste=round(event.xdata,3)
            print event.xdata,"\n"
        else: pass
    
    
    def plot(self):
        """ycalcurves"""
        for key,value in self.plotting_set.iteritems():
            if value[1]:
                self.curves[key], = self.figsub.plot(
                                 getattr(estra_output, value[0]),      
                                 getattr(estra_output, key), value[2])
        self.canvas.draw()
        
        
    def check_command(self):
        """ycalcurves"""
        for key,value in self.check_variable.iteritems():
            self.plotting_set[key][1]=value.get()
        self.refresh()
        
        
        
        
        
    def refresh(self):
        """
        refesch graph
        """

        for key,value in self.plotting_set.iteritems():
            if value[1]:
                if self.curves.has_key(key):
                    self.curves[key].set_ydata(getattr(estra_output, key))
                    self.curves[key].set_xdata(getattr(estra_output, value[0])) 
                else:
                    self.curves[key], = self.figsub.plot(
                                 getattr(estra_output, value[0]),      
                                 getattr(estra_output, key), value[2])                    
            else:
                if self.curves.has_key(key):
                    self.figsub.lines.remove(self.curves[key])
                    del self.curves[key]
                else:pass
            self.canvas.draw()                         





class EstraGraph_kMu(EstraGraph):                
    def plot(self):
        """ycalcurves"""
        e2k=lambda x:0.512 * np_sqrt(abs(x-estra_output.edge))*np_sign(x-estra_output.edge)
        for key,value in self.plotting_set.iteritems():
            if value[1]:
                self.curves[key], = self.figsub.plot(
                                 map(e2k,getattr(estra_output, value[0])),      
                                 getattr(estra_output, key), value[2])  
        self.canvas.draw()        



    def refresh(self):
        """
        refesch graph
        """ 
        e2k=lambda x:0.512 * np_sqrt(abs(x-estra_output.edge))*np_sign(x-estra_output.edge)        
        for key,value in self.plotting_set.iteritems():
            if value[1]:
                if self.curves.has_key(key):
                    self.curves[key].set_ydata(getattr(estra_output, key))
                    self.curves[key].set_xdata(map(e2k,getattr(estra_output, value[0]))) 
                else:
                    self.curves[key], = self.figsub.plot(
                                 map(e2k,getattr(estra_output, value[0])),      
                                 getattr(estra_output, key), value[2])                    
            else:
                if self.curves.has_key(key):
                    self.figsub.lines.remove(self.curves[key])
                    del self.curves[key]
                else:pass
            self.canvas.draw()                         


class EstraGraph_exa(EstraGraph): 
    def __init__(self, genitore, plotting_set):
        EstraGraph.__init__(self,genitore, plotting_set)
        #super(EstraGraph_exao, self).__init__(genitore, plotting_set)
        self.figsub2 = self.figsub.twinx()
    
        
        
    def plot(self):
        """ycalcurves"""
        for key,value in self.plotting_set.iteritems():
            if value[1]:
                if key!="wind":
                    self.curves[key], = self.figsub.plot(
                                 getattr(estra_output, value[0]),      
                                 getattr(estra_output, key), value[2])  
                else:
                    self.curves[key], = self.figsub2.plot(
                                 getattr(estra_output, value[0]),      
                                 getattr(estra_output, key), value[2]) 
        self.figsub2.set_ylim(ymax = 1.05, ymin = -0.0)            
        self.canvas.draw() 
    def refresh(self):
        """
        refesch graph
        """
        for key,value in self.plotting_set.iteritems():
            if value[1]:
                if self.curves.has_key(key):
                    self.curves[key].set_ydata(getattr(estra_output, key))
                    self.curves[key].set_xdata(getattr(estra_output, value[0])) 
                else:
                    if key!="wind":
                        self.curves[key], = self.figsub.plot(
                                 getattr(estra_output, value[0]),      
                                 getattr(estra_output, key), value[2])  
                    else:
                        self.curves[key], = self.figsub2.plot(
                                 getattr(estra_output, value[0]),      
                                 getattr(estra_output, key), value[2])                    
            else:
                if self.curves.has_key(key):
                    if key!="wind":
                        self.figsub.lines.remove(self.curves[key])
                    else:
                        self.figsub2.lines.remove(self.curves[key])                        
                    del self.curves[key]
                else:pass
        self.figsub2.set_ylim(ymax = 1.05, ymin = 0)         
        self.canvas.draw()       
    
class EstraGraph_fou(EstraGraph):
    def plot(self):
        EstraGraph.plot(self)
        self.figsub.axhline(0, c="k")


#######################################################################################################
class mymenu():
    def __init__(self,genitore):
        global statusmsg
        self.genitore=genitore
        self.genitore.option_add('*tearOff', FALSE)
        self.menubar = Menu(genitore)
        # create a pulldown menu, and add it to the menu bar
        self.filemenu = Menu(self.menubar, tearoff=0)

        ##########################  files   #########################################################
        self.filemenu.add_command(label="Open Estra input", command=self.openEstra)
        self.filemenu.add_command(label="Save Estra input", command=self.saveEstra)        
        self.filemenu.add_command(label="Open Datafiles to average", command=self.openData) 
        
        self.filemenu.add_command(label="Clear Ini and Exit", command=self.quit2)  
        self.filemenu.add_command(label="Exit", command=self.quit1)
              
        
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        
        self.helpmenu = Menu(self.menubar, tearoff=0)
        self.helpmenu.add_command(label="Shorcuts and easter eggs", command=self.help)
        self.helpmenu.add_command(label="About", command=self.version)        
        self.menubar.add_cascade(label="Help", menu=self.helpmenu)
        
        # display the menu
        genitore.config(menu=self.menubar)
        
        
    def help(self):
       Top=Toplevel()
       #-----------------------------    Help Final Label      -----------------------------------------------' 
       Label(Top, text= "When Open Datafile menu select more than one file, the program will use the average").pack(side = TOP, anchor= W, pady=0)          
       Label(Top, text= "Use negative Edge Value for fixed Edge position").pack(side = TOP, anchor= W, pady=0)       
       Label(Top, text= "Press <Enter> to run the program").pack(side = TOP, anchor= W, pady=0)
       Label(Top, text= "Double click on the graph  to copy and <right button> to paste on the GUI"
              ).pack(side = TOP, anchor= W,pady=1)        


       
    def version(self):
        print "\n version %s \n" %__version__
        
    def quit1(self):
            print "\n\nfinally someones reads the menu ; )\nthank you again\n Carmelo\n"
            writeini()
            self.genitore.quit()
            
    def quit2(self):
            print "\n\nfinally someones reads the menu ; )\nthank you again\n Carmelo\n"
            clearini()
            self.genitore.quit()            
            
    def openEstra(self):
        global estra_input
        filename=browse_single()
        os.chdir(os.path.dirname(filename))
        filename =os.path.basename(filename)
        estra_input=class_estra_input(filename)
        estra_input.read_input()
        if hasattr(estra_input, "filedata"):
            pass
        else:
            pass
        return

    def saveEstra(self):
        global estra_input
        namesave= tkFileDialog.asksaveasfilename(title="Save Estra input file")
        estra_input.write_input(namesave)



    def openData(self):
        global data_file
        filenames = browse_multiple()
        if filenames!="":
            colwin=Column_Window(filenames)
        return


#######################################################################################################
class Corr_Line:
    def __init__(self, genitore, Limits_default=None, Limits=None): 
      #-----------------------------      Declare      --------------------------------------------------
        self._check= IntVar()
        self._compo=IntVar()
        self._lim=[]
        for item in range(len(Limits_default)):
            self._lim.append(StringVar())
      #-----------------------------      Define       --------------------------------------------------  
        self._check.set(False)
        if Limits:
            for item in range(len(Limits_default)):            
                self._lim[item].set(str(Limits[item])) 
        else:
            for item in range(len(Limits_default)):            
                self._lim[item].set(str(Limits_default[item]))             

        
      #-----------------------------      Structure    --------------------------------------------------
        self.check = Checkbutton(genitore, variable=self._check, width=2)#, command=self.underline )
        self.check.grid(row= 0, column=0)
        self.Entrylist=[]
        for item in range(len(Limits_default)):         
            self.Entrylist.append(pEntry(genitore, textvariable =  self._lim[item],
                                        width=14, justify=CENTER ))
            self.Entrylist[-1].grid(row= 0, column=item+1, padx=10)#.pack(side=TOP, expand=Y, fill=BOTH)
      #-----------------------------      Functios    --------------------------------------------------
    def destroy(self):
       self.check.destroy()
       for item in range(len(self._lim)):       
            self.Entrylist[item].destroy()
       
    def return_value(self):
        return [float(self._lim[item].get()) for item in  range(len(self._lim))]
      
      
class Corr_Glitch:
    def __init__(self, genitore, Title): 
        global estra_input
        self.padre=genitore
        self.genitore=Frame(self.padre)
        self.genitore.pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        Label(self.genitore, text=Title).pack(side=TOP, expand=Y, fill=BOTH)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        self.build()
        self.Quad_But=Frame(self.genitore)
        self.Quad_But.pack(side = BOTTOM, anchor= W, expand = Y, fill = BOTH)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = BOTTOM, anchor= W, expand = Y, fill = BOTH)
        self.addButton=Button(self.Quad_But,
               command = self.add,
               text = "Add line",
               width = 7)
        self.addButton.pack(side = LEFT, anchor= W, expand = Y, fill = BOTH)
        self.removeButton=Button(self.Quad_But,
               command = self.remove,
               text = "Remove selected lines",
               width = 7)
        self.removeButton.pack(side = LEFT, anchor= W, expand = Y, fill = BOTH)
        self.padre.protocol("WM_DELETE_WINDOW", self.topcallback)
    #######################    functions
    def build(self):
        global estra_input
        self.magic_string="glitches"
        self.default=[12.5, 12.6]
        self.frames=[]
        self.lines=[]
        for i in  range(estra_input.inp_exa.key_dict[self.magic_string+"_n"]):
            self.frames.append(Frame(self.genitore))
            self.frames[-1].pack(side=TOP, fill=BOTH, expand=Y)
            self.lines.append(Corr_Line(self.frames[-1],
                                        self.default,
                                        estra_input.inp_exa.key_dict[self.magic_string+"_lim"][i][:len(self.default)]))
            
    def remove(self):
        for i,item  in  enumerate(self.lines):
            if item._check.get():
                self.lines[i].destroy()
                del self.lines[i]
                self.frames[i].destroy()
                del self.frames[i]
                estra_input.inp_exa.key_dict[self.magic_string+"_n"] -=1 
                if estra_input.inp_exa.key_dict[self.magic_string+"_n"]<0: estra_input.inp_exa.key_dict[self.magic_string+"_n"]=0 
                break
        else:
            print estra_input.inp_exa.key_dict[self.magic_string+"_n"]
            return
        self.remove()    

        
    def add(self):
            self.frames.append(Frame(self.genitore))
            self.frames[-1].pack(side=TOP, fill=BOTH, expand=Y)
            self.lines.append(Corr_Line(self.frames[-1],self.default))
            estra_input.inp_exa.key_dict[self.magic_string+"_n"] +=1
            
    def read(self):
        #print "pippo00000000000000000000000000000000000000000000"
        estra_input.inp_exa.key_dict[self.magic_string+"_lim"]=[]
        for i,item  in  enumerate(self.lines):
            estra_input.inp_exa.key_dict[self.magic_string+"_lim"].append(item.return_value())
    def topcallback(self):
        self.read()
        self.padre.destroy()                     


class Corr_Sec_Edge(Corr_Glitch):       
    def build(self):
        self.magic_string="sec_edge"
        self.default=[12.5, 0.02, 0.001]
        self.frames=[]
        self.lines=[]
        for i in  range(estra_input.inp_exa.key_dict[self.magic_string+"_n"]):
            self.frames.append(Frame(self.genitore))
            self.frames[-1].pack(side=TOP, fill=BOTH, expand=Y)
            self.lines.append(Corr_Line(self.frames[-1],
                                        self.default,
                                        estra_input.inp_exa.key_dict[self.magic_string+"_lim"][i][:len(self.default)]))
                
        




#------------------------------------------------------------------------------------------------------



class pEntry(Entry):
    def __init__(self, master, **kw):
            apply(Entry.__init__, (self, master), kw)
            self.bind("<Button-3>", self.printo)

            
    def printo(self,event):
        self.delete(0,100)
        self.insert(0,paste)

        
        
class EXA():
    def __init__(self, genitore):
        global estra_input
      #-----------------------------      Declare      --------------------------------------------------
        self._label = StringVar()
        self._estra_input = StringVar()
        self._filedata = StringVar()
        #--------------------------      Declare  Exa    -----------------------------------------  
        self._pre_edge1 = StringVar()
        self._pre_edge2 = StringVar()
        self._edge1     = StringVar()
        self._edge2     = StringVar()   
        self._splines_n = IntVar()
        self._glitches_n =IntVar()
        self._sec_edge_n= IntVar()
        self._lim_Jump_eval1= StringVar()
        self._lim_Jump_eval2= StringVar()
        self._Const_Edge_eval=[StringVar(),StringVar(),IntVar()]
        self._k_lim1    = StringVar()    
        self._k_lim2    = StringVar()
        self._k_w       = IntVar()
        self._R_par1    = StringVar()
        self._R_par2    = StringVar()
        self._apod      = StringVar()
        self._win       = StringVar()
        self._BF_lim1   = StringVar()
        self._BF_lim2   = StringVar()
        self._Rexp      = StringVar()
        self._Jump_eval = StringVar()
        
        
        

      #-----------------------------      Define      --------------------------------------------------
        self._pre_edge1.set("")
        self._pre_edge2.set("")
        self._edge1.set("")
        self._edge2.set("") 
        self._lim_Jump_eval1.set(3)
        self._lim_Jump_eval2.set(3.5)
        self._splines_n.set(5)
        self._glitches_n.set(0)
        self._sec_edge_n.set(0)  
        self._Const_Edge_eval[0].set("n")
        self._Const_Edge_eval[1].set("k")      
        self._Const_Edge_eval[2].set(3)
        self._k_lim1.set("")
        self._k_lim2.set("")
        self._k_w.set(3)
        self._R_par1.set("7")
        self._R_par2.set("0.01") 
        self._apod.set(0.1)
        self._win.set("G")
        self._BF_lim1.set("1.0")
        self._BF_lim2.set("3.0")
        self._Rexp.set("1.80")
        self._Jump_eval.set("")
        
        
        self._label.set("Not defined")
        self._filedata.set("Not Defined")  
        self._estra_input.set("Not Defined")
        

      #-----------------------------      geometry      --------------------------------------------------'
      
        
        self.genitore=Frame(genitore)
        self.genitore.pack(side = LEFT, anchor= W, expand = Y, fill = BOTH)
        self.Quadro_labes=Frame(self.genitore)
        self.Quadro_labes.pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        self.QLab_estra=LabelFrame(self.Quadro_labes, text="Input file")
        self.QLab_estra.grid(row=0, column=0,pady=5, padx=5)
        self.QLab_data=LabelFrame(self.Quadro_labes, text="Data  file")
        self.QLab_data.grid(row=0, column=1,pady=5, padx=5)
        self.QLab_label=LabelFrame(self.Quadro_labes, text="Sample")
        self.QLab_label.grid(row=0, column=2,pady=5, padx=5)        
        Entry(self.QLab_estra, text= self._estra_input, state="disabled", width=18).pack() #
        Entry(self.QLab_label, text= self._label, width=18).pack()
        Entry(self.QLab_data, text= self._filedata, state="disabled", width=18).pack()
        #-----------------------------      start param      -----------------------------------------------'        
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH) 
        #-----------------------------      start param      -----------------------------------------------'
        self.Quadro_EXA1=Frame(self.genitore)
        self.Quadro_EXA1.pack(side = TOP, anchor= W, expand = Y, fill = BOTH)   
        self.QLab_pre_edge=LabelFrame(self.Quadro_EXA1, text="Pre_edge")
        self.QLab_pre_edge.grid(row=0, column=0,pady=5, padx=5) 
        pEntry(self.QLab_pre_edge, text= self._pre_edge1, width=9).grid(row=0, column=0,pady=5, padx=5)      
        pEntry(self.QLab_pre_edge, text= self._pre_edge2, width=9).grid(row=0, column=1,pady=5, padx=5)              
        
        self.QLab_edge=LabelFrame(self.Quadro_EXA1, text="Edge Value and Tol.")
        self.QLab_edge.grid(row=0, column=1,pady=5, padx=5) 
        pEntry(self.QLab_edge, text= self._edge1, width=9).grid(row=0, column=0,pady=5, padx=5)              
        pEntry(self.QLab_edge, text= self._edge2, width=9).grid(row=0, column=1,pady=5, padx=5)
        
        self.QLab_Jumpl=LabelFrame(self.Quadro_EXA1, text= "Jump evaluation limit")
        self.QLab_Jumpl.grid(row=0, column=2,pady=5, padx=5)                         
        pEntry(self.QLab_Jumpl, text= self._lim_Jump_eval1, width=9).grid(row=2, column=1,pady=5, padx=5)      
        pEntry(self.QLab_Jumpl, text= self._lim_Jump_eval2, width=9).grid(row=2, column=2,pady=5, padx=5)
        
        #Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH)   
      #-----------------------------      spline param      -----------------------------------------------'  
        self.Quadro_EXA2=LabelFrame(self.genitore,text="Spline evaluation")
        self.Quadro_EXA2.pack(side = TOP, anchor= W, expand = N, fill = BOTH)    
        Q_Spline_n=LabelFrame(self.Quadro_EXA2, text= " Splines Number  ")
        Q_Spline_n.pack(side = TOP, anchor= W, expand = N, fill = Y, padx=10)
        Spinbox(Q_Spline_n, from_ = 1, to = 10, 
                                 command= self.cr_spl_entry,  #lambda x=self.num_der:  self.panor2_der(x), 
                                 textvariable= self._splines_n,
                                 state= "readonly",
                                 width = 3).pack(side = TOP, anchor= W, expand = N, fill = Y)

        
        self.QLab_spline_int=LabelFrame(self.Quadro_EXA2, text="Splines Intervals")
        self.QLab_spline_int.pack(side = TOP, anchor= W, expand = Y, fill = Y, padx=5)   
        self.QLab_spline_deg=LabelFrame(self.Quadro_EXA2, text="Splines degrees")
        self.QLab_spline_deg.pack(side = TOP, anchor= W, expand = Y, fill = Y, padx=5)           
        Label(self.QLab_spline_deg, text="", width=5).grid(row=0, column=0 ,pady=10)
        
        
        self.cr_spl_entry()
        
                                 
        self.Quadro_EXA5=LabelFrame(self.genitore, text="Edge Constrain")
        self.Quadro_EXA5.pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        
        

        QLab_Jump=LabelFrame(self.Quadro_EXA5, text="Jump")
        QLab_Jump.grid(row=0, column=0)
        Jump=Combobox(QLab_Jump, textvariable= self._Const_Edge_eval[0],
                                 state ="readonly",values=("n", "J"), width=3)
        Jump.pack(side=LEFT, padx=3)
        Jump.bind("<<ComboboxSelected>>",self.Jump_value)
        self.Jump_Entry=pEntry(QLab_Jump, text= self._Jump_eval, width=7)
        self.Jump_Entry.pack(side=LEFT, padx=3)
        self.Jump_value("x")        
        
                                 
        QLab_Weigh=LabelFrame(self.Quadro_EXA5, text="Estraction Weigh ")
        QLab_Weigh.grid(row=0, column=1,padx=35)        
        Combobox(QLab_Weigh, textvariable= self._Const_Edge_eval[1],
                                 state ="readonly",values=("E", "k"), width=3
                                 ).pack(side=LEFT, padx=3)

        Label(QLab_Weigh, text= "      ").pack(side=LEFT, padx=3)         
        Spinbox(QLab_Weigh, from_ = 0, to = 5, 
                                 textvariable= self._Const_Edge_eval[2],
                                 state= "readonly",
                                 width = 3).pack(side=LEFT, padx=3)
        Label(QLab_Weigh, text= "wgt").pack(side=LEFT, padx=3)                                 
                                 
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH) 
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH) 

      #-----------------------------      correction param      -----------------------------------------------' 
        self.Quadro_EXA3=LabelFrame(self.genitore, text="Correction param.")
        self.Quadro_EXA3.pack(side = TOP, anchor= W, expand = N, fill = BOTH)
        
        self.Quadro_Glit=LabelFrame(self.Quadro_EXA3, text="Glitches param.")
        self.Quadro_Glit.pack(side = LEFT, anchor= E, expand = N, fill = BOTH) 
        pEntry(self.Quadro_Glit, text= self._glitches_n, width=3, state="readonly").grid(row=0, column=0)        
        Label(self.Quadro_Glit, text= "Glitches  ").grid(row=0, column=1)      
        Button(self.Quadro_Glit,
               command = self.Set_glitches,
               text = "Set",
               width = 7).grid(row=0, column=2,pady=5, padx=5)
        ##Spinbox(self.Quadro_Glit, from_ = 0, to = 5, #command= lambda x=self.num_der:  self.panor2_der(x), 
        #                         textvariable= self._glitches_n,
        #                         state= "readonly",
        #                         width = 3).grid(row=0, column=1,pady=5, padx=5)
        

        Frame(self.Quadro_EXA3).pack(side = LEFT, anchor= E, expand = Y, fill = BOTH)
        self.Quadro_secEd=LabelFrame(self.Quadro_EXA3, text="Secondary edges")
        self.Quadro_secEd.pack(side = LEFT, anchor= E, expand = N, fill = BOTH)   
        pEntry(self.Quadro_secEd, text= self._sec_edge_n, width=3, state="readonly").grid(row=0, column=0)        
        Label(self.Quadro_secEd, text= "Secondary edges  ").grid(row=0, column=1)      
        Button(self.Quadro_secEd,
               command = self.Set_sec_edge,
               text = "Set",
               width = 7).grid(row=0, column=2,pady=5, padx=5)
        #Spinbox(self.Quadro_secEd, from_ = 0, to = 5, #command= lambda x=self.num_der:  self.panor2_der(x), 
        #                         textvariable= self._sec_edge_n,
        #                         state= "readonly",
        #                         width = 3).grid(row=0, column=1,pady=5, padx=5)                                 
                                 
                                 
                                 
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH) 
      #-----------------------------      Fourier param      -----------------------------------------------'        
        self.Quadro_FOU1=LabelFrame(self.genitore, text="Fourier Parameter")
        self.Quadro_FOU1.pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        self.QLab_Fou_klim=LabelFrame(self.Quadro_FOU1, text="k limits")
        self.QLab_Fou_klim.grid(row=0, column=0,pady=5, padx=5) 
        pEntry(self.QLab_Fou_klim, text= self._k_lim1, width=5).grid(row=0, column=0,pady=5, padx=5)      
        pEntry(self.QLab_Fou_klim, text= self._k_lim2, width=5).grid(row=0, column=1,pady=5, padx=5)
        
        Label(self.Quadro_FOU1, text="       ").grid(row=0, column=1,pady=5, padx=5)
        self.QLab_Fou_k_w=LabelFrame(self.Quadro_FOU1, text="k weight")
        self.QLab_Fou_k_w.grid(row=0, column=1,pady=5, padx=5) 
        Combobox(self.QLab_Fou_k_w, textvariable= self._k_w,
                                 state ="readonly",values=(1, 2, 3), width=3
                                 ).grid(row=0, column=1,pady=5, padx=5) 
        
        
        QLab_Fou_Rpar1=LabelFrame(self.Quadro_FOU1, text="R max")
        QLab_Fou_Rpar1.grid(row=0, column=2,pady=5, padx=5)    
        QLab_Fou_Rpar2=LabelFrame(self.Quadro_FOU1, text="DR")
        QLab_Fou_Rpar2.grid(row=0, column=3,pady=5, padx=5)    
        pEntry(QLab_Fou_Rpar1, text= self._R_par1, width=3).grid(row=0, column=0,pady=5, padx=5)      
        pEntry(QLab_Fou_Rpar2, text= self._R_par2, width=5).grid(row=0, column=1,pady=5, padx=5) 
        
        QLab_Fou_win=LabelFrame(self.Quadro_FOU1, text="Win. type")
        QLab_Fou_win.grid(row=0, column=4,pady=5, padx=5)    
        QLab_Fou_apod=LabelFrame(self.Quadro_FOU1, text="G. width")
        QLab_Fou_apod.grid(row=0, column=5,pady=5, padx=5)
        Combobox(QLab_Fou_win, textvariable= self._win,
                                 state ="readonly",values=("G", "H", "N"), width=2
                                 ).grid(row=0, column=1,pady=5, padx=5)      
        pEntry(QLab_Fou_apod, text= self._apod, width=5).grid(row=0, column=1,pady=5, padx=5)         
        

        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = Y, fill = BOTH) 
      #-----------------------------    back fourier  start param      -----------------------------------------------'        
        self.Quadro_BFOU1=LabelFrame(self.genitore, text="Back Fourier Parameter")
        self.Quadro_BFOU1.pack(side = TOP, anchor= W, expand = Y, fill = BOTH)
        QLab_BFlim=LabelFrame(self.Quadro_BFOU1, text="BF limits")
        QLab_BFlim.grid(row=0, column=0,pady=5, padx=5) 
        pEntry(QLab_BFlim, text= self._BF_lim1, width=5).grid(row=0, column=0,pady=5, padx=5)      
        pEntry(QLab_BFlim, text= self._BF_lim2, width=5).grid(row=0, column=1,pady=5, padx=5)
        Label(self.Quadro_BFOU1, text= "  ").grid(row=0, column=1,pady=5, padx=5)        
        QLab_Rexp=LabelFrame(self.Quadro_BFOU1, text="R exp.")
        QLab_Rexp.grid(row=0, column=2,pady=5, padx=5)    
        pEntry(QLab_Rexp, text= self._Rexp, width=4).grid(row=0, column=0,pady=5, padx=5)         
        
        
        
  #-----------------------------      FUNCTION      --------------------------------------------------'    


    def Set_glitches(self):
        self.Top_glitch=Toplevel(takefocus=True)
        self.Top_glitch.title("Define glitches ranges")
        self.Top_glitch.wm_attributes("-topmost", 1)
        self.glitch_cla=Corr_Glitch(self.Top_glitch, "Selec       lim1        lim2                 ")
        self.glitch_cla.addButton.config(command= lambda : self.Correct("Add_G"))
        self.glitch_cla.removeButton.config(command= lambda : self.Correct("Remove_G"))      
        

    def Set_sec_edge(self):
        self.Top_glitch=Toplevel(takefocus=True)
        self.Top_glitch.title("Define glitches ranges")
        self.Top_glitch.wm_attributes("-topmost", 1)
        self.sec_edge_cla=Corr_Sec_Edge(self.Top_glitch, "Selec                Ede            Ade             Wde               ")
        self.sec_edge_cla.addButton.config(command= lambda : self.Correct("Add_S"))
        self.sec_edge_cla.removeButton.config(command= lambda : self.Correct("Remove_S"))           


    def Correct(self,event):
        if event=="Add_G":
           self.glitch_cla.add()
        if event=="Remove_G":
           self.glitch_cla.remove()           
        if event=="Add_S":
           self.sec_edge_cla.add()
        if event=="Remove_S":
           self.sec_edge_cla.remove()                
        #print "ddddddddddddddddddddddddddd", estra_input.inp_exa.key_dict["glitches_n"]
        self.set_info()
        
    def cr_spl_entry(self):
        if hasattr(self, "list_entry_inteval"):
            for item in self.list_entry_inteval:
               item.destroy()
            for item in self.list_entry_degree:
               item.destroy()
        self.list_entry_inteval =[]
        self.list_entry_degree  =[]
        self._splines_int=[StringVar() for item in range(self._splines_n.get()+1)]
        self._splines_degree=[StringVar()for item in range(self._splines_n.get())]  
        for item in range(self._splines_n.get()+1):
            self.list_entry_inteval.append(
                     pEntry(self.QLab_spline_int, text= self._splines_int[item], width=8))
            self.list_entry_inteval[-1].grid(row=0, column=item ,pady=5, padx=5)
 
        for item in range(self._splines_n.get()):
            self.list_entry_degree.append(
                       pEntry(self.QLab_spline_deg, text= self._splines_degree[item], width=8))
            self.list_entry_degree[-1].grid(row=0, column=item+1 ,pady=5, padx=5)
        
        #mette i numeri equispaziati
        if  globals().has_key("estra_input"):
            if self._splines_n.get()!=estra_input.inp_exa.key_dict["splines_n"]:
                estra_input.default_spline_intervals(self._splines_n.get())  
                self.set_info()
        return        
            
 
            
            
    def Jump_value(self,event):
        if self._Const_Edge_eval[0].get()=="J":
            self.Jump_Entry['state']='normal'
        elif self._Const_Edge_eval[0].get()=="j":
            self.Jump_Entry['state']='normal'            
        elif self._Const_Edge_eval[0].get()=="n":
            self.Jump_Entry['state']='disabled'
        pass            
        
        
    def read_info(self):
        global estra_input
        estra_input.label =self._label.get()
        estra_input.inp_exa.key_dict["pre_edge"]=map(
                         float,[self._pre_edge1.get(),self._pre_edge2.get()])
        estra_input.inp_exa.key_dict["edge"]=  map(                          
                            float,[self._edge1.get(),self._edge2.get()])
        #estra_input.inp_exa.key_dict["glitches_n"]=int(self._glitches_n.get())
        estra_input.inp_exa.key_dict["splines_n"]=int(self._splines_n.get())  
        estra_input.inp_exa.key_dict["splines_int"] =map( 
                            float,[item.get() for item in self._splines_int])
        estra_input.inp_exa.key_dict["splines_degree"] =map( 
                            int,[item.get() for item in self._splines_degree])
        #estra_input.inp_exa.key_dict["sec_edge_n"]=int(self._sec_edge_n.get())
        estra_input.inp_exa.key_dict["lim_Jump_eval"]=map(
                            float,[self._lim_Jump_eval1.get(),self._lim_Jump_eval2.get()])                            
        estra_input.inp_exa.key_dict["Const_Edge_eval"]=map(
                            str,[item.get() for item in self._Const_Edge_eval])
        estra_input.inp_exa.key_dict["Const_Edge_eval"][2]=int(
                            estra_input.inp_exa.key_dict["Const_Edge_eval"][2])
        try:
            estra_input.inp_exa.key_dict["Edge_jump"]=float(self._Jump_eval.get())
        except: 
            pass
        # Correction read ----------------------------------
        try:
            self.glitch_cla.read()  
        except:pass    
        try:    
            self.sec_edge_cla.read()
        except:pass    
        estra_input.inp_fou.key_dict["k_lim"]=map(
                            float,[self._k_lim1.get(),self._k_lim2.get()])
        estra_input.inp_fou.key_dict["k_w"]=int(self._k_w.get())
        estra_input.inp_fou.key_dict["window"]=str(self._win.get()) 
        estra_input.inp_fou.key_dict["window_apodization"]=float(self._apod.get())                             
                            
        estra_input.inp_bfou.key_dict["R_lim"]=map(                            
                            float,[self._BF_lim1.get(),self._BF_lim2.get()])

        estra_input.inp_bfou.key_dict["R_exp"]=float(self._Rexp.get())

    def set_info(self):
        """set information on the gui after readed estra.inp 
        """
        global estra_input
        self._label.set(estra_input.label)
        self._filedata.set(estra_input.inp_exa.key_dict["filename"])
        self._estra_input.set(estra_input.input_filename)
        self._pre_edge1.set(estra_input.inp_exa.key_dict["pre_edge"][0])
        self._pre_edge2.set(estra_input.inp_exa.key_dict["pre_edge"][1])
        self._edge1.set(estra_input.inp_exa.key_dict["edge"][0])
        self._edge2.set(estra_input.inp_exa.key_dict["edge"][1])
        self._glitches_n.set(estra_input.inp_exa.key_dict["glitches_n"])
        self._splines_n.set(estra_input.inp_exa.key_dict["splines_n"]) 
        self.cr_spl_entry()
        for i,item in enumerate(estra_input.inp_exa.key_dict["splines_int"]):
            self._splines_int[i].set(item)
        for i,item in enumerate(estra_input.inp_exa.key_dict["splines_degree"]):
            self._splines_degree[i].set(item)
        self._sec_edge_n.set(estra_input.inp_exa.key_dict["sec_edge_n"])
        self._lim_Jump_eval1.set(estra_input.inp_exa.key_dict["lim_Jump_eval"][0])
        self._lim_Jump_eval2.set(estra_input.inp_exa.key_dict["lim_Jump_eval"][1])
        self._Const_Edge_eval[0].set(estra_input.inp_exa.key_dict["Const_Edge_eval"][0])
        self._Const_Edge_eval[1].set(estra_input.inp_exa.key_dict["Const_Edge_eval"][1])        
        self._Const_Edge_eval[2].set(estra_input.inp_exa.key_dict["Const_Edge_eval"][2])
        if estra_input.inp_exa.key_dict["Const_Edge_eval"][0]!="n":
                 self._Jump_eval.set(estra_input.inp_exa.key_dict["Edge_jump"])
                 self.Jump_value("")
        self._k_lim1.set(estra_input.inp_fou.key_dict["k_lim"][0])
        self._k_lim2.set(estra_input.inp_fou.key_dict["k_lim"][1])
        self._k_w.set(estra_input.inp_fou.key_dict["k_w"])
        self._win.set(estra_input.inp_fou.key_dict["window"])
        self._apod.set(estra_input.inp_fou.key_dict["window_apodization"])
        self._R_par1.set(estra_input.inp_bfou.key_dict["R_lim"][0])
        self._R_par2.set(estra_input.inp_bfou.key_dict["R_lim"][1])
        self._Rexp.set(estra_input.inp_bfou.key_dict["R_exp"])
        


        

#######################################################################################################
class GRAPH():
    def __init__(self, genitore):
        global estra_input
      #-----------------------------      Declare      --------------------------------------------------
        self._label = StringVar()
        self._estra_input = StringVar()
        self._filedata = StringVar()
        #--------------------------      Declare  Exa    -----------------------------------------  
        self._pre_edge1 = StringVar()
        self._pre_edge2 = StringVar()
        self._edge1     = StringVar()
        self._edge2     = StringVar()   
        self._splines_n = IntVar()
        

      #-----------------------------      Define      --------------------------------------------------
        self._pre_edge1.set("")
        self._pre_edge2.set("")
        self._edge1.set("")
        self._edge2.set("") 

      #-----------------------------      geometry      --------------------------------------------------'
      
        self.genitore=Frame(genitore)
        self.genitore.pack(side = LEFT, anchor= W, expand = Y, fill = BOTH)
        self.Quadrobutton=LabelFrame( self.genitore, text="Plot")
        self.Quadrobutton.grid(row=0, column=0,pady=5, padx=5)
        
        Button(self.Quadrobutton,
                       command = self.Preplot,
                       text = "plot Pre",
                       width = 13).pack(side = TOP, anchor = W,pady = 10, padx=5)
        
        Button(self.Quadrobutton, 
                       command = self.Jumpplot,
                       text = "plot Jump",
                       width = 13).pack(side = TOP, anchor = W,pady = 10, padx=5)
                       
        Button(self.Quadrobutton, 
                       command = self.kJumpplot,
                       text = "plot k-Jump",
                       width = 13).pack(side = TOP, anchor = W,pady = 10, padx=5)
                       
        Button(self.Quadrobutton,
                       command = self.Statplot,
                       text = "plot Sta",
                       width = 13).pack(side = TOP, anchor = W,pady = 10, padx=5)
        
        Button(self.Quadrobutton,
                       command = self.Exaplot,
                       text = "plot Exa",
                       width = 13).pack(side = TOP, anchor = W,pady = 10, padx=5)        

        Button(self.Quadrobutton,
                       command = self.Fouplot,
                       text = "plot Fourier",
                       width = 13).pack(side = TOP, anchor = W,pady = 10, padx=5)

        Button(self.Quadrobutton,
                       command = self.Bfouplot,
                       text = "plot Bf",
                       width = 13).pack(side = TOP, anchor = W,pady = 10, padx=5)


    def Preplot(self):
        global estra_input
        global estra_output
        self.top_Pre = Toplevel()
        self.top_Pre.title("Pre plot")        
        #--------------------------   Graphic win  --------------------------------------------------
        plotting_set={"ABSpre":["Ep",True, 'b-'], "der1":["Ep",True, 'g-'],
                      "der2":["Ep",False, 'r-'], "noiseE":["Ep",False, 'k-']}
        self.grap_pre=EstraGraph(self.top_Pre, plotting_set)
        self.grap_pre.plot()


    def Jumpplot(self):
        global estra_input
        global estra_output
        self.top_Jump = Toplevel()
        self.top_Jump.title("Jump plot")        
        #--------------------------   Graphic win  --------------------------------------------------
        plotting_set={"A_nor":["Ej",False, 'b-'], "A_0":["Ej",True, 'b-'],
                      "SPL_nor":["Ejsp",False, 'g-'], "SPL_ME_nor":["Ejsp",True, 'g-'],
                      "EknotsY":["Eknots",True, 'r.']}
        self.grap_Jump=EstraGraph(self.top_Jump,  plotting_set)
        self.grap_Jump.plot()


    def kJumpplot(self):
        global estra_input
        global estra_output
        self.top_Jump = Toplevel()
        self.top_Jump.title("Jump plot")        
        #--------------------------   Graphic win  --------------------------------------------------
        plotting_set={"A_nor":["Ej",False, 'b-'], "A_0":["Ej",True, 'r-'],
                      "SPL_nor":["Ejsp",False, 'g-'], "SPL_ME_nor":["Ejsp",True, 'y-'],
                      "EknotsY":["Eknots",True, 'r.']}
        self.grap_kJump=EstraGraph_kMu(self.top_Jump,  plotting_set)
        self.grap_kJump.plot() 
        
    def Statplot(self):
        global estra_input
        global estra_output
        self.top_Stat = Toplevel()
        self.top_Stat.title("Stat plot")        
        #--------------------------   Graphic win  --------------------------------------------------
        plotting_set={"der":["Ks",False, 'b-'],"Discont":["Ks",True, 'k-'],
                        "Spectrum_noise":["Ks",True, 'g-'],
                        "Average_noise":["Ks",True, 'b-']}
        self.grap_Stat=EstraGraph(self.top_Stat,  plotting_set)
        self.grap_Stat.plot()        
        

    def Exaplot(self):
        global estra_input
        global estra_output
        self.top_Exa = Toplevel()
        self.top_Exa.title("Exa plot")        
        #--------------------------   Graphic win  --------------------------------------------------
        plotting_set={"kchi":["k",True, 'b-'],"knchi":["kf",False, 'g-'],
                      "wind":["kf",True, 'r-']}
        self.grap_Exa=EstraGraph_exa(self.top_Exa,  plotting_set)
        self.grap_Exa.plot()  

    def Fouplot(self):
        global estra_input
        global estra_output
        self.top_Fou = Toplevel()
        self.top_Fou.title("Fourier plot")        
        #--------------------------   Graphic win  --------------------------------------------------
        plotting_set={"Amp":["R",True, 'b-'],"Re":["R",False, 'g-'],
                        "Im":["R",True, 'b-']}
        self.grap_Fou=EstraGraph_fou(self.top_Fou,  plotting_set)
        self.grap_Fou.plot()   

    def Bfouplot(self):
        global estra_input
        global estra_output
        self.top_Fou = Toplevel()
        self.top_Fou.title("Back Fourier plot")        
        #--------------------------   Graphic win  --------------------------------------------------
        plotting_set={"qchi":["q",True, 'b-'],"qnchi":["q",False, 'b-'],
                        "kchi":["k",True, 'g-'],"knchi":["kf",False, 'g-']}
        self.grap_BFou=EstraGraph(self.top_Fou,  plotting_set)
        self.grap_BFou.plot()   

#######################################################################################################
############################            Main    class     #############################################
#######################################################################################################
class EstraGui:
    def __init__(self, genitore):
        global estra_input
        global estra_output
        global spectra
        global directory
        #global xlabel
        global inivar
        
        readini()
        #menu
        self.menu=mymenu(genitore)
        self.menu.filemenu.entryconfig(index=0, command=  self.set_info)
        self.menu.filemenu.entryconfig(index=1, command=  self.write_info)        
        self.menu.filemenu.entryconfig(index=2, command=  self.read_data)         
        
        self.EXA=EXA(genitore)
        self.GRAPH=GRAPH(genitore)  
        genitore.bind_all("<Return>", self.Return) 
        #genitore.bind_all("<space>", self.Return)   



  #-----------------------------      FUNCTION  GLOBAL +-    --------------------------------------------------'
    def Return(self,event):
        global estra_output
        global estra_input  
        self.EXA.read_info()
        estra_input.write_input()
        os.system("\""+os.path.join(inivar.get("Estra", "Estra_Dir"),"Estra.exe").join("\"\"")+"\"")
        estra_output=class_estra_output(os.path.splitext(estra_input.inp_exa.key_dict["filename"])[0])
        estra_output.read_output()
        if self.EXA._Const_Edge_eval[0].get() in ["n","N"]:
            try:
                self.EXA._Jump_eval.set(str(estra_output.Jump))
            except:
                pass
        #####----------------------------------Complement estra outup ----------------------------------
        estra_output.Kknots=estra_input.inp_exa. key_dict["splines_int"]
        #e2k=lambda x: 0.512317 * np_sqrt(x - estra_output.edge)
        k2e=lambda x: (x/0.512317)**2+ estra_output.edge
        estra_output.Eknots=map(k2e,estra_output.Kknots)
        estra_output.EknotsY=np_interp(estra_output.Eknots, 
                                          estra_output.Ejsp.squeeze(), estra_output.SPL_ME_nor.squeeze())
        ####----------------------------------Complement estra outup ----------------------------------        
        try: self.GRAPH.grap_pre.refresh()
        except AttributeError:  pass 
        if hasattr(self.GRAPH,"grap_Jump"):
           self.GRAPH.grap_Jump.refresh()
        if hasattr(self.GRAPH,"grap_kJump"):
           self.GRAPH.grap_kJump.refresh()
        if hasattr(self.GRAPH,"grap_Exa"):
           self.GRAPH.grap_Exa.refresh()
        if hasattr(self.GRAPH,"grap_Stat"):
           self.GRAPH.grap_Stat.refresh()   
        if hasattr(self.GRAPH,"grap_Fou"):
           self.GRAPH.grap_Fou.refresh()
        if hasattr(self.GRAPH,"grap_BFou"):
           self.GRAPH.grap_BFou.refresh()  
        #except AttributeError:  pass
        print "\n\ncall done\n\n"
        
    def set_info(self):
        self.menu.openEstra()
        self.EXA.set_info()
        
    def write_info(self):
        self.EXA.read_info()   
        self.menu.saveEstra()
 
    def read_data(self):
        self.menu.openData()        
        self.EXA.set_info()   
 
        


























##############   Inizialization   ############################################################  
def readini():
    global inivar
    if os.name =="nt":
         path_local_data=os.path.join(os.environ['APPDATA'],"EstraFitexa")
    elif os.name =="posix":
         path_local_data="~/.local/bin"
    else :
        print os.name, "ERROR--"*5+"\n  sistem not defined\n" +"ERROR--"*5
        return
    inifile=os.path.join(path_local_data,"EstraFitexa.ini")
    if __verbose__:  print inifile
    inivar.read(inifile)
    if __verbose__ : print os.getcwd()
    inivar.sections()
    if inivar.has_section("Estra"):
        if not(inivar.has_option("Estra", "Estra_Dir") and 
                os.access(inivar.get("Estra", "Estra_Dir"), os.F_OK)):
            inivar.set("Estra", "Estra_Dir", os.getcwd())    
        
        if os.access(inivar.get("Estra", "Start_Dir"), os.F_OK):
            os.chdir(inivar.get("Estra", "Start_Dir"))
        else:
            os.chdir(os.path.join(os.environ['HOMEDRIVE'],os.environ['HOMEPATH']))
   
    else:
       inivar.add_section("Estra")
       inivar.set("Estra", "Estra_Dir", os.getcwd())
       os.chdir(os.path.join(os.environ['HOMEDRIVE'],os.environ['HOMEPATH']))
    # put the file name in the section inifile
    try:
        inivar.add_section("Inifile")
    except ConfigParser.DuplicateSectionError: pass
    inivar.set("Inifile", "Inipath", path_local_data)
    inivar.set("Inifile", "Ininame", inifile)    
    return           
   
 
#################   Inizialization   ############################################################ 
def writeini():
    global inivar
    inivar.set("Estra", "Start_Dir", os.getcwd())
    path_local_data=inivar.get("Inifile", "Inipath")
    if not(os.access(path_local_data, os.F_OK)):
         os.mkdir(path_local_data)   
         
    inifile=inivar.get("Inifile", "Ininame")
    with open(inifile, 'w') as configfile:
        inivar.write(configfile)     
        configfile.close
    return           
    
def clearini():
    global inivar
    inivar.remove_section("Estra")
    path_local_data=inivar.get("Inifile", "Inipath")
    if not(os.access(path_local_data, os.F_OK)):
         os.mkdir(path_local_data)   
         
    inifile=inivar.get("Inifile", "Ininame")
    with open(inifile, 'w') as configfile:
        inivar.write(configfile)     
        configfile.close
    radice.quit()    
    return      

    

    
def destroy():
    print "\n\n\n\nhave a nice day.....  ;-) \n\n"
    writeini()
    radice.quit()    
    

if __name__ == "__main__":
   radice = Tk()
   radice.title("Estra GUI")
   pippo = EstraGui(radice)
   radice.protocol("WM_DELETE_WINDOW", destroy)
   radice.mainloop()










  
#        os.chdir(directory)
#        os.system(os.path.join(inivar["PrestoPronto_Dir"],"feff6l.exe").join("\"\""))
#        os.chdir(Start_Dir)
