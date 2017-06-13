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

#import text
import sys
import os
import ConfigParser
import re
from Tkinter import *
from ttk import *
import tkFileDialog
import tkMessageBox
import x_text as text
from collections import OrderedDict

global paste
import matplotlib
matplotlib.interactive(False)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, cursors
from matplotlib import rc

import numpy as np
global __verbose__                                                                    
__verbose__=False#True#
global __version__
__version__= "a0.0.0"
global inivar
inivar=ConfigParser.ConfigParser()

import warnings
global fitexa_output
global fitexa_input 

global paste








#######################################################################################################
########################          fitexa_output         ################################################
#######################################################################################################
class label_outread(object):
    def __init__(self,filename):
        with open(filename,'r') as filename:
            
            readlinex =lambda :(lambda x=filename.readline(): x.strip() if 
                                x.strip()  else readlinex())()  
                                
            self.INFO=dict()
            i=0
            for i in range(300) :
                line=readlinex()
                if "Number of shells" in line: self.INFO["shell_n"]=int(line.split(":")[1])
                if "Number of points" in line: self.INFO["points_n"]=int(line.split(":")[1])
                if "Number of parameters" in line: self.INFO["param_n"]=int(line.split(":")[1])                
                if "average std.dev on Chi(k)" in line: self.INFO["stdev_Chi"]=float(line.split("=")[1]) 
                if "wgt.std.dev  (k^wt Chi)" in line: self.INFO["stdev_Chiwt"]=float(line.split("=")[1]) 
                if "S/N ratio" in line: self.INFO["S/N"]=float(line.split("=")[1])                                 
                if "S/N_w ratio" in line: self.INFO["S/N_w"]=float(line.split("=")[1])                                 
                if "Least square found" in line: self.INFO["F"]=float(line.split("=")[1])                                   
                if "R^2_w" in line: self.INFO["R^2_w"]=float(line.split("=")[1])                   
                if "R^2_o" in line: self.INFO["R^2_o"]=float(line.split("=")[1])
                if "Chi^2" in line: self.INFO["Chi^2"]=float(line.split("=")[1])                   
                if "red._Chi^2" in line: self.INFO["red_Chi^2"]=float(line.split("=")[1])                
                if "AIC_m" in line: self.INFO["AIC_m"]=float(line.split("=")[1])                                     
                if "if Filtered data" in line: self.INFO["Nind"]=float(line.split("=")[2])     
                if "Suggested errdef" in line: self.INFO["Err65"]=map(float,line.split("=")[1].split())     
                if "red._Chi^2" in line: self.INFO["red_Chi^2"]=float(line.split("=")[1])  
                if "General parameters" in line: 
                    break
            else: 
                raise IOError("Output file not readable")
                print "\nOutput file not readable\n"

                    
                
            self.GEN_PARAMETER=dict()

            while True:
                line=readlinex()
                if "structural parameters" in line: break
                else: i,j=line.split("=")  
                self.GEN_PARAMETER[i.strip()]=float(j)
   
            while True:         
                if readlinex()[1:3]=="__":break
            
            self.shells=list()
            letto =lambda x: np.NAN if "*" in x else float(x)
            for i in range(self.INFO["shell_n"]):
                line=[item.strip() for item in filename.readline().split()]
                self.shells.append(dict())
                self.shells[-1]["type"]=line[0]
                self.shells[-1]["N"]    =letto(line[1])
                self.shells[-1]["R"]    =letto(line[2])
                self.shells[-1]["sigma"]=letto(line[3])
                self.shells[-1]["DE"]   =letto(line[4])   
                self.shells[-1]["gamma"]=letto(line[5])
                self.shells[-1]["C2"]   =letto(line[6])
                self.shells[-1]["C3"]   =letto(line[7]) 
                self.shells[-1]["one_string"]=' '.join([u'Fit results:',
                    u'N={0:1.2f},'.format(self.shells[-1]["N"]),
                    u'R= {0:1.2f}\u212B,'.format(self.shells[-1]["R"]),
                    u's2={0:0.4f}\u212B\u00B2,'.format(self.shells[-1]["sigma"]),
                    u'E={0:2.1f}eV'.format(self.shells[-1]["DE"])]) 
                if self.shells[-1]["gamma"]!=0:
                    self.shells[-1]["one_string"]+= " gamma="+str(self.shells[-1]["gamma"])
                if self.shells[-1]["C2"]!=0:
                    self.shells[-1]["one_string"]+= " C2=%0.4f" %(self.shells[-1]["C2"])
                    self.shells[-1]["one_string"]+= " C2=%0.4f" %(self.shells[-1]["C3"])
        return 


class output_fitexa(object):
    def __init__(self,label="A101"):
        self.label=label
        
    def read_output(self,label=None):
        if label !=None:
            self.label=label
        filename=self.label+".out"        
        self.out=label_outread(filename)
        #------------------------------#
        filename=self.label+".fit"
        self.curves=dict()
        with open(filename,'r') as filename:
            array=np.loadtxt(filename)
            self.curves["k"]=array[:,0]
            self.curves["kchi_e"]=array[:,1]
            self.curves["kchi_t"]=array[:,2]        
            self.curves["kchi_res"]=array[:,3] 
            self.curves["kwchi_e"]=array[:,4]
            self.curves["kwchi_t"]=array[:,5]        
            self.curves["kwchi_res"]=array[:,6] 
            self.curves["sigma_ratio"]=array[:,7]            
            filename.close()
        #------------------------------#    
        filename=self.label+".par"
        self.curves["chi_par"]=OrderedDict()
        with open(filename,'r') as filename:
            array=np.loadtxt(filename)
            for item in range(array.shape[1]):
                if item>0:
                    index="shell_"+str(item)
                    self.curves["chi_par"][index]=array[:,item]
            filename.close()        
        #------------------------------#        
        filename=self.label+".fou"
        with open(filename,'r') as filename:
            array=np.loadtxt(filename)
            self.curves["R"]=array[:,0]
            self.curves["MaFT_e"]=array[:,1]
            self.curves["ImFT_e"]=array[:,2]        
            self.curves["MaFT_t"]=array[:,3] 
            self.curves["ImFT_t"]=array[:,4]  
            filename.close()
            filename.close()           
        #------------------------------#        
        filename=self.label+"_Abs.fou"
        self.curves["MaFT_par"]=OrderedDict()
        with open(filename,'r') as filename:
            array=np.loadtxt(filename)
            for item in range(array.shape[1]):
                if item>0:
                    index="shell_"+str(item)
                    self.curves["MaFT_par"][index]=array[:,item]
            filename.close()           
        #------------------------------#        
        filename=self.label+"_Imm.fou"
        self.curves["ImFT_par"]=OrderedDict()
        with open(filename,'r') as filename:
            array=np.loadtxt(filename)
            for item in range(array.shape[1]):
                if item>0:
                    index="shell_"+str(item)
                    self.curves["ImFT_par"][index]=array[:,item]
            filename.close()   

        return 
#######################################################################################################
########################          fitexa_input         ################################################
#######################################################################################################
class variabile(object):
    def __init__(self, stringa):    
        """define a class with a variable for minuit
           with attribure N,param,value,error,min,max
           input = string la liea letta nell input
        """
        if __verbose__:                            
                print stringa

        if "'" in stringa:
           stringa=stringa.replace("'","") 
            
        stringa= stringa.split()
        if len(stringa)==0:
            raise TypeError
        if len(stringa)<6:
            messager=('\nBad variable format in MINUIT table:\n'
                      '{}\nsome field are missing'.format(stringa))
            print messager
            tkMessageBox.showinfo("Bad format in MINUIT table",
                                   messager)
            raise IndexError    
            
        stringa= stringa[:6]
        self.N= int(stringa[0])
        self.param=stringa[1].strip()
        self.value=float(stringa[2])
        self.error=float(stringa[3])
        self.min=float(stringa[4])        
        self.max=float(stringa[5])
        
    def make_string(self):
        self.stringa=repr(self.N).rjust(10)
        self.stringa+=" "+self.param.ljust(9)
        self.stringa+=repr(self.value).rjust(10)
        self.stringa+=repr(self.error).rjust(10)    
        
        self.stringa+=repr(self.min).rjust(10)            
        self.stringa+=repr(self.max).rjust(10)   
        #print self.stringa
        
        


class shell_path(object):
    """class that define a path and it property
    """
    def __init__(self):
        pass
    def read(self, filename):
        if __verbose__:
            print "\nreading a shell"
        readlinex =lambda :(lambda x=filename.readline().strip(): x if x!="" else readlinex())() 
        exp_var =lambda x: x.split("&")[1].strip() if x[0]=="&" else int(x.split()[0])
        exp_varCoor =lambda x: x.split("&")[1].strip() if x[0]=="&" else [int(x.split(",")[0]), float(x.split(",")[1].split()[0])] 
        self.tipo= readlinex()[0]    
        if __verbose__:
            print "\nself.tipo",self.tipo
            
        Multi=readlinex()
        if Multi[0]=="&":
            self.Coor=exp_var(Multi)    
        else:
            self.Coor, Multi=exp_varCoor(Multi)
            if Multi!=1.0:
                self.Coor= "#{0:d} * {1:.2f}".format(self.Coor, Multi)
        if __verbose__:
            print "self.Coor=" ,self.Coor
        self.Distance= exp_var(readlinex())
        self.Deb_Wal= exp_var(readlinex())
        self.n_de, self.n_gamma =map(int,readlinex().split("  ")[0].split(","))
        if self.tipo=="C":
           self.C2,self.C3 =map(int,readlinex().split("  ")[0].split(",")) 
        self.phase_type= readlinex().split("  ")[0].strip()
        self.phase=readlinex().strip("\'").strip("\"")



    def write(self):  
        spac1=30
        spac2=35
        exp_var_w= lambda x: repr(x) if type(x)==int else "& "+x+" &"
        
        stringa= self.tipo.ljust(spac1)+ "# Shell type".ljust(spac2)+"(ch1)\n"
        if type(self.Coor)==int:
            stringa+= repr(self.Coor)+",1".ljust(spac1)+\
                  "# or: &string& # coord.num".ljust(spac2)+"(int real)\n"
        else:
            stringa+= exp_var_w(self.Coor).ljust(spac1)+\
                  "# or: &string& # coord.num".ljust(spac2)+"(int real)\n"
        stringa+= exp_var_w(self.Distance).ljust(spac1)+\
                  "# or: &string& # distance".ljust(spac2)+"(int real)\n"
                  
        stringa+= exp_var_w(self.Deb_Wal).ljust(spac1)+\
                  "# or: &string& # s^2".ljust(spac2)+"(int real)\n"

        stringa+="{0:d},{1:d}".format(self.n_de, self.n_gamma).ljust(spac1)+\
                "# n_de,n_gamma".ljust(spac2)+"(2 int)\n"
        if self.tipo=="C":
            stringa+="{0:d},{1:d}".format(self.C2, self.C3).ljust(spac1)+\
                "# second, third cumulant".ljust(spac2)+"(2 int)\n"
        stringa+= self.phase_type.ljust(spac1)+"# A/F file  11: feff****.dat)".ljust(spac2)+\
                                       "(ch4)\n"
        stringa+="'"+self.phase+"'\n"       
        return stringa 
    
    def check_integrity(self):
        ###controlla che ci siano tutti i parametri
        #
        if self.Coor=="":
            raise ValueError ("ValueError: missing coordination number")
        if   self.Distance=="":
            raise ValueError ("ValueError: missing distance")
        if self.Deb_Wal=="":
            raise ValueError ("ValueError: missing Deb_Wal")  
        if self.n_de =="":
            raise ValueError ("ValueError: missing n_del")  
        if self.n_gamma=="":
            raise ValueError ("ValueError: missing gamma")  
        if self.tipo=="C":
            if self.C2=="":
                raise ValueError ("ValueError: missing C1")  
            if self.C3=="":
                raise ValueError ("ValueError: missing C2") 
        if self.phase=="":
                raise ValueError ("ValueError: missing phase file path")         
        

class input_fitexa(object):
    """Should be the class that keep information of variable
    """
    def __init__(self,label="A101"):
        self.label=label
        
    def read_inp(self,filename):
        self.filename=filename
        with open(filename,'r') as filename:
            self.variable_list=list()
            filename.readline()
            while True:
                try:
                    self.variable_list.append(variabile(filename.readline()))
                except TypeError:
                    break
                except IndexError:
                    return
                    
                    
            if __verbose__:                            
                print "reading post variable input\n"    
                
                
            readlinex =lambda :(lambda x=filename.readline(): x.strip() if 
                                x.strip() and x.strip()[0]!="#" else readlinex())()      
                
            self.chi=readlinex().strip("\'").strip("\"")
            self.label=readlinex().split("  ")[0]
            self.col =map(int,readlinex().split("  ")[0].split(","))
            self.xranges =map(float,readlinex().split("  ")[0].split(",")) 
            self.kw =float(readlinex().split("   ")[0])
            self.so2, self.eta =map(int,readlinex().split("  ")[0].split(","))  
            
            self.n_shell= int(readlinex().split("  ")[0])
            if __verbose__:
                print "self.n_shell=",self.n_shell
            
            self.FF=readlinex()[0]
            
            FF_param=["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]
            self.FF_param={}
            for i,j in zip(FF_param,readlinex().split()):
                self.FF_param[i]=float(j)
            if __verbose__:
                print "FF_param=",self.FF_param
                
                
            self.shell=list()
            #### read the path
            for item in range(self.n_shell): 
                self.shell.append(shell_path())
                self.shell[-1].read(filename)
            
            self.stat=float(readlinex().split("  ")[0])
            
            self.minuit=str()
            while True:
                try:
                    self.minuit+=readlinex()+"\n"
                except:
                    break
        return    

    def Name2N(self, expression):
        """given an expression with names give back the same expression
               with variable #number
           input expression without &&     
        """      
        def find_number(name):
            for item in self.variable_list:
                if name== item.param:
                    return "#"+str(item.N)
            else:
                raise ValueError(5*"\n"+20*"#"+"\n"+
                                 "# %s not coorespond to any variable \n"%(name) +
                                 20*"#"+"\n") 
        var= re.compile('[a-zA-Z_][a-zA-Z0-9_]*')  
        lista_name=re.findall(var,expression)
        lista_num= map(find_number, lista_name)
        for i,j in zip(lista_name,lista_num):
            expression = expression.replace(i,j)
        return expression        
       
    
    def N2Name(self, expression):
        """given an expression with #number give back the same expression
               with variable names
           input expression without &&     
        """
        def find_name(number):
            number= int(number.strip("#"))
            for item in self.variable_list:
                if number== item.N:
                    return item.param
            else:
                raise ValueError(20*"#"+"\n"+
                                 "# %i not coorespond to any variable \n"%(number) +
                                 20*"#"+"\n") 
        var= re.compile('(#\d*)')
        lista_var=re.findall(var,expression)
        lista_name= map(find_name, lista_var)
        for i,j in zip(lista_var,lista_name):
            expression=expression.replace(i,j)            
        return expression



    def write_inp(self,filename="fitexa.inp"):
        with open(filename,'w') as filename:
            filename.write("***param-N|*param***|****value|****error|******min|******max\n")
            for var in self.variable_list:
                var.make_string()
                filename.write(var.stringa+"\n")
            filename.write("\n")
            filename.write("\'{0:s}\'\n".format(self.chi)) 
            filename.write(self.label.ljust(25)+"# 4-char estension                       (ch4)\n")
            filename.write("{0:d},{1:d}".format(*self.col).ljust(25)+
                           "# colonne k, k*chi(k)\n")
            filename.write("{0:2.2f},{1:2.2f}".format(*self.xranges).ljust(25)+
                           "# range:  xmin xmax (xmin<=0 full range)(2 real)\n" )
            filename.write("{0:1.0f}".format(self.kw).ljust(25)+
                           "# data weight k^w                       (1 real)\n")
            filename.write("{0:d},{1:d}".format(self.so2, self.eta).ljust(25)+
                           "# Number so2, Number eta if 0: eta=3.   (2 int)\n")
            filename.write("\n")
            filename.write("{0:d}".format(self.n_shell).ljust(25)+
                           "# number of shells                      (1 int)\n") 
            filename.write("\n")
            filename.write(self.FF.format(self.n_shell).ljust(25)+
                           "#N/Y  use fourier filtering       (ch1)\n") 
            filename.write("#    k1     k2     wt  T(1H-2G) apo    Rup     dr     RL     RR\n")   
            stringa=str()
            for i in ["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]:
                stringa+=repr(self.FF_param[i]).rjust(7) 
            filename.write(stringa+"\n")             
            filename.write("\n")             
            for item in self.shell:
                filename.write(item.write())
                filename.write("\n")
                
            filename.write("{0:2.2f}".format(self.stat).ljust(25) +"#  max R/sig^2 for error analysis\n\n")
            filename.write("\n\n"+self.minuit)    
        return
            

#######################################################################################################
########################             MENU()            ################################################
#######################################################################################################
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
        #self._check.set(False)
      #-----------------------------      Structure    --------------------------------------------------
        try:
                self.check = Label(genitore,   text=self._label)#, command=self.underline )
                self.check.grid(row= row, column=0)
                
 
                self.column= Spinbox(genitore, from_ = 1, to = 18, 
                                     textvariable= self._position, width = 3)  
                self.column.grid(row= row, column=2)
                self.column.configure(value=row)
                
                self.pulsante_Plot = Button(genitore ,
                                              command = self.plot,
                                              text = "Plot",
                                              width = 7)
                self.pulsante_Plot.grid(row= row, column=3)   
        except Exception as exep:
          print type(exep)
          print exep.args
        
      #-----------------------------      Functios    --------------------------------------------------
    def plot(self): 
        title= "{0:s}  column {1:d}".format(self._label, int(self._position.get()))
        x_array=[np.arange(self.array.shape[0])]
        y_array=[self.array[:,int(self._position.get())-1]]
        graph = Graphcol()
        graph.plot(x_array, y_array, title= title)



class Column_Window:
    def __init__(self,filenames):
        global estra_input
      #--------------------------   Declare-------------------------------------------------
        self.filenames=filenames
        self._ChaCom=StringVar()
        self._mode=StringVar()
      #--------------------------   Define--------------------------------------------------  
        self.column_names=["k", "kchi"]
        self.column_list=[]
        self._ChaCom.set("#")
        self._mode.set("transmission")
      #--------------------------  Top level + reload--------------------------------------------------        
            

        
        s = Style()
        s.configure('Green.TButton',  background="green")

        
        self.top = Toplevel(takefocus=True)
        self.top.title("define column of interest")    

        


        self.win_text = Frame(self.top) 
        self.win_text.pack(side=LEFT,expand=YES, fill=BOTH)
        
        text.ScrolledText( parent=self.win_text, filex=self.filenames,hor=True, active=False)

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
                                      width = 13 )
        self.Reload_B.grid(row= 0, column=2)
      #--------------------------   Params  Entries--------------------------------------------------
        self.win_column = Frame(self.top_con)
        self.win_column.pack(side=TOP,expand=YES)
      #--------------------------   Mode array --------------------------------------------------        

        self.quadro_mode = Frame(self.win_column)
        self.quadro_mode.pack(side = TOP, anchor=W, fill = X, pady= 0, ipadx = 0, ipady = 0, expand = N)
        
        self.quadro_column = Frame(self.win_column)
        self.quadro_column.pack(side = TOP,  fill = X)
        
      #--------------------------   Header-------------------------------------------------- 
      
        #Label(self.quadro_column, text="Use       ").grid(row=0, column=0) 
        #Label(self.quadro_column, text="   ").grid(row=0, column=1)        
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
                                      text = "plot kchi(k)",
                                      width = 13)

        self.buttonMu.pack(side = LEFT, anchor = W,pady = 10, padx=5)
         
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
                self.array=np.loadtxt(fname=self.filenames,comments=self._ChaCom.get(),usecols=(0,1))
            except Exception as inst:
                print  "\n"*5+"ERROR-----"*10+"\n\nERROR  change comment character"
                print  "or redefine columns and press reload\n\n"+"ERROR-----"*10+"\n"*5 
                self.array=np.zeros((18,4))
                #prin inst
                #prin type(inst)     # the exception instance
                #prin inst.args      # arguments stored in .args
                
                
            for i,item in enumerate(self.column_names):
                    self.column_list.append(Col_line_Gen(self.quadro_column, label=item, array=self.array,
                                         row=i+1))  
        else:
            k=int(self.column_list[0]._position.get())-1
            chi=int(self.column_list[1]._position.get())-1
            self.array=np.loadtxt(fname=self.filenames,comments=self._ChaCom.get(),usecols=(k,chi))
            
        return

    def Muplot(self):
        E=self.array[:,int(self.column_list[0]._position.get())-1]
        Mu=self.array[:,int(self.column_list[1]._position.get())-1]
        graph = Graphcol()
        graph.plot( [E],[Mu], title= "graph")
        

    def opens(self):
        global fitexa_input
        if globals().has_key("fitexa_input"):
            pass
        else:        
            fitexa_input=input_fitexa("fitexa.inp")
            fitexa_input.filename="fitexa.inp"
        fitexa_input.col=[int(self.column_list[0]._position.get()),
                          int(self.column_list[1]._position.get())]
        E=self.array[:,int(self.column_list[0]._position.get())-1]                  
        fitexa_input.xranges=[min(E),max(E)]  
        os.chdir(os.path.dirname(self.filenames))
        filename =os.path.basename(self.filenames)                  
        fitexa_input.chi =os.path.basename(filename)
        fitexa_input.label=filename[:4]
        self.top.destroy()
        del self
      

            
            
                
                
      
      


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, cursors
class Graphcol:
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
       self.figsub.set_xlim(xmin=x_array[0][0]-step, xmax=x_array[0][-1]+step)
       self.canvas.draw()
       self.figsub.set_autoscale_on(False)

    def plot(self, x_array, y_array, comment= None,title=None, ylabel= "", xlabel=""):
       """ycalcurves = array, calcurves = line!!!!!!!"""
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
   


#------------------------------------------------------------------------------
def browse_single(**kwargs):
        filename = tkFileDialog.askopenfilename(**kwargs)   ###-defaultextension extension
        return filename
#------------------------------------------------------------------------------
class mymenu():
    def __init__(self,genitore):
        global statusmsg
        self.genitore=genitore
        self.genitore.option_add('*tearOff', FALSE)
        self.menubar = Menu(genitore)
        # create a pulldown menu, and add it to the menu bar
        self.filemenu = Menu(self.menubar, tearoff=0)

        ##########################  files   #########################################################
        self.filemenu.add_command(label="Open FitEXA input", command=self.openFitEXA)
        self.filemenu.add_command(label="Save FitEXA input", command=self.saveFitEXA)        
        self.filemenu.add_command(label="Open Datafile", command=self.openData) 
        
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
       Label(Top, text= "press + to access to variables").pack(side = TOP, anchor= W, pady=0)
       Label(Top, text= "no more need to keep the variable number use the name ").pack(side = TOP, anchor= W, pady=0)
       Label(Top, text= "Reff is defined as efficent radius of the path").pack(side = TOP, anchor= W, pady=0)       
       Label(Top, text= "<right button> on the path to destroy it or move").pack(side = TOP, anchor= W, pady=0)
       Label(Top, text= "<right button> on r s e definition to options").pack(side = TOP, anchor= W,pady=1)        


       
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
    def openFitEXA(self):
        global fitexa_input
        filename=browse_single(filetypes = [('input files', '.inp'), ('all files', '.*')])
        os.chdir(os.path.dirname(filename))
        filename =os.path.basename(filename)
        fitexa_input=input_fitexa(filename)
        fitexa_input.read_inp(filename)
        return

    def saveFitEXA(self):
        global fitexa_input
        namesave= tkFileDialog.asksaveasfilename(title="Save Estra input file")
        fitexa_input.write_inp(namesave)



    def openData(self):
        global data_file
        filename = browse_single(filetypes = [('exa files', '.exa'), ('all files', '.*')])
        if filename!="":
            colwin=Column_Window(filename)
        return


#######################################################################################################
########################             MINUIT()          ################################################
#######################################################################################################
class Variabile_line(object):
    def __init__(self, genitore,n,Var_istance):
        self._N=IntVar()
        self._param=StringVar()
        self._value=StringVar()
        self._error=StringVar()
        self._min  =StringVar()        
        self._max  =StringVar()
      #-------------------------------------------------------------  
        self._N.set(n)
        self._param.set(Var_istance.param)
        self._value.set(Var_istance.value)
        self._error.set(Var_istance.error)
        self._min.set(Var_istance.min)
        self._max.set(Var_istance.max)
      #-------------------------------------------------------------  
        
        
        self.frameline=Frame(genitore)
        self.frameline.pack(side=TOP, fill=X, expand=True)
        self.N_Lab=Label(self.frameline, textvariable= self._N,  width=3)
        self.N_Lab.grid(row=0,column=0, sticky=(N,E,W))
        for i, item in enumerate([self._param, self._value, self._error, self._min, self._max ]): 
            Entry(self.frameline, text= item,  width=10).grid(row=0,column=i+1, sticky=(N)) #         
        self.But_des=Button(self.frameline,#command = self.destroy,
                                      text = "del",
                                      width = 6)
        self.But_des.grid(row=0,column=6, sticky=(N,E,W))
        
    def make_string(self):
        strippo = lambda x: str(x.get()) 
        stringa=strippo(self._N).rjust(10)
        stringa+=" "+strippo(self._param).ljust(9)
        stringa+=strippo(self._value).rjust(10)
        stringa+=strippo(self._error).rjust(10)        
        stringa+=strippo(self._min).rjust(10)            
        stringa+=strippo(self._max).rjust(10)
        return stringa
            
    def destroy(self):
        for item in self.frameline.winfo_children():
            item.destroy()
        self.frameline.destroy()    
        del self
        
        

#------------------------------------------------------------------------------
class Minuit_Table(object):
    def __init__(self, genitore): 
      #-----------------------------      Declare      --------------------------------------------------
        global fitexa_input
      #-----------------------------      Define       --------------------------------------------------  
      #-----------------------------      Structure    --------------------------------------------------
        Header=Frame(genitore)
        Header.pack(side=TOP, fill=X)
        Label(Header, text= "N", width=3).grid(row=0,column=0, sticky=(N,E)) #
        for i,tex in enumerate(["param", "value", "error", "min", "max",  "destroy"]): 
            Label(Header, text= tex, width=10).grid(row=0,column=i+1, sticky=(N,E), padx=1)
        #Label(Header, text= "                   ", width=3).grid(row=0,column=7, sticky=(N,E))    
            
        parent = Frame(genitore)
        parent.pack(side=TOP, expand=N, fill=BOTH )
        self.TableFrame = Canvas(parent)
        self.TableFrame.pack(side=LEFT, expand=Y, fill=BOTH )
        s = Scrollbar(parent, orient=VERTICAL, command=self.TableFrame.yview)
        s.pack(side=LEFT, expand=N, fill=Y)
        self.TableFrame['yscrollcommand'] = s.set
        self.TableFrame2 = Frame(self.TableFrame)
        self.TableFrame2.pack(side=LEFT, expand=Y, fill=Y, anchor=W )
        self.Variable_list=list()
        self.build_line(stringa="  so2    .950      .000      0.00      1.0")
        self.build_line(stringa="  gamma  0.000      .000      0.00      5.00")        

            
        self.TableFrame.create_window(0, 0, anchor=NW, window=self.TableFrame2)
        self.TableFrame.update_idletasks()        
        self.TableFrame.config(scrollregion=self.TableFrame.bbox("all"))
        
        Button_F=Frame(genitore)
        Button_F.pack(side=TOP, fill=X, anchor=N)
        Separator(Button_F, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = N, fill = BOTH)
        Button(Button_F,
                      command = self.build_line,
                      text = "Add Line",
                      width = 10).pack(side=LEFT)
                      
        for i,item in enumerate(self.Variable_list):
            self.Variable_list[i].But_des.configure(command=lambda x=self.Variable_list[i]: self.set_Del(x))                       
                      
        
                      
        
      #-----------------------------      Functios    --------------------------------------------------
    def set_Del(self,variable):
        variable.destroy()
        self.Variable_list.remove(variable)
        for i,item in enumerate(self.Variable_list):
            item._N.set(i+1)

            
            
    def set_info(self):
        global fitexa_input
        # remove previous variables
        for item in self.Variable_list:
            item.destroy()
        self.Variable_list=[]    
        for i, item in enumerate(fitexa_input.variable_list): 
            self.Variable_list.append(Variabile_line(self.TableFrame2, i+1,item))
        self.TableFrame2.update_idletasks()    
        self.TableFrame.config(scrollregion=self.TableFrame.bbox("all"))
        self.TableFrame.update_idletasks() 
        for i,item in enumerate(self.Variable_list):
            self.Variable_list[i].But_des.configure(command=lambda x=self.Variable_list[i]: self.set_Del(x)) 
        print "MINUIT TABLE setted"    

   
            
    def read_info(self):
        print "M_T read info"
        fitexa_input.variable_list=[]
        for i, item in enumerate(self.Variable_list): 
            fitexa_input.variable_list.append(variabile(item.make_string()))               
      



            
            
            
    def build_line(self, stringa="  new    0.0     0.0   0.0    0.0   "):
        global fitexa_input
        N=str(len(self.Variable_list)+1)
        self.Variable_list.append(Variabile_line(self.TableFrame2,len(self.Variable_list) +1,
                                                 variabile(N+stringa)))   
        #link the delete button of each line with the command 
        for i,item in enumerate(self.Variable_list):
            self.Variable_list[i].But_des.configure(command=lambda x=self.Variable_list[i]: self.set_Del(x))         
        self.TableFrame2.update_idletasks()
        self.TableFrame.config(scrollregion=self.TableFrame.bbox("all"))
        self.TableFrame.update_idletasks()         
      
      
class Minuit(object):
    def __init__(self, genitore):
        global fitexa_input
      #-----------------------------      Declare      --------------------------------------------------
        self._label = StringVar()
        self._estra_input = StringVar()
        self._filedata = StringVar()
      #-----------------------------      Define      --------------------------------------------------
      #-----------------------------      geometry      --------------------------------------------------'
              
        self.genitore=Frame(genitore)
        self.genitore.pack(side = LEFT, anchor= W, expand = Y, fill = BOTH)
        self.Table_minuit=Minuit_Table(self.genitore)

        #-----------------------------      start param      -----------------------------------------------'        
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = N, fill = BOTH)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= W, expand = N, fill = BOTH) 
        #-----------------------------      start param      -----------------------------------------------'
        
        self.win_text = Frame(self.genitore) 
        self.win_text.pack(side=LEFT,expand=YES, fill=BOTH)
        Label(self.win_text, text="MINUIT COMMAND").pack(side=TOP, fill=X, ipady=4) 
        testo1="\nminimize\nexit"
        self.text_Minuit=text.ScrolledText(parent=self.win_text,text=testo1 ,hor=True, active=True)
        self.text_Minuit.text.configure(height=7, width=15)
        
        
    def set_info(self):
        self.Table_minuit.set_info()
        self.text_Minuit.settext(fitexa_input.minuit)
        
    def read_info(self):     
        self.Table_minuit.read_info()
        fitexa_input.minuit=self.text_Minuit.gettext()


        

#######################################################################################################
class pEntry(Entry):
    def __init__(self, master, **kw):
            apply(Entry.__init__, (self, master), kw)
            self.bind("<Button-3>", self.printo)

            
    def printo(self,event):
        self.delete(0,100)
        self.insert(0,paste)

#######################################################################################################
########################             SHELLS()          ################################################
#######################################################################################################
class feffpath:    
    #class param(self, )
    def __init__(self,feff):                              
        self.feff_file = str(feff)
        with open(feff) as feffpath:
            feffpath=open(feff)
            regular_espression=re.compile('feff', re.IGNORECASE)
            line=feffpath.readline()
            if re.search(regular_espression, line):pass
            else:
                raise ValueError("\n ERROR Not a regular Feff phase file")
            
            while (True):
                line=feffpath.readline()
                if line[1]=="-": break
            line=feffpath.readline()           
            line = line.split()
            self.nlegs, self.deg, self.reff = line[0:3]
            line=feffpath.readline()         
            geom=[]
            while (True):
                line=feffpath.readline()
                if line[4]=="k": break
                geom.append(line.split()[5])
            self.geom= "<->".join(geom)            
        feffpath.close()





class PATH():
    def __init__(self, genitore, n, path=None):
        #-------------------------------    declare    ----------------------------------------------
        self.genitore=genitore
        self.path=path
        self._pathlabel= StringVar()    
        self._type   = StringVar()
        #-------------------------------    declare    ----------------------------------------------
        self._pathlabel.set("path info")
        self._type.set("G")
        #-------------------------------    geometry    ----------------------------------------------   
        self.quadro_Path = LabelFrame(self.genitore, text = "Path "+str(n))
        row_pos=len(self.genitore.grid_slaves())
        self.quadro_Path.grid(column=0,row=row_pos,sticky=N+W+E, pady= 3, ipadx = 0, ipady = 0)
        #self.quadro_Path.pack(side = TOP, expand= 0, fill = X, pady= 3, ipadx = 0, ipady = 0 )
        quadro_Path1 = Frame(self.quadro_Path)  
        quadro_Path1.pack(side = TOP,  fill = X, pady= 1, ipadx = 0, ipady = 0,expand=0)      
        SType=Combobox(quadro_Path1, textvariable= self._type,
                                 state ="readonly",
                                 values=("G", "M", "C", "skip"), width=4)
        SType.bind('<<ComboboxSelected>>', self.def_Tree)
      
        SType.pack(side=LEFT, padx=3, pady=1)                         
        Button(quadro_Path1, command = self.feff_path,
                      text = "Browse",
                      width = 7).pack(side=LEFT)
        self.pathlabel=Label(quadro_Path1, textvariable =self._pathlabel, width=50)
        self.pathlabel.pack(side=LEFT, expand=N, fill=X, padx=3, pady=1)
        #-------------------------------    Set info    ------------------------------------------       
        if self.path:
            self.set_info()
        else:
            self.path_param=None
        #-------------------------------    Tree    ----------------------------------------------
        self.quadro_Path2 = Frame(self.quadro_Path)    
        self.quadro_Path2.pack(side = TOP,  fill = X, pady= 1, ipadx = 0, ipady = 0) 
        self.def_Tree()
        #-------------------------------    Pop Up menu    ----------------------------------------------        
        # create a menu
        self.popup = Menu(self.quadro_Path, tearoff=0)
        self.popup.add_command(label="Option")
        self.popup.add_command(label="Move to position:")
        self.popup.add_separator()
        self.popup.add_command(label="destroy")
        # bind it
        def bindall_pop(widget):
            widget.bind("<Button-3>", self.do_popup)
            for child in  widget.winfo_children():
                child.bind("<Button-3>", self.do_popup)
        bindall_pop(quadro_Path1)    
        bindall_pop(self.Tree.Frame1)
        
        #-------------------------------    Pop Up menu    ----------------------------------------------  
      

    def do_popup(self,event):
        # display the popup menu
        try:
            self.popup.post(event.x_root, event.y_root)
        finally:
            # make sure to release the grab (Tk 8.0a1 only)
            self.popup.grab_release()

    
    def def_Tree(self, event=None):
        if hasattr(self, "Tree"): self.Tree.destroy(); del self.Tree
        if self._type.get()=="G":
            self.Tree= Btree(self.quadro_Path2,path_param=self.path_param)
            self.Tree.pack(side=TOP, expand=True, fill=BOTH)
        elif self._type.get()=="C":
            self.Tree= BtreeCum(self.quadro_Path2,path_param=self.path_param )
            self.Tree.pack(side=TOP, expand=True, fill=BOTH)        
        
        
    def feff_path(self, filename=None):
        if not(filename):
            filename = os.path.relpath(tkFileDialog.askopenfilename())
        try:    
            self.PA_path =feffpath(filename)
        except ValueError as err:
            print err.args[0]
            return
            
        self.PA_type="feff"
        label_path = os.path.basename(self.PA_path.feff_file) + "  reff =" + str(self.PA_path.reff)+"\n"
        
        label_path+= self.PA_path.geom +" mult="+self.PA_path.deg.split(".")[0]+"  nleg="+str(self.PA_path.nlegs)
        self._pathlabel.set(label_path)
        #except:
        #    self.PA_type="pippo"
        return True
        
    def set_info(self):
        global fitexa_input
        print "\n"*5
        self._type.set(self.path.tipo)
        self.path_param={}
        #--------small labda----------
        int_or_not= lambda string: "#"+str(string) if type(string)==int else string
        #--------small labda----------
        
        self.path_param["n"]=fitexa_input.N2Name(int_or_not(self.path.Coor))
        self.path_param["r"]=fitexa_input.N2Name(int_or_not(self.path.Distance))
        self.path_param["s^2"]=fitexa_input.N2Name(int_or_not(self.path.Deb_Wal))     
        self.path_param["d_E"]=fitexa_input.N2Name(int_or_not(self.path.n_de)) 
        self.path_param["gamma"]=fitexa_input.N2Name(int_or_not(self.path.n_gamma))
        if self.path.tipo=="C" :      
            self.path_param["C2"]=fitexa_input.N2Name(int_or_not(self.path.C2))     
            self.path_param["C3"]=fitexa_input.N2Name(int_or_not(self.path.C3))  
        self.feff_path(filename=self.path.phase)

               
    def read_info(self):
        self.path= shell_path()
        self.path_param= self.Tree.get()
        self.path.tipo=self._type.get()
        #--------small labda no more lambda----------
        def int_or_not(string):
            string=string.replace("Reff", str(self.PA_path.reff))
            string=string.replace("reff", str(self.PA_path.reff))
            string=string.replace("REFF", str(self.PA_path.reff))
            string=fitexa_input.Name2N(string)
            try:
                string=int(string.replace("#",""))
            except  ValueError:
                string=string
            return string
        #--------small labda no more lambda----------    
        self.path.Coor    =int_or_not(self.path_param["n"]    )
        self.path.Distance=int_or_not(self.path_param["r"]    )
        self.path.Deb_Wal =int_or_not(self.path_param["s^2"]  )     
        self.path.n_de    =int_or_not(self.path_param["d_E"]  )  
        self.path.n_gamma =int_or_not(self.path_param["gamma"])
        if self.path.tipo=="C" :      
            self.path.C2=int_or_not(self.path_param["C2"])     
            self.path.C3=int_or_not(self.path_param["C3"]) 
        self.path.phase_type= "feff"  
        self.path.phase= self.PA_path.feff_file 
        return self.path
        
    def destroy(self):   
        for item in self.quadro_Path.winfo_children():
            item.destroy()
        self.quadro_Path.destroy()    
        del self  
       
        
        

#------------------------------------------------------------------------------
class Btree(Frame):  
    def __init__(self, master, path_param=None,  **kw):
        apply(Frame.__init__, (self, master), kw)
        #------------------------------------------------------------------------------
        #self.config(text="Btree")
        self._header= StringVar()
        self.path_param=path_param
        self._header.set("            feff6l")
        #------------------------------------------------------------------------------
        Style().layout("Flat.TButton",[("Button.label", {"sticky": "W",  "expand": 1})])
        Style().configure("Flat.TButton", font='courier 10')
        self.Frame1=Frame(self)
        self.Frame1.pack(side=TOP, padx=0, pady=0, expand=True, fill=X)
        self.Btree = Button(self.Frame1, text="[+]" ,command = self.Btree_expand,    
                                      width = 3 ,style='Flat.TButton')
        self.Btree.pack(side=LEFT, pady=0, anchor=W)
        Label(self.Frame1, textvariable=self._header).pack(side=LEFT, padx=3, pady=0, expand=False, fill=X)
        ##------------------------------------------------------------------------------


       
    def Frame_expand(self):
        if not(self.path_param):
            self._vars={}
            for item in ["n","r","s^2", "d_E","gamma","C2", "C3"]:
                self._vars[item]=StringVar()
                self._vars[item].set("")
        else:        
            self._vars={}
            for item in self.path_param:
                self._vars[item]=StringVar()
                self._vars[item].set(self.path_param[item])


        self.FrameC=Frame(self)
        self.FrameC.pack(side=TOP, anchor=N, padx=3, pady=3, expand=True, fill=BOTH)
        self.box=Shel_box(self.FrameC, self._vars)
        self.box.pack(side=TOP, expand=True, fill=BOTH)
        self.master.master.master.master.update_idletasks()    
        self.master.master.master.master.config(
                    scrollregion=self.master.master.master.master.bbox("all"))
        
    def Btree_expand(self):
        self.Btree.configure(text="[-]")
        self.Btree.configure(command=self.Btree_compress)
        if hasattr(self, "Frame_expand"):
            getattr(self, "Frame_expand")()
        
    
    def Btree_compress(self):
        self.Btree.configure(text="[+]")
        self._vars=self.box.get()
        self.box.destroy2()
        self.FrameC.destroy()
        del self.FrameC
        self.Btree.configure(command=self.Btree_expand)
    def get(self):
        try:
            self._vars=self.box.get()                    #check what changd in box
            self.path_param={}
            for item in self._vars:
                self.path_param[item]=self._vars[item].get()
        except AttributeError as pi:                     #if bbox was already closed
            #print pi, "     #########################" 
            pass
        return self.path_param    
        
            
        
class BtreeCum(Btree):
    def Frame_expand(self):
        print "\n"*5
        if not(self.path_param):
            self._vars={}
            for item in ["n","r","s^2", "d_E","gamma"]:
                self._vars[item]=StringVar()
                self._vars[item].set("")
        else:        
            self._vars={}
            for item in self.path_param:
                self._vars[item]=StringVar()
                self._vars[item].set(self.path_param[item])

                
        self.FrameC=Frame(self)
        self.FrameC.pack(side=TOP, anchor=N, padx=3, pady=3, expand=True, fill=BOTH)
        self.box=Shel_boxCum(self.FrameC, self._vars)
        self.box.pack(side=TOP, expand=True, fill=BOTH)
        self.master.master.master.master.update_idletasks()    
        self.master.master.master.master.config(
                    scrollregion=self.master.master.master.master.bbox("all"))
        


class Shel_box(Frame):
    def __init__(self, master, variab, **kw):
        self.master=master
        apply(Frame.__init__, (self, master), kw)
        #------------------------------------------------------------------------------
        self._variab=variab
        self.frame_var={}
        Separator(self, orient=VERTICAL).pack(
                               side = LEFT, anchor= W, expand = N, fill = Y, padx=8, ipadx=5) 
        self.Frame1=Frame(self)
        self.Frame1.pack(side=LEFT,  fill=BOTH, expand=Y)
        for item in  ["n","r","s^2", "duo"]:
            self.frame_var[item]=Frame(self.Frame1)
            self.frame_var[item].pack(side=TOP, expand=Y, fill=X)
        for item in  ["n","r","s^2"]:
            Label(self.frame_var[item], text=item, width=5).pack(side=LEFT, fill=BOTH)
            Entry(self.frame_var[item], textvariable=self._variab[item], width=20).pack(side=LEFT, fill=BOTH, expand=True)
        Label(self.frame_var["duo"], text="d_E", width=5).pack(side=LEFT, fill=BOTH)
        Entry(self.frame_var["duo"], textvariable=self._variab["d_E"], width=20).pack(side=LEFT, fill=BOTH, expand=True) 

        Label(self.frame_var["duo"], text="    gamma", width=10).pack(side=LEFT, fill=BOTH)
        Entry(self.frame_var["duo"], textvariable=self._variab["gamma"], width=20).pack(side=LEFT, fill=BOTH, expand=True)
        
    def destroy2(self):
        for item in self.winfo_children():
            item.destroy()   
        self.destroy()
        del self
        
    def get(self):
        for key in self._variab:
            self._variab[key].get()
        return  self._variab
        
        
class Shel_boxCum(Shel_box):            
    def __init__(self, master, variab, **kw):
        apply(Shel_box.__init__, (self, master,variab), kw)
        if self._variab.has_key("C2"):
            pass
        else:
            self._variab["C3"]=StringVar()
            self._variab["C3"].set("" )
            self._variab["C2"]=StringVar()
            self._variab["C2"].set("")   
            
        self.frame_var["duo2"]=Frame(self.Frame1)
        self.frame_var["duo2"].pack(side=TOP, expand=Y, fill=X)
        Label(self.frame_var["duo2"], text="C2", width=5).pack(side=LEFT, fill=BOTH)
        Entry(self.frame_var["duo2"], textvariable=self._variab["C2"], width=20).pack(side=LEFT, fill=BOTH, expand=True) 

        Label(self.frame_var["duo2"], text="    C3", width=10).pack(side=LEFT, fill=BOTH)
        Entry(self.frame_var["duo2"], textvariable=self._variab["C3"], width=20).pack(side=LEFT, fill=BOTH, expand=True)
    
        
            
            
                

    

#------------------------------------------------------------------------------
class SHELLS():
    def __init__(self, genitore):
        global fitexa_input
      #-----------------------------      Declare      --------------------------------------------------
        self._fit_input=StringVar()
        self._label    =StringVar()
        self._filedata =StringVar()
        self._xranges1 =StringVar()
        self._xranges2 =StringVar() 
        self._kw       =StringVar()  
        self._so2      =StringVar()
        self._eta      =StringVar()
        self._n_shell  =IntVar()   
        self._Rsig     =DoubleVar()
        self._FF       =StringVar()
        self._FF_par   ={}
        for item in ["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]:
            self._FF_par[item]=StringVar()  
      #-----------------------------      Define      --------------------------------------------------
        self._fit_input.set("")
        self._label.set("")
        self._filedata.set("")
        self._kw.set("1")  
        self._so2.set("so2")
        self._eta.set("0") 
        self._n_shell.set(0)
        self._Rsig.set(7.0)
        self._FF.set("Y")
        for item in ["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]:
            self._FF_par[item].set("0.0")

      #-----------------------------      geometry      --------------------------------------------------'
        self.genitore=Frame(genitore)
        self.genitore.pack(side = LEFT, anchor= N, expand =True, fill = BOTH)
        self.Quadro_labes=Frame(self.genitore)
        self.Quadro_labes.pack(side = TOP, anchor= W, expand = False, fill = X)
        self.QLab_fit=LabelFrame(self.Quadro_labes, text="Input file")
        self.QLab_fit.grid(row=0, column=0,pady=5, padx=5)
        self.QLab_data=LabelFrame(self.Quadro_labes, text="Data  file")
        self.QLab_data.grid(row=0, column=1,pady=5, padx=5)
        self.QLab_label=LabelFrame(self.Quadro_labes, text="Label")
        self.QLab_label.grid(row=0, column=2,pady=5, padx=5)        
        Entry(self.QLab_fit, text= self._fit_input, state="disabled", width=18).pack() #
        Entry(self.QLab_label, text= self._label, width=18).pack()
        Entry(self.QLab_data, text= self._filedata, state="disabled", width=18).pack()
        #-----------------------------      start param      -----------------------------------------------'        
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= N, expand = 0, fill = X)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= N, expand = 0, fill = X) 
        Quadro_Par1=Frame(self.genitore)
        Quadro_Par1.pack(side = TOP, anchor= W, expand = False, fill = X)  
        
        QLab_kw=LabelFrame(Quadro_Par1, text="k^w")
        QLab_kw.grid(row=0, column=0,pady=5, padx=3) 
        pEntry(QLab_kw, text= self._kw, width=4).grid(row=0, column=0,pady=5, padx=3)
        
        QLab_xranges=LabelFrame(Quadro_Par1, text="k-ranges")
        QLab_xranges.grid(row=0, column=1,pady=5, padx=3) 
        pEntry(QLab_xranges, text= self._xranges1, width=5).grid(row=0, column=0,pady=5, padx=5)      
        pEntry(QLab_xranges, text= self._xranges2, width=5).grid(row=0, column=1,pady=5, padx=5)              
        
        QLab_so2=LabelFrame(Quadro_Par1, text="So2")
        QLab_so2.grid(row=0, column=2,pady=5, padx=3) 
        pEntry(QLab_so2, text= self._so2, width=5).grid(row=0, column=0,pady=5, padx=5)  
        QLab_eta=LabelFrame(Quadro_Par1, text="eta")
        QLab_eta.grid(row=0, column=3,pady=5, padx=3) 
        pEntry(QLab_eta, text= self._eta, width=3).grid(row=0, column=1,pady=5, padx=5)
        QLab_Rsig=LabelFrame(Quadro_Par1, text="Rsig")
        QLab_Rsig.grid(row=0, column=4,pady=5, padx=3) 
        pEntry(QLab_Rsig, text= self._Rsig, width=5).grid(row=0, column=2,pady=5, padx=5)
        
        QLab_FF=LabelFrame(Quadro_Par1, text="Four.filt.")
        QLab_FF.grid(row=0, column=5, pady=5, padx=3, ipady=4)
        FF=Combobox(QLab_FF, textvariable= self._FF,
                                 state ="readonly",values=("Y", "N"), width=4)
        FF.pack(side=LEFT, padx=3, pady=1)
        Button(QLab_FF, text="Set" ,command = self.set_ft,    
                                      width = 5 
                                      ).pack(side=LEFT, padx=3, pady=1) 
      #------------------------------      shell param      -----------------------------------------------' 
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= N, expand = N, fill = X)
        Separator(self.genitore, orient=HORIZONTAL).pack(side = TOP, anchor= N, expand = N, fill = X) 

        Quadro_Par2=Frame(self.genitore)
        Quadro_Par2.pack(side = TOP, anchor= W, expand =1, fill = BOTH, pady=0)    

        self.CanvasFrame = Canvas(Quadro_Par2)
        self.CanvasFrame.pack(side=LEFT, expand=True, fill=BOTH )
        slider = Scrollbar(Quadro_Par2, orient=VERTICAL, command=self.CanvasFrame.yview)
        slider.pack(side=LEFT, expand=0, fill=Y, anchor=W, padx=0)
        self.CanvasFrame['yscrollcommand'] = slider.set
        
        self.TableFrame2 =Frame(self.CanvasFrame)
        self.TableFrame2.pack(side=LEFT, expand=0, fill=BOTH, anchor=W )
        self.TableFrame2.grid(column=0,row=0, sticky=N+E+W+S )
        self.CanvasFrame.bind('<Configure>', lambda evt: wippo(evt))
        self.width=0
        def wippo(evt):
            self.CanvasFrame.itemconfig(self.canvasframeid, width=evt.width)
        #####################--------define path list
        self.Path_list=list()

        self.canvasframeid=self.CanvasFrame.create_window(0, 0, anchor=NW, window=self.TableFrame2)
        self.CanvasFrame.update_idletasks()        
        self.CanvasFrame.config(scrollregion=self.CanvasFrame.bbox("all"))
        
        Button_F=Frame(self.genitore)
        Button_F.pack(side=TOP, fill=X, anchor=N)
        
        
        self.Path_list.append(PATH(self.TableFrame2, n=len(self.Path_list)+1))
        self.path_renumber_conf()
        ###################--------button add path
        Quadro_Par3=Frame(self.genitore)
        Quadro_Par3.pack(side = TOP, anchor= W, expand =N, fill = BOTH, pady=15)         
        Button(Quadro_Par3, text="Add one path" ,command = self.add_path,    
                                      width = 13 
                                      ).pack(side=LEFT, padx=3, pady=1,anchor =N)         
        


        
        #-------------------------------    geometry    ----------------------------------------------   









    #-----------------------------      function      -----------------------------------------------'
    def set_results(self):
        global fitexa_output    
    def set_info(self):
        self._fit_input.set(fitexa_input.filename)
        self._label.set(fitexa_input.label)
        self._filedata.set(fitexa_input.chi)
        self._xranges1.set(fitexa_input.xranges[0])
        self._xranges2.set(fitexa_input.xranges[1])
        self._kw.set(fitexa_input.kw)
        self._so2.set(fitexa_input.N2Name("#"+str(fitexa_input.so2)))
        self._eta.set(fitexa_input.eta)
        self._Rsig.set(fitexa_input.stat)        
        self._FF.set(fitexa_input.FF)
        for item in ["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]:
            self._FF_par[item].set(fitexa_input.FF_param[item])
        
        for item in self.Path_list:
            try: item.destroy()
            except TclError: pass
        self.Path_list=[]    
        for item in fitexa_input.shell:
            self.Path_list.append(PATH(self.TableFrame2, n=len(self.Path_list)+1, path=item)) 
        self.path_renumber_conf()
            
    def read_info(self):
        fitexa_input.label=self._label.get()
        fitexa_input.chi  = self._filedata.get()
        fitexa_input.xranges[0]= float(self._xranges1.get())
        fitexa_input.xranges[1]= float(self._xranges2.get())
        fitexa_input.kw  = float(self._kw.get())
        fitexa_input.so2 = int(fitexa_input.Name2N(self._so2.get())[1:])
        fitexa_input.eta = int(self._eta.get())
        fitexa_input.stat= self._Rsig.get()
        fitexa_input.FF  = self._FF.get()
        if hasattr(fitexa_input, "fitexa_input"):
            pass
        else:
            fitexa_input.FF_param=dict()
        for item in ["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]:
            fitexa_input.FF_param[item] = float(self._FF_par[item].get())        
        fitexa_input.shell=[]
        for item in self.Path_list:
            if item._type.get()!="skip": 
                fitexa_input.shell.append(item.read_info())
                try:
                    fitexa_input.shell[-1].check_integrity()
                except ValueError as pi:
                    print "\n"*5
                    raise ValueError(pi.message+" for "+item.quadro_Path["text"])
        fitexa_input.n_shell=len(fitexa_input.shell)
            
    
        
    def path_destroy(self,i):
        self.Path_list[i].destroy()
        del self.Path_list[i]
        self.path_renumber_conf()

        
    def path_renumber_conf(self):
        for i,item in enumerate(self.Path_list):
            item.quadro_Path.config(text="Path "+str(i+1))
            item.popup.entryconfig(index=3, command=  lambda x=i: self.path_destroy(x))   #define the command for popup
            item.popup.entryconfig(index=1, command=  lambda x=i: self.path_move(x))      #define the command for popup
            item.quadro_Path.grid(row=i, sticky=N+W+E, pady= 3, ipadx = 0, ipady = 0)
        self.TableFrame2.grid_columnconfigure(0,weight=1)
        self.CanvasFrame.update_idletasks()    
        self.TableFrame2.update_idletasks()
        self.CanvasFrame.config(scrollregion=self.CanvasFrame.bbox("all"))        

    def path_move(self,initial):
        final=StringVar()
        info=Toplevel()
        Label(info,text="Move from position "+str(initial+1)+" to ").grid(column=0,row=0)
        Entry(info, textvar=final,width =4).grid(column=1,row=0)
        Button(info, text="Move" ,command = info.destroy).grid(column=2,row=0)
        info.wait_window()
        try:
           final=int(final.get()) 
        except: 
           return
        final= final-1 if final<= len(self.Path_list) else len(self.Path_list)-1
        final= final if final>=0 else 0
        self.Path_list.insert(final,self.Path_list[initial])
        del self.Path_list[initial+1]
        self.read_info()
        self.set_info()
        self.path_renumber_conf()         
         
    def add_path(self):
        self.Path_list.append(PATH(self.TableFrame2, n=len(self.Path_list)+1))
        self.path_renumber_conf()
        self.Path_list[-1].Tree.Btree_expand()
            
    def set_ft(self):
        if hasattr(self, "param_win"):
            self.param_win.focus()
            return
            
        
        self.param_win = Toplevel()
        self.param_win.title("FT PARAMETER")
        dict_q_FF_par={}
        
        param_ff=LabelFrame(self.param_win,text="FT PARAMETER")
        param_ff.pack(side=TOP)
        for item in ["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]:
            dict_q_FF_par[item]=LabelFrame(param_ff, text = item)
            dict_q_FF_par[item].pack(side = LEFT,  fill = X)
            Entry(dict_q_FF_par[item], width = 5, textvariable= self._FF_par[item]).pack(
                                    side = LEFT, padx = 5, ipady = 3, fill = X)

        Button(self.param_win,
                    command = self.save_ft,
                    text = "Save FT param.",
                    width = 20).pack(side = TOP, anchor = W, padx = 5, pady = 5)
        self.param_win.protocol("WM_DELETE_WINDOW", self.save_ft)
        
    def save_ft(self):
        for item in ["k1","k2","wt","typ","apo","Rup","dr","RL","RR"]:
            self._FF_par[item].get()
        for item in self.param_win.winfo_children():    
            item.destroy()
        self.param_win.destroy()   
        del self.param_win
        
        
       
         
    
    
    

#######################################################################################################
########################             MINUITLOG()       ################################################
#######################################################################################################
class MINUIT_LOG():
    def __init__(self, genitore): 
      #-----------------------------      Declare      --------------------------------------------------
        global fitexa_input
      #-----------------------------      Define       --------------------------------------------------  
      #-----------------------------      Structure    --------------------------------------------------
        Header=Frame(genitore)
        Header.pack(side=TOP, fill=X)
        Label(Header, text= "fitexa.log", width=50).grid(row=0,column=0, sticky=(N,E),ipady=10,ipadx=5) #
   
            
        parent = Frame(genitore)
        parent.pack(side=TOP, expand=Y, fill=BOTH )
        
        self.logtext=text.ScrolledText(parent=parent, text='' ,hor=True, active=False)
        self.logtext.text.configure( width=18)
       #----------------------------- --------------------------------------------------
#######################################################################################################
########################             FITEXA_OUT()       ################################################
#######################################################################################################
class FITEXA_OUT():
    def __init__(self, genitore): 
      #-----------------------------      Declare      --------------------------------------------------
        global fitexa_input
      #-----------------------------      Define       --------------------------------------------------  
      #-----------------------------      Structure    --------------------------------------------------
        Header=Frame(genitore)
        Header.pack(side=TOP, fill=X)
        Label(Header, text= "fitexa.aut", width=50).grid(row=0,column=0, sticky=(N,E),ipady=10,ipadx=5) #
   
            
        parent = Frame(genitore)
        parent.pack(side=TOP, expand=Y, fill=BOTH )
        
        self.logtext=text.ScrolledText(parent=parent, text='' ,hor=True, active=False)
        self.logtext.text.configure( width=18)
       #----------------------------- --------------------------------------------------      
    


#######################################################################################################
#######################################################################################################
class pEntry(Entry):
    def __init__(self, master, **kw):
            apply(Entry.__init__, (self, master), kw)
            self.bind("<Button-3>", self.printo)

            
    def printo(self,event):
        self.delete(0,100)
        self.insert(0,paste)
    
#######################################################################################################
########################            GRAPH()            ################################################
#######################################################################################################
class Graph:
    def __init__(self, genitore, plotting_set):
        self.top = Tk.Toplevel()
        self.top.title(title)
        self.fig = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.top)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,  self.top )
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        self.figsub = self.fig.add_subplot(111)

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
                    self.slider = Tk.Scale(self.top, from_= 0, to=1,       #
                                                     command= self.scale,   #variable= self.sh, 
                                                     orient=Tk.VERTICAL,
                                                     label= "Shift"
                                                     )
                    self.slider.pack(side = Tk.LEFT,fill = Tk.BOTH, anchor = Tk.W,pady = 15, ipady = 0)
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
               item.set_ydata(array(self.ycalcurves[i]) + i*float(event))
       if hasattr(self, "curves"):
           for i,item in enumerate(self.curves):
               item.set_ydata(array(self.ycurves[i]) + i*float(event))
               
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
        self.link = self.canvas.mpl_connect('key_press_event', self.onspace)
        pass    
            
            
    def onspace(self,event):
        global paste
        #print 'key=\'%s\', x=%d, y=%d, xdata=%f, ydata=%f'%(
        #    event.key, event.x, event.y, event.xdata, event.ydata)
        if event.key==" ":
            paste=round(event.xdata,3)
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
            self.figsub.autoscale()                
            self.canvas.draw()           



class EstraGraph_kchi(EstraGraph):
    def __init__(self, genitore, dw):    
        """Plotting set is a dictionary with key(yattrib) link to a 
        string xattribute
        boolean   plot active
        stype line style
        """
        
        global estra_output
        self.curves={}
        self.check={}
        self.check_variable={}
        self.dw=dw            


        self.genitore=Frame(genitore)
        self.genitore.grid(column=0, row=0,sticky=N+W+SE)
        genitore.rowconfigure(0,weigh=1)
        genitore.columnconfigure(0,weigh=1)
        self.fig = matplotlib.figure.Figure(figsize=(6,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.genitore)
        toolbar = NavigationToolbar2TkAgg(self.canvas,  self.genitore)
        toolbar.update()
        toolbar.pack(side=TOP, fill=X, expand=0)        
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        self.figsub = self.fig.add_subplot(111) 
        
        self.Quadro_sel=Frame(genitore)
        self.Quadro_sel.grid(column=1, row=0,sticky=W+S+E)
        genitore.columnconfigure(1,weigh=0, minsize=28)
        s = Style()
        s.configure('Green.TCheckbutton',  foreground="#3C872B")
        s.configure('Blue.TCheckbutton',  foreground="#3F4ECD")
        s.configure('Red.TCheckbutton',  foreground="#F31D1D")
        for key in ["exp.","sim."]:
            self.check_variable[key]=IntVar() 
            self.check_variable[key].set(1)   
            self.check[key] = Checkbutton(self.Quadro_sel, 
                                          variable= self.check_variable[key],
                                          text=key,
                                          width=7,
                                          command=self.check_command)
            self.check[key].pack(side=TOP,expand=N,fill=BOTH, anchor=S)
        self.check["sim."].configure(style='Red.TCheckbutton')
        
        for key in fitexa_output.curves["chi_par"]:
            self.check_variable[key]=IntVar() 
            self.check_variable[key].set(1)   
            self.check[key] = Checkbutton(self.Quadro_sel, 
                                          variable= self.check_variable[key],
                                          text=key,
                                          width=7,style='Blue.TCheckbutton',
                                          command=self.check_command)
            self.check[key].pack(side=TOP,expand=N,fill=BOTH, anchor=S)            
            

        self.check_variable["error"]=IntVar() 
        self.check_variable["error"].set(1)
        self.check["error"] = Checkbutton(self.Quadro_sel, 
                                      variable= self.check_variable["error"],
                                      text="residual",
                                      width=7, style='Green.TCheckbutton',
                                      command=self.check_command)
        self.check["error"].pack(side=TOP,expand=N,fill=BOTH, anchor=S)
        #self.link = self.canvas.mpl_connect('key_press_event', self.onspace)
        pass
    
    def check_command(self):
        #print dir(self.figsub)
        for key in self.curves:
            self.figsub.lines.remove(self.curves[key])
        self.curves={}
        self.plot()
    
    def plot(self):
        """ycalcurves"""
        plotset={"exp.":["kchi_e", 'k+'],"sim.":["kchi_t", 'r-']}
        for key in ["exp.","sim."]:
            if self.check_variable[key].get():
                self.curves[key], = self.figsub.plot(fitexa_output.curves["k"],
                                     fitexa_output.curves[plotset[key][0]]*(fitexa_output.curves["k"]**(self.dw-1)),      
                                     plotset[key][1], label=key)

        self.dk=max(abs(fitexa_output.curves["kchi_e"])*fitexa_output.curves["k"]**(self.dw-1))
        
        i=0
        for key in fitexa_output.curves["chi_par"]:
            if self.check_variable[key].get():
                i+=1
                self.curves[key], = self.figsub.plot(fitexa_output.curves["k"],
                                     fitexa_output.curves["chi_par"][key]*(fitexa_output.curves["k"]**(self.dw-1))-i*self.dk,
                                          "b", label=key) 
                
        if self.check_variable["error"].get():
                i+=1
                self.curves["error1"], = self.figsub.plot(fitexa_output.curves["k"],
                                      fitexa_output.curves["kchi_res"]*\
                                     (fitexa_output.curves["k"]**(self.dw-1))-i*self.dk,
                                     "gx", ms=2, label="residual") 
                
        self.set_label()
        self.fig.subplots_adjust(left=0.15,  bottom=0.14)   
        self.canvas.draw()
        
    def  set_label(self): 
        self.figsub.tick_params(axis='both', labelsize='small')
        self.figsub.set_xlabel(r'$k(\AA^{-1})$')
        if self.dw==1:
            self.figsub.set_ylabel(r'$k^1\chi(k)$', position=(.20,.5))  
        if self.dw==2:
            self.figsub.set_ylabel(r'$k^2\chi(k)$', position=(.20,.5))  
        if self.dw==3:
            self.figsub.set_ylabel(r'$k^3\chi(k)$', position=(.20,.5))  

class EstraGraph_FT(EstraGraph):
    def __init__(self, genitore, dw):    
        """Plotting set is a dictionary with key(yattrib) link to a 
        string xattribute
        boolean   plot active
        stype line style
        """
        
        global estra_output
        self.curves={}
        self.check={}
        self.check_variable={}
        self.dw=dw            


        self.genitore=Frame(genitore)
        self.genitore.grid(column=0, row=0,sticky=N+W+S+E)
        genitore.rowconfigure(0,weigh=1)
        genitore.columnconfigure(0,weigh=1)
        #self.genitore.pack(side=LEFT, expand =True , fill= BOTH)
        self.fig = matplotlib.figure.Figure(figsize=(5,5), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.genitore)
        toolbar = NavigationToolbar2TkAgg(self.canvas,  self.genitore)
        toolbar.update()
        toolbar.pack(side=TOP, fill=X, expand=0)        
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        self.figsub = self.fig.add_subplot(111) 
        
        self.Quadro_sel=Frame(genitore,width=13)
        #self.Quadro_sel.pack(side=LEFT, expand =0 , anchor=S, fill=None)
        self.Quadro_sel.grid(column=1, row=0,sticky=W+S+E)
        genitore.columnconfigure(1,weigh=0, minsize=29)
        self.set_check()

        
        
    def set_check(self):
        s = Style()
        s.configure('Green.TCheckbutton',  foreground="#3C872B")
        s.configure('Blue.TCheckbutton',  foreground="#3F4ECD")
        s.configure('Red.TCheckbutton',  foreground="#F31D1D")
        for key in ["exp.","sim."]:
            self.check_variable[key]=IntVar() 
            self.check_variable[key].set(1)   
            self.check[key] = Checkbutton(self.Quadro_sel, 
                                          variable= self.check_variable[key],
                                          text=key,
                                          width=7,
                                          command=self.check_command)
            self.check[key].pack(side=TOP,expand=N,fill=BOTH, anchor=S)
        self.check["sim."].configure(style='Red.TCheckbutton')

        self.check_variable["error"]=IntVar() 
        self.check_variable["error"].set(1)
        self.check["error"] = Checkbutton(self.Quadro_sel, 
                                      variable= self.check_variable["error"],
                                      text="residual",
                                      width=7, style='Green.TCheckbutton',
                                      command=self.check_command)
        self.check["error"].pack(side=TOP,expand=N,fill=BOTH, anchor=S)
        pass
 

 
 
 
    def check_command(self):
        #print dir(self.figsub)
        for key in self.curves:
            self.figsub.lines.remove(self.curves[key])
        self.curves={}
        self.plot()
    
    def plot(self):
        """ycalcurves"""
        plotset={"exp.":["MaFT_e", "ImFT_e",'k+'],"sim.":["MaFT_t", "ImFT_t", 'r-']}
        
        
        
        for key in ["exp.","sim."]:
            if self.check_variable[key].get():
                self.curves[key+"M"], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves[plotset[key][0]],      
                                     plotset[key][2], label=key)
                self.curves[key+"I"], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves[plotset[key][1]],      
                                     plotset[key][2], label=key)

        self.df=max(fitexa_output.curves["MaFT_e"])*1.5
        
        i=0
               
        if self.check_variable["error"].get():
                i+=1
                self.curves["error1_M"],  = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves["MaFT_e"]-fitexa_output.curves["MaFT_t"]-self.df,      
                                     "x", ms=2,color= '#2BAB0D',  label="residual")
                self.curves["error1_I"], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves["ImFT_e"]-fitexa_output.curves["ImFT_t"]-self.df,      
                                     "x", ms=2, color='#36F70A',  label="residual") 
                
        self.set_label()
        self.fig.subplots_adjust(left=0.15,  bottom=0.14)        
        self.canvas.draw()
        
    def  set_label(self): 
        self.figsub.tick_params(axis='both', labelsize='small')
        self.figsub.set_xlabel(r'$k(\AA^{-1})$')
        if self.dw==1:
            self.figsub.set_ylabel(r'$k^1|\chi(R)|(\AA^{-2})$', position=(.20,.5))  
        if self.dw==2:
            self.figsub.set_ylabel(r'$k^2|\chi(R)|(\AA^{-3})$', position=(.20,.5))  
        if self.dw==3:
            self.figsub.set_ylabel(r'$k^3|\chi(R)|(\AA^{-4})$', position=(.20,.5))  



class EstraGraph_FT_M(EstraGraph_FT):
    def set_check(self):
        s = Style()
        s.configure('Green.TCheckbutton',  foreground="#3C872B")
        s.configure('Blue.TCheckbutton',  foreground="#3F4ECD")
        s.configure('Red.TCheckbutton',  foreground="#F31D1D")
        for key in ["exp.","sim."]:
            self.check_variable[key]=IntVar() 
            self.check_variable[key].set(1)   
            self.check[key] = Checkbutton(self.Quadro_sel, 
                                          variable= self.check_variable[key],
                                          text=key,
                                          width=7,
                                          command=self.check_command)
            self.check[key].pack(side=TOP,expand=N,fill=BOTH, anchor=S)
        self.check["sim."].configure(style='Red.TCheckbutton')
        
        for key in fitexa_output.curves["MaFT_par"]:
            self.check_variable[key]=IntVar() 
            self.check_variable[key].set(1)   
            self.check[key] = Checkbutton(self.Quadro_sel, 
                                          variable= self.check_variable[key],
                                          text=key,
                                          width=7,style='Blue.TCheckbutton',
                                          command=self.check_command)
            self.check[key].pack(side=TOP,expand=N,fill=BOTH, anchor=S)            
            

        self.check_variable["error"]=IntVar() 
        self.check_variable["error"].set(1)
        self.check["error"] = Checkbutton(self.Quadro_sel, 
                                      variable= self.check_variable["error"],
                                      text="residual",
                                      width=7, style='Green.TCheckbutton',
                                      command=self.check_command)
        self.check["error"].pack(side=TOP,expand=N,fill=BOTH, anchor=S)
        pass
 

 
 
 

    
    def plot(self):
        """ycalcurves"""
        plotset={"exp.":["MaFT_e", "ImFT_e",'k+'],"sim.":["MaFT_t", "ImFT_t", 'r-']}
        
        for key in ["exp.","sim."]:
            if self.check_variable[key].get():
                self.curves[key], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves[plotset[key][0]],      
                                     plotset[key][2], label=key)

        self.df=max(fitexa_output.curves["MaFT_e"])/2
        
        i=0
        for key in fitexa_output.curves["MaFT_par"]:
            if self.check_variable[key].get():
                i+=1
                self.curves[key], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves["MaFT_par"][key]-i*self.df,
                                         "b", label=key) 
                
        if self.check_variable["error"].get():
                i+=1
                self.curves["error1"],  = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves["MaFT_e"]-fitexa_output.curves["MaFT_t"]-i*self.df,      
                                     "x", ms=2, color= '#2BAB0D',  label="residual")
 
                
        self.set_label()
        self.fig.subplots_adjust(left=0.15,  bottom=0.14)        
        self.canvas.draw()
 

class EstraGraph_FT_I(EstraGraph_FT_M):
    def plot(self):
        """ycalcurves"""
        plotset={"exp.":["MaFT_e", "ImFT_e",'k+'],"sim.":["MaFT_t", "ImFT_t", 'r-']}
        
        for key in ["exp.","sim."]:
            if self.check_variable[key].get():
                self.curves[key+"I"], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves[plotset[key][1]],      
                                     plotset[key][2], label=key)
                self.curves[key+"M"], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves[plotset[key][0]],      
                                     plotset[key][2], label=key)

        self.df=max(fitexa_output.curves["MaFT_e"])
        
        i=0
        for key in fitexa_output.curves["MaFT_par"]:
            if self.check_variable[key].get():
                i+=1
                self.curves[key], = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves["ImFT_par"][key]-i*self.df,
                                         "b", label=key) 
                
        if self.check_variable["error"].get():
                i+=1
                self.curves["error1"],  = self.figsub.plot(fitexa_output.curves["R"],
                                     fitexa_output.curves["ImFT_e"]-fitexa_output.curves["ImFT_t"]-i*self.df,      
                                     "x", ms=2, color= '#2BAB0D',  label="residual")
 
                
        self.set_label()
        self.fig.subplots_adjust(left=0.15,  bottom=0.14)  
        #print dir(self.figsub)
        self.canvas.draw()

#------------------------------------------------------------------------------
class FitGRAPH():
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
      
        self.genitore=Frame(genitore, width=14)
        self.genitore.pack(side = LEFT, anchor= W, expand = N, fill = BOTH)
        self.Quadrobutton=LabelFrame( self.genitore, text="Plot")
        self.Quadrobutton.pack(side=TOP, anchor=N)
        
        
        Button(self.Quadrobutton,         command = lambda x=1: self.kchiplot(x),
                       text = u"k\u03C7(k)", 
                       width = 13).pack(side = TOP, anchor = W,pady = 5, padx=5)

        Button(self.Quadrobutton,         command = lambda x=2: self.kchiplot(x),
                       text = u"k\u00B2\u03C7(k)", 
                       width = 13).pack(side = TOP, anchor = W,pady = 5, padx=5)
                       
        Button(self.Quadrobutton,         command = lambda x=3: self.kchiplot(x),
                       text = u"k\u00B3\u03C7(k)", 
                       width = 13).pack(side = TOP, anchor = W,pady = 5, padx=5)
        Separator(self.Quadrobutton).pack(fill=X)               
                       
        Button(self.Quadrobutton,         command =  self.FT1plot,
                       text = u"FT[k\u02B7\u03C7(k)]", 
                       width = 13).pack(side = TOP, anchor = W,pady = 5, padx=5) 
        Button(self.Quadrobutton,         command =  self.FT2plot,
                       text = u"|FT[k\u02B7\u03C7(k)]|", 
                       width = 13).pack(side = TOP, anchor = W,pady = 5, padx=5)
        Button(self.Quadrobutton,         command =  self.FT3plot,
                       text = u"\u2111(FT[k\u02B7\u03C7(k)])", 
                       width = 13).pack(side = TOP, anchor = W,pady = 5, padx=5) 
                       
        Frame(self.genitore).pack(side=TOP, expand=True, fill=Y)
        
        s = Style()
        s.configure('Green.TButton')
        s.map("Green.TButton",
            foreground=[('pressed', 'red'), ('active', 'green')],
            background=[('pressed', '!disabled', 'black'), ('active', 'white')]
            )     
        

        ## custom ttk styles
        #style = Style()
        #arrow_layout = lambda dir: (
        #    [('Button.focus', {'children': [('Button.%sarrow' % dir, None)]})]
        #)
        #style.layout('L.TButton', arrow_layout('left'))
        #style.layout('R.TButton', arrow_layout('right'))        
        
        #        Style().layout("Flat.TButton",[("Button.label", {"sticky": "W",  "expand": 1})])
        #        Style().configure("Flat.TButton", font='courier 10')
        #[("Button.border", {"children": [("Button.focus", {"children": [("Button.spacing",
        #{"children": [("Button.label", {"sticky": "nswe"})], "sticky": "nswe"})], 
        #"sticky": "nswe"})], "sticky": "nswe", "border": "1"})]	


        Style().configure('Green.TButton', background="green", foreground="green")
        
        
        
        self.EditButton= Button(self.genitore,#                       command = self.Preplot,
                         text = u"Edit inp",
                         width = 13)
        self.EditButton.pack(side=TOP,anchor=S, ipady=15, pady=15)        
        self.FitButton= Button(self.genitore,#                       command = self.Preplot,
                       text = u"Fit",
                       style='Green.TButton',
                       width = 13)
        self.FitButton.pack(side=TOP,anchor=S, ipady=15)
    #-----------------------------     Functions       -----------------------------------------------------
    def kchiplot(self, w):
        global estra_input
        global estra_output
        name='k'+str(w)+'chi'
        setattr(self,name,Toplevel()) 
        getattr(self,name).title(u"plot k\u03C7(k)")        
        #--------------------------   Graphic win  --------------------------------------------------
        setattr(self,'graph_'+name,EstraGraph_kchi(getattr(self,name),  w))
        getattr(self,'graph_'+name).plot()
        
    def FT1plot(self):
        global estra_input
        global estra_output
        self.FT1=Toplevel()
        self.FT1.title(u"FT[k\u02B7\u03C7(k)]")        
        #--------------------------   Graphic win  --------------------------------------------------
        self.graph_FT1=EstraGraph_FT(self.FT1, fitexa_input.kw)
        self.graph_FT1.plot()       
        
    def FT2plot(self):
        global estra_input
        global estra_output
        self.FT2=Toplevel()
        self.FT2.title(u"FT|k\u02B7\u03C7(k)|")        
        #--------------------------   Graphic win  --------------------------------------------------
        self.graph_FT2=EstraGraph_FT_M(self.FT2, fitexa_input.kw)
        self.graph_FT2.plot() 
        
    def FT3plot(self):
        global estra_input
        global estra_output
        self.FT3=Toplevel()
        self.FT3.title(u"FT|k\u02B7\u03C7(k)|")        
        #--------------------------   Graphic win  --------------------------------------------------
        self.graph_FT3=EstraGraph_FT_I(self.FT3, fitexa_input.kw)
        self.graph_FT3.plot()           
        
#######################################################################################################
#######################################################################################################
########################                               ################################################
########################        FitEXAGui              ################################################
########################                               ################################################
#######################################################################################################
#######################################################################################################
class FitEXAGui:
    def __init__(self, genitore):
        global fitexa_output
        global fitexa_input 
        global inivar
        
        inivar=ConfigParser.ConfigParser()
        readini()
        #---------------------menu--------------------------------------
        self.menu=mymenu(genitore)
        self.menu.filemenu.entryconfig(index=0, command=  self.set_info)
        self.menu.filemenu.entryconfig(index=1, command=  self.write_info)        
        self.menu.filemenu.entryconfig(index=2, command=  self.read_data) 
        #---------------------menu--------------------------------------
        
        
        
        #---------------------part1--------------------------------------        
        Part1=Frame(genitore)
        Part1.grid(sticky=N+E+S+W)
        genitore.rowconfigure(0,weigh=1)
        genitore.columnconfigure(0,weigh=1)
        
        notebook = Notebook(Part1)
        notebook.pack(side=LEFT, expand=1 ,fill=BOTH)
        page1 = Frame(notebook)
        page2 = Frame(notebook)
        page3 = Frame(notebook)       
        page4 = Frame(notebook) 
        notebook.add(page1, text=' Shell')
        notebook.add(page2, text='Minuit Var.&Comm.')
        notebook.add(page3, text='Minuit log')
        notebook.add(page4, text='FitEXA out')    
        
        self.M_P=Minuit(page2)
        self.SHEL_FIT=SHELLS(page1)
        self.LOG=MINUIT_LOG(page3)
        self.OUT=FITEXA_OUT(page4)        
        #---------------------part1--------------------------------------         


        #---------------------part2-------------------------------------- 
        Part2=Frame(genitore)
        Part2.grid(row=0,column=1,sticky=NE+SW)
        genitore.rowconfigure(0,weigh=1)
        genitore.columnconfigure(1,weigh=0,minsize=30)
        self.GRAPH_Fit=FitGRAPH(Part2) 
        #---------------------part2--------------------------------------
        
        self.GRAPH_Fit.FitButton.config(command=self.Return)
        self.GRAPH_Fit.EditButton.config(command=self.Edit)
        #genitore.bind_all("<Return>", self.Return) 
        #genitore.bind_all("<space>", self.Return)   


  #-----------------------------      FUNCTION  GLOBAL +-    --------------------------------------------------'
        
    def set_info(self):
        self.menu.openFitEXA()
        try:
           self.SHEL_FIT.set_info()
        except Exception as inst:
            if isinstance(inst,IOError): 
                print inst.args[1], "(",inst.filename,")"
                tkMessageBox.showinfo("File ERROR",
                                  '{} ({})'.format(inst.args[1],inst.filename))
            print "\nimpossible to set SHELL information\n"
        try:    
            self.M_P.set_info()
        except:
            print "\nimpossible to set MINUIT information\n"            
        
    def write_info(self):
        self.SHEL_FIT.read_info()
        self.M_P.read_info()        
        self.menu.saveFitEXA()


    def read_data(self):
        self.menu.openData()     
        self.SHEL_FIT._fit_input.set(fitexa_input.filename)
        self.SHEL_FIT._label.set(fitexa_input.label)
        self.SHEL_FIT._filedata.set(fitexa_input.chi)
        self.SHEL_FIT._xranges1.set(fitexa_input.xranges[0])
        self.SHEL_FIT._xranges2.set(fitexa_input.xranges[1])




    def Edit(self):
        global fitexa_output
        global fitexa_input
        
        ###Check if fitexa.inp is defined
        try: fitexa_input
        except NameError:
             tkMessageBox.showinfo("Error","fitexa.inp not defined")
             return
             
        top=Toplevel(takefocus=True)
        text.InpEditor(parent=top,filex=fitexa_input.filename).pack()
        top.wait_window()
        fitexa_input.read_inp(fitexa_input.filename)
        try:
           self.SHEL_FIT.set_info()
        except Exception as inst:
            if isinstance(inst,IOError): 
                print inst.args[1], "(",inst.filename,")"
                tkMessageBox.showinfo("File ERROR",
                                  '{} ({})'.format(inst.args[1],inst.filename))
            print "\nimpossible to set SHELL information\n"
        try:    
            self.M_P.set_info()
        except:
            print "\nimpossible to set MINUIT information\n"    
        return

    def Return(self,event=None):
        global fitexa_output
        global fitexa_input
        
        ###Check if fitexa.inp is defined
        try: fitexa_input
        except NameError:
             tkMessageBox.showinfo("Error","fitexa.inp not defined")
             return
            
        self.M_P.read_info()
        try:
            self.SHEL_FIT.read_info()
        except ValueError as x:
            print x
            return  
            
        fitexa_input.write_inp()
        print 5*"\n"
        print 10*"$"+"  NEW RUN "+10*"$"
        command=os.path.join(inivar.get("Fitexa", "Fitexa_Dir"),
                             "FitEXA.exe")
        #print command ,"attenzione comando disabilitato"

        import subprocess
        p = subprocess.Popen([command], stdout=subprocess.PIPE, 
                                           stderr=subprocess.PIPE,
                                           shell=True)
        out, err = p.communicate()
        print out
        print '**********************************'
        print err

        #check good running
        if err: 
            tkMessageBox.showinfo("FitEXAFS ERROR",err)
            return
        if ('ERROR' in out) or ('error' in out):
            Estart=out.find('ERROR')or out.find('error')   
            tkMessageBox.showinfo("FitEXAFS ERROR",out[Estart:])
            return
            
        #print     os.path.splitext(fitexa_input.label)[0]
        
        fitexa_output=output_fitexa(os.path.splitext(fitexa_input.label)[0])
        if os.path.getsize('{}.out'.format(fitexa_output.label))==0:
           tkMessageBox.showinfo('FitEXAFS ERROR', 
                                 'Empty {}.out file'.format(fitexa_output.label)) 
           return
           
        fitexa_output.read_output()
        i=0
        for item in self.SHEL_FIT.Path_list:
            if item._type.get()=="skip":
                item.Tree._header.set("    not evaluated")
            else:
                item.Tree._header.set(fitexa_output.out.shells[i]["one_string"])
                fitexa_output.out.shells[i]
                i+=1      
        self.LOG.logtext.settext(filex="fitexa.log")
        self.LOG.logtext.text.see(END) 
        self.OUT.logtext.settext(filex='{}.out'.format(fitexa_output.label))
        self.OUT.logtext.text.see(END)         
        if hasattr(self.GRAPH_Fit,"graph_k1chi"):
            self.GRAPH_Fit.graph_k1chi.check_command()
        if hasattr(self.GRAPH_Fit,"graph_k2chi"):
            self.GRAPH_Fit.graph_k2chi.check_command() 
        if hasattr(self.GRAPH_Fit,"graph_k3chi"):
            self.GRAPH_Fit.graph_k1chi.check_command()    
        if hasattr(self.GRAPH_Fit,"graph_FT1"):
            self.GRAPH_Fit.graph_FT1.check_command()
        if hasattr(self.GRAPH_Fit,"graph_FT2"):
            self.GRAPH_Fit.graph_FT2.check_command()
        if hasattr(self.GRAPH_Fit,"graph_FT3"):
            self.GRAPH_Fit.graph_FT3.check_command()             
        return       
##############   Inizialization   ############################################################  
def readini():
    global inivar
    if os.name =="nt":
         path_local_data=os.path.join(os.environ['APPDATA'],"EstraFitexa")
    elif os.name =="posix":
         path_local_data="~/.local/bin"
    else :
        #if __verbose__:
        print os.name, "ERROR--"*5+"\n  sistem not defined\n" +"ERROR--"*5
        return
    inifile=os.path.join(path_local_data,"EstraFitexa.ini")
    if __verbose__:  print inifile
    inivar.read(inifile)
    if __verbose__ : print os.getcwd()
    #except :
    #   #if __verbose__ : print "no ini file found"
    #   #inivar=ConfigParser.ConfigParser()
    #   #inivar.add_section("Fitexa")
    #   #inivar.set("Fitexa", "Fitexa_Dir", os.getcwd())
    #   #writeini() 
    #   #return
    inivar.sections()
    if inivar.has_section("Fitexa"):
        if os.access(inivar.get("Fitexa", "Start_Dir"), os.F_OK):
            os.chdir(inivar.get("Fitexa", "Start_Dir"))
        else:
            os.chdir(os.path.join(os.environ['HOMEDRIVE'],os.environ['HOMEPATH']))
   
        if not(inivar.has_option("Fitexa", "Fitexa_Dir") and 
            os.access(inivar.get("Fitexa", "Fitexa_Dir"), os.F_OK)):
            inivar.set("Fitexa", "Fitexa_Dir", os.getcwd())
            
    else:
       inivar.add_section("Fitexa")
       inivar.set("Fitexa", "Fitexa_Dir", os.getcwd())
       os.chdir(os.path.join(os.environ['HOMEDRIVE'],os.environ['HOMEPATH']))
    return       
    

#################   Inizialization   ############################################################ 
def writeini():
    global inivar
    inivar.set("Fitexa", "Start_Dir", os.getcwd())
    if os.name =="nt":
         path_local_data=os.path.join(os.environ['APPDATA'],"EstraFitexa")
    elif os.name =="posix":
         path_local_data="~/.local/bin"
    if not(os.access(path_local_data, os.F_OK)):
         os.mkdir(path_local_data)   
         
    inifile=os.path.join(path_local_data,"EstraFitexa.ini")     
    with open(inifile, 'w') as configfile:
        inivar.write(configfile)     
        configfile.close
    return             
    
   
def clearini():
    global inivar
    inivar.remove_section("Fitexa")
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
   radice.title("FitEXA GUI")
   #fitexa_input=input_fitexa()
   #fitexa_input.read_inp("fitexa.inp")
   pippo = FitEXAGui(radice)
   radice.protocol("WM_DELETE_WINDOW", destroy)
   radice.mainloop()










  

        



























            
            
            
            
  
