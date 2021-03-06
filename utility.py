#Name: rebin-algbm29.py
# Purpose: Some class for other script
# Author: C. Prestipino based on rebin-alg.py by Ken McIvor
#
#
# Copyright 2007 ESRF
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL ILLINOIS INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#






import Tkinter as Tk
import tkFileDialog
import bm29
import os
from numpy import loadtxt,savetxt, array, column_stack,log
__singlefile__=True
from scipy import interpolate



def save_singlefile(x_array,y_array,comment=None, tit="Save  as..."):
    radix = tkFileDialog.asksaveasfilename(title=tit) 
    data=datalize(x_array[0],*y_array)
    filewrite(radix,  data, comment)   

def datalize(*arg):
    #print arg
    return column_stack((arg))
    #return transpose((vstack((arg))))

def filewrite(filename,  newdata, comment=None, footers=None,fmt='%1.10f '):
    	"""function to write a bm29 file, define a series of array with the same name
    	of column defined in the file
        """
        if os.path.exists(filename):
        	filename += ".1"
        outFile = open(filename, 'w')
        if comment is None:
		   savetxt(outFile, newdata, fmt=fmt)
        else:
           outFile.writelines(comment)
           savetxt(outFile, newdata, fmt=fmt)#%1.10f %1.8f %1.8f %d  %d %d %d %d 
        if 	not(footers is None):
            outFile.writelines(footers) 
        outFile.close
        return
        
def browse_single():
        filenames = tkFileDialog.askopenfilename()   ###-defaultextension extension
        #print filenames
        return filenames
###################################################################################################
class Browse_Directory:
    def __init__(self, genitore, title_text="open all file in directory", singlefile=0):
         self.singlefile=singlefile
         self.mioGenitore = genitore
         self.labelfiletext = Tk.StringVar()
         self.filenames = list()
         # 'quadro_fileselezione'
         self.quadro_selezione = Tk.LabelFrame(genitore, text = title_text)
         self.quadro_selezione.pack(side = Tk.TOP, expand = Tk.YES, fill = Tk.X , anchor = Tk.W,
                                  ipadx = 5, ipady = 5)
         self.pulsanteA = Tk.Button(self.quadro_selezione, command = self.browse_command,
                                  text = "BROWSE....")
         self.pulsanteA.pack(side = Tk.LEFT, expand = Tk.NO, ipadx = 5, ipady = 8 )
         self._label_filename= Tk.Label(self.quadro_selezione, textvariable=self.labelfiletext)
         self._label_filename.pack(side = Tk.LEFT,ipadx = "3m", ipady = 8)

    def browse_command(self):
        self.filenames = []
        dirname = tkFileDialog.askdirectory()   ###-defaultextension extension
        os.chdir(dirname)                                               ##TAKE CARE 
        root,directory, self.filenames =os.walk(dirname).next()
        self.labelfiletext.set(os.path.basename(self.filenames[0])+" ...... "+os.path.basename(self.filenames[-1]))



#####################################################################################################
class Browse_file:
    def __init__(self, genitore, title_text="open file", singlefile=0, All_Column=True):
         self.All_Column = All_Column
         self.singlefile=singlefile
         self.mioGenitore = genitore
         self.labelfiletext = Tk.StringVar()
         self.filenames = list()
         # 'quadro_fileselezione'
         self.quadro_selezione = Tk.LabelFrame(genitore, text = title_text)
         self.quadro_selezione.pack(side = Tk.TOP, expand = Tk.YES, fill = Tk.X , anchor = Tk.W,
                                  ipadx = 5, ipady = 5)
         self.pulsanteA = Tk.Button(self.quadro_selezione, command = self.browse_command,
                                  text = "BROWSE....")
         self.pulsanteA.pack(side = Tk.LEFT, expand = Tk.NO, ipadx = 5, ipady = 8 )
         self._label_filename= Tk.Label(self.quadro_selezione, textvariable=self.labelfiletext)
         self._label_filename.pack(side = Tk.LEFT,ipadx = "3m", ipady = 8)


    def browse_command(self):
        self.filenames = []
        if self.singlefile:
            filenames = tkFileDialog.askopenfilename()   ###
            print filenames, type(filenames)
            self.filenames.append(filenames)
        else:
            filenames = tkFileDialog.askopenfilenames()
            filenames = Tk._default_root.tk.splitlist(filenames)
            self.filenames = sorted(filenames)
        self.spectra =[]
        try:
            for i in self.filenames:
                self.spectra.append(bm29.bm29file(i, All_Column=self.All_Column))
        except  bm29.FileFormatError, e:    
            for i in self.filenames:
                self.spectra.append(bm29.disperati(i))
        if self.singlefile or (len(self.spectra)==1):
            self.labelfiletext.set(self.spectra[0].name)
        else:
            self.labelfiletext.set(self.spectra[0].name+" ...... "+self.spectra[-1].name)    
        os.chdir(os.path.dirname(self.filenames[0]))





#####################################################################################################
class Browse_filename:
    def __init__(self, genitore, title_text="open file", singlefile=0):
         self.singlefile=singlefile
         self.mioGenitore = genitore
         self.labelfiletext = Tk.StringVar()
         self.filenames = list()
         # 'quadro_fileselezione'
         self.quadro_selezione = Tk.LabelFrame(genitore, text = title_text)
         self.quadro_selezione.pack(side = Tk.TOP, expand = Tk.YES, fill = Tk.X , anchor = Tk.W,
                                  ipadx = 5, ipady = 5)
         self.pulsanteA = Tk.Button(self.quadro_selezione, command = self.browse_command,
                                  text = "BROWSE....")
         self.pulsanteA.pack(side = Tk.LEFT, expand = Tk.NO, ipadx = 5, ipady = 8 )
         self._label_filename= Tk.Label(self.quadro_selezione, textvariable=self.labelfiletext)
         self._label_filename.pack(side = Tk.LEFT,ipadx = "3m", ipady = 8)

    def browse_command(self):
        self.filenames = []
        if self.singlefile:
            filenames = tkFileDialog.askopenfilename()   ###-defaultextension extension
            #print filenames, type(filenames)
            self.filenames.append(filenames)
        else:
            filenames = tkFileDialog.askopenfilenames()
            filenames = Tk._default_root.tk.splitlist(filenames)
            self.filenames = sorted(filenames)
        if self.singlefile:
            self.name=os.path.basename(self.filenames[0])
            self.labelfiletext.set(self.name)
        else:
            self.labelfiletext.set(os.path.basename(filenames[0])+" ...... "+os.path.basename(filenames[-1]))
        os.chdir(os.path.dirname(self.filenames[0]))



#####################################################################################################


######################################################################################################
class PloteSaveB():
    """a class that allows to create two buttton that create a plot or save a set of data
    the set of sata are three list or that contain
    genitore is the frame in wich will be contained
    x_array ia a list of abscissa each abscissa could be a list or an array
    y_array ia a list of ordinate each ordinate could be a list or an array 
    z_array ia a list of error if card error is True or an another set of ordinate  each one
       could be a list or an array   
    title title at the top of the plot
    comment = list of legend for each file
    legend: list with the string for the legend
    error define if z is a list of error or an another set of y
    """
    def __init__(self, genitore, x_array=[], y_array=[], z_array=[], 
                 ext=".dat" ,comment= None,  title=None
                 ,error=False, legend=None):
        self.error=error
        self.x_array= x_array
        self.y_array= y_array
        self.z_array= z_array
        self.comments= comment
        self.title  = title
        self.ext = ext
        self.legend=None
        if (comment) is None:
            self.comments = [None for i in x_array]
        self.Button_plot = Tk.Button(genitore,
                              command = self.plot,
                              text = "Plots",
                              background = "violet",
                              width = 7,
                              padx = "3m",
                              pady = "2m")
        self.Button_plot.pack(side = Tk.LEFT, padx = 5, pady =5, anchor = Tk.W)
        self.Button_save = Tk.Button(genitore,
                              command = self.save,
                              text = "Save",
                              background = "violet",
                              width = 7,
                              padx = "3m",
                              pady = "2m")
        self.Button_save.pack(side = Tk.LEFT, padx = 5, pady =5, anchor = Tk.W)
        

    def save(self):
        if len(self.x_array)>1:
            if not(__singlefile__):
                tit = "Save  as...namefile000n."+self.ext
                radix = tkFileDialog.asksaveasfilename(title=tit)
                for i in range(len(self.x_array)):
                    if self.z_array:   data=datalize(self.x_array[i],self.y_array[i],
                                                     self.z_array[i] )
                    else: data=datalize(self.x_array[i],self.y_array[i] )
                    filename= "%s%.4d.%s" % (radix, i, self.ext)
                    filewrite(filename,  data, self.comments[i])
                    return
            else:
                z=self.x_array[0]!=self.x_array
                pos=lambda x: x if type(x)==bool else x.all()
                if pos(z):
                    x0=self.x_array[0]
                    for i_x, j_y in enumerate(self.x_array):
                        splinexy = interpolate.splrep(j_y, self.y_array[i_x])
                        self.y_array[i_x]=interpolate.splev(x0,splinexy)
        self.save_singlefile()           
                        
                        
                    
    def save_singlefile(self):
        if len(self.x_array)>1:
            tit = "Save  as..."
            radix = tkFileDialog.asksaveasfilename(title=tit) 
            data=datalize(self.x_array[0],*self.y_array)
            filewrite(radix,  data, self.comments[0])   
            if self.z_array: 
                data=datalize(self.x_array[0],*self.z_array)
                filewrite("fit"+radix,  data, self.comments[0])          
        if len(self.x_array)==1:
            tit = "Save  as..."
            radix = tkFileDialog.asksaveasfilename(title=tit) #parent=root,filetypes=myFormats ,
            if self.z_array: data=datalize(self.x_array[0],self.y_array[0],self.z_array[0])
            else:            data=datalize(self.x_array[0],self.y_array[0])
            filewrite(radix,  data, self.comments[0])
            
    def plot(self, title=None):
        comment=self.legend
        if title==None:title= self.title 
        if self.x_array==[]:
            raise ValueError("\n \n array not defined, press perform \n\n")
        if self.y_array==[]:
            raise ValueError("\n \n array not defined, press perform\n\n")  
        if self.error:
                    self.graph = Graph(self.title)
                    self.graph.errorbar(self.x_array, self.y_array,
                                        self.z_array,title= self.title)    
        else: 
            self.graph = Graph(self.title)
            if self.z_array:
                self.graph.plot(self.x_array, self.z_array, title= self.title)
            if len(self.x_array)==1 and len(self.y_array)>1:
                self.x_array=self.x_array*len(self.x_array)
            self.graph.plot(self.x_array, self.y_array, comment=comment, title=title)
            

#######################################################################################################
class LabelEntry():
    def __init__(self, genitore, Ltext = " ", EtextVar= " ", Ewith = 8, SLtext ="", 
                  labelframepack ={"side" : Tk.LEFT, "expand" : Tk.YES, 
                                   "fill" : Tk.BOTH, "anchor" : Tk.N, "ipadx" : 5, "ipady" : 5},
                  entrypack ={"side" : Tk.LEFT, "padx" : 5, "ipady" : 3}, 
                  labelpack = {"side" : Tk.LEFT, "pady":0}):
        self.Ltext = Ltext
        self.EtextVar = EtextVar
        self.Ewith = Ewith
        self.LabFr = Tk.LabelFrame(genitore, text = self.Ltext)
        self.LabFr.pack(**labelframepack)
        self._entry= Tk.Entry(self.LabFr, width = self.Ewith, textvariable= self.EtextVar)    
        self._entry.pack(**entrypack)
        self.Label = Tk.Label(self.LabFr, text= SLtext)
        self.Label.pack(**labelpack)
        
        

#######################################################################################################
class Browsefile_plot(PloteSaveB,Browse_filename):
    """
    define a method in wich browse and open a single file, adding a button to plot 
    """
    def __init__(self, genitore, title_text="open file", singlefile=0,title="",legend=None):
        self.title=title
        self.singlefile=singlefile
        self.mioGenitore = genitore
        self.labelfiletext = Tk.StringVar()
        self.filenames = list()
        self.legend=legend
        # 'quadro_totale'
        self.quadro1=  Tk.LabelFrame(genitore)
        self.quadro1.pack(side = Tk.TOP, expand = Tk.YES, fill = Tk.X , anchor = Tk.N,
                              ipadx = 5, ipady = 5)
        #-----------------------pulsant browse ------------------------------------
        self.quadro_Br= Tk.Frame(self.quadro1)
        self.quadro_Br.pack(side = Tk.LEFT, expand = Tk.YES, fill = Tk.X , anchor = Tk.W,
                                  ipadx = 5, ipady = 5)
        self.fse= Browse_filename(self.quadro_Br, title_text="open file", singlefile=1)
        self.fse.pulsanteA.configure(command= self.browse_command2)    

        #-----------------------pulsant plot-----------------------------------------
        self.quadro_Bu= Tk.Frame(self.quadro1 )
        self.quadro_Bu.pack(side = Tk.RIGHT, expand = Tk.NO, fill = Tk.X , anchor = Tk.W,
                              ipadx = 5, ipady = 7)
        self.Button_plot = Tk.Button(self.quadro_Bu,
                              command = self.plot,
                              text = "Plots",
                              background = "violet",
                              width = 7,
                              padx = "3m",
                              pady = "2m")
        self.Button_plot.pack(side = Tk.RIGHT, padx = 5, pady =5, anchor = Tk.W)
        self.error=False
        self.z_array=False
       
        
   

    def browse_command2(self):
        self.fse.browse_command()
        self.spectra=bm29.openSinglefile(self.fse.filenames[0])
        self.x_array=[item.x for item in self.spectra]
        self.y_array=[item.y for item in self.spectra]        






#######################################################################################################

import matplotlib
matplotlib.interactive(False)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg, cursors

class Graph:
    def __init__(self,title=None):
        #self.sh = Tk.StringVar()
        self.top = Tk.Toplevel()
        self.top.title(title)
        #self.top.protocol("WM_DELETE_WINDOW", self.topcallback)
        self.fig = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = self.top)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,  self.top )
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
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
        
        
class ParamGraph:
    """
    class to have a  graph windows interface with some line to obtain some values
    genitore = quadro in wich insert the graph, (often is a top windows)
    plotting_list = list with all the object containg the data
    xattr= string that defining the attribute in wich is contained the abscissa
    yattr= list of strin defining the attributes for ordinates
    """
    def __init__(self, genitore, plotting_list, xattr, yattr):    
        """yxattr  list of string rappresenting attributes
           xattr just a string
        """
        self.plotting_list = plotting_list
        self.xattr, self.yattr= xattr, yattr
        self.fig = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = genitore)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,  genitore)
        self.toolbar.update()
        self.toolbar.pack(side=Tk.TOP, fill=Tk.X, expand=0)        
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.figsub = self.fig.add_subplot(111)
        if (len(plotting_list)>1):
            self.slider = Tk.Scale(genitore, from_= 0, to=len(plotting_list)-1,
                                         command =self.panor, 
                                         orient=Tk.HORIZONTAL,
                                         label= "Shift"
                                         )
            self.slider.pack(side = Tk.TOP,fill = Tk.X, anchor = Tk.N,pady = 5, ipady = 0)
        pass
    
    

    

    def clear(self):
        self.figsub.clear()
        if hasattr(self, "curves"): del self.curves
        self.canvas.draw() 
    
    
    
                            
                                                
    def plot(self,num):
        """ycalcurves"""
        self.num=num
        self.curves=[] 
        for item in self.yattr:
            self.curves += self.figsub.plot(
                                 getattr(self.plotting_list[num], self.xattr),
                                 getattr(self.plotting_list[num], item))
        self.figsub.set_xlim(xmax = max(getattr(self.plotting_list[num], self.xattr)+15),
                             xmin=min(getattr(self.plotting_list[num], self.xattr)-15))
        try:
            maxy=max(map(max,[getattr(item, self.yattr[0]) for item in self.plotting_list]))
            miny=min(map(min,[getattr(item, self.yattr[0]) for item in self.plotting_list]))
        except:
            maxy=max(getattr(self.plotting_list[num], self.yattr[0]))
            miny=min(getattr(self.plotting_list[num], self.yattr[0])) 
        step= abs(maxy-miny)*.1
        self.figsub.set_ylim(ymax = maxy+step,ymin=miny-step)  
                 
    def paramplot(self,param,color=None, keys=None):
        self.param=param
        self.keys=keys
        self.parmlines=[]
        if keys:self.keystxt=[]
        if color==None:
            color=["r"]*len(param)
        elif len(color)!=len(param):
            color=["r"]*len(param) 
        for i,item in enumerate(self.param):
            self.parmlines.append(self.figsub.axvline(item, c=color[i],  picker=5))
            if keys:
                self.keystxt.append(self.figsub.text(x=item, y=max(self.curves[0]._y),
                                                     s=" "+self.keys[i], fontsize=12, color=color[i]))
        #self.pick=self.canvas.mpl_connect('pick_event', self.onpick)
        #self.release=self.canvas.mpl_connect('button_release_event', self.onrelease)        
    
    
    
    #def onpick(self,event_p):
    #    if not(event_p.artist in self.parmlines): 
    #        return True
    #    else:
    #        cursors.POINTER=0
    #        self.param_num= self.parmlines.index(event_p.artist)
    #        self.mov_link=self.canvas.mpl_connect('motion_notify_event', self.onmoving)
    #        print self.onmoving
    #        return True
    #
    #
    def onmoving(self,event_p):
        print "ppppppppppppppppppppppp", event_p.xdata
        self.param[self.param_num]=event_p.xdata
        self.panor(self.num)
        
    #def onrelease(self,event_r):
    #    cursors.POINTER=1
    #    try:self.canvas.mpl_disconnect(self.mov_link)
    #    except: pass   
    
    
                
    def panor(self,event):
        """
        refresh for event of slider
        """
        self.num=int(event)
        if hasattr(self, "curves"):
            for i,item in enumerate(self.curves):        
                item.set_xdata(getattr(self.plotting_list[self.num], self.xattr))
                item.set_ydata(getattr(self.plotting_list[self.num], self.yattr[i]))      
        if hasattr(self, "parmlines"):     
            for i,item in enumerate(self.parmlines):        
                item.set_xdata(self.param[i])
                if self.keys:
                    self.keystxt[i].set_position((float(self.param[i]), max(self.curves[0]._y)))
        self.canvas.draw()
            
           
           
       #self.figsub.set_ylabel(ylabel, fontsize = 8)
       #self.figsub.set_xlabel(xlabel, fontsize = 8)
       #if any(comment): 
       #     self.figsub.legend()
       #if (title): self.figsub.set_title(title)
       #self.toolbar.update()
       #self.step=max(y_array[0])-min(y_array[0])
       #self.figsub.set_ylim(ymin=(min(self.curves[0]._y)-self.step/5))
       #self.canvas.draw()
       #
       #if len(self.curves)>1:
       #    self.slider.configure(to = self.step, resolution =self.step/100)    
    
class ParamGraph_multi(Graph,ParamGraph):
    """more curve on the same paramGraph
    """
    def __init__(self,genitore,title=None):
        #self.sh = Tk.StringVar()
        #self.top.protocol("WM_DELETE_WINDOW", self.topcallback)
        self.fig = matplotlib.figure.Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master = genitore)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas,  genitore)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=Tk.LEFT, fill=Tk.BOTH, expand=1)
        self.figsub = self.fig.add_subplot(111)
        self.num=0
        
        
    def plot(self, x_array, y_array, comment= None,title=None, ylabel= "", xlabel=""):
       """ycalcurves = array, calcurves = line!!!!!!!"""
       #print len(x_array)
       #print "*******************"
       #print len(y_array[0])
       if (comment) is None:
            comment = [None for i in x_array]
       self.curves = self.figsub.plot( x_array[0], y_array[0], label= comment[0] )  #,
       for i in  range(len(x_array)-1):
           self.curves += self.figsub.plot(x_array[i+1], y_array[i+1], label = comment[i+1])
       self.figsub.set_ylabel(ylabel, fontsize = 8)
       self.figsub.set_xlabel(xlabel, fontsize = 8)
       if any(comment): 
            self.figsub.legend()
       if (title): self.figsub.set_title(title)
       self.toolbar.update()
       self.figsub.set_xlim(xmax = max(x_array[0])+5,
                               xmin=min(x_array[0])-5) 
       self.canvas.draw()
       self.figsub.set_autoscale_on(True)
    
    def panor(self,event):
        """
        refresh for event of slider
        """
        self.num=int(event)
        if hasattr(self, "parmlines"):     
            for i,item in enumerate(self.parmlines):        
                item.set_xdata(self.param[i])
                if self.keys:
                    self.keystxt[i].set_position((float(self.param[i]), max(self.curves[0]._y)))
        self.canvas.draw()    
       
     
      
    
#######################################################################################################
from scipy.optimize import leastsq
from scipy import interpolate

def fitcalibration(x1=None, y1=None,x2=None,y2=None, param=None, esterror=None, splinex1y1=None):
    """ input:
        x1,y1 reference data
        x2,y2 data to calibrate
        param list of coefficent to use for calibration
              the dimension of param give the kind of calibration max len 3
              param[0] =constant, param[1] =linear, param[2]=quadratic
        do fit using Levenberg-Marquardt
        function =function, variable to 
        spliney1 = interpolate.splrep(self.E[L1:L2], self.Mu[L1:L2], s=s)
    """
        
    if (x1==None)&(y1==None)&(splinex1y1==None):
        raise TypeError('neither x and y, not splinex1y1 for reference \
                                                   has been definited')
    if (x2==None)&(y2==None):
        raise TypeError('x or y for sample has been definited')    
    if len(param)==3:
        calib= lambda  para, x: para[2]*x**2+para[1]*x+para[0]
    elif len((param))==2:
        calib= lambda  para, x: para[1]*x+para[0]
    elif len((param))==1:
        calib= lambda  para, x: para[0]+x    
    else:
        raise ValueError("at list constan correctio len(param)=0")
    
    if not(splinex1y1):
        splinex1y1 = interpolate.splrep(x1,y1)
    fitfunction = lambda para,x: interpolate.splev(calib(para,x),splinex1y1)
    #print fitfunction(param,x2)
    errorfunction  =lambda  param, x, y: fitfunction(param,x)-y 
    p2,cov=leastsq(errorfunction,  param, args=(x2[10:-30],y2[10:-30]))
    return p2            
#######################################################################################################
