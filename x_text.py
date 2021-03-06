from Tkinter import * 
from tkSimpleDialog import askstring
from tkFileDialog   import asksaveasfilename

from tkMessageBox import askokcancel          

class Quitter(Frame):                        
    def __init__(self, parent=None):          
        Frame.__init__(self, parent)
        self.pack()
        widget = Button(self, text='Quit', command=self.quit)
        widget.pack(expand=YES, fill=BOTH, side=LEFT)
    def quit(self):
        ans = askokcancel('Verify exit', "Really quit?")
        if ans: Frame.quit(self)


class ScrolledText(Frame):
    def __init__(self, parent=None, text='', filex=None, active=True, hor=True):
        Frame.__init__(self, parent)
        self.active=active
        self.pack(expand=YES, fill=BOTH)               
        self.makewidgets(hor)
        self.settext(text, filex)
    def makewidgets(self,hor):
        sbar = Scrollbar(self)
        text = Text(self, relief=SUNKEN, font=("courier",9))
        sbar.config(command=text.yview)                  
        text.config(yscrollcommand=sbar.set) 
        if hor:
            hbar = Scrollbar(self,orient=HORIZONTAL)
            hbar.config(command=text.xview)                  
            text.config(xscrollcommand=hbar.set, wrap=NONE)  
            hbar.pack(side=BOTTOM, fill=X)   
        sbar.pack(side=RIGHT, fill=Y)                   
        text.pack(side=LEFT, expand=YES, fill=BOTH)     
        self.text = text
    def settext(self, text='', filex=None):
        if filex: 
            text = open(filex, 'r').read()
        if not(self.active):
            self.text.config(state=NORMAL)            
        self.text.delete('1.0', END)                   
        self.text.insert('1.0', text)                  
        self.text.mark_set(INSERT, '1.0')
        if not(self.active):
            self.text.config(state=DISABLED)
        self.text.focus()                                
    def gettext(self):                               
        return self.text.get('1.0', END+'-1c')         



class SimpleEditor(ScrolledText):                        
    def __init__(self, parent=None, filex=None): 
        frm = Frame(parent)
        frm.pack(fill=X)
        Button(frm, text='Save',  command=self.onSave).pack(side=LEFT)
        Button(frm, text='Cut',   command=self.onCut).pack(side=LEFT)
        Button(frm, text='Paste', command=self.onPaste).pack(side=LEFT)
        Button(frm, text='Find',  command=self.onFind).pack(side=LEFT)
        Quitter(frm).pack(side=LEFT)
        ScrolledText.__init__(self, parent, filex=filex) 
        self.text.config(font=('courier', 9, 'normal'))
        self.filex=filex
    def onSave(self):
        filename = asksaveasfilename()
        if filename:
            alltext = self.gettext()                      
            open(filename, 'w').write(alltext)          
    def onCut(self):
        text = self.text.get(SEL_FIRST, SEL_LAST)        
        self.text.delete(SEL_FIRST, SEL_LAST)           
        self.clipboard_clear()              
        self.clipboard_append(text)
    def onPaste(self):                                    
        try:
            text = self.selection_get(selection='CLIPBOARD')
            self.text.insert(INSERT, text)
        except TclError:
            pass                                      
    def onFind(self):
        target = askstring('SimpleEditor', 'Search String?')
        if target:
            where = self.text.search(target, INSERT, END)  
            if where:                                    
                print where
                pastit = where + ('+%dc' % len(target))   
               #self.text.tag_remove(SEL, '1.0', END)     
                self.text.tag_add(SEL, where, pastit)     
                self.text.mark_set(INSERT, pastit)         
                self.text.see(INSERT)                    
                self.text.focus() 
                
class InpEditor(SimpleEditor):
    def __init__(self, parent=None, filex=None): 
        frm = Frame(parent)
        frm.pack(fill=X)
        Button(frm, text='Save',  command=self.onSave).pack(side=LEFT)
        Button(frm, text='Cut',   command=self.onCut).pack(side=LEFT)
        Button(frm, text='Paste', command=self.onPaste).pack(side=LEFT)
        ScrolledText.__init__(self, parent, filex=filex) 
        self.text.config(font=('courier', 9, 'normal'))
        self.filex=filex        
    def onSave(self):
        alltext = self.gettext()                      
        open(self.filex, 'w').write(alltext)
        

if __name__ == '__main__':
    try:
        SimpleEditor(file=sys.argv[1]).mainloop()   
    except IndexError:
        SimpleEditor().mainloop()       