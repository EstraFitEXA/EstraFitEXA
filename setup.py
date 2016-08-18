from cx_Freeze import setup, Executable

#includes = ["scipy.linalg","numpy"]
#copyDependentFiles=True

com_file=['Tutorial_XAFS',
          'MANUALS',
          'FitEXA.exe',
          'ESTRA.exe',
          'ESTRA.iss',
          'FitEXA_Gui.ico',
          'Estra_Gui.ico',
          'FitEXA_Gui.py',
          'Estra_Gui.py',  
          'text.py',        
          'Readme.txt']





excluded_mod=["collections.abc","tzdata","PyQt4","PyQt4.QtGui","win32gui","pywin", "tcl", 
              "pywin.debugger", "pywin.debugger.dbgcon","pywin.dialogs",
              "pywin.dialogs.list", "win32com.server","email","wx"]
included_mod=['FileDialog']#"Tix"]

setup(
        name = "ESTRA % FITEXA",
        version = "0.7",
        author='Meneghini Prestipino',
        author_email='carmelo.prestipino@univ-rennes1.fr',
        url='https://github.com/EstraFitEXA/EstraFitEXA',
        description = "EXAFS data analysis software",
        options = {"build_exe": {"icon" : "./Estra_Gui.ico",
                                 "includes" : included_mod,
                                 "excludes" : excluded_mod,                                 
                                 "optimize" : 2,
                                 "compressed" : True,
                                 "include_files":com_file}
                                 },
        executables = [Executable("Estra_Gui.py"),
                       Executable("FitEXA_Gui.py"),
                       Executable("Estra_FitEXA_commandline.py")
                       ]
        )        
        
        
        
        
        
        
#"includes": includes,   #"excludes": excludes,#"packages": packages,#"path": path
#"includes": includes,   #"excludes": excludes,#"packages": packages,#"path": path