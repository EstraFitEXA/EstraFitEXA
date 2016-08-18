; Script generated by the Inno Setup Script Wizard.
; SEE THE DOCUMENTATION FOR DETAILS ON CREATING INNO SETUP SCRIPT FILES!

[Setup]
; NOTE: The value of AppId uniquely identifies this application.
; Do not use the same AppId value in installers for other applications.
; (To generate a new GUID, click Tools | Generate GUID inside the IDE.)
AppId={{2FFC85D8-675A-4E97-89EF-D2E556BE6815}
AppName=ESTRA_FIT
AppVersion=0.5
;AppVerName=ESTRA_FIT 05
AppPublisher= CARLO&LELLO, Inc.
AppPublisherURL=https://estrafitexa.github.io/EstraFitEXA/
AppSupportURL=http://www.example.com/
AppUpdatesURL=http://www.example.com/
DefaultDirName={pf}\ESTRA_FIT
DefaultGroupName=ESTRA_FIT
OutputBaseFilename=ESTRA_FIT_setup
Compression=lzma
SolidCompression=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Dirs]
Name: "{app}\"; Permissions: everyone-modify

[Files]
Source: "Estra_GUI.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "FitEXA_GUI.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "Estra_FitEXA_commandline.exe"; DestDir: "{app}"; Flags: ignoreversion
Source: "*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs


[INI]

Filename: {userappdata}\Roaming\EstraFitexa\EstraFitexa.ini; Section: Fitexa; Flags: uninsdeletesection
Filename: {userappdata}\Roaming\EstraFitexa\EstraFitexa.ini; Section: Estra; Flags: uninsdeletesection



; NOTE: Don't use "Flags: ignoreversion" on any shared system files
[Icons]
Name: "{group}\Estra"; Filename: "{app}\Estra_GUI.exe"
Name: "{commondesktop}\Estra"; Filename: "{app}\Estra_GUI.exe"; Tasks: desktopicon
Name: "{group}\Fit_EXA"; Filename: "{app}\FitEXA_GUI.exe"
Name: "{commondesktop}\Fit_EXA"; Filename: "{app}\FitEXA_GUI.exe"; Tasks: desktopicon
Name: "{group}\Documentation"; Filename: "{app}\Manuals\manual_Estra.pdf"
Name: "{group}\Documentation"; Filename: "{app}\Manuals\manual_Fitexa.pdf"