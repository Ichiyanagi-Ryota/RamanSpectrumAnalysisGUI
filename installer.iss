[Setup]
AppId={{RamanSpectrumAnalysisGUI}}
AppName=RamanSpectrumAnalysisGUI
AppVersion=1.0.0
AppPublisher=RamanSpectrumAnalysisGUI

DefaultDirName={autopf}\RamanSpectrumAnalysisGUI
DefaultGroupName=RamanSpectrumAnalysisGUI

OutputDir=installer
OutputBaseFilename=RamanSpectrumAnalysisGUI_Setup

Compression=lzma
SolidCompression=yes

SetupIconFile=RamanSpectrumAnalysis.ico
UninstallDisplayIcon={app}\RamanSpectrumAnalysisGUI.exe


[Files]
Source: "dist\RamanSpectrumAnalysisGUI\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs


[Icons]
Name: "{group}\RamanSpectrumAnalysisGUI"; Filename: "{app}\RamanSpectrumAnalysisGUI.exe"; IconFilename: "{app}\RamanSpectrumAnalysisGUI.exe"

Name: "{autodesktop}\RamanSpectrumAnalysisGUI"; Filename: "{app}\RamanSpectrumAnalysisGUI.exe"; IconFilename: "{app}\RamanSpectrumAnalysisGUI.exe"