%% compileStandaloneApp.m
% Builds TMSAnalysisToolBox as a standalone Windows desktop app (no
% MATLAB license required to run it) and packages a web installer.
%
% Run this script from MATLAB with the current folder set to this
% project ("TMS EMG KIT"), or just F5 it — paths are resolved relative
% to this script's own location either way.
%
% Requires: MATLAB Compiler (and MATLAB R2022b or later for the
% ExecutableSplashScreen option below).

projectRoot = fileparts(mfilename('fullpath'));

appFile   = fullfile(projectRoot, 'TMSAnalysisToolBox_v2_1_4.mlapp');
fcnsDir   = fullfile(projectRoot, 'Fcns');
logoDir   = fullfile(projectRoot, 'Logo');
outputDir = fullfile(projectRoot, 'WindowsOutput');

if ~isfile(appFile)
    error('compileStandaloneApp:missingApp', ...
        'Could not find %s. Update appFile at the top of this script if the .mlapp has been renamed/versioned.', appFile);
end

opts = compiler.build.StandaloneApplicationOptions(appFile, ...
    'ExecutableName',         'TMSEMGKit_v2_1_4', ...
    'AdditionalFiles',        {fcnsDir, logoDir}, ...   % include Fcns AND Logo (the app loads Logo/Group_13.png at runtime)
    'ExecutableIcon',         fullfile(logoDir, 'Group 13.png'), ...   % .exe / taskbar icon
    'ExecutableSplashScreen', fullfile(logoDir, 'Group 12.png'), ...  % shown while the app is starting up
    'OutputDir',              outputDir, ...
    'Verbose',                'on');

results = compiler.build.standaloneWindowsApplication(opts);

compiler.package.installer(results, ...
    'InstallerName', 'TMSAnalysisToolBox_Installer', ...
    'OutputDir',     outputDir);
