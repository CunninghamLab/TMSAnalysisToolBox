%% compileStandaloneApp.m
% Builds TMSAnalysisToolBox as a standalone desktop app (no MATLAB
% license required to run it) and packages an installer.
%
% Cross-platform: run this script ON Windows to build the Windows app,
% and ON a Mac to build the Mac app. MATLAB Compiler cannot cross-
% compile -- there is no way to produce a Mac build from Windows.
%
% This script lives in <project root>\Build\, deliberately OUTSIDE the
% Fcns folder so it is not bundled into the compiled app. It resolves
% every path from its OWN location, so you can F5 it from the editor or
% run it with any current folder -- you do not have to cd anywhere first.
%
% Requires: MATLAB Compiler (R2022b+ for ExecutableSplashScreen).

% This file sits one level down (in Build), so go up one to reach the
% project root that holds the .mlapp, Fcns and Logo.
scriptFolder = fileparts(mfilename('fullpath'));
projectRoot  = fileparts(scriptFolder);

appFile   = fullfile(projectRoot, 'TMSAnalysisToolBox_v2_1_4.mlapp');
fcnsDir   = fullfile(projectRoot, 'Fcns');
logoDir   = fullfile(projectRoot, 'Logo');

% PNG is accepted for both the icon and the splash screen, on both
% platforms. Keep these as .png -- do NOT convert to .ico/.icns.
iconFile   = fullfile(logoDir, 'Group 13.png');   % app / taskbar / dock icon
splashFile = fullfile(logoDir, 'Group 12.png');   % Windows-only startup splash

if ~isfile(appFile)
    error('compileStandaloneApp:missingApp', ...
        ['Could not find %s\n' ...
         'This script expects to live in <project root>\\Build\\. If you moved it, ' ...
         'fix projectRoot at the top, or update appFile if the .mlapp was renamed.'], appFile);
end

if ispc
    % ---- Windows ----------------------------------------------------
    % standaloneWindowsApplication (not standaloneApplication) so the
    % .exe launches straight into the GUI with no console window.
    % ExecutableSplashScreen is ONLY honored on this path.
    outputDir = fullfile(projectRoot, 'WindowsCompiled');

    opts = compiler.build.StandaloneApplicationOptions(appFile, ...
        'ExecutableName',         'Windows_TMSEMGKit_v2_1_4', ...
        'AdditionalFiles',        {fcnsDir, logoDir}, ...
        'ExecutableIcon',         iconFile, ...
        'ExecutableSplashScreen', splashFile, ...
        'OutputDir',              outputDir, ...
        'Verbose',                'on');

    results = compiler.build.standaloneWindowsApplication(opts);

else
    % ---- macOS (and Linux) ------------------------------------------
    % No standaloneWindowsApplication here, and no splash screen: the
    % splash is a Windows-only feature of MATLAB Compiler. The icon
    % still works and becomes the .app bundle icon.
    outputDir = fullfile(projectRoot, 'MacCompiled');

    opts = compiler.build.StandaloneApplicationOptions(appFile, ...
        'ExecutableName',  'MAC_TMSEMGKit_v2_1_4', ...
        'AdditionalFiles', {fcnsDir, logoDir}, ...
        'ExecutableIcon',  iconFile, ...
        'OutputDir',       outputDir, ...
        'Verbose',         'on');

    results = compiler.build.standaloneApplication(opts);
end

compiler.package.installer(results, ...
    'InstallerName', 'TMSAnalysisToolBox_Installer', ...
    'OutputDir',     outputDir);
