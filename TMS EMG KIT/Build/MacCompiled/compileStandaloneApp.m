%% compileStandaloneApp.m
% Builds TMSAnalysisToolBox as a standalone desktop app (no MATLAB
% license required to run it) and packages an installer.
%
% Cross-platform: run this script ON Windows to build the Windows app,
% and ON a Mac to build the Mac app. MATLAB Compiler cannot cross-
% compile -- there is no way to produce a Mac build from Windows.
%
% Run from MATLAB with the current folder set to this project
% ("TMS EMG KIT"), or just F5 it -- paths resolve relative to this
% script's own location either way.
%
% Requires: MATLAB Compiler (R2022b+ for ExecutableSplashScreen).

projectRoot = fileparts(mfilename('fullpath'));

appFile   = fullfile(projectRoot, 'TMSAnalysisToolBox_v2_1_4.mlapp');
fcnsDir   = fullfile(projectRoot, 'Fcns');
logoDir   = fullfile(projectRoot, 'Logo');

% PNG is accepted for both the icon and the splash screen, on both
% platforms. Keep these as .png -- do NOT convert to .ico/.icns.
iconFile   = fullfile(logoDir, 'Group 13.png');   % app / taskbar / dock icon
splashFile = fullfile(logoDir, 'Group 12.png');   % Windows-only startup splash

if ~isfile(appFile)
    error('compileStandaloneApp:missingApp', ...
        'Could not find %s. Update appFile at the top of this script if the .mlapp has been renamed/versioned.', appFile);
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
