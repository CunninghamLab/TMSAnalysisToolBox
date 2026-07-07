%{
AnalysisPluginTemplate - starting point for a custom onset/offset detection method.

HOW TO USE THIS TEMPLATE
  Edit only the three ZONE blocks below:
    ZONE 1  set numVar, your parameter labels, and optional default values
    ZONE 2  unpack your parameters from UserVar
    ZONE 3  write your detection method and fill the outputs
  Everything else is handled for you by createPluginFigure.
  Optional: add a diagnostic plot to see how inputs change detection -
  see ZScoreOnsetOffsetDetect.m for a worked example.

INPUTS (provided by the app - do not change)
  app                - handle to the main app (e.g. app.Time, in seconds)
  existFig           - whether the parameter pop-up is already open
  PluginsFolderName  - folder for this plugin's settings file
  AnalyzeSampleRate  - sample rate after processing, Hz
  PreStimData        - baseline window, numSamples x numTrials
  SelectedTrialsData - trial signals in VOLTS, numTrials x 1 cell

OUTPUTS (you must return these shapes and units)
  AllOnOffsetTime - numTrials x 2: onset (col 1) & offset (col 2), SECONDS
  OnsetLimit      - 1 x numTrials, onset threshold,  VOLTS   (or [])
  OffsetLimit     - 1 x numTrials, offset threshold, VOLTS   (or [])
  meanPreStimData - 1 x numTrials baseline mean,     VOLTS   (or [])
  CustomAnalysisOpts - the pop-up object, returned untouched
  Return [] for any output your method does not compute.
%}
function [AllOnOffsetTime, OnsetLimit, OffsetLimit, meanPreStimData, CustomAnalysisOpts]=AnalysisPluginTemplate(app, existFig, PluginsFolderName, AnalyzeSampleRate, PreStimData, SelectedTrialsData)

% ======================= ZONE 1: your parameters =======================
numVar = 3;                                   % how many parameters you need
ListofVariableLabels = {'Variable 1 Label', 'Variable 2 Label', 'Variable 3 Label'};
DefaultValues        = [0, 0, 0];             % first-run defaults, same order/units as labels ([] for none)
% =======================================================================
assert(numel(ListofVariableLabels)==numVar, 'numVar must equal the number of labels.');

% --------------------------- DO NOT EDIT -------------------------------
% Builds the parameter pop-up and loads/saves this plugin's settings.
% mfilename tells the helper which settings file belongs to this plugin.
% DefaultValues pre-fills the pop-up the first time (before settings exist).
[CustomAnalysisOpts, UserVar] = createPluginFigure(app, existFig, ...
    app.CustomAnalysisOpts, PluginsFolderName, numVar, ListofVariableLabels, mfilename, DefaultValues);
% -----------------------------------------------------------------------

% ======================= ZONE 2: unpack parameters =====================
% UserVar is (numVar+1) x 2; column 2 holds the values.
% Index 1 is ALWAYS the auto-added Start Time. YOUR parameters start at 2.
StartTime = UserVar{1,2}*0.001;   % ms -> s, to match app.Time
Var1      = UserVar{2,2};
Var2      = UserVar{3,2};
Var3      = UserVar{4,2};
% =======================================================================

% ======================= ZONE 3: detection method ======================
% Loop over SelectedTrialsData (each cell = one trial, in volts) and set,
% for every trial i:
%     AllOnOffsetTime(i,1) = onset time   (seconds)
%     AllOnOffsetTime(i,2) = offset time  (seconds)
% Preallocate as NaN so any trial with no detection stays NaN.
nTrials = numel(SelectedTrialsData);
AllOnOffsetTime = nan(nTrials, 2);
OnsetLimit = []; OffsetLimit = []; meanPreStimData = [];   % set these if your method uses them

% --- PLACEHOLDER so the plugin runs as-is. Replace with your method. ---
% This dummy just marks each trial active from Start Time to the trial end.
for i = 1:nTrials
    % your detection for trial i goes here, using SelectedTrialsData{i,1}
    AllOnOffsetTime(i,1) = StartTime;        % onset  (seconds) - placeholder
    AllOnOffsetTime(i,2) = app.Time(end);    % offset (seconds) - placeholder
end
% -----------------------------------------------------------------------

% OPTIONAL: echo the snapped Start Time back into the pop-up field.
% Finds the field by grid row (row 1), so it never depends on child order.
% ef = CustomAnalysisOpts.Children.Children;
% ef = ef(arrayfun(@(c) string(class(c))=="matlab.ui.control.NumericEditField", ef));
% ef(arrayfun(@(c) c.Layout.Row==1, ef)).Value = StartTime*1000;

% ----- OPTIONAL diagnostic plot (see ZScoreOnsetOffsetDetect.m) ---------
% To let users visualize detection:
%   1) add a numeric parameter in ZONE 1, e.g. 'Plot trial (0=off)', and
%      unpack it in ZONE 2:  PlotTrial = round(UserVar{N,2});
%   2) inside the loop, stash one trial's data when i == PlotTrial;
%   3) after the loop, if PlotTrial >= 1, call your own local plot function
%      (define it below the main function, under a DO NOT EDIT banner).
% -----------------------------------------------------------------------
% =======================================================================

end %end function
