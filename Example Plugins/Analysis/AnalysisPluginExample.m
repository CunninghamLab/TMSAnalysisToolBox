%{
AnalysisPluginExample - simple example script to show users how to make analysis plugins

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
function [AllOnOffsetTime, OnsetLimit, OffsetLimit, meanPreStimData, CustomAnalysisOpts]=AnalysisPluginExample(app, existFig, PluginsFolderName, AnalyzeSampleRate, PreStimData, SelectedTrialsData)

% ======================= ZONE 1: your parameters =======================
numVar = 2;                                   % how many parameters you need
ListofVariableLabels = {'Onset Threshold (mV)','Offset Threshold (mV)'};
DefaultValues        = [0.2, 0.3];             % first-run defaults, same order/units as labels ([] for none)
% =======================================================================
assert(numel(ListofVariableLabels)==numVar, 'numVar must equal the number of labels.');

% --------------------------- DO NOT EDIT -------------------------------
% Builds the parameter pop-up and loads/saves this plugin's settings.
% mfilename tells the helper which settings file belongs to this plugin.
% DefaultValues pre-fills the pop-up the first time (before settings exist).
[CustomAnalysisOpts, UserVar] = createAnalysisPluginFigure(app, existFig, ...
    app.CustomAnalysisOpts, PluginsFolderName, numVar, ListofVariableLabels, mfilename, DefaultValues);
% -----------------------------------------------------------------------

% ======================= ZONE 2: unpack parameters =====================
% UserVar is (numVar+1) x 2; column 2 holds the values.
% Index 1 is ALWAYS the auto-added Start Time. YOUR parameters start at 2.
StartTime = UserVar{1,2}*0.001;   % ms -> s, to match app.Time
OnsetThreshold = UserVar{2,2}*0.001; %convert to mV, data is in volts, so convert to mV
OffsetThreshold = UserVar{3,2}*0.001; %convert to mV, data is in volts, so convert to mV

% =======================================================================

% ======================= ZONE 3: detection method ======================
% Loop over SelectedTrialsData (each cell = one trial, in volts) and set,
% for every trial i:
%     AllOnOffsetTime(i,1) = onset time   (seconds)
%     AllOnOffsetTime(i,2) = offset time  (seconds)
% Preallocate as NaN so any trial with no detection stays NaN.
nTrials = numel(SelectedTrialsData);
AllOnOffsetTime = nan(nTrials, 2);
OnsetLimit = nan(1,nTrials); OffsetLimit = nan(1,nTrials); meanPreStimData = [];   % set these if your method uses them

%--------------------------------------------------------------------------
%Determine Start time, round to nearest time value
Tol=eps("double");
DiffTimeStart=abs(app.Time-(StartTime));
minDiffTimeStart=find(DiffTimeStart == min(DiffTimeStart));
StartTime=app.Time(minDiffTimeStart); %seconds
CustomAnalysisOpts.Children.Children(2).Value=StartTime*1000; %display new start time in ms
%CustomAnalysisOpts.Children.Children = the edit fields for the pop-up figure

SelectedTrialsData=cellfun(@abs,SelectedTrialsData,'UniformOutput',false);
TrialTime=app.Time;
for i=1:length(SelectedTrialsData) %for each trial
    TrialData=SelectedTrialsData{i,1};
    StartIndx=find(abs(app.Time - StartTime) < Tol);
    TrialDataOn=TrialData(StartIndx:end);
    TrialTimeOn=TrialTime(StartIndx:end);

    %Find Onset
    for p=2:length(TrialDataOn)
        CurrPoint=TrialDataOn(p);
        PrevPoint=TrialDataOn(p-1);

        if CurrPoint > OnsetThreshold && PrevPoint < OnsetThreshold %if the data crosses the threshold, mark as offset
            OnsetTime=TrialTimeOn(p);
            AllOnOffsetTime(i,1)=OnsetTime;
            break;
        end

        if p==length(TrialDataOn) %reach the end without finding an onset
            AllOnOffsetTime(i,1)=nan;
        end
    end


    %Find Offset
    if isnan(AllOnOffsetTime(i,1))%if no onset was found, don't look for an offset
        AllOnOffsetTime(i,2)=nan;
    else
        TrialDataOff=TrialData(TrialTime >= OnsetTime);
        TrialTimeOff=TrialTime(TrialTime >= OnsetTime);
        for p=2:length(TrialDataOff)
            CurrPoint=TrialDataOff(p);
            PrevPoint=TrialDataOff(p-1);

            if CurrPoint < OffsetThreshold && PrevPoint > OffsetThreshold %if the data crosses the threshold, mark as offset
                AllOnOffsetTime(i,2)=TrialTimeOff(p);
                break;
            end

            if p==length(TrialDataOff) %reach the end without finding an onset
                AllOnOffsetTime(i,2)=nan;
            end
        end

    end


OnsetLimit(i)=OnsetThreshold; 
OffsetLimit(i)=OffsetThreshold;

end %end for each trial

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
