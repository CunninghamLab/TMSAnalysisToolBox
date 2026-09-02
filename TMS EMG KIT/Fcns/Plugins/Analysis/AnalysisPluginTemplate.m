%{
AnalysisPluginTemplate - template for creating analysis plugins

HOW TO USE THIS TEMPLATE
  Edit only the three ZONE blocks below:
    ZONE 1  set numVar, your parameter labels, and optional default values
    ZONE 2  unpack your parameters from UserVar
    ZONE 3  initialize custom outputs variable
    ZONE 4  write your Analysis method and fill the outputs
  Everything else is handled for you by createPluginFigure.
  Optional: add a diagnostic plot to see how inputs change detection

INPUTS (provided by the app - do not change)
  app                - handle to the main app (e.g. app.Time, in seconds)
  existFig           - whether the parameter pop-up is already open
  PluginsFolderName  - folder for this plugin's settings file
  AnalyzeSampleRate  - sample rate after processing, Hz
  PreStimData        - baseline window, numSamples x numTrials
  SelectedTrialsData - trial signals in VOLTS, numTrials x 1 cell
  MissingAnalyze     - number of trials that are not analyzed due to the onset or offset not being found
  Start              - the starting index of of the trial (index in app.Time that equals the Onset time)
  End                - the ending index of of the trial (index in app.Time that equals the Onset time)

OUTPUTS (you must return these shapes and units)
  MissingAnalyze     - number of trials that are not analyzed due to the onset or offset not being found, double scalar
  CustomOutputs      - add the custom outputs to this scruct (e.g. CustomOutputs.Latency), column vector, each row is the result from a trial
  CustomAnalysisOpts - the pop-up object, returned untouched
%}

%{
Other MEP metrics
Lat is the time interval between the pulse delivery time and the MEP onset 
Thickness is the ratio of the area under curve (AUC) to Amp
The number of turns (NT) is counted as the significant peaks occurring during the MEP Dur
The number of phased (NP) is counted by the zero-crossing points between the MEP onset and endpoint
Source: https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2024.1415257/full#sec3 (MEPFeatX paper)
%}

function [MissingAnalyze, CustomOutputs, CustomAnalysisOpts]=AnalysisPluginTemplate(app,existFig,PluginsFolderName,AnalyzeSampleRate,PreStimData,SelectedTrialsData,MissingAnalyze,Start,End)

% ======================= ZONE 1: your parameters =======================
numVar = 2;                                   % how many parameters you need
ListofVariableLabels = {'Variable 1 Label','Variable 2 Label'};
DefaultValues        = [0, 0];             % first-run defaults, same order/units as labels ([] for none)

% =======================================================================
assert(numel(ListofVariableLabels)==numVar, 'numVar must equal the number of labels.');

% --------------------------- DO NOT EDIT -------------------------------
% Builds the parameter pop-up and loads/saves this plugin's settings.
% mfilename tells the helper which settings file belongs to this plugin.
% DefaultValues pre-fills the pop-up the first time (before settings exist).
[CustomAnalysisOpts, UserVar] = createAnalysisPluginFigure(existFig, ...
    app.CustomAnalysisOpts, PluginsFolderName, numVar, ListofVariableLabels, mfilename, DefaultValues);
% -----------------------------------------------------------------------

% ======================= ZONE 2: unpack parameters =====================
% UserVar is (numVar) x 2; column 2 holds the values.
Variable1 = UserVar{1,2}; 
Variable2 = UserVar{2,2};

% =======================================================================
nTrials=length(SelectedTrialsData(:,1));

% ======================= ZONE 3: Initialize custom outputs =============
%Initialize custom outputs with nans
%EX: app.CustomOutputs.Latency=nan(nTrials,1);

% =======================================================================

% ======================== vv DO NOT EDIT vv ==============================
for i=1:length(SelectedTrialsData) %for each trial
    [Start,End,MissingAnalyze,Analyze]=checkOverride(i,app,MissingAnalyze,Start,End); %check if override is used or if the trial doesn't have an onset/offset time

    if Analyze == 1
% ======================== ^^ DO NOT EDIT ^^ ==============================

        % ==============================================================================================================================
        % ======================= ZONE 4: Analysis ======================
    
        %Data to be analyzed
        AnalyzeData=TrialDataNR(Start:End);

        %Optional Auto calculate MEP or SP default metrics
        %MEP = amplitude and area under the curve
        %SP = percent decrease and normalized area under the curve
        PluginAutoCalcMEPandSP(i,app,PreStimData,AnalyzeData,"MEP"); %type in "MEP" or "SP"

        %If the user would like to use their own method for calculating these metrics,
        %the variables app.MEP_Amp and app.MEP_Area or app.SP_PercDecrease and app.SP_Area need to be filled in within this function
        %and the above PluginAutoCalcMEPandSP() function should be removed

        %Other Analyses
       
    end %end if Analyze, don't edit the overall structure of this if statement

end %end for each trial

% =======================================================================

end %end function






