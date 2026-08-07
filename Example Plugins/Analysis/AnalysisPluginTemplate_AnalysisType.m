%{
AnalysisPluginTemplate_AnalysisType - simple example script to show users how to make plugins for analysis
%==================== How to name this file ====================
filename_AnalysisType.m 
filename - can be whatever the user wants
AnalysisType - 'MEP' or 'SP', to indicate if this is an MEP or SP analysis

!! AnalysisType 'Custom' has not been tested yet, Custom example coming soon
AnalysisType can also be 'Custom' if the analysis is neither MEP or SP, 
but the built-in Garvey and SD onset/offset detection methods won't work with this 
(The built-in detection methods work differently depending on if the the analysis is MEP or SP)
Using 'Custom' with an onset/offset plugin works, since the user writes the plugin based on what data they expect it to intake

%===============================================================
HOW TO USE THIS TEMPLATE
  Edit only the three ZONE blocks below:
    ZONE 1  set numVar, your parameter labels, and optional default values
    ZONE 2  unpack your parameters from UserVar
    ZONE 3  write your Analysis method and fill the outputs
    ZONE 4  fill in the app.CustomOutputs variable to trials not analyzed
  Everything else is handled for you by createPluginFigure.
  Optional: add a diagnostic plot to see how inputs change detection -
  see ZScoreOnsetOffsetDetect.m for a worked example.

INPUTS (provided by the app - do not change)
  app                - handle to the main app (e.g. app.Time, in seconds)
  existFig           - whether the parameter pop-up is already open
  PluginsFolderName  - folder for this plugin's settings file
  AnalyzeSampleRate  - sample rate after processing, Hz
  PreStimData        - baseline window, numSamples x numTrials
  SelectedTrials     - trial signals in VOLTS, numTrials x 1 cell
  MissingAnalyze     - number of trials that are not analyzed due to the onset or offset not being found
  Start              - the starting index of of the trial (index in app.Time that equals the Onset time)
  End                - the ending index of of the trial (index in app.Time that equals the Onset time)
    Note: If the user uses the override button, the Start and End values will be calculated in the main app. 
          The user does not need to define them here. 
          It would be best to add an IF statement that checks if override was used or not using the variable app.OverrideUsed

OUTPUTS (you must return these shapes and units)
  MissingAnalyze     - number of trials that are not analyzed due to the onset or offset not being found, double scalar
  CustomOutputs      - add the custom outputs to this scruct (e.g. CustomOutputs.Latency), double column vector, each row is the result from a trial
  CustomAnalysisOpts - the pop-up object, returned untouched
%}

function [MissingAnalyze, CustomOutputs, CustomAnalysisOpts]=AnalysisPluginTemplate_AnalysisType(app,existFig,PluginsFolderName,AnalysisType,AnalyzeSampleRate,PreStimData,SelectedTrialsData,MissingAnalyze,Start,End)

% ======================= ZONE 1: your parameters =======================
numVar = 3;                                   % how many parameters you need
ListofVariableLabels = {'Variable 1','Variable 2', 'Variable 3'};
DefaultValues        = [0, 0, 0];             % first-run defaults, same order/units as labels ([] for none)

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
% UserVar is (numVar+1) x 2; column 2 holds the values.
% Index 1 is ALWAYS the auto-added Start Time. YOUR parameters start at 2.
Variable1 = UserVar{1,2};
Variable2 = UserVar{2,2};

% =======================================================================

% ======================== vv DO NOT EDIT vv ==============================
Analyze=1; %Boolean, whether or not to analyze the trial, if override is used, this is always 1,
%if override isn't used this is handled in the if statement

for i=1:length(SelectedTrialsData) %for each trial
    if app.OverrideUsed == 0 %Override is not used
        %Determine the index for the onset and offset
        OnsetTime=app.AllOnOffsetTime(i,1);
        OffsetTime=app.AllOnOffsetTime(i,2);
        Tol=eps("double");
        Start=find(abs(app.Time - OnsetTime) < Tol);
        End=find(abs(app.Time - OffsetTime) < Tol);
        %if either are not found, skip this trial
        if isempty(Start) || isempty(End)
            MissingAnalyze=MissingAnalyze+1; %add to the Missing Analyze counter, MissingAnalyze is initialized at 0 in the main app
            Analyze=0; %don't analyze this trial, fill with nans
        else
            Analyze=1; %analyze this trial
        end
    else


    end
    if Analyze == 1
% ======================== ^^ DO NOT EDIT ^^ ==============================

        % ==============================================================================================================================
        % ======================= ZONE 3: Analysis ======================
        %put your analyses here
        %Define AnalyzeData=1xnumSamples vector of trial data 

        %Auto calculate MEP or SP default metrics
        %MEP = amplitude and area under the curve
        %SP = percent decrease and normalized area under the curve
        %If the user would like to use their own method for calculating these metrics,
        %the variables app.MEP_Amp, app.MEP_Area, app.SP_PercDecrease, app.SP_Area need to be filled in in this function
        PluginAutoCalcMEPandSP(i,app,PreStimData,AnalyzeData,AnalysisType);

      
        %Define app.CustomOutputs= struct of all the custom outputs


        % ==============================================================================================================================
        % ==============================================================================================================================


    else %don't analyze becuase an onset/offset wasn't found
        % ==============================================================================================================================
        % ======================= ZONE 4: Missed Analysis Fill In ======================
        %fill in app.CustomOutputs with nan
        %CustomOutputs.Variable(i,:)=nan; %fill in the custom outputs with NaNs
        % ==============================================================================================================================
        % ==============================================================================================================================

    end %end if Analyze, don't edit the overall structure of this if statement

    

end %end for each trial

% ----- OPTIONAL diagnostic plot (see ZScoreOnsetOffsetDetect.m) ---------
% To let users visualize detection:
%   1) add a numeric parameter in ZONE 1, e.g. 'Plot trial (0=off)', and
%      unpack it in ZONE 2:  PlotTrial = round(UserVar{N,2});
%   2) inside the loop, stash one trial's data when i == PlotTrial;
%   3) after the loop, if PlotTrial >= 1, call your own local plot function
%      (define it below the main function, under a DO NOT EDIT banner).


end %end function
