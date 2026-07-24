%{
AnalysisPluginExample_SP - simple example script to show users how to make plugins for silent period analysis

HOW TO USE THIS TEMPLATE
  Edit only the three ZONE blocks below:
    ZONE 1  set numVar, your parameter labels, and optional default values
    ZONE 2  unpack your parameters from UserVar
    ZONE 3  write your Analysis method and fill the outputs
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
  MissingAnalyze     - number of trials that are not analyzed due to the onset or offset not being found
  CustomOutputs      - add the custom outputs to this scruct (e.g. CustomOutputs.Latency)
  CustomAnalysisOpts - the pop-up object, returned untouched
  Return [] for any output your method does not compute.
%}

%{
Other SP metrics
iSP Latency: The time elapsed between the TMS pulse and iSP onset.
Average iSP Depth: Calculation of average iSP depth involves taking the mean EMG signal for the entire iSP duration and normalizing this depth to the average pre-stimulus EMG level.
Maximum iSP Depth: The maximum iSP depth is indicated by the pink point. Maximum iSP depth is typically normalized to the average pre-stimulus EMG level
Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC8276277/
%}

function [MissingAnalyze, CustomOutputs, CustomAnalysisOpts]=AnalysisPluginExample_SP(app,existFig,PluginsFolderName,AnalysisType,AnalyzeSampleRate,PreStimData,SelectedTrialsData,MissingAnalyze,Start,End)

% ======================= ZONE 1: your parameters =======================
numVar = 2;                                   % how many parameters you need
ListofVariableLabels = {'Pulse Time (s)','Plot TriaL (0=off)'};
DefaultValues        = [0, 0];             % first-run defaults, same order/units as labels ([] for none)

useDefaultMEPorSPAnalysis = 1;  %1-yes, 0-no, if the user would like to use their own method for calculating these metrics
%if the user would like to use their own method for calculating these metrics,
%the variables app.MEP_Amp, app.MEP_Area, app.SP_PercDecrease, app.SP_Area need to be filled in under Zone 3

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
PulseTime = UserVar{1,2}; %in seconds
PlotTrial = round(UserVar{2,2}); %plot

% =======================================================================

% ======================= ZONE 3: Analysis ======================
Analyze=1; %Boolean, whether or not to analyze the trial, if override is used, this is always 1, if override isn't used this is handled in the if statement

% Diagnostic-plot setup: only active when PlotTrial points at a real trial.
nTrials = length(SelectedTrialsData(:,1));
doPlot   = PlotTrial >= 1 && PlotTrial <= nTrials;
plotData = [];

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

        %Data to be analyzed
        AnalyzeData=SelectedTrialsData{i,1}(Start:End);

         %Auto calculate MEP or SP default metrics
        %MEP = amplitude and area under the curve
        %SP = percent decrease and normalized area under the curve
        %If the user would like to use their own method for calculating these metrics,
        %the variables app.MEP_Amp, app.MEP_Area, app.SP_PercDecrease, app.SP_Area need to be filled in in this function
        PluginAutoCalcMEPandSP(i,app,PreStimData,AnalyzeData,AnalysisType);
        %% Analyses 

        %The time elapsed between the TMS pulse and iSP onset.
        OnsetTime=app.Time(Start); %app.AllOnOffsetTime(i,1)
        CustomOutputs.Latency(i,1)=OnsetTime-PulseTime;

        %Calculation of average iSP depth involves taking the mean EMG signal for the entire iSP duration and normalizing this depth to the average pre-stimulus EMG level.
        meanEMG=mean(AnalyzeData);
        meanPreStim=mean(PreStimData(:,i));
        CustomOutputs.NormDepth(i,1)=meanEMG./meanPreStim;

        %The maximum iSP depth is indicated by the pink point. Maximum iSP depth is typically normalized to the average pre-stimulus EMG level
        [MaxDepth, MaxDepthLoc]=min(AnalyzeData);
        CustomOutputs.NormMaxDepth(i,1)=MaxDepth./meanPreStim;

        % Stash this trial's data for the optional diagnostic plot.
        if doPlot && i == PlotTrial
            plotData=struct('AnalyzeData',AnalyzeData,'AvgDepth',meanEMG,'MaxDepth',MaxDepth,'MaxDepthLoc',MaxDepthLoc,'meanPreStim',meanPreStim);

        end

    

    else
        % Stash this trial's data for the optional diagnostic plot
        if doPlot && i == PlotTrial
            plotData=struct('AnalyzeData',[],'AvgDepth',0,'MaxDepth',0,'MaxDepthLoc',0,'meanPreStim',0);

        end
        CustomOutputs.Latency(i,:)=nan;
        CustomOutputs.NormDepth(i,:)=nan;
        CustomOutputs.NormMaxDepth(i,:)=nan;
    end
end %end for each trial

% Optional diagnostic plot: shows the envelope, threshold, and detected
% onset/offset for one trial, so you can see how the inputs change detection.
if doPlot
    plotTrialFunction(plotData);
elseif PlotTrial >= 1
    warning('Plot trial %d exceeds the number of trials (%d).', PlotTrial, nTrials);
end
% =======================================================================

end %end function

% ===================== DO NOT EDIT BELOW =====================
%function for plotting trial data
function plotTrialFunction(plotData)
AnalyzeData=plotData.AnalyzeData;
avgDepth=plotData.AvgDepth;
MaxDepth=plotData.MaxDepth;
MaxDepthLoc=plotData.MaxDepthLoc;
meanPreStim=plotData.meanPreStim;
figure(); grid on;
plot(AnalyzeData,'b'); hold on;
yline(meanPreStim,'k--');
yline(avgDepth,'c--');
plot(MaxDepthLoc,MaxDepth,'ro');
legend('Data','meanPreStim','Average Depth','Max Depth');

end
