%{
AnalysisPluginExample_MEP - simple example script to show users how to make plugins for MEP analysis

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

%{
Other MEP metrics
Lat is the time interval between the pulse delivery time and the MEP onset 
Thickness is the ratio of the area under curve (AUC) to Amp
The number of turns (NT) is counted as the significant peaks occurring during the MEP Dur
The number of phased (NP) is counted by the zero-crossing points between the MEP onset and endpoint
Source: https://www.frontiersin.org/journals/neuroscience/articles/10.3389/fnins.2024.1415257/full#sec3 (MEPFeatX paper)
%}

function [MissingAnalyze, CustomOutputs, CustomAnalysisOpts]=AnalysisPluginExample_MEP(app,existFig,PluginsFolderName,AnalysisType,AnalyzeSampleRate,PreStimData,SelectedTrialsData,MissingAnalyze,Start,End)

% ======================= ZONE 1: your parameters =======================
numVar = 2;                                   % how many parameters you need
ListofVariableLabels = {'Pulse Time (s)','Plot TriaL (0=off)'};
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
% UserVar is (numVar+1) x 2; column 2 holds the values.
% Index 1 is ALWAYS the auto-added Start Time. YOUR parameters start at 2.
PulseTime = UserVar{1,2}; %in seconds
PlotTrial = round(UserVar{2,2}); %plot

% =======================================================================

% ======================== vv DO NOT EDIT vv ==============================
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
        % ======================== ^^ DO NOT EDIT ^^ ==============================

        % ==============================================================================================================================
        % ======================= ZONE 3: Analysis ======================
        %For MEP data use non-rectified data if needed
        if ~isempty(app.Processed_Conditions_DataAll{3}) && AnalysisType=="MEP" %if the third cell isn't empty then the non-rectified data was saved and if MEP is done, use the non-rectified data
            if app.AverageCheckBox.Value == 1
                TrialDataNR=app.Processed_Conditions_DataAll{3}(1,:);
            else
                TrialDataNR=app.Processed_Conditions_DataAll{3}(app.Analyze_TrialsUsed(i),:);
            end
        else
            TrialDataNR=SelectedTrialsData{i,1};
        end

        %Data to be analyzed
        AnalyzeData=TrialDataNR(Start:End);

        %Auto calculate MEP or SP default metrics
        %MEP = amplitude and area under the curve
        %SP = percent decrease and normalized area under the curve
        %If the user would like to use their own method for calculating these metrics,
        %the variables app.MEP_Amp, app.MEP_Area, app.SP_PercDecrease, app.SP_Area need to be filled in in this function
        PluginAutoCalcMEPandSP(i,app,PreStimData,AnalyzeData,AnalysisType);

        %Analyses=====
        %The number of turns (NT) is counted as the significant peaks occurring during the MEP Dur
        [Pks,Locs]=findpeaks(AnalyzeData);
        [PksN,LocsN]=findpeaks(-AnalyzeData);
        Pks=[Pks PksN];  Locs=[Locs LocsN];
        CustomOutputs.NT(i,:)=length(Pks);

        %The number of phases (NP) is counted by the zero-crossing points between the MEP onset and endpoint
        LessThanzero=find(AnalyzeData <=0);
        GreaterThanzero=find(AnalyzeData > 0);
        if length(LessThanzero) < length(GreaterThanzero)
            Diff=length(GreaterThanzero) - length(LessThanzero);
            LessThanzero(end+1:end+Diff)=0;
            PosCross=ismember(LessThanzero+1,GreaterThanzero);
            PosCross(end-Diff+1:end)=[];
            NegCross=ismember(GreaterThanzero+1,LessThanzero);
        elseif length(LessThanzero) > length(GreaterThanzero)
            Diff=length(LessThanzero) - length(GreaterThanzero);
            GreaterThanzero(end+1:end+Diff)=0;
            PosCross=ismember(LessThanzero+1,GreaterThanzero); %Indx in Analyze data before Neg cross happend = GreaterThanzero(NegCross)
            NegCross=ismember(GreaterThanzero+1,LessThanzero); %Indx in Analyze data before Pos cross happend = LessThanzero(PosCross)
            NegCross(end-Diff+1:end)=[];
        end
        CustomOutputs.NP(i,1)=sum([PosCross NegCross]);

        %Latency is the time interval between the pulse delivery time and the MEP onset
        MEPOnsetTime=app.Time(Start); %app.AllOnOffsetTime(i,1)
        CustomOutputs.Latency(i,1)=MEPOnsetTime-PulseTime;

        %Thickness is the ratio of the area under curve (AUC) to Amp
        CustomOutputs.Thickness(i,1)=app.MEP_Area(i,1)./app.MEP_Amp(i,1);

        % Stash this trial's data for the optional diagnostic plot.
        if doPlot && i == PlotTrial
            NTData=[Locs' Pks'];
            NPData=[GreaterThanzero(NegCross) LessThanzero(PosCross); GreaterThanzero(NegCross)+1 LessThanzero(PosCross)+1];
            plotData=struct('AnalyzeData',AnalyzeData,'NT',NTData,'NP',NPData');


        end

        % ==============================================================================================================================
        % ==============================================================================================================================

    else %don't analyze becuase an onset/offset wasn't found
        % ==============================================================================================================================
        % ======================= ZONE 4: Missed Analysis Fill In ======================
        % Stash this trial's data for the optional diagnostic plot
        if doPlot && i == PlotTrial
            NTData=[0 0];
            NPData=[0 0; 0 0];
            plotData=struct('AnalyzeData',[],'NT',NTData,'NP',NPData');


        end
        CustomOutputs.NT(i,:)=nan;
        CustomOutputs.NP(i,:)=nan;
        CustomOutputs.Latency(i,:)=nan;
        CustomOutputs.Thickness(i,:)=nan;
        % ==============================================================================================================================
        % ==============================================================================================================================

    end %end if Analyze, don't edit the overall structure of this if statement

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
NTData=plotData.NT;
NPData=plotData.NP;

c=colororder;
figure(); grid on;
plot(AnalyzeData,'b'); hold on;
plot(-AnalyzeData,'c--');
yline(0,'k--')
plot(NTData(:,1),NTData(:,2),'ro','LineWidth',2);
e=1;
for i3=1:length(NPData(:,1))
    plot(NPData(i3,1),AnalyzeData(NPData(i3,1)),'color',[c(e,1) c(e,2) c(e,3)],'marker','x','MarkerSize',12,'Linewidth',0.1);
    plot(NPData(i3,2),AnalyzeData(NPData(i3,2)),'color',[c(e,1) c(e,2) c(e,3)],'marker','x','MarkerSize',12,'Linewidth',0.1);
    e=e+1;
    if e > length(c(:,1))
        e=1;
    end
    if i3 == 1
        lgd=legend('Data','Negated Data','','Peaks','Zero Cross');
        lgd.AutoUpdate='off';
    end
end
end
