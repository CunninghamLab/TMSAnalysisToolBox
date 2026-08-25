%{
AnalysisPluginExample - simple example script to show users how to make analysis plugins

%8/25/26: Need to update this How-To section
%Note: this function is run in a loop within the main app for each trial, so this function is analyzing one trial's data, this is different than the onset/offset analysis plugin

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
  app.CustomOutputs      - add the custom outputs to this scruct (e.g. app.CustomOutputs.Latency), double column vector, each row is the result from a trial
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

function [UserVar]=AnalysisPluginExample2(whichTrial,app,existFig,PluginsFolderName,UserVar,NumTrials,AnalyzeSampleRate,PreStimData,TrialData,Start,End)

if whichTrial == 0 %Run this at the beginning, do not edit this line
    % ======================= ZONE 1: your parameters =====================
    numVar = 2;                                   % how many parameters you need
    ListofVariableLabels = {'Pulse Time (s)','Plot TriaL (0=off)'};
    DefaultValues        = [0, 0];             % first-run defaults, same order/units as labels ([] for none)

    % =====================================================================
    assert(numel(ListofVariableLabels)==numVar, 'numVar must equal the number of labels.');

    % --------------------------- DO NOT EDIT -----------------------------
    % Builds the parameter pop-up and loads/saves this plugin's settings.
    % mfilename tells the helper which settings file belongs to this plugin.
    % DefaultValues pre-fills the pop-up the first time (before settings exist).
    [app.CustomAnalysisOpts, UserVar] = createAnalysisPluginFigure(existFig, ...
        app.CustomAnalysisOpts, PluginsFolderName, numVar, ListofVariableLabels, mfilename, DefaultValues);
    % ---------------------------------------------------------------------

    % ======================= ZONE 2: Initialize custom outputs ===========
    %Initialize custom outputs
    app.CustomOutputs.MEPAmp=nan(NumTrials,1);
    app.CustomOutputs.MEPArea=nan(NumTrials,1);
    app.CustomOutputs.NT=nan(NumTrials,1);
    app.CustomOutputs.NP=nan(NumTrials,1);
    app.CustomOutputs.Latency=nan(NumTrials,1);
    app.CustomOutputs.Thickness=nan(NumTrials,1);
    % =====================================================================
    return;
end %end if whichTrial==0, do not edit this line

% ======================= ZONE 3: unpack parameters =======================
% UserVar is (numVar+1) x 2; column 2 holds the values.
% Index 1 is ALWAYS the auto-added Start Time. YOUR parameters start at 2.
PulseTime = UserVar{1,2}; %in seconds
PlotTrial = round(UserVar{2,2}); %plot

% =========================================================================

% ======================= ZONE 4: Analysis ================================

%---------------------------- Optional ------------------------------------
%Auto calculate MEP or SP default metrics
%MEP = amplitude and area under the curve
%SP = percent decrease and normalized area under the curve
%[MEPAmpValue, MEPAreaValue,SPPercDecreaseValue,SPAreaValue]=PluginAutoCalcMEPandSP(whichTrial,app,PreStimData,TrialData,"MEP"); %type in "MEP" or "SP"
%Note that these results are in volts
[MEPAmpValue, MEPAreaValue,~,~]=PluginAutoCalcMEPandSP(whichTrial,app,PreStimData,TrialData,"MEP"); %type in "MEP" or "SP"
app.CustomOutputs.MEPAmp(whichTrial,1)=MEPAmpValue*1000; %convert to mV
 app.CustomOutputs.MEPArea(whichTrial,1)=MEPAreaValue*1000; %convert to mV 
%--------------------------------------------------------------------------

%The number of turns (NT) is counted as the significant peaks occurring during the MEP Dur
[Pks,Locs]=findpeaks(TrialData);
[PksN,LocsN]=findpeaks(-TrialData);
Pks=[Pks PksN];  Locs=[Locs LocsN];
app.CustomOutputs.NT(whichTrial,:)=length(Pks);

%The number of phases (NP) is counted by the zero-crossing points between the MEP onset and endpoint
LessThanzero=find(TrialData <=0);
GreaterThanzero=find(TrialData > 0);
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
app.CustomOutputs.NP(whichTrial,1)=sum([PosCross NegCross]);

%Latency is the time interval between the pulse delivery time and the MEP onset
MEPOnsetTime=app.Time(Start); 
app.CustomOutputs.Latency(whichTrial,1)=MEPOnsetTime-PulseTime;

%Thickness is the ratio of the area under curve (AUC) to Amp
app.CustomOutputs.Thickness(whichTrial,1)=app.CustomOutputs.MEPArea(whichTrial,1)./app.CustomOutputs.MEPAmp(whichTrial,1);

%if this is the selected plot trial, plot the trial
if whichTrial == PlotTrial
    NTData=[Locs' Pks'];
    NPData=[GreaterThanzero(NegCross) LessThanzero(PosCross); GreaterThanzero(NegCross)+1 LessThanzero(PosCross)+1];
    plotData=struct('TrialData',TrialData,'NT',NTData,'NP',NPData');
    plotTrialFunction(plotData); %run ploting function

end

end %end function

% ===================== Plotting function =================================
%function for plotting trial data
function plotTrialFunction(plotData)
TrialData=plotData.TrialData;
NTData=plotData.NT;
NPData=plotData.NP;

c=colororder;
figure(); grid on;
plot(TrialData,'b'); hold on;
plot(-TrialData,'c--');
yline(0,'k--')
plot(NTData(:,1),NTData(:,2),'ro','LineWidth',2);
e=1;
for i3=1:length(NPData(:,1))
    plot(NPData(i3,1),TrialData(NPData(i3,1)),'color',[c(e,1) c(e,2) c(e,3)],'marker','x','MarkerSize',12,'Linewidth',0.1);
    plot(NPData(i3,2),TrialData(NPData(i3,2)),'color',[c(e,1) c(e,2) c(e,3)],'marker','x','MarkerSize',12,'Linewidth',0.1);
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
