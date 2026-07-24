%{
AnalysisPluginExample_Custom - simple example script to show users how to make plugins for a custom analysis

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



function [MissingAnalyze, CustomOutputs, CustomAnalysisOpts]=AnalysisPluginExample_Custom(app,existFig,PluginsFolderName,AnalysisType,AnalyzeSampleRate,PreStimData,SelectedTrials,MissingAnalyze,Start,End)

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

%Auto calculate MEP or SP default metrics
%MEP = amplitude and area under the curve
%SP = percent decrease and normalized area under the curve
if useDefaultMEPorSPAnalysis == 1
    PluginAutoCalcMEPandSP(app,PreStimData,SelectedTrials,AnalysisType);
end


% ======================= ZONE 2: unpack parameters =====================
% UserVar is (numVar+1) x 2; column 2 holds the values.
% Index 1 is ALWAYS the auto-added Start Time. YOUR parameters start at 2.
PulseTime = UserVar{1,2}; %in seconds
PlotTrial = round(UserVar{2,2}); %plot

% =======================================================================

% ======================= ZONE 3: Analysis ======================
Analyze=1; %Boolean, whether or not to analyze the trial, if override is used, this is always 1, if override isn't used this is handled in the if statement
for i=1:length(SelectedTrials) %for each trial
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
        %For MEP data use non-rectified data if needed
        if ~isempty(app.Processed_Conditions_DataAll{3}) && AnalyzeMethod=="MEP" %if the third cell isn't empty then the non-rectified data was saved and if MEP is done, use the non-rectified data
            if app.AverageCheckBox.Value == 1
                TrialDataNR=app.Processed_Conditions_DataAll{3}(1,:);
            else
                TrialDataNR=app.Processed_Conditions_DataAll{3}(app.Analyze_TrialsUsed(i),:);
            end
        else
            TrialDataNR=SelectedTrials{i,1};
        end

        %Data to be analyzed
        AnalyzeData=TrialDataNR(Start:End);

        %% Analyses 
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

        if i==PlotTrial
            plotData=AnalyzeData;
            NTData=[Locs' Pks'];
            NPData=[GreaterThanzero(NegCross) LessThanzero(PosCross); GreaterThanzero(NegCross)+1 LessThanzero(PosCross)+1];

        end

    else
        CustomOutputs.NT(i,:)=nan;
        CustomOutputs.NP(i,:)=nan;
        CustomOutputs.Latency(i,:)=nan;
        CustomOutputs.Thickness(i,:)=nan;
    end
end

%plot data to visualize the peaks
if PlotTrial >= 1
    c=colororder;
    figure(); grid on;
    plot(plotData,'b'); hold on;
    plot(-plotData,'c--');
    yline(0,'k--')
    plot(NTData(:,1),NTData(:,2),'ro','LineWidth',2);
    e=1;
    for i3=1:length(NPData(1,:))
        plot(NPData(1,i3),plotData(NPData(1,i3)),'color',[c(e,1) c(e,2) c(e,3)],'marker','x','MarkerSize',12,'Linewidth',0.1);
        plot(NPData(2,i3),plotData(NPData(2,i3)),'color',[c(e,1) c(e,2) c(e,3)],'marker','x','MarkerSize',12,'Linewidth',0.1);
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

% =======================================================================

end %end function
