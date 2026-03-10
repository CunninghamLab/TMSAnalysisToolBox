

%%
%
clear; clc;

load('SignalDataTest_CFS_OLD');
app.SignalData=appSignalData;
fname{1}='FileName';

frameinfo =[];
values = [];
StateInfo = {};
app.totEvent = 0;
app.SignalData.values = [];

values = cat(3,values,app.SignalData.CfsFile.D.data);

app.FileName.Value = fname{1};
app.totChannel = size(app.SignalData.CfsFile.D.data,2); %need to determine if data chan or other
app.totEvent = size(app.SignalData.CfsFile.Chans(1).DS,2) + app.totEvent;
app.timeUnit = app.SignalData.CfsFile.Chans(1).xUnits;
app.time = (1:size(app.SignalData.CfsFile.Chans(1).DS(1).Data,1)).* (app.SignalData.CfsFile.Chans(1).DS(1).xScale);
app.origSampleRate = 1/(app.SignalData.CfsFile.Chans(1).DS(1).xScale);
app.SignalData.frameinfo = [];
app.reloadInd = 0;                                                   % the first time find the events, save all data into one variable, so
app.ProcessInd = 0;
app.SignalData.start = app.SignalData.CfsFile.Chans(1).DS(1).xOffset;
app.PreTriggerTime.Value = abs(double(app.SignalData.start)) *1000;

StateInfoDesc = {app.SignalData.CfsFile.DSVars.Desc};
State = find(contains(StateInfoDesc,'State'));
StateL = find(contains(StateInfoDesc,'StateL'));
%StateInfo = [StateInfo,CfsFile.DSVars(State(1)).Values];
StateInfo = [StateInfo,app.SignalData.CfsFile.DSVars(State(1)).Values];

app.SignalData.values = values;
%app.SignalData.frameinfo = frameinfo; %add concat eventinfo
app.PreDurs.Value =  num2str(abs(app.SignalData.CfsFile.Chans(1).DS(1).xOffset));
app.PostDur.Value = num2str(app.time(end)-abs(app.SignalData.CfsFile.Chans(1).DS(1).xOffset));


%StateInfoDesc = {app.SignalData.CfsFile.DSVars.Desc};


%State = find(contains(StateInfoDesc,'State'));
%StateL = find(contains(StateInfoDesc,'StateL'));
%StateInfo = CfsFile.DSVars(State(1)).Values;

%if isempty(CfsFile.DSVars(StateL)) == 0
if isempty(app.SignalData.CfsFile.DSVars(StateL)) == 0
    %StateInfoL = CfsFile.DSVars(StateL(1)).Values;
    StateInfoL = app.SignalData.CfsFile.DSVars(StateL(1)).Values;
    StateInfoLCheck = cellfun(@isempty,StateInfoL);
else
    StateInfoLCheck = 1;
end



if sum(StateInfoLCheck) == size(StateInfoLCheck,2)
    app.SignalData.StateInfo = StateInfo;
    cell2mat(app.SignalData.StateInfo)
    app.SignalData.StateInfo = strsplit(num2str(ans));
else
    app.SignalData.StateInfo = StateInfo;
end

%}

%% .cfs conversion NEW current
%
clear; clc;

load('SignalDataTest_CFS');
app.Signal_Data=appSignal_Data;
ConvertV=1;

%extract from loaded data
app.AllSampleRate=round(1/app.Signal_Data.D.param.xScale(1));
app.SampleRateLabel.Text=['Sample Rate: ' num2str(app.AllSampleRate) 'Hz'];
app.PreTriggerTime=abs(round(app.Signal_Data.D.param.xOffset(1)));
PreTriggerTime_Units=app.Signal_Data.D.param.xUnits{1};        %check PreTriggerTime units

%app.Signal_Data.values=app.Signal_Data.D.data;
app.numChannels=size(app.Signal_Data.values,2);
IsWave=[app.Signal_Data.Chans.IsWaveData];
ChTitles1={app.Signal_Data.Chans.Name}';
ChTitles=ChTitles1(IsWave);
DataUnits1={app.Signal_Data.Chans.yUnits}';
DataUnits=DataUnits1(IsWave);

app.numBlocks=size(app.Signal_Data.values,3);
StateInfoLoc=find(string({app.Signal_Data.DSVars.Desc})=="State"); %change this to take whatever is in the Desc %!!!
StateA=app.Signal_Data.DSVars(StateInfoLoc).Values;
if StateA{1} == 0 %zero index
    app.numStates=StateA{end}+1;                                   %assumes largest state value to be at the end
elseif StateA{1} == 1 %one index
    app.numStates=StateA{end};
end
LabelsLoc=find({app.Signal_Data.DSVars.Desc}=="StateL");
labelsAll=app.Signal_Data.DSVars(LabelsLoc).Values;



Block1=cell(app.numBlocks,1); StateB=cell(app.numStates,1); e=1;
for i=1:app.numBlocks
    %Block1{i,1}=app.Signal_Data.values(:,:,i); %!!!!!! 
    DataValues=app.Signal_Data.values(:,:,i); %numPointsxnumChannels
    for i2=1:app.numChannels
        Block1{i2,e}=DataValues(:,i2)'*ConvertV;
    end
    label=labelsAll{i};

    %check for zero vs one indexed %!!!
    if StateA{1}==0
        State=StateA{i}+1;
        if (i~=1) && (State ~= StateA{i-1}+1)  %reset for when state changes, assumes states occur in chronological order
            e=1;
        end
    else
        State=StateA{i};
        if (i~=1) && (State ~= StateA{i-1})    %reset for when state changes, assumes states occur in chronological order
            e=1;
        end
    end

    StateB{State,1}{e,1}=Block1; %organize by state, StateB contains cells with block data for each state, convert to volts
    StateB{State,2}=label;
    e=e+1;
end
%There should be 10 channels x 40 blocks 
%StateB=1x2cell, {1}=10x40cell {1}{1}=1x30000 data
%                {2}=[], there are no labels for the states

%}
%% .mat import
%

clear; clc;

load('SignalDataTest');
app.Signal_Data=appSignal_Data;
app.Block=appBlock;
ConvertV=1;

app.numChannels=double(app.Signal_Data.chans);
StateA=[app.Signal_Data.frameinfo.state]';
labelsAll={app.Signal_Data.frameinfo.label}';
app.numBlocks=app.Signal_Data.frames; 
app.numStates=StateA(end); %assumes largest state value to be at the end

Block1=cell(app.numChannels,1); StateB=cell(app.numStates,1); e=1;
for i=1:app.numBlocks %separate Blocks into their "States"
    DataValues=app.Signal_Data.values(:,:,i); %numPointsxnumChannels
    for i2=1:app.numChannels
        Block1{i2,e}=DataValues(:,i2)'*ConvertV; 
    end
    label=labelsAll{i};
    State=StateA(i);

    StateB{State,1}=Block1;
    StateB{State,2}=label;

    %reset for when state changes, assumes states occur in chronological order
    if (i~=app.numBlocks) && (State ~= StateA(i+1)) 
        e=1;
        Block1=cell(app.numChannels,1);
    else
        e=e+1;
    end
   

end
%}

%%
% StateS=6;
% Channel=2;
% SelectedTrial=10;
%
% figure()
% plot(StateB{StateS,1}{Channel,SelectedTrial});









