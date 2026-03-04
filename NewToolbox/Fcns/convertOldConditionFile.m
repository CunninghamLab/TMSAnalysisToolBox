%Convert old Condition files
%Inputs: structure from old files
%Outputs: structure with data in new format

function New=convertOldConditionFile(Old)

%% ConditionsMatrix
%convert pages (Trial) to commam delimited
New=struct();
for i=1:length(Old.ConditionInfo(:,1))
    Trials=Old.ConditionInfo{i,1};
    TrialsSplit=strsplit(Trials);
    NewTrials=convert2commadel(TrialsSplit,',');
    New.ConditionsMatrix{i,1}=NewTrials;
    New.ConditionsMatrix{i,2}=Old.ConditionInfo{i,2};
    New.ConditionsMatrix{i,3}=Old.ConditionInfo{i,3};
end

%% ConditionsData

%Column 1: Data
for i=1:length(Old.ConditionData(:,1)) %for each condition
    numTrialsCurrCondition=length(Old.ConditionData{i,1}{1,1}(:,1));
    NewConditionData=cell(numTrialsCurrCondition,1);
    Data=Old.ConditionData{i, 1}{1,1};
    for e=1:numTrialsCurrCondition %for each trial in current condition
        NewConditionData{e,1}=Data(e,:)*0.001; %Old condition files saved data in mV, convert to V;
        NewConditionData{e,2}=Old.time;
    end
    New.ConditionsData{i,1}=NewConditionData;
end

%Column 2: Event Type and Column 3: Event Name or Events parameters
%Extra data only saved if there is an Event Name (comments)
if ~isempty(Old.EventsToConditions)
    for i=1:size(Old.EventsToConditions,1)
        New.ConditionsData{i,2}=2;
        New.ConditionsData{i,3}=string(Old.EventsToConditions(i,:));
    end
else
    New.ConditionsData(:,2)=cell(size(Old.ConditionInfo,1),1);
    New.ConditionsData(:,3)=cell(size(Old.ConditionInfo,1),1);
end

%Check for empty data cells 
New.ConditionsData(cellfun(@isempty,New.ConditionsData(:,1)),:)=[];

%% Channel Titles
%channel titles are not saved in old condition files
New.ChTitles=[];

%% Sample Rate
New.SampleRate=Old.sampleRate;

%% File Name
New.FileName{1,1}=Old.FileName.Value;

%% All Data and numTrialsperCondition

%Make all data same length
numConditions=length(New.ConditionsData(:,1));
TimeRanges=zeros(numConditions,2);
for i=1:numConditions
    TimeRanges(i,1)=New.ConditionsData{i}{1,2}(1);
    TimeRanges(i,2)=New.ConditionsData{i}{1,2}(end);
end

%Determine max time length, and which condition has this max length
minTimeRange=min(TimeRanges(:,1));
maxTimeRange=max(TimeRanges(:,2));

%create time base on the smallest and largest time
newTime=minTimeRange:1/Old.sampleRate:maxTimeRange;

newConditions_Data=cell(numConditions,length(New.ConditionsData(1,:))+1);
newConditions_Data(:,3:end)=New.ConditionsData(:,2:end);
newConditions_Data{1,2}=newTime;
for i=1:numConditions
    Data=New.ConditionsData{i}(:,1);
    DataTime=New.ConditionsData{i}{1,2};
    %check if time already matches new Time
    if (length(DataTime) ~= length(newTime)) || all(DataTime ~= newTime) %if time doesn't match newTime
        %find start index
        [~,StartIndx]=min(abs(newTime - DataTime(1)));
        [~,EndIndx]=min(abs(newTime - DataTime(end)));
        DataMat=cell2mat(Data);
        newData=[nan(length(DataMat(:,1)),StartIndx-1) DataMat nan(length(DataMat(:,1)),length(newTime)-EndIndx)];
        newConditions_Data{i,1}=newData;

    else
        newConditions_Data{i,1}=cell2mat(Data);


    end

end

New.AllData{1,1}=cell2mat(newConditions_Data(:,1));
New.AllData{1,2}=newTime;



%Create list of all Conditions Trials
[numTrialsperCondition1, ~]=cellfun(@size, New.ConditionsData(:,1), 'UniformOutput', false);
New.numTrialsperCondition=zeros(length(numTrialsperCondition1),3);
New.numTrialsperCondition(:,1)=cell2mat(numTrialsperCondition1);
%New.AllData=cell(sum(New.numTrialsperCondition(:,1)),3);
CStart=1;
for i=1:length(New.numTrialsperCondition(:,1)) %for each trial in each condition
    CEnd=CStart+New.numTrialsperCondition(i,1)-1;
    New.numTrialsperCondition(i,2)=CStart;
    New.numTrialsperCondition(i,3)=CEnd;
    %New.AllData(CStart:CEnd,1:2)=New.ConditionsData{i,1};
    CStart=CEnd+1;
end
%New.AllData(:,3)=mat2cell(num2str([1:sum(New.numTrialsperCondition)]'),ones(1,sum(New.numTrialsperCondition(:,1))));



end




