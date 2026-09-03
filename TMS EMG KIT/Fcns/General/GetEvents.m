%Get data for Data Type Events
%Inputs: app properties
%Outputs: numTrials and Events for Data Type Events

function [numTrials, Events] = GetEvents(app)
    %Determine which Detections Setting is being used, currently there is only one option
    app.DetectionSetting = app.DetectionSettingsDropDown.Value;
    switch app.DetectionSetting
        case 'Threshold Crossing'                 %simple threshold
            TriggerLoc = findTriggers_ThreshCross(app.numBlocks,app.Block,app.EventThresh,app.EventChannel,app.numChannels);
    
    end %end of detectionsetting switch
    
    %If no events are found, end action
    if TriggerLoc(1,:) == zeros(1,2)
        numTrials=0; Events=[];
        msgbox('No Events Found')
        return
    elseif TriggerLoc(1,:) == [0 1] %for error "Select Valid Event Channel"
        numTrials=0; Events=[];
        msgbox('Select Valid Event Channel')
        return
    end
    
    %Extract Events data
    [numTrials,Events] = Extract_EventsData(TriggerLoc,app.SearchISI,app.AllSampleRate,app.PreDur,app.PostDur,app.numChannels,app.Block);

end


%Find events trigger locations using Event parameters
%Inputs: Number of blocks, Block data, event threshold
%Outputs: Trigger locations (numTriggersx2) [Block trigger is located in Frame within block trigger is located in]        

function [TriggerLoc] = findTriggers_ThreshCross(numBlocks,Block,EventThresh,EventChannel,numChannels)
%check for valid event channel
if EventChannel==0 || EventChannel > numChannels
    TriggerLoc=[0 1;0 1];
        return
end

%Find trigger locations
w=1; TriggerLoc=zeros(10,2);
for i=1:numBlocks
    Data=Block.Data{EventChannel,i}*1000;              %apply mulitplier to put data in mV
    for e=2:length(Data)
        if Data(e) > EventThresh && Data(e-1) < EventThresh     %find points that cross the set threshold (Triggers)
            TriggerLoc(w,1)=i;                                          %Block trigger is located in
            TriggerLoc(w,2)=e;                                          %Frame within block trigger is located in
            w=w+1;
        end
    end
end


end

%Parse out events data
%Inputs: Trigger locations, search interval, sample rate, pre-duration, post-duration, number of channels, Block data
%Outputs: Number of trials, Events Data

function [numTrials,Events] = Extract_EventsData(TriggerLocOri,SearchISI,AllSampleRate,PreDur,PostDur,numChannels,Block)

%Apply search interval
TriggerLoc(1,:)=TriggerLocOri(1,:);
e=1;
for i=2:length(TriggerLocOri) %For each trigger
    if TriggerLocOri(i,1) == TriggerLoc(e,1)
        if TriggerLocOri(i,2)-TriggerLoc(e,2) >= SearchISI*AllSampleRate %compare length since last trigger
            TriggerLoc(e+1,:)=TriggerLocOri(i,:);                        %if length is larger than the interval update the triggerlist
            e=e+1;
        end
    else %New block
        TriggerLoc(e+1,:)=TriggerLocOri(i,:); %add first of each block
        e=e+1;
    end

end

numEvents=length(TriggerLoc(:,2)); %total number of events
numTrials=numEvents;

%Parse out data around events
PreDurValue=PreDur*AllSampleRate;
PostDurValue=PostDur*AllSampleRate;
Data2=zeros(numChannels,length(-PreDurValue:PostDurValue));
for i=1:length(TriggerLoc) %for each trigger location, create data array
    StartDiff=0; EndDiff=0;
    Start=TriggerLoc(i,2)-PreDurValue; End=TriggerLoc(i,2)+PostDurValue;
    if Start <= 0 %if start is less than zero, fill in with zeros
        StartDiff=-Start+1;
        Start=1;
    end
    DataLength=length(Block.Data{1,TriggerLoc(i,1)});
    if End > DataLength %if end is greater than the length of the data, fill with NaNs
        EndDiff=End-DataLength;
        End=DataLength;
    end
    Data=Block.Data(:,TriggerLoc(i,1));
    %reformat and cut data
    for i3=1:numChannels
        if isempty(Data{i3})
            Data2(i3,:)=nan(1,length(-PreDurValue:PostDurValue));
        else
            Data2(i3,:)=[nan(1,StartDiff) Data{i3}(Start:End) nan(1,EndDiff)]; %fill beginning and end with NaNs if needed
        end
    end
    Events.Data{i,1}=Data2;
end

end