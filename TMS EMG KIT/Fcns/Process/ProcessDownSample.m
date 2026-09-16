%Down samples data and time by a user selected factor
%Inputs: data time, data, factor, sample rate
%Outputs: new time, data, and sample rate

function [newTime,newData,NewSampleRate]=ProcessDownSample(app,DataTime,DataAll,Factor,SampleRate)

NewSampleRate=SampleRate/Factor;
if rem(NewSampleRate,Factor) ~= 0
    NewSampleRate=round(SampleRate/Factor);
    Factor=SampleRate/NewSampleRate;
    if string(class(app)) ~= "double"
        app.FactorEditField.Value=Factor;
    end

end

newData=[];

if rem(Factor,1) == 0 %integer
    %use matlab function downsample
    newTime=downsample(DataTime,Factor);
    newData=downsample(DataAll',Factor)';

else %factor is not an integer

    for i=1:length(DataAll(:,1))
        Data=DataAll(i,:);
        %downsample using resample function
        %resample assumes the values before and after the data given are zeros,
        %added padding to help with the affects of this on the filter
        RPadding=(SampleRate/Factor)*4;
        DataP=[repmat(Data(1),1,RPadding) Data repmat(Data(end),1,RPadding)];
        DataB=resample(DataP,SampleRate/Factor,SampleRate);
        DataB(1:floor((RPadding/Factor))-1)=[]; DataB(end-ceil((RPadding/Factor)):end)=[];
        newData(i,:)=DataB;
    end

    %downsample time
    % TimeP=[repmat(DataTime(1),1,RPadding) DataTime repmat(DataTime(end),1,RPadding)];
    % TimeB=resample(TimeP,SampleRate/Factor,SampleRate);
    % if rem((RPadding/Factor),1) ~= 0
    %     disp('RPadding is not an integer');
    %     TimeB(1:floor((RPadding/Factor))-1)=[]; TimeB(end-ceil((RPadding/Factor)):end)=[];
    % else
    %     TimeB(1:(RPadding/Factor)-1)=[]; TimeB(end-(RPadding/Factor):end)=[];
    % 
    % end
    % newTime=TimeB;

    newTimeStep=1/NewSampleRate;
    TimeStart=DataTime(1);
    newTime=TimeStart:newTimeStep:((length(newData(1,:))-1)*newTimeStep)+TimeStart; %use starting time then fill to the length of the newData

end


end

