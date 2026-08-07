%Down samples data and time by a user selected factor
%Inputs: data time, data, factor, sample rate
%Outputs: new time, data, and sample rate

function [newTime,newData,newSampleRate]=ProcessDownSample(DataTime,Data,Factor,SampleRate)

newData=Data(:,1:Factor:end);
newTime=DataTime(1:Factor:end);
newSampleRate=SampleRate/Factor;



end