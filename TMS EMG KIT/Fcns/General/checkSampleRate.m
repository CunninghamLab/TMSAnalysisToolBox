%checkSampleRate
%checks if the sample rates of condition files loaded in are different 
%Downsamples everything to the smallest sample rate%checkSampleRate
%Inputs: Conditions_DataAll before from the loaded conditions before they are combined, all the file sample rates
%Outputs: the new sample rate (minimum sample rate of the files), new Conditions_DataAll values

function [newAllSampleRate,Conditions_DataAll]=checkSampleRate(Conditions_DataAll, AllSampleRates)

%lowest sample rate
MinSampleRate=min(AllSampleRates);

DataWithDiffSR=find(AllSampleRates ~= MinSampleRate); %data with different sample rates

for i=DataWithDiffSR
    %update in Conditions_DataAll
    %Downsample to min sample rate
    CurrSampleRate=AllSampleRates(i);
    Factor=CurrSampleRate/MinSampleRate;

    DataTime=Conditions_DataAll{i+1,2};
    DataAll=Conditions_DataAll{i+1,1};
    [newTime,newData,newSampleRate]=ProcessDownSample(1,DataTime,DataAll,Factor,CurrSampleRate);
    if newSampleRate ~= MinSampleRate
        error('Downsampled to the wrong sample rate?');
    end

    Conditions_DataAll{i+1,1}=newData;
    Conditions_DataAll{i+1,2}=newTime;

end
newAllSampleRate=MinSampleRate;


end


