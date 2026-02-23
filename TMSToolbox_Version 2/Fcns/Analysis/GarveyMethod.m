%Uses either the Garvey et all method or the SD limit method
%specify which method and if there are 50% algorithm checks
%Inputs:app data, PreStimData, SelectedTrialsData, PreStim stdv, MCD or SD ('MCD' or 'SD'), varargin(OnsetFifty, OffsetFifty (boolean))
%Outputs: All Onset and Offset time for each trial

function [AllOnOffsetTime, OnsetLimit,OffsetLimit,meanPreStimData]=GarveyMethod(app,AnalyzeSampleRate,PreStimData,SelectedTrialsData,AnalyzeMethod,MCDorSD,OnsetFifty,OffsetFifty)


switch MCDorSD
    case 'MCD'
        % find limits using MCD
        PreStimData=abs(PreStimData); %rectify
        Difference=abs(PreStimData(1:end-1,:)-PreStimData(2:end,:));
        MCD=mean(Difference);
        meanPreStimData=mean(PreStimData);
        MCDConstant=app.MCDConstantEditField.Value;

        switch app.AnalysisOptionDropDown.Value
            case app.AnalysisOptionDropDown.ItemsData{1} %MEP
                %LowEndLimit=meanPreStimData - MCD*MCDConstant;
                HighEndLimit=meanPreStimData + MCD*MCDConstant;

            case app.AnalysisOptionDropDown.ItemsData{2} %SP
                HighEndLimit=meanPreStimData - MCD*MCDConstant;
                
            case app.AnalysisOptionDropDown.ItemsData{3} %MEP

        end
        

        OnsetLimit=HighEndLimit;
        OffsetLimit=HighEndLimit;

    case 'SD'
        %Determine Limits
        meanPreStimData=mean(abs(PreStimData)); %mean of rectified data
        PreStimSTDV=std(abs(PreStimData));

        switch app.AnalysisOptionDropDown.Value
            case app.AnalysisOptionDropDown.ItemsData{1} %MEP
                OnsetLimit=meanPreStimData + app.OnsetSDEditField.Value*PreStimSTDV;
                OffsetLimit=meanPreStimData + app.OffsetSDEditField.Value*PreStimSTDV;

            case app.AnalysisOptionDropDown.ItemsData{2} %SP
                OnsetLimit=meanPreStimData - app.OnsetSDEditField.Value*PreStimSTDV;
                OffsetLimit=meanPreStimData - app.OffsetSDEditField.Value*PreStimSTDV;
            case app.AnalysisOptionDropDown.ItemsData{3} %EMG

        end


end


% find onset and offset
AllOnOffsetTime=zeros(length(SelectedTrialsData(:,1)),2);
NoOnset=0; %NoOffset=0;
TrialTime=app.Time;
for i=1:length(SelectedTrialsData(:,1)) %for each trial
    TrialDataNR=SelectedTrialsData{i,1};
    %TrialTime=SelectedTrialsData{i,2};
    TrialData=abs(TrialDataNR); %rectify

    %find onset time
    Tol=eps("double");
    %Start=find(abs(appTime - OnsetTime) < Tol);
    Start=app.StartTimemsEditField.Value*0.001;
    StartSamples=Start*AnalyzeSampleRate;
    StartSamples=round(StartSamples);           %round if start time isn't within the sample rate
    if StartSamples ~= Start*AnalyzeSampleRate  %if start time was changed update edit field
        app.StartTimemsEditField.Value=StartSamples*(1/AnalyzeSampleRate)*1000; %update start time to new value in milliseconds
    end


    StartIndx=find(abs(TrialTime-app.StartTimemsEditField.Value*0.001) < Tol);
    TrialTimeO=TrialTime(StartIndx:end);
    TrialDataO=TrialData(StartIndx:end);
    MinDur=app.OnsetMinDurmsEditField.Value*0.001; %convert ms to seconds
    MinSamplesDur=MinDur*AnalyzeSampleRate; %how many data point in a row needs to be above threshold (LowEndLimit)
    MinSamplesDur=round(MinSamplesDur);  %round if start time isn't within the sample rate
    if MinSamplesDur < 1  %if start time was changed update edit field
        MinSamplesDur=1;
        app.OnsetMinDurmsEditField.Value=MinSamplesDur*(1/AnalyzeSampleRate)*1000; %update start time to new value in milliseconds
    elseif MinSamplesDur ~= MinDur*AnalyzeSampleRate
        app.OnsetMinDurmsEditField.Value=MinSamplesDur*(1/AnalyzeSampleRate)*1000; %update start time to new value in milliseconds

    end

    %Find Onset
    for i2=2:length(TrialDataO) %for each point in trial data
        DataPt_Prev=TrialDataO(i2-1);
        DataPt=TrialDataO(i2);

        %if thresold is crossed
        if AnalyzeMethod == "MEP" && DataPt > OnsetLimit(i) && DataPt_Prev < OnsetLimit(i)
            DurationEnd=i2+MinSamplesDur-1;
            if DurationEnd > length(TrialDataO) %Onset was found too close to the end of the trial data
                OnsetTime=nan;
                OffsetTime=nan;
                NoOnset=1;
                break;
            end
            DurationData=TrialDataO(i2:DurationEnd);
            %If min duration is met
            %No 50% rule-100% of the duration data must be above the limit, 50% rule-at least 50% of the duration data needs to be above the limit
            if (OnsetFifty == 0 && sum(DurationData > OnsetLimit(i)) == length(DurationData)) || (OnsetFifty == 1 && sum(DurationData > OnsetLimit(i))/length(DurationData)*100 >= 50)
                OnsetTime=TrialTimeO(i2);
                NoOnset=0;
                break;
            end

        elseif AnalyzeMethod == "SP" && DataPt < OnsetLimit(i) && DataPt_Prev > OnsetLimit(i)
            DurationData=TrialDataO(i2:i2+MinSamplesDur-1);
            %If min duration is met
            %No 50% rule-100% of the duration data must be above the limit, 50% rule-at least 50% of the duration data needs to be above the limit
            if (OnsetFifty == 0 && sum(DurationData < OnsetLimit(i)) == length(DurationData)) || (OnsetFifty == 1 && sum(DurationData < OnsetLimit(i))/length(DurationData)*100 >= 50)
                OnsetTime=TrialTimeO(i2);
                NoOnset=0;
                break;
            end

        end

        %if the end of the data is reached and no onset if found
        if i2 == length(TrialDataO)
            OnsetTime=nan;
            OffsetTime=nan;
            NoOnset=1;
        end
    end %end for each point, end of searching for onset


    if NoOnset==0
        TrialDataO=TrialData(TrialTime >= OnsetTime);
        TrialTimeO=TrialTime(TrialTime >= OnsetTime);
        MinDur=app.OffsetMinDurmsEditField.Value*0.001; %convert ms to seconds
        MinSamplesDur=MinDur*AnalyzeSampleRate; %how many data point in a row needs to be above threshold (LowEndLimit)
        MinSamplesDur=round(MinSamplesDur);  %round if start time isn't within the sample rate
        if MinSamplesDur < 1 %if start time was changed update edit field
            MinSamplesDur=1;
            app.OffsetMinDurmsEditField.Value=MinSamplesDur*(1/AnalyzeSampleRate)*1000; %update start time to new value in milliseconds
        elseif MinSamplesDur ~= MinDur*AnalyzeSampleRate
            app.OffsetMinDurmsEditField.Value=MinSamplesDur*(1/AnalyzeSampleRate)*1000; %update start time to new value in milliseconds
        end

        %Find Offset
        for i2=2:length(TrialDataO)
            DataPt_Prev=TrialDataO(i2-1);
            DataPt=TrialDataO(i2);
            %if thresold is crossed
            if AnalyzeMethod == "MEP" && DataPt < OffsetLimit(i) && DataPt_Prev > OffsetLimit(i)
                DurationEnd=i2+MinSamplesDur-1;
                if DurationEnd > length(TrialDataO) %Onset was found too close to the end of the trial data
                    OffsetTime=nan;
                    break;
                end
                DurationData=TrialDataO(i2:i2+MinSamplesDur-1);
                %If min duration is met
                %No 50% rule-100% of the duration data must be above the limit, 50% rule-at least 50% of the duration data needs to be above the limit
                if (OffsetFifty == 0 && sum(DurationData < OffsetLimit(i)) == length(DurationData)) || (OffsetFifty == 1 && sum(DurationData < OffsetLimit(i))/length(DurationData)*100 >= 50)
                    OffsetTime=TrialTimeO(i2);
                    break;
                end

            elseif AnalyzeMethod == "SP" && DataPt > OffsetLimit(i) && DataPt_Prev < OffsetLimit(i)
                DurationData=TrialDataO(i2:i2+MinSamplesDur-1);
                %If min duration is met
                %No 50% rule-100% of the duration data must be above the limit, 50% rule-at least 50% of the duration data needs to be above the limit
                if (OffsetFifty == 0 && sum(DurationData > OffsetLimit(i)) == length(DurationData)) || (OffsetFifty == 1 && sum(DurationData > OffsetLimit(i))/length(DurationData)*100 >= 50)
                    OffsetTime=TrialTimeO(i2);
                    break;
                end

            end
            
            if i2 == length(TrialDataO)
                OffsetTime=nan;
                break;
            end
        end %end for each frame

    end %end if an onset was found

    % figure() %for troublshooting
    % hold on;
    % plot(app.Time,abs(TrialDataNR),'r')
    % % plot(TrialTimeO,TrialDataO,'r')
    % % plot(TrialTimeO(i2:i2+MinSamplesDur-1),DurationData,'b--')
    % % plot(PreStimData,'b--')
    % yline(OnsetLimit(i),'c')
    % yline(OffsetLimit(i),'m--');
    % xline(OnsetTime,'k','LineWidth',2)
    % xline(OffsetTime,'g','LineWidth',2)
    % legend('Data','Onset Limit','Offset Limit','Onset Time','Offset Time')

    AllOnOffsetTime(i,1)=OnsetTime;
    AllOnOffsetTime(i,2)=OffsetTime;



end %end for each trial



end %end function

