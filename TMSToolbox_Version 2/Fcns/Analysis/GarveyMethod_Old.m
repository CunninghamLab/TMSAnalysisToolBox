%Uses either the Garvey et all method or the SD limit method
%specify which method and if there are 50% algorithm checks
%Inputs:app data, PreStimData, SelectedTrialsData, PreStim stdv, MCD or SD ('MCD' or 'SD'), varargin(OnsetFifty, OffsetFifty (boolean))
%Outputs: All Onset and Offset time for each trial


function [AllOnOffsetTime, OnsetLimit,OffsetLimit,meanPreStimData]=GarveyMethod_Old(app,AnalyzeSampleRate,PreStimData,SelectedTrialsData,AnalyzeMethod,MCDorSD,OnsetFifty,OffsetFifty)


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
    %StartIndx=find(single(TrialTime) == single(app.StartTimemsEditField.Value*0.001)); %convert ms to seconds, convert to single for comparison
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

    Counter=0; e=1; %For troubleshooting: 
    CounterAll=[]; BAll=zeros(1,MinSamplesDur); 
    i2=2;
    while i2 <= length(TrialDataO) %for each point in trial data
        DataPt_Prev=TrialDataO(i2-1);
        DataPt=TrialDataO(i2);
        if (AnalyzeMethod == "MEP" && DataPt > OnsetLimit(i)) || (AnalyzeMethod == "SP" && DataPt < OnsetLimit(i)) %if the data point is greater than the High end limit
            Counter=Counter+1;      %add to counter
            if OnsetFifty==1
                B(e)=1;
            end
            if Counter == 1         %if this is the first count,
                if (AnalyzeMethod == "MEP" && DataPt_Prev < OnsetLimit(i)) || (AnalyzeMethod == "SP" && DataPt_Prev > OnsetLimit(i))  %check if the previous point is below the High end limit
                    OnsetTimeTemp=TrialTimeO(i2); %if so, save the time
                    OnsetFrameTemp=i2;
                else
                    Counter=0;                    %if not, reset counter
                    if OnsetFifty==1
                        B=[];
                    end
                end
            end

        elseif  OnsetFifty==1 && Counter > 0
            B(e)=0;        %add to the counter and vector
            Counter=Counter+1;
        else                                      %if the data point is not greater than the High end limit, reset the counter
            Counter=0;
            if OnsetFifty==1
                B=[]; e=1;
            end
        end

        %as long as the counter is not reset, make sure to increment for 50% check vector
        if OnsetFifty == 1 && Counter > 0
            e=e+1;
        end

        if Counter >= MinSamplesDur    %if the counter makes it to the minimum duration set by the user
            if OnsetFifty == 0
                OnsetTime=OnsetTimeTemp;   %define the onset time
                break;                     %stop search
            elseif OnsetFifty == 1
                if sum(B)/length(B)*100 >= 50      %check if 50% of the data was below the Low end limit
                    OnsetTime=OnsetTimeTemp;    %if so, set offset time
                    break;                        %end loop
                else                              %if 50% of the data was not below the Low end limit
                    B=[]; e=1;                    %reset the counter and 50% check vector
                    Counter=0;
                    i2=OnsetFrameTemp;
                end
            end
        end
        if i2 == length(TrialDataO)    %if the end of the data is reached and no onset if found
            OnsetTime=nan;
            OffsetTime=nan;
            NoOnset=1;
        end

        %For troubleshooting============================ !!!
        CounterAll(i2,1)=Counter;
        % if OnsetFifty ==1
        %     BAll(i2,1)=B;
        % end
        %====================================
        i2=i2+1;
    end %end for each point, end of searching for onset
    if NoOnset ==1
        msgbox('No Onset found. Please try again.');
    elseif NoOnset==0
        %find offset time
        MinDur=app.OffsetMinDurmsEditField.Value*0.001; %convert ms to seconds
        MinSamplesDur=MinDur*AnalyzeSampleRate; %how many data point in a row needs to be above threshold (LowEndLimit)
        MinSamplesDur=round(MinSamplesDur);  %round if start time isn't within the sample rate
        if MinSamplesDur < 1 %if start time was changed update edit field
            MinSamplesDur=1;
            app.OffsetMinDurmsEditField.Value=MinSamplesDur*(1/AnalyzeSampleRate)*1000; %update start time to new value in milliseconds
        elseif MinSamplesDur ~= MinDur*AnalyzeSampleRate
            app.OffsetMinDurmsEditField.Value=MinSamplesDur*(1/AnalyzeSampleRate)*1000; %update start time to new value in milliseconds
        end

        TrialDataO=TrialData(TrialTime >= OnsetTime);
        TrialTimeO=TrialTime(TrialTime >= OnsetTime);

        Counter=0; e=1; B=[];
        i2=2;
        while i2 <= length(TrialDataO)
            DataPt_Prev=TrialDataO(i2-1);
            DataPt=TrialDataO(i2);
            %check if the data has fallen below the Low end limit
            if (AnalyzeMethod == "MEP" && DataPt < OffsetLimit(i)) || (AnalyzeMethod == "SP" && DataPt > OffsetLimit(i)) %if the data has fallen below the Low end limit
                Counter=Counter+1;     %add to counter
                if OffsetFifty == 1
                    B(e)=1;                %add to 50% check vector
                end

                if Counter == 1 %if this is the first encounter,
                    if (AnalyzeMethod == "MEP" && DataPt_Prev > OffsetLimit(i)) || (AnalyzeMethod == "SP" && DataPt_Prev < OffsetLimit(i))  %check that the previous value was above the Low end limit
                        OffsetTimeTemp=TrialTimeO(i2); %if so, save the offset time
                        OnsetFrameTemp=i2;
                    else                               %if not, reset counter and vector
                        Counter=0;
                        if OffsetFifty == 1
                            B=[];
                        end
                    end
                end
            elseif OffsetFifty == 1 && Counter > 0 %if the data if not below the limit and the counter is above zero
                B(e)=0;        %add to the counter and vector
                Counter=Counter+1;
            else               %if the data is not below the limit and the counter is at zero
                Counter=0;     %reset the counter and vector
                if OffsetFifty == 1
                    B=[]; e=1;
                end
            end

            %as long as the counter is not reset, make sure to increment for 50% check vector
            if OffsetFifty == 1 && Counter > 0
                e=e+1;
            end

            %check if the counter has reached the minimum duration
            if Counter >= MinSamplesDur           %if it has
                if OffsetFifty == 0
                    OffsetTime=OffsetTimeTemp;   %define the onset time
                    break;                     %stop search
                elseif OffsetFifty == 1
                    if sum(B)/length(B)*100 >= 50      %check if 50% of the data was below the Low end limit
                        OffsetTime=OffsetTimeTemp;    %if so, set offset time
                        break;                        %end loop
                    else                              %if 50% of the data was not below the Low end limit
                        B=[]; e=1;                    %reset the counter and 50% check vector
                        Counter=0;
                        i2=OnsetFrameTemp;
                    end
                end
            end
            %send warning if there is no offset found
            if i2 == length(TrialDataO)
                msgbox('No Offset found. Please try again.');
                OffsetTime=nan;
                break;
            end
            i2=i2+1;
        end %end for each frame
    end %end if an onset was found

    % figure() %for troublshooting
    % plot(app.Time,abs(TrialDataNR),'r')
    % hold on
    % %plot(PreStimData,'b--')
    % yline(OnsetLimit(i),'c')
    % yline(OffsetLimit(i),'m--');
    % xline(OnsetTime,'k','LineWidth',2)
    % xline(OffsetTime,'g','LineWidth',2)
    % legend('Data','Onset Limit','Offset Limit','Onset Time','Offset Time')


    AllOnOffsetTime(i,1)=OnsetTime;
    AllOnOffsetTime(i,2)=OffsetTime;



end %end for each trial



end %end function

