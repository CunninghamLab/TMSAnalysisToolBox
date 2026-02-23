%Process Data
%Inputs: app, Data, Data time
%Outputs:updated Data and Data time


function [DataTime,Data]=ProcessingData(app,Data,DataTime)

%%
%app.ProcessOrderListBox.ItemsData=[1 2 5 8 7 ]; %DS,Delay,rect,lowpass

Order=app.ProcessOrderListBox.ItemsData;
%ItemData 8 == Rectify
RectifyPres= Order == 8;
if any(RectifyPres) %if rectify is present in the process order
    if length(Order) == 1 %if there is only one process
        SaveNonRectData{1}=1; %0-don't save non-rectified data, 1-save non-rectified data
    else 
        RectifyLoc2=find(RectifyPres);
        if RectifyLoc2 ~= length(RectifyPres) %if Rectify isn't at the end
            if any(Order == 7) && Order(end) == 7 && Order(end-1) == 8  %if Rectify is second to last and the last is delay
                SaveNonRectData{1}=1; %need to save data after the delay
                %change so delay is before rectify 
                app.ProcessOrderListBox.Items{end}='Rectify';
                app.ProcessOrderListBox.Items{end-1}='Delay';
                app.ProcessOrderListBox.ItemsData(end)=8;
                app.ProcessOrderListBox.ItemsData(end-1)=7;
            else
                SaveNonRectData{1}=0;
            end %end specific check for Delay
            
        else
            SaveNonRectData{1}=1;
        end %end if rectify is at the end
    end %end if there is only one process
else
    SaveNonRectData{1}=0;
end %end if rectify is present

%%
%Process data
D=0; x=1; app.ProcessOrder=cell(length(app.ProcessOrderListBox.ItemsData),2);
for i=app.ProcessOrderListBox.ItemsData

    %if i is the second to last item and the data needs to be saved do that, then run last process which should be Rectify
    if  (SaveNonRectData{1} == 1) && (length(app.ProcessOrderListBox.ItemsData) == 1 || i == app.ProcessOrderListBox.ItemsData(end-1)) 
        app.Processed_Conditions_DataAll{1,3}=Data;
        app.Processed_Conditions_DataAll{1,4}=DataTime;
    end

    switch i

        case 7 %Delay
            DataTime=ProcessDelay(DataTime,app.DelayEditField.Value);
            app.ProcessOrder{x,2}=app.DelayEditField.Value;

        case 8 %Rectify
            Data=ProcessRectify(Data);

        case 1 %Down sample
            [DataTime,Data,newSampleRate]=ProcessDownSample(DataTime,Data,app.FactorEditField.Value,app.AllSampleRate(1));
            app.SampleRateLabel.Text=['Sample Rate: ' num2str(newSampleRate) 'Hz']; 
            %app.AllSampleRate
            app.ProcessOrder{x,2}=app.FactorEditField.Value;
            D=1;

        case 2 %Notch
            if D==1
                SampleRate=newSampleRate;
            else
                SampleRate=app.AllSampleRate(1);
            end

            Data=ProcessNotch(app,Data, SampleRate,app.NotchEditField.Value);
            app.ProcessOrder{x,2}=app.NotchEditField.Value;

        case 3 %Band pass
            if D==1
                SampleRate=newSampleRate;
            else
                SampleRate=app.AllSampleRate(1);
            end

            Data=ProcessBandPass(app,Data, SampleRate,app.BandPassStartEditField.Value,app.BandPassEndEditField.Value);
            app.ProcessOrder{x,2}=[app.BandPassStartEditField.Value app.BandPassEndEditField.Value];

        case 4 %Low pass
            if D==1
                SampleRate=newSampleRate;
            else
                SampleRate=app.AllSampleRate(1);
            end

            Data=ProcessLowPass(app,Data, SampleRate,app.LowPassEditField.Value);
            app.ProcessOrder{x,2}=app.LowPassEditField.Value;

        case 5 %High pass
            if D==1
                SampleRate=newSampleRate;
            else
                SampleRate=app.AllSampleRate(1);
            end

            Data=ProcessHighPass(app,Data, SampleRate,app.HighPassEditField.Value);
            app.ProcessOrder{x,2}=app.HighPassEditField.Value;

        case 6 %Derivative
            if D==1
                SampleRate=newSampleRate;
            else
                SampleRate=app.AllSampleRate(1);
            end

            [DataTime,Data]=ProcessDerivative(DataTime,Data,SampleRate,app.DerivativeEditField.Value);
            app.ProcessOrder{x,2}=app.DerivativeEditField.Value;


        case 9 %Custom
            disp('Run Custom'); %!!! custom feature is not implemented

    end

    x=x+1;

end %end process loop

if D==1
    app.DownSampledRate=newSampleRate;
end

end