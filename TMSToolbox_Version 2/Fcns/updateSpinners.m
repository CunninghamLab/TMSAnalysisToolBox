%Update plots based on spinner values
%Inputs:public app properties, Trial selected
%Outputs:Data to be plotted, current channel title, time

function [Plotted_Data, CurrChTitle, Time] = updateSpinners(app,TrialSelected)

%Determine Data Type
switch app.DataType 
    case {1,2,4,5,6,8} 
        switch app.EventType
            case 1 %Block
                Plotted_Data=app.Block.Data{app.Ch,TrialSelected};
                CurrChTitle=app.Block.Titles(app.Ch);
                Time=0:1/app.AllSampleRate:length(Plotted_Data)/app.AllSampleRate-(1/app.AllSampleRate);
                if app.PreTriggerTime ~= 0
                    Time=Time-app.PreTriggerTime/1000;                                       %Adjust for Pre trigger setting
                end
                if app.PostTriggerTime ~=0                                                             %Adjust for Post trigger setting
                    FramesAfter=app.PostTriggerTime/1000*app.AllSampleRate+app.PreTriggerTime/1000*app.AllSampleRate;
                    %check for if post trigger time is beyond the dataset
                    if FramesAfter > length(Plotted_Data)
                        FramesAfter=length(Plotted_Data);
                        app.PostTriggermsEditField.Value=FramesAfter/app.AllSampleRate*1000;
                    end
                    Time=Time(1:FramesAfter);
                    Plotted_Data=Plotted_Data(1:FramesAfter);
                end
            case 2 %Comments
                if app.autoImportEvents == 1
                    comCol=app.ImportEventsNameL;                       
                else
                    Com=string(app.SelectedComtxt);
                    comCol=app.EventNameListBox.ItemsData == Com; %ver 11 and up
                end
                Plotted_Data=app.Comments.Data{TrialSelected,comCol}(app.Ch,:);
                CurrChTitle=app.Block.Titles(app.Ch);
                Time=-app.PreDur:1/app.AllSampleRate:app.PostDur;
            case 3 %Events
                Plotted_Data=app.Events.Data{TrialSelected}(app.Ch,:);
                CurrChTitle=app.Block.Titles(app.Ch);
                Time=-app.PreDur:1/app.AllSampleRate:app.PostDur;
        end
        
    case 3 %Signal 
        %Block with states
        Plotted_Data=app.Block.Data{app.StateV}{TrialSelected}(:,app.Ch);
        Time=0:1/app.AllSampleRate:length(Plotted_Data)/app.AllSampleRate-(1/app.AllSampleRate); %Time defined by sample rate and pre/post duration settings
        if app.PreTriggerTime ~= 0
            Time=Time-app.PreTriggerTime/1000;                                                 %Adjust for Pre trigger setting
        end
        if app.PostTriggerTime ~=0                                                             %Adjust for Post trigger setting
            FramesAfter=app.PostTriggerTime/1000*app.AllSampleRate+app.PreTriggerTime/1000*app.AllSampleRate;
            %check for if post trigger time is beyond the dataset
            if FramesAfter > length(Plotted_Data)
                FramesAfter=length(Plotted_Data);
                app.PostTriggermsEditField.Value=FramesAfter/app.AllSampleRate*1000;
            end
            Time=Time(1:FramesAfter);
            Plotted_Data=Plotted_Data(1:FramesAfter);
        end
        CurrChTitle=app.Block.Titles{app.Ch};

    

end %end DataType switch

end





