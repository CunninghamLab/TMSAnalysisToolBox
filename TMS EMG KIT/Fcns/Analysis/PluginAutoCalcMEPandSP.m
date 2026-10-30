%Auto calculate MEP or SP default metrics
%MEP = amplitude and area under the curve
%SP = percent decrease and normalized area under the curve

function PluginAutoCalcMEPandSP(i,app,PreStimData,AnalyzeData,AnalysisType)

    switch AnalysisType
        case "MEP"
            %Calculate amplitude and area
            if ~isempty(app.ProcessOrder) && any(string(app.ProcessOrder(:,1)) == "Rectify") && isempty(app.Processed_Conditions_DataAll{3})
                app.MEP_Amp=nan;
            else
                Peak=max(AnalyzeData);
                Trough=min(AnalyzeData);
                app.MEP_Amp(i,:)=Peak-Trough;
            end
            
            %remove nans %!!!!!
            AnalyzeData(isnan(AnalyzeData))=[];
            app.MEP_Area(i,:)=trapz(abs(AnalyzeData)); %trapezoidal integration method, for MEP area should be of absolute value of data

        case "SP"
            app.SP_PercDecrease(i,:)=100-(mean(AnalyzeData,'omitnan')/mean(PreStimData(:,i),'omitnan')*100); % percent of silent period mean value decreased compared with pre stimulation period

            %Normalized area of SP
            %remove nans %!!!!!
            AnalyzeData(isnan(AnalyzeData))=[];
            SPArea=trapz(AnalyzeData); %trapezoidal integration method
            %mean preStimulation times length of silent period duration equals the area
            areaMeanSP = abs(mean(PreStimData(:,i),'omitnan')*length(AnalyzeData));
            app.SP_Area(i,:) = (1-SPArea/areaMeanSP)*100;      %normalized area of slient period unit
    
    end

end