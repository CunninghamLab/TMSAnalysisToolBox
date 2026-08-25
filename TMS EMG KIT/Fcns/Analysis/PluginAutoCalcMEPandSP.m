%Auto calculate MEP or SP default metrics
%MEP = amplitude and area under the curve
%SP = percent decrease and normalized area under the curve

function [MEPAmp,MEPArea,SPPercDecrease,SPArea ]=PluginAutoCalcMEPandSP(i,app,PreStimData,AnalyzeData,AnalysisType)

    switch AnalysisType
        case "MEP"
            %Calculate amplitude and area
            if ~isempty(app.ProcessOrder) && any(string(app.ProcessOrder(:,1)) == "Rectify") && isempty(app.Processed_Conditions_DataAll{3})
                MEPAmp=nan;
            else
                Peak=max(AnalyzeData);
                Trough=min(AnalyzeData);
                MEPAmp=Peak-Trough;
            end
            MEPArea=trapz(abs(AnalyzeData)); %trapezoidal integration method, for MEP area should be of absolute value of data
            
            SPPercDecrease=nan; SPArea=nan;
    
        case "SP"
            SPPercDecrease=100-(mean(AnalyzeData)/mean(PreStimData(:,i))*100); % percent of silent period mean value decreased compared with pre stimulation period
    
            %Normalized area of SP
            SPArea=trapz(AnalyzeData); %trapezoidal integration method
            %mean preStimulation times length of silent period duration equals the area
            areaMeanSP = abs(mean(PreStimData(:,i))*length(AnalyzeData));
            SPArea = (1-SPArea/areaMeanSP)*100;      %normalized area of slient period unit

            MEPAmp=nan; MEPArea=nan;
    
    end

end