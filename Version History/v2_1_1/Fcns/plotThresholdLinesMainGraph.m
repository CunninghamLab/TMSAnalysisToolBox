%Plots thresholds onto main graph
%Inputs:app, Analyze Method ("MEP","SP"), median of the onset thresholds, median of the offset thresholds
%Outputs: adds y lines to UIAxes

function plotThresholdLinesMainGraph(app,AnalyzeMethod,MedianOnsetThresh,MedianOffsetThresh)

%There is processing and the data is rectified
if ~isempty(app.ProcessOrder) && any(string(app.ProcessOrder(:,1)) == "Rectify") 
    switch AnalyzeMethod
        case "MEP"
            app.OnsetLimitLines=yline(app.UIAxes,MedianOnsetThresh,'r:','LineWidth',2);
            app.OffsetLimitLines=yline(app.UIAxes,MedianOffsetThresh,'b--','LineWidth',2);
        case "SP"
            app.OnsetLimitLines=yline(app.UIAxes,MedianOnsetThresh,'r:','LineWidth',2);  
            app.OffsetLimitLines=yline(app.UIAxes,MedianOffsetThresh,'b--','LineWidth',2);
    end
end
%}


end