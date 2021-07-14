function [CentralizedData, CentralizedName] = Centralization(Data,Name,Percentile)
%CENTRALIZATION Summary of this function goes here
%   Detailed explanation goes here
for iCondition = 1: size(Data)
    iData = Data{iCondition};
    iName = Name{iCondition};
    DeleteNum = floor(size(iData,2)*Percentile/100);  % number of trials to be deleted in one side 
    if DeleteNum>=1
        [xs, index] = sort(iData);
        resultLow = iData(sort(index(1:DeleteNum)));
        resultHigh = iData(sort(index(end-DeleteNum+1:end)));
        allResult = [resultLow resultHigh];
        for kResult = 1: length(allResult)
            kOrder = find(iData==allResult(kResult));
            iData(kOrder) = [];
            iName(kOrder) = [];
            Data{iCondition} = iData;
            Name{iCondition} = iName;
        end
    end
    
end

CentralizedData = Data;
CentralizedName = Name;

end

