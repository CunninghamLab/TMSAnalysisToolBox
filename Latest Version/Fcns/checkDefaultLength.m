%check default value lengths and resize if needed
%Inputs: the user's DefaultValues variable from the plugin function, the expected number of variables
%Outputs: DefaultValues, either adjusted to match the numVar or is the same if the lengths already match

function DefaultValues=checkDefaultLength(DefaultValues,numVar)


DefaultValuesLength=length(DefaultValues);
if  DefaultValuesLength ~= numVar
    if DefaultValuesLength > numVar 
        DiffVarLength=DefaultValuesLength-numVar;
        DefaultValues(end-DiffVarLength+1:end)=[];

    elseif DefaultValuesLength < numVar 
        DiffVarLength=numVar-DefaultValuesLength;
        DefaultValues(end+1:end+DiffVarLength)=zeros(1,DiffVarLength);
    end
end



end