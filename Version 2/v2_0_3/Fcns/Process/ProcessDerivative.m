%Take the derivative of the data
%Inputs: time, Data, Order
%Outputs: derivative of data

function [newTime,DerivedData]=ProcessDerivative(DataTime,Data,SampleRate,Order)

%take derivative
DiffData=diff(Data,Order,2);
DerivedData=DiffData./(1/SampleRate)^Order;

newTime=DataTime(Order+1:end);

end

