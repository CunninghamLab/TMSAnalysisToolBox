%Shift time by the amount inputted by the user
%Inputs:Time,Delay amount in milliseconds
%Outputs:time shifted by delay amount


function ShiftedTime=ProcessDelay(Time,Delay)

%Move time back 
ShiftedTime=Time-Delay/1000; %delay is in ms while time is saved in seconds



end