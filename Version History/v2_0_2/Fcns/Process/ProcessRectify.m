%Rectify data
%Inputs: data
%Outputs: rectified data

function RecData=ProcessRectify(app, Data)

%Rectify data
if app.FullWaveButton.Value == 1
    app.HalfWaveButton.Value = 0;
    RecData=abs(Data);
elseif app.HalfWaveButton.Value == 1
    app.FullWaveButton.Value = 0;
    RecData = Data;
    RecData(RecData < 0) = 0;
end

end