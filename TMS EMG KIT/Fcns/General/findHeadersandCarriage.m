

function [Headers,CarriageReturn,FileData]=findHeadersandCarriage(ResultFile)

ResultFile_FID=fopen(ResultFile,'r');
FileData=fread(ResultFile_FID); %double dec
fclose(ResultFile_FID);

CarriageReturn=find(FileData == 13); %13 dec == carriage return
if isempty(CarriageReturn) %files written on macOS/Linux use LF only (10 dec)
    CarriageReturn=find(FileData == 10);
end
CarriageReturn=CarriageReturn(1);
HeaderLine=char(FileData(1:CarriageReturn-1))';
Headers=string(strsplit(HeaderLine,','));

end