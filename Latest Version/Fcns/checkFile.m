%checkFile
%Check csv file headers to see if any need to be added
%Inputs: FileData,Headers,expected HeaderLineText to look for,CarriageReturn location,ResultFile name
%Outputs: updates the csv file if needed

function checkFile(FileData,Headers,HeaderLineText,CarriageReturn,ResultFile)
%see if the header has the expected header line text
ContainsHeaderLine=find(contains(Headers,HeaderLineText));

if isempty(ContainsHeaderLine) %none the the header texts are present
    updateFile(FileData,HeaderLineText,CarriageReturn,ResultFile);
elseif length(ContainsHeaderLine) ~= length(HeaderLineText) %only some were found
    %find which ones are present
    Present=Headers(ContainsHeaderLine);
    if any(contains(HeaderLineText,Present)) %if any of the present headers are within the header text line
        HeadersToAdd=HeaderLineText(~contains(HeaderLineText,Present));
        if ~isempty(HeadersToAdd) %if everything that is present isn't within the header line text update the file
            updateFile(FileData,HeadersToAdd,CarriageReturn,ResultFile);
        end
    end
end

end

%updates the csv file
function updateFile(FileData,HeaderLineText,CarriageReturn,ResultFile)

HeaderLineTextC=char(join(HeaderLineText,','))'; %make the header text
NewFileData=[FileData(1:CarriageReturn-1);double(','); double(HeaderLineTextC); FileData(CarriageReturn:end)]; %insert new headers into file
ResultFile_FID=fopen(ResultFile,'w+');
fwrite(ResultFile_FID,NewFileData);
fclose(ResultFile_FID);




end