%Convert cell array to comma or underscore delimited character array
%Input: Cell array containing the trial numbers in character format {'1','2','3'}
%Output: Comma delimited character array '1,2,3'

function ConvertedArray=convert2commadel(CellArray,Del)



CellArrayFormat=cell(1,length(CellArray)*2-1); %Create new cell array 
CellArrayFormat(1:2:end)=CellArray;            %add original cell array data
CellArrayFormat(2:2:end)={Del};                %insert commas between


if Del == ','
    CellArray=char(CellArrayFormat)';              %convert to character
    if length(CellArray(:,1)) > 1     %If the number is double digits and up
        CellArrayFormat=[];
        for i=1:length(CellArray(1,:))
            CellArrayFormat=[CellArrayFormat CellArray(:,i)'];
        end
        CellArrayFormat(CellArrayFormat==' ')=[];
        CellArray=CellArrayFormat;
    end

else
    CellArray=CellArrayFormat;
end

ConvertedArray=CellArray;

end
