%Extract comments data from LabChart data
%Creates data for each comment based on pre and post durations
%Inputs: LabChart data, LabChart Block data, number of channels, number of comments, pre duration, post duration, samples rate
%Outputs: Comments data, number of trials for each comment type

function [Comments,numTrialsC] = ExtractLB_CommentsData(labChart_Data,Block,numChannels,numComments,PreDur,PostDur,AllSampleRate)


Comments.CommentTxtOptions=string(labChart_Data.comtext); %text of comments
numComTxtOpt=length(Comments.CommentTxtOptions);              %number of comment texts

%Extract and organize data
PreDurValue=PreDur*AllSampleRate;
PostDurValue=PostDur*AllSampleRate;
e=1; Data2=zeros(numChannels,length(-PreDurValue:PostDurValue)); numTrialsC=zeros(1,numComTxtOpt);
for i2=1:numComTxtOpt %for each comment text
    for i=1:numComments %extract and save the data
        CommentInfo=labChart_Data.com(i,:);  %channel, block, frame, type of comment, comment indx in comtext
        CommentBlk=CommentInfo(2); CommentTxt=CommentInfo(5); CommentLoc=CommentInfo(3);
        if Comments.CommentTxtOptions(CommentTxt) == Comments.CommentTxtOptions(i2) %only extract data for the current comment text (numComTxtOpt)
            %Parse out desired data
            StartDiff=0; EndDiff=0;
            Start=CommentLoc-PreDurValue; End=CommentLoc+PostDurValue;
            if Start <= 0 %if start is less than zero, fill in with NaNs
                StartDiff=-Start+1;
                Start=1;
            end
            DataLength= length(Block.Data{1,CommentBlk});
            if End > DataLength %if end is greater than the length of the data, fill with NaNs
                EndDiff=End-DataLength;
                End=DataLength;
            end
            Data=Block.Data(:,CommentBlk);
            %reformat and cut data
            for i3=1:numChannels
                Data2(i3,:)=[nan(1,StartDiff) Data{i3}(Start:End) nan(1,EndDiff)]; %fill beginning and end with NaNs if needed
            end
            Comments.Data{e,i2}=Data2; %each column is for a comment text, order of comment texts is the same
            e=e+1;
        end

    end
    numEmpty=cellfun(@isempty,Comments.Data);          %Find number of empty vectors to determine
    numTrialsC(i2)=length(find(numEmpty(:,i2) == 0));  %number of comments for each comment text
    e=1;
end %end extract and organize data

end