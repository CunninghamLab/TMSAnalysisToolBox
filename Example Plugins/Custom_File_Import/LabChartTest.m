%Custom File Upload Example LabChart
%{
Use this script as a guide for creating a custom file upload pluggin.
Data must be put in the "Block" form. This is the base for the other event 
types.
Block form is created by separating the data using "natural" breaks in the 
data.

Output Variables:
Variable Name       DataType                            Description
numChannels         1x1 double                          number of data channels  
numBlocks           1x1 double                          number of Blocks
Block               See below                           Block: data, channel titles, coordinates (if applicable)
AllSampleRate       1x1 double                          sample rate of data collected
EventTypeOptions    1xnumEventTypes double vector       Event Types to be included (must include Blocks)
Data                variable                            original load in of the data
CommentLocAll       See below                           Locations of all comments
numTrialsC          1xnumCommentTxtOpts double vector   number of comments for each comment text option
CommentTxtOptions   numCommentTxtOptsx1 string vector   comment text options

%=======================================================================================================================
Block variable 
Block=Struct()
Block.Data=numChannels x numBlocks cell (each cell 1 x numFrames double vector)
Block.Titles=numChannels x num string vector 
Block.Coordinates=numBlocks x 3 (xyz)    

If there is a layer above the blocks (i.e. Signal Data)
Block.Data:
	First column=numblocks x 1 cell for each state, contains numFrames x numChannels cell
	Second column=1x1 cells of the state labels
%=======================================================================================================================
CommentLocAll Variable
CommentLocAll=numComments x numCommentTextOptions cell array
			[Comment Block   CommentLoc]
			First column = locations for first comment text option
			Second column = locations for second comment text option

%}


function [numChannels,numBlocks,Block,AllSampleRate,EventTypeOptions,Data,CommentLocAll,numTrialsC,CommentTxtOptions]=LabChartTest(app,i)

%% Event Type Options
%{
    Which event types would you like to be available for this data type
    1 - Block
    2 - Comments
    3 - Events
%}
EventTypeOptions=[1 2 3]; %fill vector with options % =========


%% Import Data
Data=load(fullfile(app.Path, app.Filename{i}));      % =========


%% Extract blocks

%Extract base data
%Channel Titles
ChannelTitles=string(Data.titles);

%Number of channels and blocks
numChannels=length(Data.datastart(:,1)); % =========
numBlocks=length(Data.datastart(1,:));   % =========


%Extract and organize data by blocks
for i=1:numBlocks
    for e=1:numChannels
        [Ch_Blk, ChTitle,SampleRate]=ExtractData(Data,ChannelTitles,e,i);
        Block.Data{e,i}=Ch_Blk;                                                     % =========
        if i==1 %for first block only
            Block.Titles(e,1)=ChTitle; %Assumed same channel titles between blocks  % =========
            AllSampleRate=SampleRate; %Assumed same sample rate between blocks      % =========
        end
    end
end


%% Extract Comments
%If no Comments
% CommentLocAll=[];
% numTrialsC=[];
% CommentTxtOptions="";

%Otherwise add code to find comment locations
numComments=length(Data.com(:,1));                            %number of comments
Comments.CommentTxtOptions=string(Data.comtext);              %comment text options
numComTxtOpt=length(Comments.CommentTxtOptions);              %number of comment texts


%Extract and organize data
e=1; e2=1;
CommentLocAll=cell(1,numComTxtOpt);
for i2=1:numComTxtOpt %for each comment text
    for i=1:numComments %look through each comment
        CommentInfo=Data.com(i,:);  %channel, block, frame, type of comment, comment indx in comtext
        CommentBlk=CommentInfo(2); CommentTxt=CommentInfo(5); CommentLoc=CommentInfo(3);

        if Comments.CommentTxtOptions(CommentTxt) == Comments.CommentTxtOptions(i2) %if this comment text matches the current commen text
            CommentLocAll{e,e2}(1,2)=CommentLoc;                % =========
            CommentLocAll{e,e2}(1,1)=CommentBlk;                % =========
            e=e+1;
        end

    end
    e=1;
    e2=e2+1;
end
numTrialsC=zeros(1,numComTxtOpt);
for i=1:numComTxtOpt
    TrialC=CommentLocAll(:,i);
    numEmpty=sum(cellfun(@isempty,TrialC));
    numTrialsC(1,i)=length(TrialC)-numEmpty;           % =========
end

CommentTxtOptions=Comments.CommentTxtOptions;          % =========

                    

end %end function

%% Functions
%Create Block data
%Inputs: LabChartData, channel titles, channel number, block number
%Outputs: Data for a single block for one channel, that channel title, the sample rate

function [Ch_Blk, ChTitle,SampleRate] = ExtractData(Data,ChannelTitles,Ch,Blk)

    %channel data in block
    Ch_Blk=Data.data(Data.datastart(Ch,Blk):Data.dataend(Ch,Blk));

    %check units=
    Units=Data.unittext;
    if length(Units) ==1
        %all channels and blocks have the same unit
        DataUnit=Units;           
    end
    
    %convert data to volts if needed
    switch DataUnit  
        case 'V' %volts
            Multiplier_V=1;
        case 'mV' %millivolts
            Multiplier_V=0.001;
    end
    Ch_Blk=Ch_Blk*Multiplier_V;

    %channel titles
    ChTitle=ChannelTitles(Ch);

    %Sample Rate
    SampleRate=Data.samplerate(Ch,Blk);

end












