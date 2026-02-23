%Custom File Upload Template
%{
Use this script as a guide for creating a custom file upload pluggin.
Change the function name to the desired file type name. (EX: LabChart)
Place this function in the "Plugins" folder.
Data must be put in the "Block" form. This is the base for the other event 
types.
Save data in Volts. Convert if needed.


Output Variables:
Variable Name       DataType (row x column)             Description
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
Block.Titles=numChannels x 1 string vector 

%=======================================================================================================================
CommentLocAll Variable
CommentLocAll=numComments x numCommentTextOptions cell array
			First column = locations for first comment text option
			Second column = locations for second comment text option
            Each cell contains  [Block_the_comment_is_located_in   Frame_in_Block_the_Comment_is_located_in]
%}


function [numChannels,numBlocks,Block,AllSampleRate,EventTypeOptions,Data,CommentLocAll,numTrialsC,CommentTxtOptions]=CustomFileTemplate(app,i)

%% Event Type Options
%{
    Which event types would you like to be available for this data type
    1 - Block (Required)
    2 - Comments
    3 - Events
%}
EventTypeOptions=[1]; %fill vector with options, must contain 1 (Blocks)


%% Import Data
File=fullfile(app.Path, app.Filename{i}); %Path and Filename will be selected by the user before this function runs

%Load in the data from "File"


%% Extract blocks
%Block form is created by separating the data using "natural" breaks in the 
%data (ex: start/stop recording times)
%Output Variables
%{
numChannels
numBlocks
Block 
AllSampleRate
%}

%Add code to separate data into Blocks-------------------------------------


%% Extract Comments
%Output Variables
%{
CommentLocAll
numTrialsC
CommentTxtOptions
%}

%Add code to find comment locations within the blocks----------------------


%Check for Comments Data, if there is no code above, create empty variables
if ~exist('CommentLocAll','var') || ~exist('numTrialsC','var') || ~exist('CommentTxtOptions','var')
    CommentLocAll=cell(1,1);
    numTrialsC=[];
    CommentTxtOptions="";
end

end %end function

%% Functions
%Add any additional functions used for extracting the data here



