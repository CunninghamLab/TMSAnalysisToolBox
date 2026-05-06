%Extract LabChart Block Data
%Inputs: LabChart data
%Outputs: Number of channels, number of blocks, Block data, sample rate for all blocks (assumed to be the same between all blocks)

function [numChannels,numBlocks,Block,AllSampleRate]=ExtractLB_BlockData(labChart_Data)

%Extract base data
%Channel Titles
ChannelTitles=string(labChart_Data.titles);

%Number of channels and blocks
numChannels=length(labChart_Data.datastart(:,1));
numBlocks=length(labChart_Data.datastart(1,:));


%Extract and organize data by blocks
for i=1:numBlocks
    for e=1:numChannels
        [Ch_Blk, ChTitle,SampleRate]=ExtractData(labChart_Data,ChannelTitles,e,i);
        Block.Data{e,i}=Ch_Blk;
        if i==1 %for first block only
            Block.Titles(e,1)=ChTitle; %Assumed same channel titles between blocks
            AllSampleRate=SampleRate; %Assumed same sample rate between blocks
        end
    end
end



end %end function

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


%% Notes on LabChart Data Structs
%{
From https://www.adinstruments.com/support/knowledge-base/how-does-matlab-open-exported-data?srsltid=AfmBOoq49cn7uXr7wo2pcaPN9wSH4k67rZlViX8Pt7HqVebdgXSnTVN7

data=[ch1_blk1, ch1_blk2, ch1_blkn,ch2_blk1...];   
       Data vector that will be parsed from

datastart=[ch1_blk1,ch1_blk2,ch1_blkn;              
            ch2_blk1, ch2_blk2,ch2_blkn;...];
            Start values in data vector

dataend=[ch1_blk1,ch1_blk2,ch1_blkn;              
          ch2_blk1, ch2_blk2,ch2_blkn;...];
          End values in data vector

          EX: Data in channel 1 block 2
              data12=data(datastart(1,2):dataend(1,2));

blocktimes=[time of day when first sample in block 1 was recorded, time of day when first sample in block 1 was recorded...];
             These are serial date numbers
             If only a selection is exported this is the first sample in the selection

tickrate=[block 1 max sample rate, block 2 max sample rate, ...];
           Max sample rate of each block

samplerate[ch1_blk1, ch1_blk2, ch1_blkn;
            ch2_blk1, ch2_blk2,ch2_blkn; ...];
            Sample rate of each channel and block
            Empty channels will have sample rate of 0

unittext=vector of possible data units
         [V;N;%CO2];

unittextmap=[unittext index for ch1_blk1, unittext index for ch1_blk2,unittext index for ch1_blkn;
              unittext index for ch2_blk1, unittext index for ch2_blk2, unittext index for ch2_blkn; ...];

titles=names of channels
        Character array, convert to string to parse out titles

rangemin=[min value in ch1_blk1, min value in ch1_blk2, min value in ch1_blkn;
           min value in ch2_blk1, min value in ch2_blk2, min value in ch2_blkn;...];

rangemax=[max value in ch1_blk1, max value in ch1_blk2, max value in ch1_blkn;
           max value in ch2_blk1, max value in ch2_blk2, max value in ch2_blkn;...];
           Empty channels will have a min/max value of 0

scaleoffset and scale units matrices
        LabChart.m converts the 16-bit data into the Chart View units by adding the "scaleoffset" 
        to the data or range and multiplying it with the appropriate "scaleunits" factor. 
        Empty channels in a block will have scaleunits = scaleoffset = 0.

comtext=vector of comment text
        Character array, convert to string to part out comments

com=[comment channel (-1 means all channels), 
      block number the comment was made in, 
      the position of the comment in the block about the tick rate of the block, frame of comment based on max sample rate
      type of comment (user comments (1) and event markers (2), 
      comtext index that holds the string of the commment];
     One row for each comment

firstsampleoffset= number between 0 and 1
       In multi-rate cases, a block can start part-way through a low-rate sample.
       This matrix contains a number between 0 and 1, which indidcates the fractional offset in units of samples.
       It shows how long after the start of the first sample in the block, the block started.
        
       EX: firstsampleoffset=0.9s
           The block start 90% of the way through the first sample stored in that channel for the block.

%}











