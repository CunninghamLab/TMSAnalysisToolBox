function out = readcfs(cfsn,plotdata)
% This function reads the CED's signal file (.CFS) and gives a structure as
% an output.
%
% cfsn: full file name of the CFS file (including directory, if not in the
% same folder)
%
% out: structure with the following sub-structures:
% out.GH: General Header
% out.FVI: File Variable Information
% out.ChI: Channel Information (if different across data segments, a vector
% will be present instead of scalar)
% out.DSVI: Data Segment Variable Information (if different across data
% segments, a vector will be present instead of scalar)
% out.GDSH: General Data Segment Header
%
% If the number of data points are the same across channels:
% out.ChData: Channel data, organized as [data points x channels]
% out.Chts: Channel timestamps, organized as [data points x channels]
% If segment sizes are same across the channels, out.Chts is a vector
%
% If the number of data points across the channels are different:
% out.ChData: Channel data, organized as matrix of [time points x channels]
% or cell array {data segments x channels} if the number of data points are
% different across channels
% ChData can be converted from cells to array (given the channels share
% the same number of elements) using cell2mat(ChData(:,Ch#));
% out.Chts: Channel timestamps for the data, organized as a vector (if all
% channels share the timestamps) or matrix of [time points x channels] (if
% number of timestamps are the same across the channels or cell array {data segments x channels}
% Chts can be converted from cells to array (given the timestamps share
% the same number of elements) using cell2mat(Chts(:,Ch#));
%
% Based on PDF document titled: CFS - The CED Filing System - October 2006
% This PDF is part of: http://ced.co.uk/files/MS54.exe
%
% By,
% Pratik Y. Chhatbar, PhD, MBBS
% Assistant Professor, Duke Neurology
% pratik.chhatbar@duke.edu
% 27 May 2020

if ~exist('plotdata','var') || isempty(plotdata)
    plotdata=true; % plotting of data by default
end

fid = fopen(cfsn);

%% File Header

% A file header with n file variables, c data channels and d data section
% variables can be shown as:
%
% GH: General header
%
% ChI:
% Channel 0 information
% Channel 1 information
% ...
% Channel c-1 information
%
% FVI:
% File variable 0 information
% File variable 1 information
% ...
% File variable n-1 information
% A system file var to calculate the size of var n-1
%
% DSVI:
% DS variable 0 information
% DS variable 1 information
% ...
% DS variable d-1 information
% A system DS var to calculate the size of var d-1
%
% FV (combined as a part of FVI as a sub-structure)
% File variable 0
% File variable 1
% ...
% File variable n-1


% GH: General Header
%
%       Description                             Bytes   Offset	Type
% 1     Marker to identify the file             8       0x00    char[8]
% 2     File name (up to 12 characters)         14      0x08    string[13]
% 3     File size (in bytes)                    4       0x16    long
% 4     Time when file was created              8       0x1a    char[8]
% 5     Date when file was created              8       0x22    char[8]
% 6     Channels in data section (0-99)         2       0x2a    short
% 7     Number of file variables (0-99)         2       0x2c    short
% 8     Number of data section variables (0-99) 2       0x2e    short
% 9     Byte size of file header                2       0x30    short
% 10    Byte size of data section header        2       0x32    short
% 11    Last data section header offset         4       0x34    long
% 12    Number of data sections                 2       0x38    WORD
% 13    Disk block size rounding (1=none)       2       0x3a    WORD
% 14    File comment (up to 72 characters)      74      0x3c    string[73]
% 15    Pointer table offset                    4       0x86    long
% 16    Reserved space for future expansion     40      0x8a    Zero filled

GH.Marker = fread(fid,[1,8],'uint8=>char'); % The last character is the
% revision level of the CFS (!=1, "=2 and so on in ASCII sequence - 33,34,35..).
GH.FileName = fread(fid,[1,14],'uint8=>char');
GH.FileSz = fread(fid,1,'ulong=>double');
GH.FileTime = fread(fid,[1,8],'uint8=>char');
GH.FileDate = fread(fid,[1,8],'uint8=>char');
NCh = fread(fid,1,'ushort=>double'); % Number of Channels
GH.NCh = NCh;
NFV = fread(fid,1,'ushort=>double'); % Number of File Variables
GH.NFV = NFV;
NDSV = fread(fid,1,'ushort=>double'); % Number of Data Section Variables
GH.NDSV = NDSV;
GH.FHSz = fread(fid,1,'ushort=>double'); % File Header Size
GH.DSHSz = fread(fid,1,'ushort=>double'); % Data Section Header Size
GH.LastDSHOffset = fread(fid,1,'ulong=>double');
NDS = fread(fid,1,'ushort=>double'); % Number of Data Sections
GH.NDS = NDS;
GH.DiscBlkSzRound = fread(fid,1,'ushort=>double');
are = fread(fid,[1,74],'uint8=>char');
GH.PtrTblOffset = fread(fid,1,'ulong=>double');
fseek(fid,40,'cof'); % skipping 40 bytes of reserved space

out.GH = GH;

% Channel Information:
% Each data section can hold up to 99 channels.
%
%       Description                                     Bytes   Offset  Type
%   1   Channel name (up to 20 characters)              22      0x00    string[21]
%   2   Y axis units (up to 8 characters)               10      0x16    string[9]
%   3   X axis units (up to 8 characters)               10      0x20    string[9]
%   4   Data type e.g. int2, rl4                        1       0x2a    TDataType
%   5   Data kind (equalSpaced, matrix or subsidiary)   1       0x2b    TCFSKind
%   6   Byte space between elements                     2       0x2c    short
%   7   Next (matrix) or master (subsidiary) channel    2       0x2e    short
%
% TDataType: A set of 8 constants are defined which are used to describe a data type. These are used in
% describing file or data section variables and channel data:
% 0 INT1  Data in the range -128 to 127 stored in a single byte of memory.
% 1 WRD1  Data in the range 0 to 255 stored in a single byte of memory.
% 2 INT2  Data in the range -32768 to 32767 stored in two bytes of memory
%       (for example 1401 analogue input data).
% 3 WRD2  Data in the range 0 to 65535 stored in two bytes of memory.
% 4 INT4  Data stored as 4 byte integer data in the range -2147483648 to 2147483647.
% 5 RL4   IEEE format 4 byte floating point numbers.
% 6 RL8   IEEE format 8 byte floating point numbers.
% 7 LSTR  Character data.
%
% TCFSKind: Three constants are defined which describe the storage method for channel data:
% #define EQUALSPACED 0
% #define MATRIX 1
% #define SUBSIDIARY 2
% "equalSpaced" data is typically analogue data where the time interval between each data
% point is fixed and the sequential position of the data in the channel is important.
% "subsidiary" data is like "equalspaced" data, but is extra data associated with an
% equalspaced channel such as errors. "matrix" data implies an N dimensional array
% (where N can be 1 for a simple list of values). This is used for a set of (x,y) positions, for
% example. The next channel description holds the channel number of the next channel in
% the matrix, or the current channel number for a one dimensional matrix. The last channel
% in a matrix points back at the first channel.

ChI(NCh,1) = struct;
for ii = 1:NCh
    ChI(ii).Name = fread(fid,[1,22],'uint8=>char');
    ChI(ii).Yunit = fread(fid,[1,10],'uint8=>char');
    ChI(ii).Xunit = fread(fid,[1,10],'uint8=>char');
    ChI(ii).DataType = fread(fid,1);
    ChI(ii).DataKind = fread(fid,1);
    ChI(ii).ByteSpace = fread(fid,1,'ushort=>double');
    ChI(ii).NextChan = fread(fid,1,'ushort=>double');
end

% File Variable Information:
% By convention the first file variable should be used to store a dummy value, with the
% description and unit fields used to store identification for where the file was produced.
%
%       Description                         Bytes   Offset  Type
% 1     Description (up to 20 characters)   22      0x00    string[21]
% 2     Variable type (e.g. int2 or rl4)    2       0x16    TDataType
% 3     Variable units (up to 8 characters) 10      0x18    string[9]
% 4     Byte offset of variable             2       0x22    short

FVI(NFV+1,1) = struct;
for ii = 1:NFV+1
    FVI(ii).Desc = fread(fid,[1,22],'uint8=>char');
    FVI(ii).FVType = fread(fid,1,'ushort=>double');
    FVI(ii).FVUnit = fread(fid,[1,10],'uint8=>char');
    FVI(ii).ByteOffset = fread(fid,1,'ushort=>double');
end

% Data Section Variable Information:
% The format is exactly the same as for the file variable information, except that
% offset to the variable is from the start of the data section variable area.

DSVI(NDSV+1,1) = struct;
for ii = 1:NDSV+1
    DSVI(ii).Desc = fread(fid,[1,22],'uint8=>char');
    DSVI(ii).DSType = fread(fid,1,'ushort=>double');
    DSVI(ii).DSUnit = fread(fid,[1,10],'uint8=>char');
    DSVI(ii).ByteOffset = fread(fid,1,'ushort=>double');
end

% File Variables:
% These are stored as a continuous list and are accessed using the offsets stored in the file
% variable information frame.
FVbo = diff([FVI(:).ByteOffset].'); % File Variable byte offset
for ii = 1:NFV
    switch FVI(ii).FVType
        case 0, conv = 'int8=>double'; curfvbo = FVbo(ii);
        case 1, conv = 'uint8=>char';  curfvbo = [1,FVbo(ii)];
        case 2, conv = 'short=>double'; curfvbo = FVbo(ii)/2;
        case 3, conv = 'short=>char';  curfvbo = [1,FVbo(ii)/2];
        case 4, conv = 'long=>double'; curfvbo = FVbo(ii)/4;
        case 5, conv = 'single';  curfvbo = FVbo(ii)/4;
        case 6, conv = 'double';  curfvbo = FVbo(ii)/8;
        case 7, conv = 'uint8=>char';  curfvbo = [1,FVbo(ii)];
    end
    FVI(ii).FV = fread(fid,curfvbo,conv);
end

out.FVI = FVI;

if ftell(fid)~=GH.FHSz
    disp('Mismatch between described header size and the file! Check the file integrity!');
end

% End File Header

%% Data Sections

% Each data section has two parts: the data section header and
% the channels of data, referred to as the data section data.

% The data sections are linked together by a list of pointers which point
% to the previous section.

DSbo = diff([DSVI(:).ByteOffset].'); % Data Section byte offset

for ds = NDS:-1:1
    
    % Data Section Header
    
    % General data section header
    %
    % Channel 0 information
    % Channel 1 information
    % ...
    % Channel c information
    %
    % Data section variable 0
    % Data section variable 1
    % ...
    % Data section variable d
    
    if ds==NDS
        fseek(fid,GH.LastDSHOffset,'bof');
    else
        % fseek(fid,DS(ds+1).GDSH.PrevDSHPtr,'bof');
        fseek(fid,GDSH.PrevDSHPtr(ds+1),'bof');
    end
    
    
    % General Data Section Header
    %       Description                                     Bytes   Offset	Type
    % 1     Pointer to previous data section header         4       0x00    long
    % 2     Pointer to start of channel data for this DS	4       0x04    long
    % 3     Size of this channel data area                  4       0x08    long
    % 4     Flags for marking data sections                 2       0x0c    TSFlags
    % 5     Reserved space for future expansion             16      0x0e    Zero filled
    GDSH.PrevDSHPtr(ds,1) = fread(fid,1,'ulong=>double');
    GDSH.StartChDataPtr(ds,1) = fread(fid,1,'ulong=>double');
    GDSH.ChDataSz(ds,1) = fread(fid,1,'ulong=>double');
    GDSH.DSFlags(ds,1) = fread(fid,1,'ushort=>double');
    fseek(fid,16,'cof'); % skipping 16 bytes of reserved space
    
    % Channel Information across data section
    % The channel parameters which can change between data sections.
    %   Description                                     Bytes   Offset  Type
    % 1	Offset in data section to first byte            4       0x00    long
    % 2 Data points (not bytes)                         4       0x04    long
    % 3 Y scale (integer/word data)                     4       0x08    float
    % 4 Y offset (integer/word data)                    4       0x0c    float
    % 5 X increment (equalSpaced and subsidiary data)	4       0x10    float
    % 6 X offset (equalSpaced and subsidiary data)      4       0x14    float
    for ch = 1:NCh
        DSChI(ch,1).DSoffset(ds,1) = fread(fid,1,'ulong=>double');
        DSChI(ch,1).DataPoints(ds,1) = fread(fid,1,'ulong=>double');
        DSChI(ch,1).Yscale(ds,1) = fread(fid,1,'single=>double');
        DSChI(ch,1).Yoffset(ds,1) = fread(fid,1,'single=>double');
        DSChI(ch,1).Xinc(ds,1) = fread(fid,1,'single=>double');
        DSChI(ch,1).Xoffset(ds,1) = fread(fid,1,'single=>double');
    end
    
    % Data Section variable across data section
    for dsv = 1:NDSV
        switch DSVI(dsv).DSType
            case 0, DSDSVI(dsv,1).DSV(ds,1) = fread(fid,DSbo(dsv),'int8=>double');
            case 1, DSDSVI(dsv,1).DSV{ds,1} = fread(fid,[1,DSbo(dsv)],'uint8=>char');
            case 2, DSDSVI(dsv,1).DSV(ds,1) = fread(fid,DSbo(dsv)/2,'short=>double');
            case 3, DSDSVI(dsv,1).DSV{ds,1} = fread(fid,[1,DSbo(dsv)/2],'ushort=>char');
            case 4, DSDSVI(dsv,1).DSV(ds,1) = fread(fid,DSbo(dsv)/4,'long=>double');
            case 5, DSDSVI(dsv,1).DSV(ds,1) = fread(fid,DSbo(dsv)/4,'single');
            case 6, DSDSVI(dsv,1).DSV(ds,1) = fread(fid,DSbo(dsv)/8,'double');
            case 7, DSDSVI(dsv,1).DSV{ds,1} = fread(fid,[1,DSbo(dsv)],'uint8=>char');
        end
    end
end

% If the headers are exact same across data segments, reduce to one value
% before appending the values to the structures that are part of the header

for ch = 1:NCh
    ChI(ch).TotalDataPoints = sum(DSChI(ch).DataPoints);
    if length(unique(DSChI(ch).DSoffset))==1
        ChI(ch).DSoffset = DSChI(ch).DSoffset(1);
    else
        ChI(ch).DSoffset = DSChI(ch).DSoffset;
    end
    if length(unique(DSChI(ch).DataPoints))==1
        ChI(ch).DSDataPoints = DSChI(ch).DataPoints(1);
    else
        ChI(ch).DSDataPoints = DSChI(ch).DataPoints;
    end
    if length(unique(DSChI(ch).Yscale))==1
        ChI(ch).DSYscale = DSChI(ch).Yscale(1);
    else
        ChI(ch).DSYscale = DSChI(ch).Yscale;
    end
    if length(unique(DSChI(ch).Yoffset))==1
        ChI(ch).DSYoffset = DSChI(ch).Yoffset(1);
    else
        ChI(ch).DSYoffset = DSChI(ch).Yoffset;
    end
    if length(unique(DSChI(ch).Xinc))==1
        ChI(ch).DSXinc = DSChI(ch).Xinc(1);
    else
        ChI(ch).DSXinc = DSChI(ch).Xinc;
    end
    if length(unique(DSChI(ch).Xoffset))==1
        ChI(ch).DSXoffset = DSChI(ch).Xoffset(1);
    else
        ChI(ch).DSXoffset = DSChI(ch).Xoffset;
    end
end

for dsv = 1:NDSV
    if length(unique(DSDSVI(dsv).DSV))==1
        DSVI(dsv).DSV = DSDSVI(dsv).DSV(1);
    else
        DSVI(dsv).DSV = DSDSVI(dsv).DSV;
    end
end

%% Fill up ChData and Chts - channel data and timestamps

% find data points per data segment per channel
emptych = [ChI.TotalDataPoints].'==0; fillch=cumsum(~emptych);
utotaldp = unique([ChI(~emptych).TotalDataPoints]);
sametotaldp = length(utotaldp)==1; 
% sametotaldp = ~sametotaldp; % debug line to check if cell format works...
if sametotaldp % all channels have the same numbers of total data points
    ChData = nan(utotaldp,fillch(end));
    samedsdp=size(unique([DSChI(~emptych).DataPoints].','rows'),1)==1;
    if samedsdp
        Chts = nan(utotaldp,1);
    else
        Chts = ChData;
    end
else
    ChData = cell(NDS,NCh);
    Chts = ChData;
end



for ch = 1:NCh
    sumdp=0;
    for ds = 1:NDS
        % Data Section Data / Channel data
        fseek(fid,GDSH.StartChDataPtr(ds)+DSChI(ch).DSoffset(ds),'bof');
        switch ChI(ch).DataType
            case 0, conv = 'int8=>double'; useb = 1;
            case 1, conv = 'uint8=>char'; useb = 1;
            case 2, conv = 'short=>double'; useb = 2;
            case 3, conv = 'ushort=>char'; useb = 2;
            case 4, conv = 'long=>double'; useb = 4;
            case 5, conv = 'single'; useb = 4;
            case 6, conv = 'double'; useb = 8;
            case 7, conv = 'uint8=>char'; useb = 1;
        end
        curdp = DSChI(ch).DataPoints(ds);
        curChData = fread(fid,curdp,conv,ChI(ch).ByteSpace-useb)*DSChI(ch).Yscale(ds)-DSChI(ch).Yoffset(ds);
        curChts = DSDSVI(4).DSV(ds)+(0:curdp-1)'*DSChI(ch).Xinc(ds)-DSChI(ch).Xoffset(ds);
        % Assuming DSDSVI(4) is start of data segment in seconds.
        if sametotaldp
            ChData(sumdp+(1:curdp),fillch(ch)) = curChData;
            if samedsdp
                if ch==1
                    Chts(sumdp+(1:curdp),1)=curChts;
                end
            else
                Chts(sumdp+(1:curdp),fillch(ch))=curChts;
            end
        else % USE THIS ONLY IF DO NOT WANT THE DATA STRUCTURE TO CHANGE FROM FILE TO FILE!!!
            ChData{ds,ch}=curChData;
            if curdp % necessary as transpose did something that prevented cell2mat function by assigning value to an empty cell
                Chts{ds,ch}=curChts;
            end
        end
        sumdp=sumdp+curdp;
    end
end
fclose(fid);


out.ChI = ChI;
out.DSVI = DSVI;
out.GDSH = GDSH;
out.ChData = ChData;
out.Chts = Chts;

%% plotting of ChData

if plotdata
    % Chts = repmat(Chts,1,size(ChData,2)); % debug line to check if matrix format works
    ufillch=unique(fillch);
    figure;
    if iscell(Chts) % cell
        for ch = 1:NCh
            for ds = 1:NDS
                if ~isempty(Chts{ds,ch})
                    plot(Chts{ds,ch},ChData{ds,ch}); hold on;
                end
            end
        end
    elseif isvector(Chts) % vector
        chtsds = reshape(Chts,ChI(ufillch(1)).DSDataPoints,[]);
        for ii=1:size(ChData,2)
            chdatads = reshape(ChData(:,ii),ChI(ufillch(1)).DSDataPoints,[]);
            plot(chtsds,chdatads); hold on;
        end
    else % matrix
        for ii = 1:size(Chts,2)
            curdsdp = ChI(ufillch(ii)).DSDataPoints;
            if length(curdsdp)>1
                csdsdp=cumsum(curdsdp);
                for jj=1:length(csdsdp)
                    if jj==1, bb=1; else, bb=csdsdp(jj-1)+1; end
                    ee=csdsdp(jj);
                    if ee>bb
                        plot(Chts(bb:ee,jj),ChData(bb:ee,jj)); hold on;
                    end
                end
            else
                plot(reshape(Chts(:,ii),ChI(ufillch(ii)).DSDataPoints,[]),reshape(ChData(:,ii),ChI(ufillch(ii)).DSDataPoints,[]));
                hold on;
            end
        end
    end
end