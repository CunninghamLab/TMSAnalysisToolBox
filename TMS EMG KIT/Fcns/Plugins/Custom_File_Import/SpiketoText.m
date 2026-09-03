%SpiketoText - Worked Example of CustomFileTemplate (Spike2 / CED text export)
%{
PURPOSE:
  Import function for Spike2 (CED) text exports - the "*_export.txt" files
  produced by exporting a .smr recording as text. Built from
  CustomFileTemplate.m. Steps specific to THIS format are labelled
  [DATA-SPECIFIC]; steps that apply to any format are labelled [GENERIC].

THIS FILE FORMAT (Spike2 text export):
  - Tab-separated, CRLF (\r\n) line endings, quoted string fields
  - Structure:
        "INFORMATION"
        <source .smr filename>
        "SEGMENT_METADATA"
            "Duration (s):"            <tab> 3667.116
            "Samples per channel:"     <tab> 18335579
        "SUMMARY"
            "1" "Waveform" "EMG 1"     <tab> 5000 <tab> "mV" ...
            "2" "Waveform" "EMG 2"     <tab> 5000 <tab> "mV" ...
            "3" "Waveform" "RawTorque" <tab> 5000 <tab> "Nm" ...
            "4" "Evt+"     "DigMark"   <tab> 150  <tab> ""   ...
        "CHANNEL" "1"  ... "START" <start> <interval>  then one value per line
        "CHANNEL" "2"  ...
        "CHANNEL" "3"  ...
        "CHANNEL" "n"  (Marker channel - no START line)
            <time> <tab> "A???" <tab> 66 <tab> 0 <tab> 0 <tab> 0

  - Waveform channels are a bare single column of numbers, one per line.
    Sample n is at absolute time n*interval (interval from the START line).
  - The recording is ONE continuous segment: START is 0.000000 and never
    resets, so there are no natural block breaks.  ->  numBlocks = 1.
  - Marker (Evt+) channels hold the stimulus codes as 4-character Spike2
    marker strings, e.g. "A???" where '?' is unused padding.
  - EMG channels are recorded in mV; torque is in Nm.

MARKER CHANNEL SELECTION [DATA-SPECIFIC, IMPORTANT]:
  A recording may contain SEVERAL Evt+ channels - e.g. Digitimer, Magstim,
  Keyboard and DigMark - and they are NOT interchangeable:
    - some are empty,
    - some carry timestamps with no text code at all (Magstim triggers),
    - Keyboard and DigMark often duplicate the same stimuli a few ms apart.
  Merging them would double-count every trial and add a blank "" comment
  type, so this function reads each marker channel SEPARATELY and then asks
  the user which ONE to use as the comment/trial source. Channels are listed
  with their marker count and code letters so the choice is informed. If only
  one usable channel exists it is taken automatically with no prompt.

  The prompt appears on EVERY import by default, because the right channel is
  a per-recording decision. It carries a "Remember this choice" tick box for
  when you are importing a batch and do not want to be asked each time:
      unticked (default) - answer applies to this file only, ask again next time
      ticked             - saved with setpref and reused silently from then on
  To go back to being asked every time:
      rmpref('SpiketoText','MarkerChannel')
  Set app.SpikeMarkerChannel = "DigMark" beforehand to skip the prompt
  entirely without saving anything (useful for scripted/batch imports).

PERFORMANCE NOTE [IMPORTANT]:
  These files are large - ~430 MB and ~55 million lines for a one-hour
  recording. Reading them line-by-line (fgetl / textscan with '%s') takes
  many minutes. This function instead seeks past each channel header and
  issues a SINGLE bulk textscan('%f', N) per waveform channel, using the
  exact sample count declared in SEGMENT_METADATA. That reduces the read to
  a few seconds. Expect roughly 8*3*numSamples bytes of memory for the
  waveform data (~440 MB for a 3-channel, 3667 s recording at 5 kHz).

OUTLIER EXCLUSION [DATA-SPECIFIC]:
  If a sibling file
      <app.Path>/results/All_stims_data_auto_summary_<basename>.csv
  exists, its Is_Outlier / Manual_Excluded columns are used to drop rejected
  stimuli at import. Its 'Segment' column is the 1-based CHRONOLOGICAL index
  of that stimulus within its StimType. If the file is absent, ALL markers
  are imported and a warning is issued.

Output Variables:  (required by the toolbox - do not change names)
Variable Name       DataType (row x column)             Description
numChannels         1x1 double                          number of data channels
numBlocks           1x1 double                          number of Blocks
Block               See below                           Block data and channel titles
AllSampleRate       1x1 double                          sample rate in Hz
EventTypeOptions    1xN double vector                   event types available (must include 1)
Data                variable                            record of what was loaded
CommentLocAll       See below                           locations of all comments
numTrialsC          1xnumCommentTxtOpts double vector   number of comments per text option
CommentTxtOptions   numCommentTxtOpts x 1 string vector unique comment text options

%=======================================================================================================================
Block variable
Block = struct()
Block.Data   = numChannels x numBlocks cell  (each cell is a 1 x numFrames double row vector)
Block.Titles = numChannels x 1 string vector (channel names)

%=======================================================================================================================
CommentLocAll variable
CommentLocAll = numComments x numCommentTextOptions cell array
    Each cell contains: [Block_number   Frame_within_Block]
    First column  = locations of the first comment text type
    Second column = locations of the second comment text type
%}


function [numChannels,numBlocks,Block,AllSampleRate,EventTypeOptions,Data,CommentLocAll,numTrialsC,CommentTxtOptions]=SpiketoText(app,i)

%% File Extension
%{
    [GENERIC] - Set this to match your data file type so the toolbox file
    picker shows only relevant files by default.

    [DATA-SPECIFIC] - Spike2 text exports are saved as .txt
%}
%FileExtension='*.txt'; %set to match your data file type


%% Event Type Options
%{
    [GENERIC] - Which event types does this format support?
        1 - Block    (Required)
        2 - Comments (Include if the file has labeled comment/trigger markers)
        3 - Events   (Include if the file has discrete event markers)

    [DATA-SPECIFIC] - The Spike2 export carries one or more Marker channels
    (DigMark, Keyboard, ...) holding the stimulus codes, so Comments are
    available. Those markers are treated as Comments rather than Events
    because they carry text labels (A, B, C, ...) that the toolbox uses to
    group trials by stimulus type. See "MARKER CHANNEL SELECTION" above.
%}
EventTypeOptions=[1,2]; %must contain 1 (Blocks)


%% Import Data
%{
    [GENERIC] - Load the raw file. app.Path and app.Filename{i} are set by
    the user in the toolbox UI.

    [DATA-SPECIFIC] - Parsed in two passes over a single open file handle:
      pass 1 (cheap, line-by-line) reads the INFORMATION / SEGMENT_METADATA /
              SUMMARY header, which declares every channel and the sample count
      pass 2 (bulk) seeks through each CHANNEL section and reads its data
    The header tells us how many samples to expect, so each waveform channel
    is read with one preallocated textscan call instead of a loop.
%}
File=fullfile(app.Path, app.Filename{i}); %Path and Filename selected by user

fid=fopen(File,'r');
if fid==-1
    error('SpiketoText:CannotOpen','SpiketoText: Cannot open file: %s', File);
end
closeFile=onCleanup(@() fclose(fid)); %guarantees the file is closed even on error

TAB=sprintf('\t');

% --- [DATA-SPECIFIC] Pass 1: read the header ---------------------------------
nSamplesDeclared=NaN;
chanNum={}; chanType={}; chanTitle={}; chanRate={}; chanUnits={};

while true
    pos=ftell(fid);
    ln=fgetl(fid);
    if ~ischar(ln)
        error('SpiketoText:NoChannels','SpiketoText: Reached end of file before any CHANNEL section: %s', File);
    end
    s=strtrim(ln);

    if beginsWith(s,'"Samples per channel:"')
        tk=strsplit(s,TAB);
        if numel(tk)>=2
            nSamplesDeclared=str2double(unquote(tk{2}));
        end

    elseif beginsWith(s,'"SUMMARY"')
        % Each following row describes one channel, until a blank line.
        while true
            p2=ftell(fid);
            l2=fgetl(fid);
            if ~ischar(l2), break; end
            s2=strtrim(l2);
            if isempty(s2), break; end
            if beginsWith(s2,'"CHANNEL"')
                fseek(fid,p2,'bof'); %no blank separator - hand it back
                break
            end
            tk=strsplit(s2,TAB);
            if numel(tk)<5, break; end
            k=numel(chanNum)+1;
            chanNum{k}   = str2double(unquote(tk{1})); %#ok<AGROW>
            chanType{k}  = unquote(tk{2});             %#ok<AGROW>
            chanTitle{k} = unquote(tk{3});             %#ok<AGROW>
            chanRate{k}  = str2double(unquote(tk{4})); %#ok<AGROW>
            chanUnits{k} = unquote(tk{5});             %#ok<AGROW>
        end

    elseif beginsWith(s,'"CHANNEL"')
        fseek(fid,pos,'bof'); %rewind - pass 2 starts here
        break
    end
end

if isempty(chanType)
    error('SpiketoText:NoSummary','SpiketoText: SUMMARY section not found in file: %s', File);
end

% [DATA-SPECIFIC NOTE] SUMMARY calls waveform channels "Waveform" and the
% marker channel "Evt+"; inside the CHANNEL section the marker channel calls
% itself "Marker". We rely on SUMMARY for titles, units and channel order.
isWave = strcmpi(chanType,'Waveform');
waveIdx = find(isWave);
numChannels = numel(waveIdx);
if numChannels==0
    error('SpiketoText:NoWaveforms','SpiketoText: No Waveform channels listed in SUMMARY: %s', File);
end

Data=struct();
Data.File=File;
Data.SourceFile=[];
Data.ChannelSummary=struct('Number',chanNum,'Type',chanType,'Title',chanTitle, ...
                           'Rate',chanRate,'Units',chanUnits);
Data.NumSamplesDeclared=nSamplesDeclared;
% [GENERIC NOTE] 'Data' is meant to be a record of what was loaded. We store
% the parsed header and marker list rather than all 55 million raw text lines,
% which would otherwise multiply memory use several times over.


%% Extract Blocks
%{
    [GENERIC] - This section must produce:
        numChannels, numBlocks, Block (.Data and .Titles), AllSampleRate

    [DATA-SPECIFIC] - The Spike2 export is one continuous recording with a
    single START at t=0 and no gaps or resets, so there is exactly one block
    spanning the whole session. Comments are then located by frame index
    within that single block.

    If you later export in segments (Spike2 writes one SEGMENT_METADATA /
    CHANNEL group per segment), you would loop this section per segment and
    make each segment a block.
%}

Block=struct();
Block.Titles=string(chanTitle(waveIdx))';   %numChannels x 1 string vector
Block.Data=cell(numChannels,1);             %numChannels x numBlocks cell
numBlocks=1;

AllSampleRate=NaN;
sampleInterval=NaN;
% [DATA-SPECIFIC] Markers are kept PER CHANNEL - see "MARKER CHANNEL
% SELECTION" in the header. Accumulating them into one flat list would merge
% Keyboard/DigMark duplicates and Magstim's blank codes into a single set.
markerTimesByChan=repmat({zeros(0,1)},1,numel(chanType));
markerCodesByChan=repmat({strings(0,1)},1,numel(chanType));
wc=0; %waveform channel counter

for c=1:numel(chanType)

    % Position the file just after the next "CHANNEL" header line. Doing this
    % explicitly (rather than counting lines) resynchronises after each bulk
    % read, which stops mid-line at the last number it consumed.
    if ~seekPastChannelHeader(fid)
        break %no further channel sections
    end

    % Skip this channel's sub-header. Waveform sections end with a START line
    % (which carries the sample interval); marker sections have no START, so
    % we stop as soon as a line looks like data and rewind onto it.
    %
    % [DATA-SPECIFIC] A marker channel can be completely EMPTY (e.g. an unused
    % Digitimer channel), in which case the next thing after its sub-header is
    % the following '"CHANNEL"' line. That must terminate the skip too,
    % otherwise the parser runs straight on and reads the NEXT channel's rows
    % under this channel's name, shifting every later channel by one.
    interval=NaN;
    while true
        pos=ftell(fid);
        ln=fgetl(fid);
        if ~ischar(ln), break; end
        s=strtrim(ln);
        if beginsWith(s,'"START"')
            tk=strsplit(s,TAB);
            if numel(tk)>=3
                interval=str2double(unquote(tk{3}));
            end
            break %data begins on the next line
        elseif beginsWith(s,'"CHANNEL"')
            fseek(fid,pos,'bof'); %this channel is empty - hand the header back
            break
        elseif ~isempty(s) && looksNumeric(s(1))
            fseek(fid,pos,'bof');
            break
        end
    end

    if isWave(c)
        % --- [DATA-SPECIFIC] Bulk read one waveform channel ---
        wc=wc+1;
        if isnan(nSamplesDeclared) || nSamplesDeclared<=0
            nWant=Inf; %fall back to "read until it stops parsing"
        else
            nWant=nSamplesDeclared;
        end
        % '\r' is listed as whitespace so CRLF line endings parse cleanly.
        raw=textscan(fid,'%f',nWant,'Whitespace',' \b\t\r\n');
        v=raw{1};
        if ~isnan(nSamplesDeclared) && numel(v)~=nSamplesDeclared
            warning('SpiketoText:SampleCount', ...
                'SpiketoText: Channel "%s" returned %d samples but the header declared %d.', ...
                chanTitle{c}, numel(v), nSamplesDeclared);
        end

        % --- [DATA-SPECIFIC] Convert to Volts where the unit is electrical ---
        % The toolbox convention is Volts. EMG is recorded in mV, so it is
        % scaled by 1e-3. Torque (Nm) has no Volt equivalent and passes
        % through unchanged - scaling is driven by the unit string in the
        % header, not by channel position, so a re-scaled export still works.
        Block.Data{wc,1}=scaleToVolts(v(:).', chanUnits{c});

        if isnan(sampleInterval) && ~isnan(interval) && interval>0
            sampleInterval=interval;
            AllSampleRate=round(1/interval);
        end

    else
        % --- [DATA-SPECIFIC] Read a marker channel line-by-line ---
        % Marker channels hold only a few hundred rows, so a simple loop is
        % fine here and avoids textscan's error-recovery bookkeeping.
        tThis=zeros(0,1);
        cThis=strings(0,1);
        while true
            pos=ftell(fid);
            ln=fgetl(fid);
            if ~ischar(ln), break; end
            s=strtrim(ln);
            if isempty(s), continue; end
            if ~looksNumeric(s(1))
                fseek(fid,pos,'bof'); %next CHANNEL header - hand it back
                break
            end
            tk=strsplit(s,TAB);
            t=str2double(tk{1});
            if isnan(t), continue; end
            tThis(end+1,1)=t; %#ok<AGROW>
            if numel(tk)>=2
                % Spike2 pads marker codes to 4 characters with '?'.
                cThis(end+1,1)=string(regexprep(unquote(tk{2}),'\?+$','')); %#ok<AGROW>
            else
                cThis(end+1,1)=""; %#ok<AGROW>
            end
        end
        markerTimesByChan{c}=tThis;
        markerCodesByChan{c}=cThis;
    end
end

if isnan(AllSampleRate)
    % Fall back to the rate column in SUMMARY if no START line was found.
    r=chanRate{waveIdx(1)};
    if ~isnan(r) && r>0
        AllSampleRate=round(r);
        sampleInterval=1/AllSampleRate;
    else
        error('SpiketoText:NoSampleRate','SpiketoText: Could not determine sample rate for file: %s', File);
    end
end

% Guard against a channel that failed to read, which would break the toolbox
% downstream with a confusing error.
for ch=1:numChannels
    if isempty(Block.Data{ch,1})
        error('SpiketoText:EmptyChannel','SpiketoText: Channel "%s" read as empty in file: %s', ...
            Block.Titles(ch), File);
    end
end
numFrames=numel(Block.Data{1,1});

Data.SampleRate=AllSampleRate;
Data.SampleInterval=sampleInterval;
Data.NumFrames=numFrames;


%% Choose the Marker Channel
%{
    [DATA-SPECIFIC] - Pick ONE marker channel to act as the comment source.
    A channel is a candidate only if it actually contains markers; channels
    whose codes are all blank (e.g. bare Magstim triggers) are listed but
    flagged, because they cannot be used to group trials by stimulus type.

    Resolution order:
      1. app.SpikeMarkerChannel, if the caller set it (scripted imports)
      2. a saved preference, but ONLY if the user ticked "Remember this
         choice" in a previous prompt - otherwise we ask every time
      3. automatic, if there is exactly one candidate
      4. a modal pop-up listing every marker channel with its counts/codes
%}
markerIdx=find(~isWave);
markerIdx=markerIdx(~cellfun(@isempty,markerTimesByChan(markerIdx)));

if isempty(markerIdx)
    warning('SpiketoText:NoMarkerChannel', ...
        'SpiketoText: No marker channel in %s contained any markers. Importing without comments.', File);
    markerTimes=zeros(0,1);
    markerCodes=strings(0,1);
    chosenTitle="";
else
    % Build one description line per candidate, e.g.
    %   Ch 7  DigMark    126 markers   7 codes: A, B, C, D, F, G, H
    nCand=numel(markerIdx);
    titles=strings(nCand,1);
    labels=cell(nCand,1);
    hasCodes=false(nCand,1);
    for k=1:nCand
        c=markerIdx(k);
        codes=markerCodesByChan{c};
        uc=unique(codes(strlength(codes)>0));
        hasCodes(k)=~isempty(uc);
        titles(k)=string(chanTitle{c});
        if hasCodes(k)
            codeStr=strjoin(cellstr(sort(uc(:))'),', ');
            labels{k}=sprintf('Ch %d  %-12s  %4d markers   %d codes: %s', ...
                chanNum{c}, chanTitle{c}, numel(codes), numel(uc), codeStr);
        else
            labels{k}=sprintf('Ch %d  %-12s  %4d markers   (no text codes - cannot group trials)', ...
                chanNum{c}, chanTitle{c}, numel(codes));
        end
    end

    % Default highlight: prefer DigMark, else the coded channel with the most
    % distinct codes, else the first candidate.
    defaultPick=find(strcmpi(titles,'DigMark') & hasCodes,1);
    if isempty(defaultPick)
        nUnique=zeros(nCand,1);
        for k=1:nCand
            cc=markerCodesByChan{markerIdx(k)};
            nUnique(k)=numel(unique(cc(strlength(cc)>0)));
        end
        [~,defaultPick]=max(nUnique);
    end

    pick=[];

    % (1) Caller-supplied override
    wanted=getOptionalProp(app,'SpikeMarkerChannel');
    if ~isempty(wanted)
        pick=find(strcmpi(titles,string(wanted)),1);
        if isempty(pick)
            warning('SpiketoText:UnknownMarkerChannel', ...
                'SpiketoText: app.SpikeMarkerChannel = "%s" is not a marker channel in %s. Falling back to the prompt.', ...
                string(wanted), File);
        end
    end

    % (2) Saved preference - only present if the user ticked "Remember this
    %     choice" at some point. Without it we deliberately ask every time.
    % NOTE: "" is a 1x1 string, so isempty() is false for it - test strlength.
    savedChannel=getMarkerPref();
    if isempty(pick) && strlength(savedChannel)>0
        pick=find(strcmpi(titles,savedChannel),1);
        if isempty(pick)
            warning('SpiketoText:SavedChannelMissing', ...
                ['SpiketoText: Saved marker channel "%s" is not present in %s. ' ...
                 'Asking again for this file.'], savedChannel, File);
        else
            fprintf(['SpiketoText: using saved marker channel "%s". ' ...
                     'Run rmpref(''SpiketoText'',''MarkerChannel'') to be asked again.\n'], savedChannel);
        end
    end

    % (3) Only one candidate - there is nothing to choose between
    if isempty(pick) && nCand==1
        pick=1;
    end

    % (4) Ask the user. The tick box is what turns a one-off answer into a
    %     saved preference; left unticked, the next import asks again.
    if isempty(pick)
        if strlength(savedChannel)>0
            preselect=find(strcmpi(titles,savedChannel),1);
        else
            preselect=[];
        end
        if isempty(preselect), preselect=defaultPick; end

        [pick,rememberChoice]=promptMarkerChannel(app,labels,preselect,app.Filename{i});
        if isempty(pick)
            error('SpiketoText:MarkerChannelCancelled', ...
                'SpiketoText: Import cancelled - no marker channel selected for %s.', File);
        end
        if rememberChoice
            setMarkerPref(titles(pick));
        end
    end

    chosenIdx=markerIdx(pick);
    chosenTitle=titles(pick);
    markerTimes=markerTimesByChan{chosenIdx};
    markerCodes=markerCodesByChan{chosenIdx};

    fprintf('SpiketoText: using marker channel "%s" (%d markers).\n', chosenTitle, numel(markerTimes));

    if ~hasCodes(pick)
        warning('SpiketoText:BlankMarkerCodes', ...
            'SpiketoText: Marker channel "%s" has no text codes, so all trials will share one blank comment type.', ...
            chosenTitle);
    end

    % [DATA-SPECIFIC] Flag codes that differ from another only by case, e.g. a
    % stray lowercase "b" from a keypress with Shift released. These create a
    % spurious extra stimulus type. Reported, not silently corrected.
    uc=unique(markerCodes(strlength(markerCodes)>0));
    [~,~,g]=unique(upper(uc));
    for gi=1:max(g)
        if sum(g==gi)>1
            warning('SpiketoText:MarkerCaseMismatch', ...
                ['SpiketoText: Marker channel "%s" contains codes differing only by case (%s). ' ...
                 'These will be treated as separate stimulus types.'], ...
                chosenTitle, strjoin(cellstr(uc(g==gi)'),', '));
        end
    end
end

Data.MarkerChannel=chosenTitle;
Data.MarkerChannelTitles=string(chanTitle(~isWave))';
Data.MarkerTimesAll=markerTimes;
Data.MarkerCodesAll=markerCodes;
Data.MarkerTimesByChannel=markerTimesByChan(~isWave);
Data.MarkerCodesByChannel=markerCodesByChan(~isWave);


%% Extract Comments
%{
    [GENERIC] - Produces:
        CommentLocAll     - cell array locating each comment within blocks/frames
        numTrialsC        - count of comments per comment text type
        CommentTxtOptions - the unique comment text labels

    [DATA-SPECIFIC] - Comments come from the single marker channel selected
    above (Data.MarkerChannel). Two steps:
      1. Optionally drop stimuli flagged as outliers by the analysis pipeline.
      2. Convert each surviving marker's absolute time to a frame index in
         block 1, and group by stimulus letter.
%}

keepMask=true(numel(markerTimes),1);

% --- [DATA-SPECIFIC] Apply outlier exclusions, if the summary file exists ---
% 'Segment' in that CSV is the 1-based chronological index of the stimulus
% within its StimType, which is exactly the ordering markerCodes is already in.
[~,baseName]=fileparts(app.Filename{i});
summaryFile=fullfile(app.Path,'results',['All_stims_data_auto_summary_' baseName '.csv']);

if exist(summaryFile,'file')==2
    try
        T=readtable(summaryFile);
        stimTypeCol = findColumn(T,'StimType');
        segmentCol  = findColumn(T,'Segment');
        outlierCol  = findColumn(T,'Is_Outlier');
        manualCol   = findColumn(T,'Manual_Excluded');

        if isempty(stimTypeCol) || isempty(segmentCol) || isempty(outlierCol)
            warning('SpiketoText:SummaryColumns', ...
                'SpiketoText: %s is missing expected columns; importing all markers.', summaryFile);
        elseif height(T)~=numel(markerTimes)
            warning('SpiketoText:SummaryMismatch', ...
                ['SpiketoText: %s has %d rows but the file contains %d markers. ' ...
                 'Skipping outlier exclusion.'], summaryFile, height(T), numel(markerTimes));
        else
            isOut = toLogicalCol(T.(outlierCol));
            if ~isempty(manualCol)
                isOut = isOut | toLogicalCol(T.(manualCol));
            end
            stimTypes = toStringCol(T.(stimTypeCol));
            segments  = toDoubleCol(T.(segmentCol));

            nDropped=0;
            for r=1:height(T)
                if ~isOut(r), continue; end
                thisCode = strtrim(stimTypes(r));
                ofType   = find(markerCodes==thisCode); %chronological order
                seg      = segments(r);
                if seg>=1 && seg<=numel(ofType)
                    keepMask(ofType(seg))=false;
                    nDropped=nDropped+1;
                end
            end
            fprintf('SpiketoText: excluded %d of %d markers flagged as outliers.\n', ...
                nDropped, numel(markerTimes));
        end
    catch ME
        warning('SpiketoText:SummaryUnreadable', ...
            'SpiketoText: Could not read %s (%s); importing all markers.', summaryFile, ME.message);
    end
else
    % [DATA-SPECIFIC NOTE] Not every recording has been through the analysis
    % pipeline (e.g. the right-limb export in this dataset has no results/
    % folder). Importing everything is the safe fallback.
    warning('SpiketoText:NoSummaryFile', ...
        ['SpiketoText: No outlier summary found at %s. ' ...
         'Importing all %d markers without exclusion.'], summaryFile, numel(markerTimes));
end

markerTimesKept=markerTimes(keepMask);
markerCodesKept=markerCodes(keepMask);
Data.MarkerKept=keepMask;

% --- [GENERIC] Build CommentTxtOptions / numTrialsC / CommentLocAll ---
if ~isempty(markerCodesKept)
    CommentTxtOptions=unique(markerCodesKept);      %unique() returns sorted
    CommentTxtOptions=sort(CommentTxtOptions(:));   %explicit: alphabetical, numCommentTxtOpts x 1
    nOpt=numel(CommentTxtOptions);

    counts=zeros(1,nOpt);
    for c=1:nOpt
        counts(c)=sum(markerCodesKept==CommentTxtOptions(c));
    end
    numTrialsC=counts;

    CommentLocAll=cell(max(counts),nOpt);
    for c=1:nOpt
        idx=find(markerCodesKept==CommentTxtOptions(c));
        for m=1:numel(idx)
            % Sample n sits at absolute time n*interval, so frame (1-based)
            % is round(t/interval)+1. Clamped so a marker at the very end of
            % the recording cannot index past the data.
            frame=round(markerTimesKept(idx(m))/sampleInterval)+1;
            frame=min(max(frame,1),numFrames);
            CommentLocAll{m,c}=[1, frame]; %[block, frame] - always block 1 here
        end
    end
end


%Check for Comments Data, if there is no code above, create empty variables
%[GENERIC] - Keep this safety check. It ensures the toolbox always receives
%            valid comment variables, even when your format has no comments.
if ~exist('CommentLocAll','var') || ~exist('numTrialsC','var') || ~exist('CommentTxtOptions','var')
    CommentLocAll=cell(1,1);
    numTrialsC=[];
    CommentTxtOptions="";
end

end %end function


%% Functions
%{
    [GENERIC] - Helper functions used above. These must go AFTER the main
    function's 'end' statement.
%}

function tf=beginsWith(s,prefix)
% Prefix test written with strncmp rather than startsWith for compatibility
% with older MATLAB releases.
tf = numel(s)>=numel(prefix) && strncmp(s,prefix,numel(prefix));
end


function tf=looksNumeric(ch)
% True if this character could start a number - used to detect where a
% channel's header ends and its data begins.
tf = (ch>='0' && ch<='9') || ch=='-' || ch=='+' || ch=='.';
end


function s=unquote(s)
% Strip the surrounding double quotes Spike2 puts around string fields.
s=strtrim(s);
if numel(s)>=2 && s(1)=='"' && s(end)=='"'
    s=s(2:end-1);
end
end


function ok=seekPastChannelHeader(fid)
% Advance to the next '"CHANNEL"' line and leave the file positioned just
% after it. Returns false at end of file. This resynchronises the parser
% after a bulk textscan, which stops mid-line on the last value it read.
ok=false;
while true
    ln=fgetl(fid);
    if ~ischar(ln), return; end
    if beginsWith(strtrim(ln),'"CHANNEL"')
        ok=true;
        return
    end
end
end


function val=getOptionalProp(obj,name)
% Read obj.(name) if it exists, for either a struct or a class/App Designer
% object. Returns [] when absent or empty, so callers can just test isempty.
val=[];
try
    if isstruct(obj)
        if isfield(obj,name), val=obj.(name); end
    elseif isobject(obj) || ishandle(obj)
        if isprop(obj,name), val=obj.(name); end
    end
catch
    val=[]; %a property that errors on read is treated as absent
end
% Treat "" and '' as "not set" too. Note isempty("") is FALSE in MATLAB - a
% zero-length string is still a 1x1 string scalar - so test strlength.
if (ischar(val) || isstring(val)) && all(strlength(string(val))==0)
    val=[];
end
if isempty(val), val=[]; end
end


function [pick,remember]=promptMarkerChannel(app,labels,defaultPick,fileName)
% Modal pop-up asking which marker channel to use as the comment source.
%   pick     - index into 'labels', or [] if the user cancelled
%   remember - true if the user ticked "Remember this choice", meaning the
%              answer should be saved and reused instead of asking again
%
% Uses a uifigure list dialog so the text stays monospaced and the columns in
% each label line up. Falls back to listdlg if uifigure is unavailable (for
% example under -nodisplay), and to the default pick if neither can run. The
% fallback path cannot offer the tick box, so it never sets 'remember'.

pick=[];
remember=false;
n=numel(labels);

try
    fig=uifigure('Name','Select marker channel', ...
                 'Position',[100 100 640 290], ...
                 'Resize','off', ...
                 'WindowStyle','modal');
    cleanupFig=onCleanup(@() delete(fig(isvalid(fig)))); %#ok<NASGU> - fires on exit

    % Centre on the calling app's window if we can find it.
    parentFig=findAppFigure(app);
    if ~isempty(parentFig)
        p=parentFig.Position;
        fig.Position(1:2)=[p(1)+(p(3)-fig.Position(3))/2, p(2)+(p(4)-fig.Position(4))/2];
    end

    uilabel(fig,'Position',[20 248 600 22], ...
        'Text',sprintf('%s contains %d marker channels.',fileName,n), ...
        'FontWeight','bold');
    uilabel(fig,'Position',[20 226 600 22], ...
        'Text','Choose the ONE channel to use for trial comments (do not combine them):');

    lb=uilistbox(fig,'Position',[20 90 600 130], ...
        'Items',labels, ...
        'ItemsData',1:n, ...
        'Value',defaultPick, ...
        'Multiselect','off', ...
        'FontName','Consolas');

    % Unticked by default: the safe behaviour is to ask again next time,
    % because the right channel can differ between recordings.
    cb=uicheckbox(fig,'Position',[20 58 600 22], ...
        'Text','Remember this choice and stop asking (undo: rmpref(''SpiketoText'',''MarkerChannel''))', ...
        'Value',false);

    chosen=[]; chosenRemember=false;
    uibutton(fig,'Position',[430 20 90 28],'Text','OK', ...
        'ButtonPushedFcn',@(~,~) okPressed());
    uibutton(fig,'Position',[530 20 90 28],'Text','Cancel', ...
        'ButtonPushedFcn',@(~,~) delete(fig));
    % Double-clicking a row is the same as pressing OK. Guarded on its own
    % because DoubleClickedFcn only exists from R2022a.
    try %#ok<TRYNC>
        lb.DoubleClickedFcn=@(~,~) okPressed();
    end

    uiwait(fig);
    pick=chosen;
    remember=chosenRemember;
    return

catch ME
    if exist('fig','var') && ~isempty(fig) && isvalid(fig)
        delete(fig);
    end
    % uifigure is not available - fall back to the classic dialog.
    try
        [sel,ok]=listdlg('PromptString',{'Choose the marker channel to use', ...
                                         'for trial comments:',''}, ...
                         'ListString',labels, ...
                         'SelectionMode','single', ...
                         'InitialValue',defaultPick, ...
                         'ListSize',[560 160], ...
                         'Name','Select marker channel');
        if ok, pick=sel; end
        return
    catch
        warning('SpiketoText:NoDialog', ...
            ['SpiketoText: Could not show the marker-channel dialog (%s). ' ...
             'Using the default channel instead.'], ME.message);
        pick=defaultPick;
        return
    end
end

    function okPressed()
        chosen=lb.Value;
        chosenRemember=logical(cb.Value);
        delete(fig);
    end
end


function val=getMarkerPref()
% Read the saved marker-channel preference, or "" if the user has never
% ticked "Remember this choice". Wrapped in try/catch so a read-only or
% missing prefdir degrades to "ask every time" rather than erroring.
val="";
try
    if ispref('SpiketoText','MarkerChannel')
        val=string(getpref('SpiketoText','MarkerChannel'));
    end
catch
    val="";
end
if strlength(val)==0, val=""; end
end


function setMarkerPref(name)
% Save the marker-channel preference. Stored in MATLAB's prefdir, so it
% survives restarts and adds no files to the study folder.
try
    setpref('SpiketoText','MarkerChannel',char(name));
    fprintf(['SpiketoText: marker channel "%s" saved - you will not be asked again. ' ...
             'Undo with rmpref(''SpiketoText'',''MarkerChannel'').\n'], name);
catch ME
    warning('SpiketoText:PrefNotSaved', ...
        'SpiketoText: Could not save the marker channel preference (%s). It will ask again next time.', ...
        ME.message);
end
end


function fig=findAppFigure(app)
% Locate the calling app's window so the dialog can be centred on it.
fig=[];
try
    if isobject(app)
        props=properties(app);
        for k=1:numel(props)
            v=app.(props{k});
            if isa(v,'matlab.ui.Figure') && isvalid(v)
                fig=v;
                return
            end
        end
    end
catch
    fig=[];
end
end


function v=scaleToVolts(v,units)
% Convert a channel to Volts based on its unit string. Non-electrical units
% (Nm, deg, etc.) are left untouched - there is no meaningful conversion.
u=lower(strtrim(units));
u=strrep(u,char(181),'u'); %micro sign
u=strrep(u,char(956),'u'); %greek mu
switch u
    case {'v','volt','volts'}
        %already in Volts
    case {'mv','millivolt','millivolts'}
        v=v/1e3;
    case {'uv','microvolt','microvolts'}
        v=v/1e6;
    case {'nv','nanovolt','nanovolts'}
        v=v/1e9;
    otherwise
        %not an electrical unit - leave as recorded
end
end


function name=findColumn(T,wanted)
% Locate a table column by name, tolerating the underscore/case mangling
% readtable applies to headers such as 'Is_Outlier'.
name='';
vars=T.Properties.VariableNames;
target=lower(strrep(wanted,'_',''));
for k=1:numel(vars)
    if strcmpi(vars{k},wanted) || strcmp(lower(strrep(vars{k},'_','')),target)
        name=vars{k};
        return
    end
end
end


function out=toLogicalCol(col)
% Normalise a column that may arrive as logical, numeric, cellstr or string
% into a logical vector. readtable renders 'True'/'False' inconsistently
% across MATLAB releases, hence the belt-and-braces handling.
if islogical(col)
    out=col(:);
elseif isnumeric(col)
    out=col(:)~=0;
else
    s=lower(strtrim(string(col(:))));
    out=(s=="true" | s=="1" | s=="yes");
end
end


function out=toStringCol(col)
if isnumeric(col)
    out=string(col(:));
else
    out=string(col(:));
end
end


function out=toDoubleCol(col)
if isnumeric(col)
    out=double(col(:));
else
    out=str2double(string(col(:)));
end
end
