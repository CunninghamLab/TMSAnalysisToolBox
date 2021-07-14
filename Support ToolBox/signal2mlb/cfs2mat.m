% cfs2mat
% =======
%
% Script for translating CED CFS files into Matlab files
%
% Written by William Collins, May 2021
% Derived from a script written by Gilad Jacobson, May 2016
%
% Prerequisites: Install the matcfs64c library developed by J.G. Colebatch
% from the CED website. Extract the contents of the ZIP file into a folder
% and add that folder to the Matlab path.
%
% To run, simply type cfs2mat.Convert() in the Matlab prompt and select
% the files that you want to convert
%
% The output data from each file is saved as a series of .mat files, with
% the samename and location as the source file. Each .mat file contains a
% single structure called 'CfsFile' which contains the date and time of
% creation, any file comment and a list of file and frame (Data Section)
% variables. Finally, there is a list of channels, each of which holds the
% channel name, units, scaling and the data itself, divided into frames.

classdef cfs2mat
    methods (Static)
        function Convert()
            
            % Choose files to convert
            [fNames, cfsDir] = uigetfile(["*.cfs", "CFS files (*.cfs)"], "Choose a file to load","multiselect", "on");
            
            if numel(fNames) == 1 && all(fNames == 0)
                error("No valid file selected");
            end
            
            % Var for file open method
            READ = 0;
            
            % Vars for distinguishing between file and frame variables
            FILEVAR = 0;
            DSVAR = 1;
            
            % Vars for channel types
            EQUALSPACED = 0;
            MATRIX = 1;
            SUBSIDIARY = 2;
            
            if isstr(fNames)
                fNames = {fNames};
            end
            
            % Loop through all files
            for iFile = 1:length(fNames)
                
                fileHandle = matcfs64c('cfsOpenFile', [cfsDir fNames{iFile}],READ,0);
                
                if (fileHandle < 0)
                    error("File opening error: " + fileHandle);
                end
                
                % Initialise output structure
                CfsFile = struct;
                
                % Read file parameters
                [CfsFile.Time, CfsFile.Date, CfsFile.Comment] =...
                    matcfs64c('cfsGetGenInfo', fileHandle);
                
                [chans, nfVars, nDSVars, DSs] = ...
                    matcfs64c('cfsGetFileInfo', fileHandle);
                
                % Fetch file variables
                for ifVar = 1:nfVars
                    [Size, Type, CfsFile.FileVars(ifVar).Units, CfsFile.FileVars(ifVar).Desc] ...
                        = matcfs64c('cfsGetVarDesc', fileHandle, ifVar - 1, FILEVAR);
                    
                    CfsFile.FileVars(ifVar).Value = matcfs64c('cfsGetVarVal', fileHandle, ifVar - 1, FILEVAR, 0, Type, Size);
                end
                
                % loop over all frames (data sections)
                for iDS = 1:DSs
                    
                    % fetch frame variables
                    for iDSVar = 1:nDSVars
                        % Only get name and units once
                        if iDS == 1
                            [Size, Type, CfsFile.DSVars(iDSVar).Units, CfsFile.DSVars(iDSVar).Desc] ...
                                = matcfs64c('cfsGetVarDesc', fileHandle, iDSVar - 1, DSVAR);
                            
                            % Remove the string 'UserN' if a name exists
                            if iDSVar > 5 && iDSVar < 22
                                tempstr = strrep(CfsFile.DSVars(iDSVar).Desc, "User" + (iDSVar - 5), "");
                                if tempstr ~= ""
                                    CfsFile.DSVars(iDSVar).Desc = tempstr;
                                end
                            end
                        else
                            [Size, Type, ~, ~] ...
                                = matcfs64c('cfsGetVarDesc', fileHandle, iDSVar - 1, DSVAR);
                        end
                        
                        CfsFile.DSVars(iDSVar).Values{iDS} = matcfs64c('cfsGetVarVal', ...
                            fileHandle, iDSVar - 1, DSVAR, iDS, Type, Size);
                    end
                    
                    iChanOrig = 1;  % The index we access the source with
                    iChan = 1;      % The index we add to the list with
                    
                    % loop over all channels
                    while iChanOrig < chans && iChan < chans
                        
                        % read trace parameters
                        [CfsFile.Chans(iChan).Name, CfsFile.Chans(iChan).yUnits, CfsFile.Chans(iChan).xUnits, ...
                            dataType, dataKind, ~, other] = matcfs64c('cfsGetFileChan', fileHandle, iChanOrig - 1);
                        CfsFile.Chans(iChan).IsWaveData = dataKind == EQUALSPACED; % If not, matrix / marker
                        
                        % Get frame parameters
                        [~, points, CfsFile.Chans(iChan).DS(iDS).yScale, CfsFile.Chans(iChan).DS(iDS).yOffset, ...
                            CfsFile.Chans(iChan).DS(iDS).xScale, CfsFile.Chans(iChan).DS(iDS).xOffset] = ...
                            matcfs64c('cfsGetDSChan', fileHandle, iChanOrig - 1, iDS);
                        
                        % Read data
                        if dataKind == EQUALSPACED
                            if points > 0
                                CfsFile.Chans(iChan).DS(iDS).Data = matcfs64c('cfsGetChanData', ...
                                    fileHandle, iChanOrig - 1, iDS, 0, points, dataType);
                                
                                % Scale data
                                CfsFile.Chans(iChan).DS(iDS).Data = (CfsFile.Chans(iChan).DS(iDS).Data ...
                                    * CfsFile.Chans(iChan).DS(iDS).yScale) + CfsFile.Chans(iChan).DS(iDS).yOffset;
                            else
                                CfsFile.Chans(iChan).DS(iDS).Data = [];
                            end
                        elseif dataKind == MATRIX
                            if CfsFile.Chans(iChan).Name(1:6) == "Marker"
                                if other ~= iChanOrig
                                    error("Marker channel " + iChanOrig + " data not sequential (" + other + ")")
                                end
                                
                                if points > 0
                                    % Get marker times
                                    tempTimes = matcfs64c('cfsGetChanData', ...
                                        fileHandle, iChanOrig - 1, iDS, 0, points, dataType);
                                    % Scale marker times
                                    tempTimes = tempTimes ...
                                        * CfsFile.Chans(iChan).DS(iDS).yScale + CfsFile.Chans(iChan).DS(iDS).yOffset;
                                else
                                    tempTimes = [];
                                end
                                
                                % Get info and real name for attached channel
                                [CfsFile.Chans(iChan).Name, ~, ~, dataType, dataKind, ~, ~] ...
                                    = matcfs64c('cfsGetFileChan', fileHandle, other);
                                if dataKind ~= MATRIX
                                    error("Marker channel " + iChanOrig + " linked to non-matrix channel " + other);
                                end
                                [~, points, ~, ~, ~, ~] = matcfs64c('cfsGetDSChan', fileHandle, other, iDS);
                                if (points == 0 && ~isempty(tempTimes)) ...
                                        || (points > 0 && (isempty(tempTimes) || length(tempTimes) ~= points))
                                    error("Marker channel " + other + " data count (" + points ...
                                        + ") doesn't match number of timestampes (" + length(tempTimes) + ")");
                                end
                                if points > 0
                                    % Combine marker values with times
                                    tempVals = matcfs64c('cfsGetChanData', ...
                                        fileHandle, iChanOrig - 1, iDS, 0, points, dataType);
                                    CfsFile.Chans(iChan).DS(iDS).Data = [tempTimes,tempVals];
                                end
                                iChanOrig = iChanOrig + 1;
                            else
                                nDim = 1;
                                if points > 0
                                    % Get this channel's data
                                    CfsFile.Chans(iChan).DS(iDS).Data(nDim,:) = matcfs64c('cfsGetChanData', ...
                                        fileHandle, iChanOrig - 1, iDS, 0, points, dataType);
                                else
                                    CfsFile.Chans(iChan).DS(iDS).Data = {};
                                end
                                % Combine pure matrix data
                                while other ~= 0
                                    nDim = nDim + 1;
                                    if other ~= iChanOrig
                                        error("Marker channel " + iChanOrig + " data not sequential (" + other + ")")
                                    end
                                    % Get info for next channel
                                    [~, ~, ~, dataType, dataKind, ~, ~] = matcfs64c('cfsGetFileChan', fileHandle, other);
                                    if dataKind ~= MATRIX
                                        error("Matrix channel " + iChanOrig + " linked to non-matrix channel " + other);
                                    end
                                    [~, points, ~, ~, ~, ~] = matcfs64c('cfsGetDSChan', fileHandle, other, iDS);
                                    if (points == 0 && ~isempty(CfsFile.Chans(iChan).DS(iDS).Data())) ||...
                                            (points ~= 0 && points ~= length(CfsFile.Chans(iChan).DS(iDS).Data(1)))
                                        error("Matrix channel " + other + " data count (" + points ...
                                            + ") doesn't match that of first matrix channel (" ...
                                            + length(CfsFile.Chans(iChan).DS(iDS).Data(:)) + ")");
                                    end
                                    if points > 0
                                        % Add data points to matrix
                                        CfsFile.Chans(iChan).DS(iDS).Data(nDim,:) = matcfs64c('cfsGetChanData', ...
                                            fileHandle, other, iDS, 0, points, dataType);
                                    end
                                    [~, ~, ~, ~, ~, ~, other] = matcfs64c('cfsGetFileChan', fileHandle, other);
                                    iChanOrig = iChanOrig + 1;
                                end
                            end
                        elseif dataKind == SUBSIDIARY
                            % TODO
                        end
                        
                        iChan = iChan + 1;
                        iChanOrig = iChanOrig + 1;
                    end
                end
                
                %%%%%%%%D Structure
                D = struct; % initialise Matlab output structure
                
                % read file parameters
                [D.param.fTime,D.param.fDate,D.param.fComment] = matcfs64c('cfsGetGenInfo',fileHandle);
                [D.param.channels,fileVars,DSVars,D.param.dataSections] = matcfs64c('cfsGetFileInfo',fileHandle);
                
                % loop over all trials (data sections)
                
                for dsCount = 1:D.param.dataSections
                    
                    flagSet = matcfs64c('cfsDSFlags', fileHandle, dsCount, READ);  % setit = 0 to read
                    dSbyteSize = matcfs64c('cfsGetDSSize',fileHandle,dsCount);
                    
                    % loop over all channels
                    
                    for chCount = 1:D.param.channels
                        
                        % read trace parameters
                        [D.param.tOffset(chCount),points, ...
                            D.param.yScale(chCount), ...
                            D.param.yOffset(chCount), ...
                            D.param.xScale(chCount), ...
                            D.param.xOffset(chCount)] = ...
                            matcfs64c('cfsGetDSChan',fileHandle,chCount-1,dsCount);
                        [D.param.channelName{chCount}, ...
                            D.param.yUnits{chCount},D.param.xUnits{chCount}, ...
                            dataType,dataKind,spacing,other] = ...
                            matcfs64c('cfsGetFileChan',fileHandle,chCount-1);
                        
                        % zero corresponds to EQUALSPACED, or normal adc data as
                        % opposed to matrix data, which in Signal usually designates
                        % markers, or subsidiary data
                        if (dataKind == 0)
                            % read actual data
                            D.data(:,dsCount,chCount) = ...
                                matcfs64c('cfsGetChanData',fileHandle, ...
                                chCount-1,dsCount,0,points,dataType);
                            D.data(:,dsCount,chCount) = ...
                                (D.data(:,dsCount,chCount) * D.param.yScale(chCount)) + ...
                                D.param.yOffset(chCount);
                        end
                    end % for chCount
                end  % for dsCount
                
                
                ret = matcfs64c('cfsCloseFile',fileHandle); % close the CFS file
                
                outName = [cfsDir fNames{iFile}(1:end-4)];
                
                D.data =  permute(D.data, [1 3 2]);
                
                CfsFile.D = D;
                %save(outName,'D') % save the Matlab structure
                
                
                % close the CFS file
                matcfs64c('cfsCloseFile', fileHandle);
                
                % save the Matlab structure
                save([cfsDir fNames{iFile}(1:end-4)], "CfsFile")
                
            end
        end
    end
end
