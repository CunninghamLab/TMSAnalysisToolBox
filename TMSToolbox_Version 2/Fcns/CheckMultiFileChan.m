%Check if there are different channels for different files
%Inputs: app data
%Outputs: none, updates variables in function

function CheckMultiFileChan(app)

if length(app.Filename) > 1

    %check if numTitles are the same for each block
    if any(any(ismissing(app.Block.Titles))) %if there are missing values
        %There are a different number of channels for each file
        [BTitlesData,UniqueTitles] = MultiFileDiffChan(app);

        app.Block.Data=BTitlesData;
        app.Block.Titles=UniqueTitles;
        app.numChannels=length(UniqueTitles);

    else %there are the same number of channels for each block

        %check if all channels are the same
        if all(app.Block.Titles(:,1)==app.Block.Titles(:,2:end)) %if all the column titles are the same and in the same order
            %Don't change Block.Data
            app.Block.Titles=app.Block.Titles(:,1);
            app.numChannels=length(app.Block.Titles);
        else %see if they are just in a different order, but the channels are the same
            %reorder for easier comparison
            Reorder=sort(app.Block.Titles,1);
            Comp=Reorder(:,1)==Reorder(:,2:end); %compare all columns with column one
            if all(Comp) %if all the columns match
                %All the column titles are the same, but in a different order
                %Use order of first file (column)
                DiffAll=app.Block.Titles(:,1)==app.Block.Titles(:,2:end);
                BTitles=app.Block.Titles; BTitlesData=app.Block.Data;
                %Reorder Titles
                for i=1:length(BTitles(1,2:end)) %for each column after the first
                    Diff=DiffAll(:,i);
                    DiffLoc=find(Diff == 0);        %determine where there are differences
                    DiffLocData=struct();
                    for i3=1:length(DiffLoc)
                        DiffLocData(i3).Data=BTitlesData(DiffLoc(i3),i+1);
                        DiffLocData(i3).Title=BTitles(DiffLoc(i3),i+1);
                    end
                    for i2=1:length(DiffLoc)
                        %find title being switched
                        DesiredTitle=BTitles(DiffLoc(i2),1);
                        for i4=1:length(fieldnames(DiffLocData)) %match titles with data
                            if DiffLocData(i4).Title == DesiredTitle
                                DesiredData=DiffLocData(i4).Data;
                            end
                        end
                        BTitlesData(DiffLoc(i2),i+1)=DesiredData;
                        BTitles(DiffLoc(i2),i+1)=DesiredTitle;     %fill in based on column one

                    end
                end

                app.Block.Titles=BTitles(:,1);
                app.Block.Data=BTitlesData;
                app.numChannels=length(app.Block.Titles);

            else %there are the same number of titles, but they are not all the same
                [BTitlesData,UniqueTitles] = MultiFileDiffChan(app);

                app.Block.Data=BTitlesData;
                app.Block.Titles=UniqueTitles;
                app.numChannels=length(UniqueTitles);

            end %end if there are the same number of channels, but they are not in the same order
        end %end if there are the same number of channels and titles are in the same order
    end %end if there are different number of channels


else %only one file
    app.Block.Titles=app.Block.Titles(:,1);
    app.numChannels=length(app.Block.Titles);

end %end if there is more than one file, check numTitles
end %end function
