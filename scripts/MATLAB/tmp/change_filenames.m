function change_filenames(sn)
    if sn >= 7  && sn<=13 %% From S07, they have long-ITI data in the first session
        folder = fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/imaging_data_raw',sprintf('S%02d',sn));

        % Get a list of all files in the folder
        fileList = dir(fullfile(folder, '*'));

        % Loop through each file in the folder
        for i = 1:length(fileList)
            % Skip directories
            if fileList(i).isdir
                continue;
            end

            % Get the current file name
            oldName = fileList(i).name;

            % Check if the file name contains 'run_09'
            if contains(oldName, 'run_09')
                % Create the new file name
                newName = strrep(oldName, 'run_09', 'run_17');        
                % Get the full paths for the old and new file names
                oldFilePath = fullfile(folder, oldName);
                newFilePath = fullfile(folder, newName);

                % Rename the file
                movefile(oldFilePath, newFilePath);
            end
        end
        disp('File renaming complete.');
    end
        
  %  directory = fullfile('/Volumes/diedrichsen_data$/data/SeqSpatialSupp_fMRI/imaging_data_raw',sprintf('R%02d',sn));
    directory = fullfile('/srv/diedrichsen/data/SeqSpatialSupp_fMRI/imaging_data_raw',sprintf('R%02d',sn));
%     if ~exist(directory,"dir")
%         mkdir(directory)
%     end
%       
    cd(directory);
% Get a list of all files in the directory
    files = dir(directory);

    % Loop through each file
    for i = 1:length(files)
        % Get the file name
        oldName = files(i).name;

        % Check if the file name contains 'R'
        if contains(oldName, 'R')
            % Replace 'R' with 'S'
            newName = strrep(oldName, 'R', 'S');

            % Find the pattern 'run_' followed by two digits
            runNumberPattern = 'run_\d{2}';
            runNumberMatch = regexp(newName, runNumberPattern, 'match', 'once');

            if ~isempty(runNumberMatch)
                % Extract the two-digit number from the match
                oldRunNumber = regexp(runNumberMatch, '\d{2}', 'match', 'once');

                % Convert to number, add 8, and convert back to two-digit string
                newRunNumber = sprintf('%02d', str2double(oldRunNumber) + 8);

                % Replace the old run number with the new one in the match
                newRunNumberMatch = strrep(runNumberMatch, oldRunNumber, newRunNumber);

                % Replace the old run number match with the new one in the filename
                newName = strrep(newName, runNumberMatch, newRunNumberMatch);
            end

            % Construct full file paths
            oldFilePath = fullfile(directory, oldName);
            newFilePath = fullfile(directory, newName);

            % Rename the file
            movefile(oldFilePath, newFilePath);

            % Display the change
            fprintf('Renamed: %s -> %s\n', oldName, newName);
        end
    end
    % elseif sn >=7 & sn<=13
    %     directory = fullefile('/Volumes/diedrichsen_data$/data/SeqSpatialSupp_fMRI/imaging_data',sprintf('S%02d',sn));
    % 

    
