function fileDateCurrent = getRightFileDate(originDate, fileDates, fileNames, location_name, pathogen_name_full)

% Start with the first file date following the origin date
iFile = find(fileDates >= originDate, 1, 'first')-1;

% Keep moving forwards through the file dates until you get to one that has an origin date on or after the current origin (or you get the last file date)
foundRightFile = false;
while ~foundRightFile & iFile < length(fileDates)
    iFile = iFile + 1;

    % Get date information for this file date
    date_info = getDateInfo(fileNames, fileDates(iFile));
    
    % Get the origin date corresponding to this file date
    checkOrigin = date_info.origin_date(date_info.location == location_name & date_info.pathogen == pathogen_name_full);

    % Check if this origin date satisfies the criteria
    foundRightFile = checkOrigin >= originDate;
end
fileDateCurrent = fileDates(iFile);

