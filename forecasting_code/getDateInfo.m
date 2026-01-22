function date_info = getDateInfo(fileNames, fileDate)

fileDateString = datestr(fileDate, 'YYYY-mm-DD');

fName = fileNames.dataFolder + fileDateString + "/" + fileNames.date_info + fileDateString + ".csv";
opts = detectImportOptions(fName);
opts = setvartype(opts, {'location', 'pathogen'}, 'categorical');
date_info = readtable(fName, opts);


