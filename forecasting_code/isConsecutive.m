function isConsFlag = isConsecutive(dates)

isConsFlag = all( days(diff(dates)) == 1  );



