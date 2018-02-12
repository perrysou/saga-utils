function status = check_gps(ABS)

% Check the validity of Data in the struct ABS
% ABS.sec
% status of 0 implies data is fine, anything else is an error

status = 0;

if (isempty(ABS.sec) == 1) %wasn't a number
    status = 1;
    disp('Absolute Time information not numeric');
    return;
end;

if (ABS.sec < 0 )
    status = 1;
    disp('Absolute Time in seconds cannot be less than 0');
    return;
end;

return;