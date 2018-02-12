function eph = readeph(file_name);
% eph = readeph(file_name);
%
% Script to read in ephemeris data from an ASCII formatted RINEX file.
% Input:
%    file_name - string containing the RINEX formatted ephemeris file 
%                to be read (nx1)
%
% Output:
%    gps_ephem - ephemeris matirx for all satellites (nx24), with columns of
%                 [prn,M0,delta_n,e,sqrt_a,long. of asc_node at GPS week epoch,
%                 i,perigee,ra_rate,i_rate,Cuc,Cus,Crc,Crs,Cic,Cis,Toe,IODE,
%                 GPS_week,Toc,Af0,Af1,Af2,perigee_rate] 
%                 Ephemeris parameters are from ICD-GPS-200 with the 
%                 exception of perigee_rate.
% 
% See propgeph for details on propagating the ephemeris matrix.

% Written by: Jimmy LaMance March 1998 
% Copyright (c) 1998 by Constell, Inc.

% Open the RINEX Navigation file to read
eph = [];
fid = fopen(file_name);
if (fid < 1)
    fprintf('Open failed on file %s.\n',file_name);
    return
end

fid_temp = fopen('temp.dat','w');

% Set the number of characters in the ephemeris data line
num_to_read = 34;

% Verify that a valid file ID was returned.
if fid < 1
    fprintf('Invalid file ID.\n')
    return
end % if fid < 1

% Start reading in the header of the RINEX data file
found_eoh = 0;    % flag indicating the end of the header
while found_eoh == 0
    s = fgetl(fid);
    if ~isempty(findstr(s,'END OF HEADER'))
        found_eoh = 1;
    end % if findstr(s,'END OF HEADER');
end % while found_eoh == 0

% Now that the file pointer is set to the first satellite ephemeris,
% read it and store the data for conversion to the standard Toolbox
% ephemeris format. Set the end of file flag (eof_flag) to zero.
eof_flag = 0;                 % end of file not found
eph_count = 1;

% Read in all of the data as characters.  Replace the Ds with Es.  Write
% the modified data out to the temp file.
[data, num_read] = fscanf(fid,'%c');
ID = findstr(data,'D');

if ~isempty(ID)
    data(ID) = ones(size(ID)) * 'E';
end % 

c = fprintf(fid_temp,'%c',data);
fclose(fid_temp);

fid_temp = fopen('temp.dat','r');

while eof_flag == 0
    % Read in the data a line at a time from the file.  We have to read 
    % this way because someRINEX writters include the spare fields and some don't.
    d = [];
    for i = 1:8
        if ~feof(fid_temp)
            line = fgetl(fid_temp);
            d = [d; sscanf(line,'%g')];
        end
    end
    
    data = d;
    
    % Make sure that the correct number of ephemeris elements are read (38)
    if size(data) < num_to_read
        eof_flag = 1; 
        break
    end % if num_read ~= num_to_read
    
    % Find the GPS Week from the UTC time tag (TOC)
    [GPS_week, GPS_sec] = utc2gps(data(2:7)');
    
    % This is now deemed a valid ephemeris, save in the nx24 GPS ephemeris
    % format...
    %   1      2       3       4        5                   6
    % [prn,   M0,   delta_n,   e,    sqrt_a,     long. of asc_node at GPS week epoch,
    %  7       8          9      10        11     12     13     14     15    16    17    18
    %  i,   perigee,  ra_rate,  i_rate,   Cuc,   Cus,   Crc,   Crs,   Cic,  Cis,  Toe,  IODE,
    %     19        20     21     22     23         24
    %  GPS_week,   Toc,   Af0,   Af1,   Af2,   perigee_rate]
    eph(eph_count,:) = [data(1)  data(14) data(13) data(16) data(18) ...
            data(21) data(23) data(25) data(26) data(27) ... 
            data(15) data(17) data(24) data(12) data(20) ...
            data(22) data(19) data(11) GPS_week GPS_sec ...
            data(8)  data(9)  data(10) 0]; 
    health(eph_count) = data(32);
    accuracy(eph_count) = data(31);
    eph_count = eph_count + 1;
end % while eof_flag == 0

% Close the data file
fclose(fid);  
fclose(fid_temp); 