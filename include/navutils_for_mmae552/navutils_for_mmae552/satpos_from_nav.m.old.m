% Function [x, y, z] = satpos_from_nav(prn, stn,GPS_TIME)
% reads in the relevant nav file and calculates
% x, y, and z position of the sv.
%
% Seebany Datta-Barua
% 28 Jan 2005 Use idlfile as template, pass to satpos_from_alm for
% computation.

function [x,y,z] = satpos_from_nav(prn, stn,gps_time)

originaldir = cd;

% cd(dirname);
gpsconst;
x = zeros(size(gps_time,1),1);
y = zeros(size(gps_time,1),1);
z = zeros(size(gps_time,1),1);

GPS.week = gps_time(1,1);
GPS.sec = gps_time(1,2);
[absolute, utc, rinex] = convert_gps(GPS);
clear absolute utc

% Find correct nav file.
% dirname = ['c:\Documents and Settings\seebany\My Documents\' ...
%     'data\nav\'];
dirname = '/Users/seebany/Desktop/';
filelist = dir(dirname);
for i = 1:length(filelist)
    if findstr(filelist(i).name, num2str(rinex(1).day)) ...
            & findstr(filelist(i).name, ['.0' num2str(rinex(1).year - 2000) 'n'])
        filename = filelist(i).name;
    end
end
fid = fopen([dirname filename],'r');

% Skip past the header
line = '';
while ~feof(fid) & ( isempty(findstr(line, 'END OF HEADER')) | ...
        ~findstr(line,'END OF HEADER'))
    line = fgetl(fid);
end

% Loop through nav file, find eph for desired times.
while ~feof(fid)
    line = fgetl(fid);%, '%d %d %d %d %d %d %f %f %f %f');
    prn_compare = str2num(line(1:2));%A(1);
    utc.year = str2num(line(3:5))+2000;%A(2);
    utc.mon = str2num(line(6:8));%A(3);
    utc.day = str2num(line(9:11));%A(4);
    utc.hour = str2num(line(12:14));%A(5);
    utc.min = str2num(line(15:17));%A(6);
    utc.sec = str2num(line(18:22));%A(7);
    bias = str2num(line(23:41));%A(8);
    drift = str2num(line(42:60));%A(9);
    driftrate = str2num(line(61:end));%A(10);

    clear GPS
    %    tic
    % Convert_utc.m runs really slow.
    %    [GPS,absolute,rinex] = convert_utc(utc);
    [GPS_week, GPS_sec] = utc2gps([utc.year utc.mon utc.day ...
        utc.hour utc.min utc.sec]);
    GPS.week = GPS_week;
    GPS.sec = GPS_sec;
    clear absolute rinex
    if prn == prn_compare & isempty(find(x ~= 0))
        rows = 1:size(gps_time,1);
    else
        rows = find(gps_time(:,1)*sec_per_wk + gps_time(:,2) >= ...
            GPS.week*sec_per_wk + GPS.sec);
    end % if prn == prn_compare & isempty(find(x ~= 0))

    if prn == prn_compare & ~isempty(rows)
        line = fgetl(fid);
        iode = str2num(line(1:22));
        crs = str2num(line(23:41));
        d_n = str2num(line(42:60));
        M_0 = str2num(line(61:end));

        line = fgetl(fid);
        cuc = str2num(line(1:22));
        ecc = str2num(line(23:41));
        cus = str2num(line(42:60));
        sqrtA = str2num(line(61:end));

        line = fgetl(fid);
        toe_tow = str2num(line(1:22));
        cic = str2num(line(23:41));
        Omg_0 = str2num(line(42:60));
        cis = str2num(line(61:end));

        line = fgetl(fid);
        i0 = str2num(line(1:22));
        crc = str2num(line(23:41));
        w = str2num(line(42:60));
        Omgdot = str2num(line(61:end));

        line = fgetl(fid);
        idot = str2num(line(1:22));
        l2codes = str2num(line(23:41));
        toe_week = str2num(line(42:60));
        l2flag = str2num(line(61:end));

        line = fgetl(fid);
        acc = str2num(line(1:22));
        health = str2num(line(23:41));
        tgd = str2num(line(42:60));
        iodc = str2num(line(61:end));

        % Extra line of unneeded info.
        line = fgetl(fid);

        % Calculate position, adapted from idlfile satpos_from_nav.
        n0 = sqrt(mu_e) / (sqrtA^3);
        t = gps_time(rows,1)*sec_per_wk+gps_time(rows,2) - ...
            (toe_week*sec_per_wk + toe_tow);
        checkrows = find(t < -sec_per_wk/2);
        t(checkrows) = t(checkrows)+sec_per_wk;
        clear checkrows
        checkrows = find(t > sec_per_wk/2);
        t(checkrows) = t(checkrows)-sec_per_wk;
        clear checkrows

        n = n0 + d_n;
        M = M_0 + n*t;
        % min(n*t)
        % max(n*t)
        % toe_tow
        % min(M)
        % max(M)
        checkrows = find(M < -pi);
        M(checkrows) = M(checkrows) + 2*pi;
        clear checkrows
        checkrows = find(M > pi);
        M(checkrows) = M(checkrows) - 2*pi;

        Eold = M;
        E = M + ecc * sin(Eold);
        test = abs(E-Eold) > 1e-4;
        tic
        while (~isempty(find(test)))
            checktime = toc;
            if checktime > 1
                keyboard
            end
            Eold = E;
            E = M + ecc * sin(Eold);
            test = abs(E-Eold) > 1e-4;
        end

        checkrows = find(E < -pi);
        E(checkrows) = E(checkrows) + 2*pi;
        clear checkrows
        checkrows = find(E > pi);
        E(checkrows) = E(checkrows) - 2*pi;
        clear checkrows

        num = sqrt(1-ecc^2)*sin(E);
        den = (cos(E) - ecc);

        nu = atan2(num,den);
        phi = nu + w;
        phi = mod(phi, 2*pi);

        du = cus*sin(2*phi) + cuc*cos(2*phi);
        dr = crs*sin(2*phi) + crc*cos(2*phi);
        di = cis*sin(2*phi) + cic*cos(2*phi);

        u = phi + du;
        r = sqrtA^2 * (1-ecc*cos(E)) + dr;
        inc = i0 + di + idot*t;

        xprime = r.*cos(u);
        yprime = r.*sin(u);

        Omg = Omg_0 + (Omgdot - w_e)*t - w_e*toe_tow;
        checkrows = find(Omg < -pi);
        Omg(checkrows) = Omg(checkrows) + 2*pi;
        clear checkrows
        checkrows = find(Omg > pi);
        Omg(checkrows) = Omg(checkrows) - 2*pi;
        clear checkrows

        x(rows,1) = xprime.*cos(Omg) - yprime.*cos(inc).*sin(Omg);
        y(rows,1) = xprime.*sin(Omg) + yprime.*cos(inc).*cos(Omg);
        z(rows,1) = yprime.*sin(inc);

    else
        % Skip past all the lines we don't need.
        for i = 1:7
            line = fgetl(fid);
        end
    end %if ~isempty(rows)
end %while ~feof(fid)
fclose(fid);
cd(originaldir)