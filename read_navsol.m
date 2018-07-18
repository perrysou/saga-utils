function [] = read_navsol(offsetyear, offsetdoy, in_path_offset, ...
    doy, year, sep, folder_path, op_path, signal)
%given time periods and receiver ID, read xyz coords of from navsol.log
%------------------
% doy = '342';
% year = '2013';
% folder_path = '/media/yang/Misc/Dropbox/navsol/grid108/';
% op_path = folder_path;
% sep = filesep;
%------------------
if strcmp(folder_path(end-7:end-1), 'ASTRArx') || strcmp(folder_path(end-16:end-10), 'ASTRArx')
    navfilestruct = dir([folder_path, 'txt', sep, 'navsol*.log']);
    navfilestruct0 = [];
else
    navfilestruct = dir([folder_path, 'txt', sep, 'navsol*', year, '_', doy, '_*.log']);
    navfilestruct0 = dir([in_path_offset, 'txt', sep, 'navsol*', offsetyear, '_', offsetdoy, '_*.log']);
end

if ~isempty(navfilestruct0)
    navfilename0 = dlmread([in_path_offset, 'txt', sep, navfilestruct0(end).name]);
else
    navfilename0 = [];
end

if ~isempty(navfilestruct)
    NAV = [];
    for jj = 1:size(navfilestruct, 1)
        navfile = navfilestruct(jj, :);
        navfilename = dlmread([folder_path, 'txt', sep, navfile.name]);
        st = datestr(gps2utc(navfilename(1, 1:2)), 'HHMM-');
        se = datestr(gps2utc(navfilename(end, 1:2)), 'HHMM UT');
        %         disp([navfile.name,' actually has data for ',st,se]);
        NAV = [NAV; navfilename];
    end
    
    NAV = [NAV; navfilename0];
    
    NAV = NAV(NAV(:, 1) <= 3640, :);
    NAV = sortrows(NAV, [1, 2, 3]);
    ORTW = NAV(:, 1);
    ORTS = NAV(:, 2) + NAV(:, 3);
    XYZECEF = NAV(:, 4:6);
    
    NAVDATA = [ORTW, ORTS, XYZECEF];
else
    NAVDATA = [];
end

outfilename = strcat('prn_files_', signal, sep, 'navdata.mat');
save([op_path, outfilename], 'NAVDATA');
end
