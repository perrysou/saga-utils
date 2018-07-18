function [] = read_channel(doy, year, sep, folder_path, op_path, signal)

if strcmp(folder_path(end-7:end-1), 'ASTRArx') || strcmp(folder_path(end-16:end-10), 'ASTRArx')
    channelfilestruct = dir([folder_path, 'txt', sep, 'channel*.log']);
else
    channelfilestruct = dir([folder_path, 'txt', sep, 'channel*', year, '_', doy, '_*.log']);
end

if ~isempty(channelfilestruct)
    CHANNEL = [];
    for jj = 1:size(channelfilestruct, 1)
        channelfile = channelfilestruct(jj, :);
        channelfilename = dlmread([folder_path, 'txt', sep, channelfile.name]);
        st = datestr(gps2utc(channelfilename(1, 1:2)), 'HHMM-');
        se = datestr(gps2utc(channelfilename(end, 1:2)), 'HHMM UT');
        %         disp([channelfile.name,' actually has data for ',st,se]);
        CHANNEL = [CHANNEL; channelfilename];
    end
    
    CHANNEL = CHANNEL(CHANNEL(:, 3) <= 3640, :);
    CHANNEL = sortrows(CHANNEL, [3, 4, 5]);
    ORTW = CHANNEL(:, 3);
    ORTS = CHANNEL(:, 4) + CHANNEL(:, 5);
    CARRIER = CHANNEL(:, 7);
    PSEUDORANGE = CHANNEL(:, 8);
    ISVALID = CHANNEL(:, 10);
    CYCLESLIPQ = CHANNEL(:, 11);
    STYPE = CHANNEL(:, 13);
    PRN = CHANNEL(:, 14);
    
    CHANNELDATA = [ORTW, ORTS, CARRIER, PSEUDORANGE, ISVALID, CYCLESLIPQ, STYPE, PRN];
    %     unique(PRN,'stable')
else
    CHANNELDATA = [];
end

outfilename = strcat('prn_files_', signal, sep, 'channeldata.mat');
save([op_path, outfilename], 'CHANNELDATA');
end
