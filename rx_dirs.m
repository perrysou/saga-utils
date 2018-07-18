function [rx_struct] = rx_dirs(cases_folder, year, doy)
sep = filesep;
if strcmp(doy, '342') && strcmp(year, '2013')
    comm = ['dir -d1 ', cases_folder, year, sep, doy, sep, '{grid*,ASTRA*}'];
elseif strcmp(doy, '233') && strcmp(year, '2017')
    comm = ['ls -d1 ', cases_folder, year, sep, doy, sep, 'grid*'];
else
    comm = ['dir -d1 ', cases_folder, year, sep, doy, sep, 'grid*'];
end
[err, dirs] = system(comm);
if err == 0
    struct = textscan(dirs, '%s', 'whitespace', '\n');
    rx_struct = char(struct{1});
    rx_struct = rx_struct(:, end-6:end);
    
    %special case for 2013/342
    %[grid108;grid163;grid162;ASTRArx;grid161]
    if strcmp(doy, '342') && strcmp(year, '2013')
        for i = size(rx_struct, 1):-1:1
            rx_struct1(size(rx_struct, 1)-i+1, :) = rx_struct(i, :);
        end
        rx_struct1 = [rx_struct1(end-1, :); rx_struct1(1:end-3, :); rx_struct1(end, :); rx_struct1(end-2, :)];
        rx_struct = rx_struct1;
    end
else
    rx_struct = [];
    disp(['No receiver structure for doy ', num2str(doy, '%03i')]);
end