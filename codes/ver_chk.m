function [prop, out_path, in_path] = ver_chk
%check matlab version for proper configurations

if verLessThan('matlab', '8.4.0')
    prop = 'clim';
else
    prop = 'limits';
end

if strcmp('GLNXA64', computer)
    %matlab on apollo linux server
    prop = 'limits';
    out_path = '/data1/home/ysu27/Dropbox/';
    in_path = '/data1/home/ysu27/Dropbox/research/';
    
elseif strcmp('PCWIN64', computer) && verLessThan('matlab','8.7')
    %on the lab windows desktop
    out_path = 'C:\Users\Yang Su\Dropbox\';
    in_path = 'C:\Users\Yang Su\Dropbox\research\';
    
elseif strcmp('MACI64', computer)
    %my own mac laptop
    out_path = '/Users/yangsu/Dropbox/';
    in_path = '/Users/yangsu/Dropbox/research/';
    
else
    %on my own windows desktop
    out_path = 'D:\Dropbox\';
    in_path = 'D:\Dropbox\research\';
end

end


