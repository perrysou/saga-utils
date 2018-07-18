function [prop, out_path, in_path] = ver_chk()
%check matlab version for proper configurations

global cases_folder home_dir mat_dir sep;
if verLessThan('matlab', '8.4.0')
    prop = 'clim';
else
    prop = 'limits';
end

if strcmp('GLNXA64', computer)
    %matlab on apollo linux server
    prop = 'limits';
    out_path = [home_dir, sep, mat_dir, sep];
    in_path = [home_dir, sep, mat_dir, sep];
    
elseif strcmp('PCWIN64', computer) && verLessThan('matlab', '8.7')
    %on the lab windows desktop
    out_path = [home_dir, sep, mat_dir, sep];
    in_path = [home_dir, sep, mat_dir, sep];
    
elseif strcmp('MACI64', computer)
    %my own mac laptop
    out_path = [home_dir, sep, mat_dir, sep];
    in_path = [home_dir, sep, mat_dir, sep];
    
else
    %on my own windows desktop
    out_path = 'D:\Dropbox\';
    in_path = 'D:\Dropbox\research\';
end

end
