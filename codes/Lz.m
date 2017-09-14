function [R_rytov,k_v,index_valid] = Lz(V,Psi_v,az,ze,thickness,height,f)
%Lz rytov approximation
% Given: az, azimuth angle; ze, zenith angle; thickness; height; 
% f, nyquist frequencies corresponding to given sampling rate 
% Return: R_rytov, Rytov ratio vector; k_par, corresponding wavenumber vector
% See also spectral

    %convert L,z to km
    thickness = thickness;
    height = height;
    
    %speed of light in vaccum
    c = 299792458;
    
    %frequency ranges from 0 to half the sampling rate 100Hz
%     ff = [linspace(0,Fs,NFFT/2+1)]';
%     ff = [0:Fs/NFFT:Fs]';
    k_v = 2 * pi * f / V;    
    
    %limit the wavenumber within [10^-3, 0.5*10^-2]
%     k_par_ind = find(k_par<=10^-1&k_par>=10^-3);
    index_valid = find(k_v >= 0);
    index_valid = find(f > 0.2);
    
    %incident wave number of the signal
    k = 2 * pi * 1575.42 * 10^6 / c;
%     az is measured from north in North East Down Coords
%     az= -az+90;
%     if az > 180
%     az = az-360;
%     end

    % Dr. Bust's way 
    k_x = k_v*cos(Psi_v);
    k_y = k_v*sin(Psi_v);
    k_squared = repmat(k_v .^2, 1, length(az)) + (k_x * cos(az) + k_y * sin(az)) .^2 .* tan(ze) .^2; 
    
    % Our way
    k_squared = k_v .^ 2 * (sin(Psi_v + az) .^ 2 .* sec(ze) .^ 2 + cos(Psi_v + az) .^ 2);
    
    %effective thickness and height along signal ray path
    thickness_effective = thickness * sec(ze);
    height_effective = height * sec(ze);
    
    ratio = k_squared ./ k;
    
    dummy = 2 ./ (ratio .* thickness_effective) ...
        .* sin(.5 * ratio .* thickness_effective) ...
        .* cos(ratio .* (height_effective - .5 * thickness_effective));
    
    if thickness == 0
        S_minus = sin(ratio .* height_effective) .^2;
        S_plus = cos(ratio .* height_effective) .^2;
    else
        S_minus = 1 - dummy;
        S_plus = 1 + dummy;
    end
    
    R_rytov = S_minus ./ S_plus;
end

