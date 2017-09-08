function [R_rytov,k_par,k_par_ind] = Lz(Fs,NFFT,V,Psi_v,az,ze,L,z)
    load lz.mat
    %convert L,z to km
    L = L;
    z = z;
    %speed of light in vaccum
    c = 299792458;
    %frequency ranges from 0 to half the sampling rate 100Hz
    ff = [linspace(0,Fs,NFFT/2+1)]';
%     ff = [0:Fs/NFFT:Fs]';
    k_par = 2*pi*ff/V;    
    %limit the wavenumber within [10^-3, 0.5*10^-2]
%     k_par_ind = find(k_par<=10^-1&k_par>=10^-3);
    k_par_ind = find(k_par>=0);
    k = 2*pi*(1575.42*10^6)/c;
    k_x = k_par*cos(Psi_v);
    k_y = k_par*sin(Psi_v);
    k_perp = k_par;
%     k_perp == k_par
%     k_perp = 0;
    
%     %az is measured from north in North East Down Coords
%     az= -az+90;
%     if az > 180
%     az = az-360;
%     end

    k_prime_sq = k_perp.^2 + (k_x*cos(az) + k_y*sin(az)).^2*tan(ze)^2; 
    k_prime_sq = k_par .^ 2 * (sin(Psi_v + az) .^ 2 * sec(ze) .^ 2 + cos(Psi_v + az) .^ 2);
    
    ratio = k_prime_sq./k;
    if L == 0
        S_minus = sin(ratio/2*z).^2;
        S_plus = cos(ratio/2*z).^2;
    else
        S_minus = 1 - 2./(L*ratio).*sin(L*ratio/2).*cos((z-L/2)*ratio);
        S_plus = 1 + 2./(L*ratio).*sin(L*ratio/2).*cos((z-L/2)*ratio);
    end
    R_rytov = S_minus./S_plus;
end

