% Used for xyz2enu and enu2xyz
%
% Seebany Datta-Barua
% 14 Nov 2003
% ----------------------------------------
% Support function
% ----------------------------------------

function Vout = quicktimes(T,Vin)

[m n] = size(T);
if m ~= 3 | n ~= 3
    error('xyz2enu: Rotation matrix must be 3 x 3');
end

[m n] = size(Vin);
if m ~= 3
    error('xyz2enu: Input *vector* must be 3 x n');
end

Vout(1,:) = T(1,1) * Vin(1,:) + T(1,2) * Vin(2,:) + T(1,3) * Vin(3,:);
Vout(2,:) = T(2,1) * Vin(1,:) + T(2,2) * Vin(2,:) + T(2,3) * Vin(3,:);
Vout(3,:) = T(3,1) * Vin(1,:) + T(3,2) * Vin(2,:) + T(3,3) * Vin(3,:);

return
