% function gENU = svxyz2genu(lat, lon, ht, vx, vy, vz)
% rotates a vector expressed in satellite body coordinates
% vx (in-track forward), vy (cross-track to the left when facing "forward"
% in the direction of motion), vz (locally upward); into geographic
% east-north-up coordinates.

function gENU = svxyz2genu(lat, lon, ht, vx, vy, vz)

%% The xyz to gENU rotation is needed for all rotations.
% Rotate the in-track, cross-track coordinates to the geographic ENU.
% At each point compute how far ENU the next point on the track is.
xyz = wgslla2xyz(lat, lon, ht);%840e3);
% Assuming time steps forward consecutively, use the next data row as
% the relative position.
enu = xyz2enu_new(xyz(2:end,:), lat(1:end-1), lon(1:end-1), ht(end));
% Azimuth angle, repeat the last one to keep number of rows the same.
theta = [atan2(enu(:,2), enu(:,1)); atan2(enu(end,2), enu(end,1))];
% Get rotation matrix.  Notice for the direction of our transform, we
% want to undo the angle so make the angle negative. R = [3 x 3 x n].
R = rot_new(-theta*180/pi, 3);
% Rotate the DMSP vectors into the ENU directions. dmsp_v = [n x 3].
for i = 1:size(R,3)
    gENU(i,:) = R(:,:,i)*[ vx(i); vy(i); vz(i)];
end

% gENU = gENU'; % Transpose to have dimension [n x 3].