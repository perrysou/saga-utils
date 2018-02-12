% Function prn = svn2prn(svn) returns
% the prn corresponding to the svn.
% Taken from http://164.214.2.59/GandG/sathtml/satinfo.html
% Nat'l Imagery and Mapping Agency Satellite Geodesy Home Page
%
% Seebany Datta-Barua
% 17.09.01
% 9 June 2005 Updated with listing at
% http://earth-info.nima.mil/GandG/sathtml/satinfo.html

function prn = svn2prn(svn)

% column 1 is prn.  column 2 is svn.
table =    [ 1 32;
             2 13;
             3 33;
             4 34;
             5 35;
             6 36;
             7 37;
             8 38;
             9 39;
            10 40;
            11 46;
            13 43;
            14 41;
            15 15;
            16 56;
            17 17;
            18 54;
            19 19;
            20 51;
            21 45;%21;
            22 22;
            23 23;
            24 24;
            25 25;
            26 26;
            27 27;
            28 44;
            29 29;
            30 30;
            31 31];

prn = table(find(table(:,2) == svn),1);

return
