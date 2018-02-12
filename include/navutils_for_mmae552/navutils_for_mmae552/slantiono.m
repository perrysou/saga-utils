% Function Is = slantiono(Iv, el) takes in values
% of Iv and elevation as found in truth___.dat files
% (produced by waas_to_su.m) and removes the obliquity
% factor.  Works for lists of Iv's and el's also.
%
% Seebany Datta-Barua
% 07.06.01

function Is = slantiono(Iv,el)

Ob = cos(asin(0.94797966*cos(el)));
Is = Iv./Ob;
% Make sure Iv has been converted from TEC to m via the TECU2l1m constant
return
