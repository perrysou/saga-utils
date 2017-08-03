function [ tspan,dt ] = find_tspan(MAX,span)
%find_tspan: given max sigmaphi data, find the corresponding time interval
%for high rate data processing
%
%INPUT
%MAX: timestamped data containing maximi of sigmaphi data 
%example of MAX: [time max sp prn# receiver ID]
%span: user-defined time interval in seconds [s]
%
%OUTPUT
%tspan: [start time;end time] in the form of datenum. its length equals span
%dt length of the zoom-in data segment

% reftime = MAX(MAX(:,2)==max(MAX(:,2)),1);
reftime = MAX(:,1);
% reftime = sortrows(reftime);
tspan = [reftime(1)-1/24/3600*span/2;reftime(end)+1/24/3600*span/2];
if reftime(end)-reftime(1)>1/24/3600*span
    dt = diff(tspan)/2;
else   
    dt = diff(tspan)/1;
end

end

