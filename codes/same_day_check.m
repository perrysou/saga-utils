function [ ind ] = same_day_check( time )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
diff_day = find(diff(time)<0);
    if isempty(diff_day)~=1
        if length(diff_day)==2
        ind = diff_day(1)+1:diff_day(2);
        elseif (length(time)-diff_day)> abs(1-diff_day)
            ind = diff_day+1:length(time);
        else ind = 1:diff_day;  
        end
    else ind = 1:length(time);
    end
end

