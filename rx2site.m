function [sitenum_op] = rx2site(rcvr_op)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% load testdata.mat;
sitenum_op = cell(size(rcvr_op, 1), 1);
for i = 1:size(rcvr_op, 1)
    if strcmp(rcvr_op(i, :), 'grid108')
        sitenum_op{i, :} = 'IIT-1';
    elseif strcmp(rcvr_op(i, :), 'grid163')
        sitenum_op{i, :} = 'IIT-3';
    elseif strcmp(rcvr_op(i, :), 'grid162')
        sitenum_op{i, :} = 'IIT-11';
    elseif strcmp(rcvr_op(i, :), 'grid161')
        sitenum_op{i, :} = 'IIT-15';
    elseif strcmp(rcvr_op(i, :), 'grid160')
        sitenum_op{i, :} = 'IIT-9';
    elseif strcmp(rcvr_op(i, :), 'grid154')
        sitenum_op{i, :} = 'IIT-16';
    elseif strcmp(rcvr_op(i, :), 'ASTRArx')
        sitenum_op{i, :} = 'IIT-13';
    end
end

end
