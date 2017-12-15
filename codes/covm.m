function [covmat] = covm(obs,hat)
%compute the covariance matrix for a column vector whose rows are
%measurements with known mean estimates
tic;
obs = [7,1,NaN,4;2,3,9,10;8,1,7,1];
hat = [3.25000000000000;6;4.25000000000000];
covmat = zeros(length(hat));
nnsbl = size(obs,2);
for i = 1:length(hat)
    nnani = length(find(isnan(obs(i,:))));
    for j = i:length(hat)       
        nnanj = length(find(isnan(obs(j,:))));
        tmp = ~isnan(obs(i,:))&~isnan(obs(j,:));
        covmat(i,j) = (obs(i,tmp) - hat(i)) * (obs(j,tmp) - hat(j))' ...
            / (nnsbl - max(nnani, nnanj) - 1); 
        covmat(j,i) = covmat(i,j);
    end
end
toc;
end