function [ ] = xcorrmatrix(s1,s2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tic
s1 = randn(3000,1);
s2 = randn(3000,1);
n = size(s1,1);
S_mega = [];
s2_mega = [];
corrv = zeros(2*n-1,1);
lag = (-n+1:n-1)';
for i = n:-1:1
%     S = triu(ones(n,n),i-1)-triu(ones(n,n),i);
    S = sparse(diag(ones(n-i+1,1),i-1));
    corrv(n-i+1,:) = s1'*S*s2;
    corrv(n+i-1,:) = s1'*S'*s2;
    i
%     S_mega = sparse(blkdiag(S_mega,S));
end

% s1_mega = repmat(s1',[1 2*n-1]);
% 
% for j = 1:2*n-1
%     s2_mega = sparse(blkdiag(s2_mega,s2));
% end
% 
% corrv_ = s1_mega*S_mega*s2_mega;

toc

corrv = corrv.';
subplot(311);
stem(lag,corrv);
% subplot(312);
% stem(lag,corrv_');
subplot(313);
stem(xcorr(s1,s2));
end

