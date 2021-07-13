% Date: March 2016      
% Coded by: Khanh T. Ngoc Luong
% This function is applied in G-GNMF method, it create a matrix R from fea
% such that R(i,j) = fea(i,j) if j within p neighbors of i and vice versa.
% this implementation is similar to constructW of DengCai
function R = constructR(fea,options)

nSmp = size(fea,1);
nFea = size(fea,2);
G = zeros(nSmp*(options.p+1),3);
smpIdx = 1:nSmp;
relatedness = fea;
relatedness_old = relatedness;
nSmpNow = length(smpIdx);
dump = zeros(nSmpNow,options.p+1);
idx = dump;
for j = 1:options.p+1
[dump(:,j),idx(:,j)] = max(relatedness,[],2);
temp = (idx(:,j)-1)*nSmpNow+[1:nSmpNow]';
relatedness(temp) = 0;
end
%   G(1:nSmp*(options.p+1),1) = repmat(smpIdx',[options.p+1,1]);    
%   G(1:i*BlockSize*(options.p+1),2) = idx(:);
G(1:nSmp*(options.p+1),1) = repmat(smpIdx',[options.p+1,1]);
G(1:nSmp*(options.p+1),2) = idx(:);

G_cosine = G;
G_binary = G;
G_cosine(1:nSmp*(options.p+1),3) = dump(:);
G_binary(1:nSmp*(options.p+1),3) = 1;

 R = sparse(G_cosine(:,1),G_cosine(:,2),G_cosine(:,3),nSmp,nFea);
%R = sparse(G_binary(:,1),G_binary(:,2),G_binary(:,3),nSmp,nFea);
