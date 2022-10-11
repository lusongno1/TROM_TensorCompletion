%% Last modified on 08/07/2022
%  lusong@lsec.cc.ac.cn
clc
clear
close all
data = load('../../data_missing');
Xk = data.X;
X = double(Xk);
X(X==0)=NaN;
rate = sum(isnan(X),'all')./numel(X);
RCP = 10;
Options(1) = 1e-3;
model1 = parafac(X,RCP,Options); % difficult data might require a lower value.
%[Factors,G] = tucker(X,[3 3 3 3 3 3],Options);% For tucker decomposition

Xkrec = ktensor(model1);
Xm = double(Xkrec);
SSX=sum(misssum(misssum(X.^2)),'all');
E = X-Xm;
W = isnan(X);
err = sqrt(sum((E(~W).^2),'all')/SSX);

viz(Xkrec)


% A = rand(4,5); %<-- First column is a_1, second is a_2.
% B = rand(3,5); %<-- Likewise for B.
% C = rand(2,5); %<-- Likewise for C.
% X = ktensor({A,B,C}); %<-- Create the ktensor.
% X = double(X);
% [Factors,G] = tucker(X,[1 1 1]);
