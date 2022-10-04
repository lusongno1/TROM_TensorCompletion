clc
clear
close all
data_path = '../../';
load([data_path 'data274D4222T100.mat']);
Phi = double(Phi_tucker);
sz = size(Phi);
S=ndSparse(Phi,sz);
KK = S-Phi
sum(S - Phi,'all')