clc
clear
close all
format rat
A = rand(1,1000)
[N,D] = numden(sym(A))
for i=1:size(A,2)
    i
    scatter(i,sum(N(1:i))./sum(D(1:i)));
    hold on;
    drawnow;
end