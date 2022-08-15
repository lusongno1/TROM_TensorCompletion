% Last modifed on 08/15/2022
% lusong@lsec.cc.ac.cn
clc
clear
close all
missingRates = [0.75 0.65 0.55];
Rcps = [12 24 48 96 150];
ERR = zeros(length(missingRates),length(Rcps));
ii = 0;
jj = 0;
options.maxIts = 50;
for missingRate = missingRates
    ii = ii+1;
    jj = 0;
    for Rcp = Rcps
        jj = jj+1;
        global Y_true Yd Wd Err rng_s Wdval Wdmiss
        rng_seed = 99;
        rng(rng_seed);
        Err = [];
        Y_true = [];
        Yd = [];
        Wd = [];
        Wdval = [];
        Wdmiss = [];
        %%
        load ../../data574.mat
        Phits = tensor(Phi);
        Phik = cp_als(Phits,Rcp);
        ndim = ndims(Phik);
        ncpt = ncomponents(Phik);
        sz  = size(Phi);
        Y_true = Phi;
        X = Phi;
        %%
        creat_missing;
        Wd = double(W);
        X = Y_true.*W;
        %checkRate = nnz(X)/prod(sz);
        %%
        options.factr = 1e7;
        tic
        [Y,~,output] = cp_wopt(X, W, Rcp, 'opt_options',options);
        cputime = toc;
        %% calculate accuracy
        Yd = double(Y);
        errW = cal_acc_avail(Y_true,Yd,Wd)
        err = cal_acc(Y_true,Yd)
        %%
        figure
        Err2 = Err;
        Err2(:,1) = sqrt(2*Err(:,1))./sqrt(sum((Wd.*Y_true).^2,'all'));
        plot(Err2(:,1),'Marker','.');
        hold on;
        plot(Err2(:,2),'Marker','.');
        xlabel('Iteraion');
        ylabel('Accuracy');
        title(['Missing:' num2str(missingRate) ',  ' 'Rank:' num2str(Rcp) ' CPU time of ' ...
            num2str(options.maxIts) ' iters:',num2str(cputime) 's']);
        [~,ind] = min(Err2(:,2));
        text_point(ind,Err2(ind,2))
        legend('$\hat \varepsilon$','$\varepsilon$','minimum','Interpreter','LaTex')
        drawnow;
        %%
        ERR(ii,jj) = Err2(ind,2);
    end
end
figure;
for i=1:size(ERR)
    plot(Rcps,ERR(i,:),'Marker','.')
    hold on;
end
for i=1:size(ERR)
    legend_cell{i} = ['missing=' num2str(missingRates(i))];
end
legend(legend_cell,'Location','best');
xlabel('Ranks');
ylabel('Accuracy');
title(['Accuracy by Rank: maxIts=' num2str(options.maxIts) ' ,rng\_seed=' num2str(rng_seed) ]);
