function err = cal_acc_avail_std(Y_true,Y,W)
    err = sum((W.*(Y_true-Y)).^2,'all')/sum((W.*Y_true).^2,'all');
    err = sqrt(err);
end