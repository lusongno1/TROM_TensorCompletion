function err = cal_acc_avail(Y_true,Y,W)
    %err = norm(Y_true-Y,2)/norm(Y_true,2);
    %err = sum((W.*(Y_true-Y)).^2,'all')/sum((W.*Y_true).^2,'all');
    err = sum((W.*(Y_true-Y)).^2,'all');
    %err = sqrt(err);
    err = err*0.5;
end