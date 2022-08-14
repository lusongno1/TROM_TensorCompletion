function err = cal_acc(Y_true,Y)
    %err = norm(Y_true-Y,2)/norm(Y_true,2);
    err = sum((Y_true-Y).^2,'all')/sum(Y_true.^2,'all');
    err = sqrt(err);
end