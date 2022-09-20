global Y_true ERR W
epsilon = cal_acc(Y_true(:),Xm(:));
tilde_epsilon = cal_acc_avail_std(Y_true(:),Xm(:),double(W(:)));
errs = [epsilon,tilde_epsilon]
ERR(end+1,:) = errs;
