global Y_true Wd Err Wdval Wdmiss
Y = ktensor(tt_cp_vec_to_fac(x, Y_true));
Y = double(Y);
errW = cal_acc_avail(Y_true,Y,Wd);
err = cal_acc(Y_true,Y);
Err(outerIter,1:2) = [errW,err];
%fprintf('\n mask_acc = %0.8e, acc = %.8e   ',errW, err);
if(exist('split_avail_data')==2&&(~isempty(Wdval)))
    errVal = cal_acc_avail(Y_true,Y,Wdval);
    Err(outerIter,3) = errVal;
end
if(exist('get_missing_inds')==2&&~isempty(Wdmiss))
    errMiss = cal_acc_avail(Y_true,Y,Wdmiss);
    Err(outerIter,4) = errMiss;
end