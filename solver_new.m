function [eig_val_critical,eig_vec_critical,flag,max_eig,all_eig]=solver_new(A,B,num,p,sigma1)
opts = struct('disp', 0, 'maxit', 500, 'UseParallel', true);
[eig_vec, D] = eigs(A,B,p,sigma1,opts);

eig_val=diag(D);
all_eig=eig_val;
% Obtain location of critical eigen values
threshold_magnitude = 1e6;
eig_val_critical = eig_val(real(eig_val) > 0 & abs(eig_val) < threshold_magnitude);
[~, eig_critical_location] = ismember(eig_val_critical, eig_val);

%whether to return eigen vector or not, num==1->return else do not return
if num==1
    eig_vec_critical=eig_vec(:,eig_critical_location);
else
    eig_vec_critical='not asked';
end

if isempty(eig_val_critical)
    flag=0;
    max_eig=[];
else
    flag=1;
    max_eig=max(eig_val_critical);
end

end