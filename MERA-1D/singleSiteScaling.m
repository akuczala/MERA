function [scaleDims, scaleOps] = singleSiteScaling(w,wdag,num)
    ds = size(w);
    din = ds(1);
    dout = ds(4);
    superOperator = ncon({w,wdag},{[1 -2 2 -4],[1 -1 2 -3]},[1 2]);
    
    superOperator = reshape(superOperator,din*din,dout*dout);
    %[evecs,evals] = eigs(superOperator);
    [evecs,D] = eig(superOperator);
    evals = diag(D);
    disp(evals(1:num));
    scaleDims = -log(evals(1:num))/log(3);
    scaleOps = evecs;
end