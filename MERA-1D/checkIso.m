function diff = checkIso(w,dimin,dimout)
    concatindices = [1:length(dimin)];
    out1 = -[1:length(dimout)];
    out2 = -length(dimout) - out1;
    legLinks = {horzcat(concatindices,out1),horzcat(concatindices,out2)};
    testid = ncon({w,conj(w)},legLinks,concatindices);
    id = eye(prod(dimout));
    diff = reshape(testid,prod(dimin),prod(dout)) - id;
    diff = max(diff(:));
    
end