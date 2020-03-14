function diff = checkwUnitarity(w)
    d = length(w);
    testid = ncon({w,conj(w)}, ...
        {[1 2 3 4 5 -1],[ 1 2 3 4 5 -2]},[1 2 3 4 5]);
    id = eye(d);
    diff = reshape(testid,d,d) - id;
    diff = max(diff(:));
end