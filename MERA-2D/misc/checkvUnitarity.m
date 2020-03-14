function diff = checkvUnitarity(v)
    d = length(v);
    testid = ncon({v,conj(v)}, ...
        {[1 2 3 4 -3 -4],[ 1 2 3 4 -1 -2]},[1 2 3 4]);
    id = eye(d^2);
    diff = reshape(testid,d^2,d^2) - id;
    diff = max(diff(:));
end