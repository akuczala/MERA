function diff = checkuUnitarity(u)
    d = length(u);
    testid = ncon({u,conj(u)}, ...
        {[1 2 3 4 -5 -6 -7 -8],[ 1 2 3 4 -1 -2 -3 -4]},[1 2 3 4]);
    id = eye(d^4);
    diff = reshape(testid,d^4,d^4) - id;
    diff = max(diff(:));
end