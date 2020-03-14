function diff = checkuUnitarity(u)
    d = length(u);
    testid = ncon({u,conj(u)},{[1 2 -3 -4],[ 1 2 -1 -2]},[1 2]);
    id = eye(d*d);
    diff = reshape(testid,d^2,d^2) - id;
    diff = max(diff(:));
end