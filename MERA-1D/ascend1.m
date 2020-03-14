function out = ascend1(w,op)
    ds = size(w);
    din = ds(1);
    dout = ds(4);
    out = ncon({w,op,conj(w)},{[1 3 2 -2],[4 3],[1 4 2 -1]},[1 2 3 4]);
end