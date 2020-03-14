function out = checkuSymmetry(u,S)
    out = ncon({S,S,u,conj(S),conj(S)},{[3 -3],[4 -4],[1 2 3 4],[1 -1],[2 -2]}) - u;
end