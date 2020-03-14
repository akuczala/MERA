function out = dirtySymmetrizeW(w,S)
out = 0.5*(ncon({S,w,conj(S),conj(S),conj(S)},{[4 -4],[1 2 3 4],[1 -1],[2 -2],[3,-3]}) + w);