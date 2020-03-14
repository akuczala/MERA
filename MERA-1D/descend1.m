function out =  descend1(w,op)
%acts single site descending superoperator on operator
d = length(w);


tensors = {conj(w),op,w};
links = {[1 -2 2 4],[3 4],[1 -1 2 3]};
out = ncon(tensors,links,[1 2 3 4]);
end