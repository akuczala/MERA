function out = ascend3c(op,w,u,wdag,udag)

tensors = {op,w,w,w,u,u,udag,udag,wdag,wdag,wdag};
links = {[13 16 19 12 15 18 ],[1 2 8 -4 ],[11 15 17 -5 ],[5 4 3 -6 ],[9 12 8 11 ],[18 6 17 5 ],[9 13 10 14 ],[19 6 20 7 ],[1 2 10 -1 ],[14 16 20 -2 ],[7 4 3 -3 ]};
seq = [1 2 3 4 11 20 16 19 13 14 6 18 9 12 15 17 5 7 8 10];
out = ncon(tensors,links,seq);
end