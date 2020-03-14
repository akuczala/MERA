function out = descend3c(op,w,u,wdag,udag)
tensors = {w,w,w,u,u,udag,udag,wdag,wdag,wdag,op};
legLinks = {[1 2 3 15 ],[6 -2 8 16 ],[10 13 14 17 ],[4 -1 3 6 ],[-3 11 8 10 ],[4 -4 5 7 ],[-6 11 9 12 ],[1 2 5 18 ],[7 -5 9 20 ],[12 13 14 19 ],[15 16 17 18 20 19 ]};
seq = [6 9 1 2 13 14 17 19 15 18 12 20 5 7 10 11 3 4 8 16];
out = ncon(tensors,legLinks,seq);
end