d = 2;
u = eye(d*d,d*d); u = reshape(u,d,d,d,d);
udag = conj(u);
randmat = rand(d*d*d,d)+ rand(d*d*d,d)*i;
[U,S,V] = svd(randmat);
w = ctranspose(U(1:d,:));
w = reshape(w,d,d,d,d);
wdag = conj(w);
rho = rand(d,d,d,d);
h = eye(d*d,d*d); h = reshape(h,d,d,d,d);

%environment tensor network for w,u (9 diagrams)
wtensors = {rho,w,u,h,udag,wdag,wdag};
utensors = {rho,w,w,h,udag,wdag,wdag};
%all contracted, h middle
%allLinks = {[3 4 1 2],[5 6 7 3],[8 9 10 4],[11 12 7 8],[13 14 11 12],[13 14 15 16],[5 6 15 1],[16 9 10 2]};
%left, right w missing
legLinks ={[-4 4 1 2], [8 9 10 4], [11 12 -3 8], [13 14 11 12], [13 14 15 16], [-1 -2 15 1], [16 9 10 2]};
sequence = [4 8 9 10 2 16 14 12 11 13 15 1];
Upw = ncon(wtensors,legLinks,sequence);
legLinks ={[3 -4 1 2], [5 6 7 3], [11 12 7 -1], [13 14 11 12], [13 14 15 16], [5 6 15 1], [16 -2 -3 2]};
sequence = [3 7 6 5 1 15 13 11 12 14 16 2];
Upw = Upw + ncon(wtensors,legLinks,sequence);
%u missing
legLinks = {[3 4 1 2], [5 6 -3 3], [-4 9 10 4], [13 14 -1 -2], [13 14 15 16], [5 6 15 1], [16 9 10 2]};
sequence = [3 5 6 1 15 13 14 16 2 9 10 4];
Upu = ncon(utensors,legLinks,sequence);
%all contracted, h left
allLinks = {[3 4 1 2],[5 6 7 3],[8 9 10 4],[11 13 7 8],[14 12 6 11],[12 13 15 16],[5 14 15 1],[16 9 10 2]};
%left, right w missing
legLinks = {[-4 4 1 2], [8 9 10 4], [11 13 -3 8], [14 12 -2 11], [12 13 15 16], [-1 14 15 1], [16 9 10 2]};
sequence = [4 8 9 10 2 16 13 11 12 14 15 1];
Upw = Upw + ncon(wtensors,legLinks,sequence);
legLinks = {[3 -4 1 2], [5 6 7 3], [11 13 7 -1], [14 12 6 11], [12 13 15 16], [5 14 15 1], [16 -2 -3 2]};
sequence = [3 7 6 11 5 1 14 15 12 13 16 2];
Upw = Upw + ncon(wtensors,legLinks,sequence);
%u missing
legLinks = {[3 4 1 2], [5 6 -3 3], [-4 9 10 4], [14 12 6 -1], [12 -2 15 16], [5 14 15 1], [16 9 10 2]};
sequence = [3 5 1 14 6 12 15 16 2 9 10 4];
Upu = Upu + ncon(utensors,legLinks,sequence);
%h right (copy above but horizontal mirror tensors)
%left, right w missing
legLinks = {[4 -3 2 1], [10 9 8 4], [13 11 8 -4], [12 14 11 -1], [13 12 16 15], [15 -2 14 1], [9 10 16 2]};
sequence = [4 8 9 10 2 16 13 11 12 14 15 1];
Upw = Upw + ncon(wtensors,legLinks,sequence);
legLinks = {[-3 3 2 1], [7 6 5 3], [13 11 -2 7], [12 14 11 6], [13 12 16 15], [15 14 5 1], [-4 -1 16 2]};
sequence = [3 7 6 11 5 1 14 15 12 13 16 2];
Upw = Upw + ncon(wtensors,legLinks,sequence);
%u missing
legLinks = {[4 3 2 1], [-4 6 5 3], [10 9 -3 4], [12 14 -2 6], [-1 12 16 15], [15 14 5 1], [10 9 16 2]};
sequence = [3 5 1 14 6 12 15 16 2 9 10 4];
Upu = Upu + ncon(utensors,legLinks,sequence);
%SVD u, w to get updated u,w
UpwMat = reshape(Upw,d^3,d);
[U,S,V] = svd(UpwMat);
S(S ~= 0) = 1 %set nonzero schmidt values to 1,get an identity operator
neww = -conj(U*S*V) %don't transpose - order of indices takes care of this
neww = reshape(neww,d,d,d,d);
UpuMat = reshape(Upw,d^2,d^2);
[U,S,V] = svd(UpuMat);
newu = -conj(U*V) %don't transpose - order of indices takes care of this
newu = reshape(neww,d,d,d,d);