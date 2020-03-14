function [neww,newu,newwdag,newudag] = updateMERAonce(w,u,wdag,udag,rho,h)

dims = size(w);din = dims(1); dout = dims(4);
%environment tensor network for w,u (9 diagrams)
%9/9/14 optimized contraction sequence with netcon
%%%%%%u environment
utensors = {rho,w,w,h,udag,wdag,wdag};
%h middle
legLinks = {[3 4 1 2], [5 6 -3 3], [-4 9 10 4], [13 14 -1 -2], [13 14 15 16], [5 6 15 1], [16 9 10 2]};
%sequence = [3 5 6 1 15 13 14 16 2 9 10 4];
sequence =  [5 6 9 10 2 4 1 3 15 16 13 14];
Upu = ncon(utensors,legLinks,sequence);
%h left
legLinks = {[3 4 1 2], [5 6 -3 3], [-4 9 10 4], [14 12 6 -1], [12 -2 15 16], [5 14 15 1], [16 9 10 2]};
%sequence = [3 5 1 6 14 15 12 16 2 9 10 4];
sequence = [9 10 2 4 1 15 16 12 14 3 5 6];
Upu = Upu + ncon(utensors,legLinks,sequence);
%h right
legLinks = {[4 3 2 1], [-4 6 5 3], [10 9 -3 4], [12 14 -2 6], [-1 12 16 15], [15 14 5 1], [10 9 16 2]};
%sequence = [3 5 1 6 14 15 12 16 2 9 10 4];
sequence = [9 10 2 4 1 15 16 12 14 3 5 6];
Upu = Upu + ncon(utensors,legLinks,sequence);

%SVD u environment to get new u

UpuMat = reshape(Upu,din^2,din^2);
[U,S,V] = svd(UpuMat,'econ');
%disp(max(max(UpuMat - U*S*ctranspose(V))));

S(S ~= 0) = 1;
newu = U*S*ctranspose(V); %don't transpose - order of indices takes care of this
newu = conj(newu); %removing this seemed to improve energy ratio of first to 0th
%layer
newu = reshape(newu,din,din,din,din);
newu = -newu;

%STEP SIZE
%newu = u*(1-stepSize) + stepSize*newu;
%newu = sqrtm(new*ctranspose(x))*x;

newudag = conj(newu);


%%%%%%%%w environment
wtensors = {rho,w,u,h,udag,wdag,wdag};
%all contracted, h middle
%allLinks = {[3 4 1 2],[5 6 7 3],[8 9 10 4],[11 12 7 8],[13 14 11 12],[13 14 15 16],[5 6 15 1],[16 9 10 2]};
%left, right w missing
legLinks ={[-4 4 1 2], [8 9 10 4], [11 12 -3 8], [13 14 11 12], [13 14 15 16], [-1 -2 15 1], [16 9 10 2]};
%sequence = [4 8 9 10 2 16 14 12 11 13 15 1];
sequence =  [9 10 13 14 11 12 8 16 2 4 1 15];
Upw = ncon(wtensors,legLinks,sequence);
legLinks ={[3 -4 1 2], [5 6 7 3], [11 12 7 -1], [13 14 11 12], [13 14 15 16], [5 6 15 1], [16 -2 -3 2]};
%sequence = [3 7 6 5 1 15 13 11 12 14 16 2];
sequence =  [5 6 13 14 11 12 7 15 1 3 2 16];
Upw = Upw + ncon(wtensors,legLinks,sequence);

%all contracted, h left
%allLinks = {[3 4 1 2],[5 6 7 3],[8 9 10 4],[11 13 7 8],[14 12 6 11],[12 13 15 16],[5 14 15 1],[16 9 10 2]};
%left, right w missing
legLinks = {[-4 4 1 2], [8 9 10 4], [11 13 -3 8], [14 12 -2 11], [12 13 15 16], [-1 14 15 1], [16 9 10 2]};
%sequence = [4 8 9 10 2 16 13 11 12 14 15 1];
sequence =  [9 10 2 4 1 15 16 12 14 8 11 13];
Upw = Upw + ncon(wtensors,legLinks,sequence);
legLinks = {[3 -4 1 2], [5 6 7 3], [11 13 7 -1], [14 12 6 11], [12 13 15 16], [5 14 15 1], [16 -2 -3 2]};
%sequence =  [3 5 1 6 14 15 12 11 13 7 16 2];
sequence = [15 12 14 11 13 5 6 7 1 3 2 16];
Upw = Upw + ncon(wtensors,legLinks,sequence);

%h right (copy above but horizontal mirror tensors)
%right, left w missing
%fixed wrongness 8/20/14
legLinks = {[4 -4 2 1],[9 10 8 4],[13 11 8 -1],[12 14 11 -2],[13 12 16 15],[15 14 -3 1],[9 10 16 2]};
%this one was WRRONG 8/20/14
%legLinks = {[4 -4 2 1], [10 9 8 4], [13 11 8 -1], [12 14 11 -2], [13 12 16 15], [15 -3 14 1], [9 10 16 2]};
%sequence = [4 8 9 10 2 16 13 11 12 14 15 1];
sequence =  [9 10 2 4 1 15 16 12 14 8 11 13];
Upw = Upw + ncon(wtensors,legLinks,sequence);
legLinks = {[-4 3 2 1], [7 6 5 3], [13 11 -3 7], [12 14 11 6], [13 12 16 15], [15 14 5 1], [-1 -2 16 2]};
%sequence =  [3 5 1 6 14 15 12 11 13 7 16 2];
sequence = [15 12 14 11 13 5 6 7 1 3 2 16];
Upw = Upw + ncon(wtensors,legLinks,sequence);

%SVD w to get updated w
UpwMat = reshape(Upw,din^3,dout);
[U,S,V] = svd(UpwMat,'econ');

%disp(max(max(UpwMat - U*S*ctranspose(V))));

S(S ~= 0) = 1; %set nonzero schmidt values to 1,get an identity operator

neww = U*S*ctranspose(V); %don't transpose - order of indices takes care of this

neww = conj(neww);%removing this seemed to improve energy ratio of first to 0th
%layer
neww = reshape(neww,din,din,din,dout);
neww = -neww;

%STEP SIZE
%neww = w*(1-stepSize) + neww*stepSize;

newwdag = conj(neww);

%Z = [1 0; 0 -1];
%newu = dirtySymmetrizeU(u,Z);
%neww = dirtySymmetrizeW(w,Z);



end
