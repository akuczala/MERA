utensors = {rho0,w0,w0,h0,conj(u0),conj(w0),conj(w0)};
din = 2;
%h middle
legLinks = {[3 4 1 2], [5 6 -3 3], [-4 9 10 4], [13 14 -1 -2], [13 14 15 16], [5 6 15 1], [16 9 10 2]};
%sequence = [3 5 6 1 15 13 14 16 2 9 10 4];
sequence = [3 5 6 1 15 13 14 16 2 9 10 4];
Upu = ncon(utensors,legLinks,sequence);

UpuMat = reshape(Upu,din^2,din^2);
[U,S,V] = svd(UpuMat);
%disp(max(max(UpuMat - U*S*ctranspose(V))));

S(S ~= 0) = 1;
newu = U*S*ctranspose(V); %don't transpose - order of indices takes care of this
newu = conj(newu); %removing this seemed to improve energy ratio of first to 0th
%layer
newu = reshape(newu,din,din,din,din);
newu = -newu;