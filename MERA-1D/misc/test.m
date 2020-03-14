d = 5;
X = rand(5); Y = rand(5);
%tr[X A(Y)]
%tr1 = trace(X*ascend1(w,Y));
%tr2 = trace(descend1(w,X)*Y);
%tr1 - tr2
%trace(descend1(w,X))/trace(X)

%compute eigenvalues of descending operator
tensors = {conj(w2),w2};
links = {[1 -4 2 -2],[1 -3 2 -1]};

%vertical flip of diagram conjugates evals
d=5;
D = ncon(tensors,links,[1 2]);
Dmat = transpose(reshape(D,d*d,d*d)); %transpose D for eigops of tensor!!!
[V,Dev] = eigs(Dmat); %first 2 evals
evals = diag(Dev);
disp(evals)
evec1 = V(:,1);evec2 = V(:,2);

eop1 = reshape(evec1,d,d); eop2 = reshape(evec2,d,d);
eop3 = reshape(V(:,3),d,d);
descend1(w2,eop2) - evals(2)*eop2
trace(descend1(w2,eop2))-evals(2)*trace(eop2)
trace(descend1(w2,eop2))-trace(eop2)
%det(descend1(w,eop2)-evals(2)*eye(4))
%disp(evals);