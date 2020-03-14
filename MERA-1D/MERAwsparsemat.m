%parameters
h = 0.5; J = 1;
d = 2; %site dimension
chi = 10; %scale invariant dimension
chi1 = 4; %transitional dimension
randDensity = 0.4;
%initialize
Z = [1 0; 0 -1]; X = [0 1; 1 0]; id = [1 0; 0 1]; 
h0 = J*(kron(Z,Z) + h*(kron(X,id)+kron(id,X)));
%inital MERA matrices
u0 = speye(2*d); %transitional disentanglers
u1 = speye(2*chi1);
u = speye(2*chi); %invariant layer disentangler
%construct random isometries w satisfying w^dag w = 1
randmat = rand(3*d,chi1); %transitional isometries
[U,S,V] = svd(randmat);
w0 = ctranspose(U(1:chi1,:)); % 3*d x chi1
randmat = rand(3*chi1,chi); 
[U,S,V] = svd(randmat);
w1 = ctranspose(U(1:chi,:)); % 3*chi1 x chi
randmat = rand(3*chi,chi);
[U,S,V] = svd(randmat);
w = ctranspose(U(1:chi,:)); %invariant layer isometry (3*chi x chi)
%initialize rho etc as empty matrices that will be determined in first
%iteration
rho0 = sparse(2*d,2*d);
rho1 = sparse(2*chi1,2*chi1);
rho = sparse(2*chi,2*chi);
h1 = sparse(2*chi1,2*chi1);
h = sparse(2*chi,2*chi);
%LOOP

%construct scale invariant descending operator
db = chi; dt = chi; %"bottom", "top" dimensions
buildD = sparse(6*db,2*dt);
nu = nnz(u); nw = nnz(w);% # of nonzero elements in u,w
for i = 1:nu*nw*nw
        buildD(i,j) = u(
end
