function [w,u,wdag,udag] = initTensors(d1,d2,makeId)

 %construct random unitary disentanglers u
randmat = rand(d1^2,d1^2)+rand(d1^2,d1^2)*i;
[U,S,V] = svd(randmat);
u = ctranspose(U(1:d1^2,:));

%DEBUG
%make random orthogonal matrix
%randmat = rand(d1^2,d1^2);
%[U,S,V] = svd(randmat);
%u = U;

if makeId
    u = eye(d1^2,d1^2); 
end

u = reshape(u,d1,d1,d1,d1);
udag = conj(u);

%construct random isometries w satisfying w^dag w = 1
randmat = rand(d1^3,d2)+rand(d1^3,d2)*i; %transitional isometries
[U,S,V] = svd(randmat);
w = ctranspose(U(1:d2,:)); % d^3 x chi1

%DEBUG
%make random orthogonal matrix
%randmat = rand(d1^3,d2);
%[U,S,V] = svd(randmat);
%w = transpose(U(1:d2,:));

if(makeId)
    w = eye(d1^3,d2); %make w identity matrix
end

w = reshape(w,d1,d1,d1,d2);
wdag = conj(w);

end