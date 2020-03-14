function [w,u,wdag,udag] = initTensors(d1,d2,makeId)


randmat = rand(d1^2,d1^2)+rand(d1^2,d1^2)*i;
[U,S,V] = svd(randmat);
u = ctranspose(U(1:d1^2,:));
if makeId
    u = eye(d1^2,d1^2); 
end
u = reshape(u,d1,d1,d1,d1);
udag = conj(u);

randmat = rand(d1^3,d2)+rand(d1^3,d2)*i; %transitional isometries
[U,S,V] = svd(randmat);
w = ctranspose(U(1:d2,:)); % d^3 x chi1
if(makeId)
    w = eye(d1^3,d2); %make w identity matrix
end
w = reshape(w,d1,d1,d1,d2);
wdag = conj(w);

end