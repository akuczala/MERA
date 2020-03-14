function out = vary(x,f)
    [U,S,V] = svd(f);
    id = S;
    id(id ~= 0) = 1;
    smallid = id;
    new = 1*ctranspose(U*(smallid)*ctranspose(V));
    %Y = new*ctranspose(x);
    %[V,D] = eig(Y);
    %sqrtY = V*sqrt(D)*ctranspose(V);
    %out = sqrtY*x;
    %out = sqrtm(new*ctranspose(x))*x;
    out = new;
    %disp('new min');
    %disp(trace(out*f));
    %disp('min');
    %disp(-trace(S));
end