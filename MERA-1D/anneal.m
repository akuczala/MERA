function [wout,uout] = anneal(w,u,t)
dims = size(w);din = dims(1); dout = dims(4);
r = rand(din*din,din*din)+ rand(din*din,din*din)*i;
r = t*r*r'/trace(r*r');
randu = expm(i*r);
%abs(trace(randu))
randu = reshape(randu,din,din,din,din);
r = rand(dout,dout)+ rand(dout,dout)*i;
r = r*r';
r = t*r;
randout = expm(i*r);
r = rand(din^3,din^3)+ rand(din^3,din^3)*i;
r = r*r';
r = t*r;
randin = expm(i*r);
randin = reshape(randin,din,din,din,din,din,din);
randout = reshape(randout,dout,dout);
wout = ncon({randin,w,randout},{[-1 -2 -3 1 2 3],[1 2 3 4],[4 -4]});
uout = ncon({randu,u},{[-1 -2 1 2],[1 2 -3 -4]},[1 2]);
end