din = 2; dout = 3;
r = rand(din*din,din*din)+ rand(din*din,din*din)*i;
r = 0.1*r*r'/trace(r*r');
randu = expm(i*r);
abs(trace(randu))
randu = reshape(randu,din,din,din,din);
r = rand(dout*dout*dout,dout*dout*dout)+ rand(dout*dout*dout,dout*dout*dout)*i;
r = r*r';
r = .1*r;
randout = expm(i*r);
r = rand(din,din)+ rand(din,din)*i;
r = r*r';
r = .1*r;
randin = expm(i*r);
randin = reshape(randin,din,din);
randout = reshape(randout,dout,dout);