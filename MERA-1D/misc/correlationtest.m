mera = critlong;
w0 = mera{6}{1};w = mera{6}{2};
u0 = mera{7}{1};u = mera{7}{2};
%top = mera{9};
layers = 2;
chiFix = 2;
%make density matrices from MERA
% rhoTop = top*ctranspose(top);
% rhoTop = rhoTop;
% rhoTop = rhoTop/trace(rhoTop);
% rhoTop = reshape(rhoTop,chiFix,chiFix,chiFix,chiFix);
% rho{layers} = updateRho(w{layers},u{layers},conj(w{layers}),conj(u{layers}),rhoTop);
% for tau = fliplr(1:layers-1)
%     rho{tau} = updateRho(w{tau},u{tau},conj(w{tau}),conj(u{tau}),rho{tau+1});
% end
% rho0 = updateRho(w0,u0,conj(w0),conj(u0),rho{1});
% 
 c = 1:3^5; c= c*0;
%take density matrices from MERA
rhoTop = mera{5}{3};
rho = mera{5}{2};
rho0 = mera{5}{1};
Z = [1 0;0 -1]; X = [0 1; 1 0];
op0 = Z;
prod = ncon({op0,op0},{[-1 -3],[-2 -4]},[]);
prod0 = prod;
c(1) =vev2(prod,rho0);
op{1} = ascend1(w0,op0);
prod = ncon({op{1},op{1}},{[-1 -3],[-2 -4]},[]);
c(3) =vev2(prod,rho{1});
for tau = 2:layers
op{tau} = ascend1(w{tau-1},op{tau-1});
prod = ncon({op{tau},op{tau}},{[-1 -3],[-2 -4]},[]);
c(3^(tau-1)) = vev2(prod,rho{tau-1});
end
c(3^(layers)) = vev2(prod,rhoTop);
%some intermediates (more complicated diagrams)
%r = 2
% tensors = {rho{1},w0,w0,u0,op0,op0,conj(u0),conj(w0),conj(w0)};
% links = {[1 2 3 4],[5 6 7 1],[8 9 10 2],[11 12 7 8], ...
%     [15 11],[17 9],[15 12 18 16],[5 6 18 3],[16 17 10 4]};
% c(4) = ncon(tensors,links,[1 5 6 7 3 4 2 8 9 10 11 12 15 18 16 17]);


tensors = {w0,w0,u0,prod0,conj(w0),conj(w0),conj(u0),rho{1}};
links = {[9 8 16 12 ],[1 10 11 13 ],[2 7 16 1 ],[3 6 2 7 ],[9 8 4 14 ],[5 10 11 15 ],[3 6 4 5 ],[12 13 14 15 ]};
ncon(tensors,links,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16])



r = [1 3 9 27 3^4 3^5];
cs = [c(1) c(3) c(9) c(27) c(81) c(3^5)];
plot(real(c));
%arrayfun(@(x) log(x),r)./arrayfun(@(x) log(x),cs)