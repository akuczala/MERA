function [neww,newu,newwdag,newudag] = updateMERA(w,u,wdag,udag,rho,h,steps,T)
    
dims = size(w);din = dims(1); dout = dims(4);
    for j = 1:steps
       %kick it
       
       %calculate initial energy
       Ei = vev2(updateH(w,u,wdag,udag,h),rho);
       wi = w; ui = u; wdagi = wdag; udagi = udag;
       %minimize
       [w,u,wdag,udag] = updateMERApar(w,u,wdag,udag,rho,h);
       if T > 0
           %calculate new energy
           Ef = vev2(updateH(w,u,wdag,udag,h),rho);
           dE = real(Ef - Ei);
           
           %probability of doing random
           p = exp(-abs(dE)/T);
           %disp(dE)
           %disp(T)
           disp(p)
           %disp(p);
           r = rand(1);
           if p>r 
               [w,u] = anneal(w,u,0.2);
               wdag = conj(w); udag = conj(u);
               disp('kick');
           end
       end
    end
    neww = w; newu = u; newwdag = wdag; newudag = udag;
    %[neww,newu,newwdag,newudag] = updateMERAonce(w,u,wdag,udag,rho,h);
end