function rounded = roundTensor(T,prec)

rounded = round(T*10^(prec))*10^(-prec); %round
%rounded = (T > 10^(-prec)).*T; %set all numbers smaller than precision to 0
end