function diff = checkHerm2(A)
d= length(A);
B = reshape(A,d^4,d^4);
diff = B - transpose(conj(B));
diff = max(diff(:));
end