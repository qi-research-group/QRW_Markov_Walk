function H = buildH(a,b,c)
% H = buildH(a,b,c)
% m = number of states
% a = oﬀ diag left
% b = diag
% c = oﬀ diag right
m = size(a,1);
H = zeros(m);
H(1,[1 2])= [ b(1) c(1)];
H(m,[m-1 m]) = [a(m) b(m)];
for k = 2:(m-1)
    K = [(k-1) k (k+1)] ;
    H(k,K) = [a(k) b(k) c(k)];
end
