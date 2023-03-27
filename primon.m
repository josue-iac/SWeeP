function mret = primon(n)
%Retorna os n primeiros numeros primos (2:n)
m = n;
nmax = 100;
ps = find(isprime(1:nmax));
while length(ps)<m
    nmax = nmax*2;
    ps = find(isprime(1:nmax));
end    
mret = ps(1:m);