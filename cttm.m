function mret = cttm(varfas)
n = length(varfas);
mret = zeros(n,1);
for i=1:n
    mret(i) = length(varfas(i).Sequence);
end

    
