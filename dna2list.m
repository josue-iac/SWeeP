%Cria matriz w com padroes (exameros) por linha para um trecho dado de dna
function mret = dna2list(sdna, q)
n = length(sdna);
xcol = uint32(1:(n-q+1))';
%
inds = repmat(xcol,1,q) + repmat(uint32(1:q),(n-q+1),1)-1;
mret = sdna(inds);
if q==1
    mret = mret';
end    
