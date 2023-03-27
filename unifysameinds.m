function [vals inds] = unifysameinds(ind0, x, func)
% retorna lista de pares de indices unificados (ind) cujo valor (x) associado seja
% transformado por @func ex:  
% [a b] = unifysameinds([1 2 2 2 3 3],[5 3 7 4 5 8],@max) 
%
% a =
% 
%      5     7     8
% 
% 
% b =
% 
%      1     2     3
%%   
if min(size(x))==1
    itn = vetincol(x)';
else
    itn = x;
end
ids = vetincol(ind0)';
%
[idsord ind] = sort(ids);
[xx ii] = size(itn);
if xx > 1
    itnord = itn(ind,:);
else
    itnord = itn(ind)';
end
[uidx ii] = unique(idsord);
ii = vetincol(ii)'; %linha adicionada em 06/2019
list = [(1+ii-[1 diff(ii)])' ii'];



list(1,1) = 1;
% Antigo
% imx = intervallist2inds(list);
% iun = mat2celllines(imx);
% iun = cellfun(@unique,iun,'UniformOutput',false);
%
iun = (cellfun(@(x) x(1):x(end),mat2celllines(list),'UniformOutput',false));
% Acima em teste
[xx cols] = size(itnord);
vals = cell2mat(cellfun(@(x) func(itnord(x,1:cols)),iun,'UniformOutput',false));
inds = uidx;
