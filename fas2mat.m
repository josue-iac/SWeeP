function W160k = fas2mat(xfas, varargin)
% Define vetor com vetores 160k para as sequencias de xfas
% opcionalmente o vetor pode ser reordenado de acordo com iord
V = varargin;
%
withPos = varginfind(V,'Withpos',0);
toSum = varginfind(V,'toSum',0);
iord = varginfind(V,'Ordem',1:160000);
%
Ps = primon(1);
nmx = max(cttm(xfas));
%
if withPos
    Ps = primon(nmx);
end
%
F = struct2cell(xfas);
S = F(2,:);
% M = cellfun(@(x) aa2mat(x,withPos),S,'UniformOutput',false);
M = cellfun(@(x) matinline(aa2mat(x, Ps, withPos, toSum)),S,'UniformOutput',false);
W160k = cell2mat(M');
W160k = W160k(:,iord);
end

function [mret xy] = aa2mat(xseq, Ps, withPos, toSum)
% xseq - dna sequence
% varargin - sampling length. Eg. 1 monopeptide, 2 dipeptide...
%
l = 2;
%
s = 20^l;
L5 = dna2list(xseq,l*2+1);
%
%xy = round(([aa2num2(L5(:,1:2)) aa2num2(L5(:,4:5))]-1) * s / 400 + 0.5);
xy = [aa2num2(L5(:,1:l)) aa2num2(L5(:,l+2:l*2+1))];
%
inot = find(~prod(1-double(xy<=0),2));
xy(inot,:) = [];
%
inds = ij2inds(xy,s);
inds(inds>(s^2)) = [];
M = zeros(s,s,'double');
if withPos
    if length(Ps)<length(inds)
        Ps = primon(length(inds));
    end
    [vals inds] = unifysameinds(inds, Ps(1:length(inds)), @(x) prod(x.^(1/length(x))));
    M(inds) = vals;
else
    if toSum
        [vals inds] = unifysameinds(inds, ones(1,length(inds)), @(x) sum(x));
        M(inds) = vals;
    else
        M(inds) = 1;
    end
end
%
mret = M;
end