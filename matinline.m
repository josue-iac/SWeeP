function mret = matinline(M)
% Lineariza matriz por colunas
mret = M(find(M==M))';