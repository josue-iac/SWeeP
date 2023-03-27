function mret = mat2celllines(mat)
% gera cell com linhas de mat
mret = mat2cell(mat,ones(1,length(mat(:,1))),length(mat(1,:)));