function mret = aa2num2(xseq)
%  comments
    [n m] = size(xseq);
    vls = double(aa2int(upper(xseq)))-1;
    pot = repmat(double(0:(m-1)),n,1);
    mret = sum((repmat(20,n,m).^pot).*vls,2)+1;
    inot = find(~prod(1-double(vls<0 | vls>19),2));
    mret(inot) = -1;
end
  