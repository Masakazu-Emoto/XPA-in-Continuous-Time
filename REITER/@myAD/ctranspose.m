function x=ctranspose(x)
  % by SeHyoun Ahn, Jan 2016

  % Note:
  %    Only allows real derivatives, so this does not conjugate derivative values

  [n,m] = size(x.values);
  try
    x.derivatives = dertransp(x.derivatives,n);
  catch
    l = size(x.derivatives,2);
    [i,j,v] = find(x.derivatives);
    aux1 = mod(i,n);
    aux1(aux1==0) = n;
    x.derivatives=sparse(((aux1-1)*m+ceil(i/n)),j,v,m*n,l);
  end
  x.values=x.values';
end
