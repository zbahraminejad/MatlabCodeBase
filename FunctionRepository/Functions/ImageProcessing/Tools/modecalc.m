function calcedmode=modecalc(values)
[nout,xout]=ksdensity(values);
[~,idx]=max(nout);
calcedmode=xout(idx);