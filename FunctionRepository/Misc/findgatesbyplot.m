function findgatesbyplot(dapilist,cebpbvalues,ppargvalues)


temp = regexp(dapilist,'_[A-Z][0][0-3]_','match');
tempind = ~cellfun(@(x) isempty(x), temp);

% Check for gating parameters

dscatter([cebpbvalues{:}]',[ppargvalues{:}]')
title('All population');
figure,scatter([cebpbvalues{tempind}]',[ppargvalues{tempind}]','b');
hold on
scatter([cebpbvalues{~tempind}]',[ppargvalues{~tempind}]','r');
hold off
title('Control (blue) versus Treatment (red)');
figure ,dscatter([cebpbvalues{tempind}]',[ppargvalues{tempind}]')
title('Control Population');
figure ,dscatter([cebpbvalues{~tempind}]',[ppargvalues{~tempind}]')
title('Treatment Population');

end