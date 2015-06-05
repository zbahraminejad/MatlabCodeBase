function p=fitlineTotalLeastSquares(x,y)
% p(1)=slope
% p(2)=intercept
% written following the method given on the wikipedia page
% (http://en.wikipedia.org/wiki/Total_least_squares#Computation)

A=[x(:) ones(size(x(:)))];
B=y(:);
C=[A B];

%Check for and handle NaN data
nanInd=isnan(sum(C,2));
if ~isempty(nanInd)
    fprintf('Removing NaN-containing data points.\n');
    C(nanInd,:)=[];
end

[U,S,V]=svd(C,0);
VAB=V(1:2,3:end);
VBB=V(3:end,3:end);
X=-1*VAB/VBB;
p=X(:);
