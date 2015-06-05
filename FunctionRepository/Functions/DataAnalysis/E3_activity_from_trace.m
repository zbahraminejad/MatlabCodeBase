function E3=E3_activity_from_trace(trace,k_synthesis)
% tr should be either the geminin or cdt1 values, k_syn should be your
% estimate of the promoter stength. To estimate, take the slope of the line
% when geminin/cdt1 increases linearly. (I use a k_syn value of 0.0126,
% which was calculated from normalized geminin traces)
if nargin<2
    k_synthesis=0.0163;
end

for i=1:size(trace,1);
    H_n1 = trace(i,2 : end); % vectorized H_{n+1}
    H_n0 = trace(i,1 : end - 1); % vectorized H_n
    E3_initial = -2 * (H_n1 - H_n0 - k_synthesis) ./ (H_n0 + H_n1);

    E3_initial(E3_initial<0)=0;
    E3_initial2=conv(E3_initial,normpdf(-3:3,0,2),'same');
    E3_min=E3_initial2-min(E3_initial2);
    E3_norm=E3_min/max(E3_min);
    E3(i,:)=E3_norm;
end