N = [81; 161; 321; 641];
%ts = [0.007; 0.013; 0.078; 0.5503]; % when using linear solve and jacobian
ts = [0.0021; 0.0038; 0.0059; 0.0116]; % when using sparseThomas

p = polyfit(log(N), log(ts), 1);

plot(N, ts);
