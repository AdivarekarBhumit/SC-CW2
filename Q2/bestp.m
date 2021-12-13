% Q2a
% N = [27; 216; 729];
% ts = [0.002587; 0.069477; 1.266354]; 
% 
% p = polyfit(log(N), log(ts), 1);
% 
% plot(N, ts);

%Q2d
N = [2197; 4913; 9261];
ts = [0.05945; 0.17599; 0.26773]; 

p = polyfit(log(N), log(ts), 1);

plot(N, ts);


