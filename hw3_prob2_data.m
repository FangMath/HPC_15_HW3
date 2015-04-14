close all; clear all;
x = 2.^[0:9];
ygs = [110.98, 55.56, 28.15, 14.59, 110.33, 11.62,11.53, 12.33, 14.78, 19.99];
yjc = [108.61, 55.26, 27.45, 14.25, 47.75, 10.90, 10.83, 10.99, 12.81, 16.21];

bar(log(x)/log(2),ygs, 'r'); hold on;
bar(log(x)/log(2),yjc,'b');
legend('GS','Jacobi');
xlabel('Number of threads is 2^k');
ylabel('timings (secends)');
title('timings with different number of threads');
