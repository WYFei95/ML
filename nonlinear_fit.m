clc;
clear;
close all;
exp_LTP_raw=[   1	1.53E-08
2	3.16E-08
3	7.11E-08
4	1.36E-07
5	2.70E-07
6	4.83E-07
7	7.76E-07
8	1.08E-06
9	1.47E-06
10	1.91E-06
11	2.37E-06
12	2.85E-06
13	3.31E-06
14	3.76E-06
15	4.18E-06
16	4.64E-06
17	5.02E-06
18	5.36E-06
19	5.69E-06
20	5.97E-06
21	6.22E-06
22	6.42E-06
23	6.62E-06
24	6.77E-06
25	6.92E-06
26	7.04E-06
27	7.15E-06
28	7.23E-06
29	7.30E-06
30	7.37E-06
31	7.44E-06
32	7.48E-06
33	7.53E-06
34	7.58E-06
35	7.62E-06
36	7.66E-06
37	7.68E-06
38	7.69E-06
39	7.72E-06
40	7.74E-06];

exp_LTD_raw=[   1	6.84E-06
2	5.83E-06
3	4.90E-06
4	4.05E-06
5	3.30E-06
6	2.63E-06
7	2.05E-06
8	1.58E-06
9	1.18E-06
10	8.72E-07
11	6.23E-07
12	4.44E-07
13	3.17E-07
14	2.17E-07
15	1.55E-07
16	1.09E-07
17	7.60E-08
18	5.49E-08
19	3.99E-08
20	2.96E-08
21	2.21E-08
22	1.73E-08
23	1.38E-08
24	1.08E-08
25	8.93E-09
26	7.40E-09
27	6.16E-09
28	5.19E-09
29	4.51E-09
30	3.94E-09
31	3.46E-09
32	3.10E-09
33	2.73E-09
34	2.49E-09
35	2.33E-09
36	2.05E-09
37	1.88E-09
38	1.80E-09
39	1.67E-09
40	1.55E-11];
xf_ltp = floor(max(exp_LTP_raw(:,1)));
xf_ltd = floor(max(exp_LTD_raw(:,1)));
x_step_ltp = 1/xf_ltp;
x_step_ltd = 1/xf_ltd;
yf_ltp = max(exp_LTP_raw(:,2));
yf_ltd = max(exp_LTD_raw(:,2));
yi_ltp = min(exp_LTP_raw(:,2));
yi_ltd = min(exp_LTD_raw(:,2));

exp_LTP(:,1) = exp_LTP_raw(:,1)/xf_ltp;
exp_LTD(:,1) = exp_LTD_raw(:,1)/xf_ltd;
exp_LTP(:,2) = (exp_LTP_raw(:,2)-yi_ltp)/(yf_ltp-yi_ltp);
exp_LTD(:,2) = (exp_LTD_raw(:,2)-yi_ltd)/(yf_ltd-yi_ltd);

plot(exp_LTP(:,1), exp_LTP(:,2), 'bo', 'LineWidth', 1);
hold on;
plot(exp_LTD(:,1), exp_LTD(:,2), 'ro', 'LineWidth', 1);

xf = 1;
A_LTP = 0.5;
B_LTP = 1./(1-exp(-1./A_LTP));
A_LTD = -0.2;
B_LTD = 1./(1-exp(-1./A_LTD));

% LTP fitting
var_amp = 0.035;    % LTP cycle-to-cycle variation
rng(103);
x_ltp(1) = 0;
y_ltp(1) = 0;
for n=1:1/x_step_ltp+1
    x_ltp(n+1) = x_ltp(n)+x_step_ltp;
    y_ltp(n+1) = B_LTP(1)*(1-exp(-x_ltp(n+1)/A_LTP(1)));
    delta_y = (y_ltp(n+1)-y_ltp(n)) + randn*var_amp;
    y_ltp(n+1) = y_ltp(n) + delta_y;   
    if y_ltp(n+1)>=1
        y_ltp(n+1)=1;
    elseif y_ltp(n+1)<=0
        y_ltp(n+1)=0;
    end
    x_ltp(n+1) = -A_LTP(1)*log(1-(y_ltp(n+1))/B_LTP(1));
end
plot((0:n-1)/(n-1), y_ltp(1:n), 'b', 'linewidth', 2);

% LTD fitting
var_amp = 0.025;    % LTD cycle-to-cycle variation
rng(898);
x_ltd(1) = 1;
y_ltd(1) = 1;
for n=1:1/x_step_ltd+1
    x_ltd(n+1) = x_ltd(n)-x_step_ltd;
    y_ltd(n+1) = B_LTD(1)*(1-exp(-x_ltd(n+1)/A_LTD(1)));
    delta_y = (y_ltd(n+1)-y_ltd(n)) + randn*var_amp;
    y_ltd(n+1) = y_ltd(n) + delta_y;
    if y_ltd(n+1)>=1
        y_ltd(n+1)=1;
    elseif y_ltd(n+1)<=0
        y_ltd(n+1)=0;
    end
    x_ltd(n+1) = -A_LTD(1)*log(1-(y_ltd(n+1))/B_LTD(1));
end
x_start = numel(x_ltd(:));
x_end = numel(x_ltd(:)) - n;
plot((n-1:-1:0)/(n-1), y_ltd(1:n), 'r', 'linewidth', 2);

xlabel('Normalized Pulse #');
ylabel('Normalized Conductance');
legend('Exp data (LTP)','Exp data (LTD)', 'Fit (LTP)', 'Fit (LTD)', 'location', 'southeast');