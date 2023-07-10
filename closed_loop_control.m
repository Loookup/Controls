%% Initialize
close all
clear
clc

%% Function info

t = 0:0.01:1;
wn = 19;        % Natural Frequnecy
zeta = 8/19;    % Damping Coefficient ( Underdamped : zeta < 1)

%% Peak Time (Overshoot)

tp = pi/sqrt(297);      % First Element among derivative of c(t) is zero

str_peak = ['Peak tp : ', num2str(tp)];

%% Rise Time

PI = atan(8/sqrt(297));

tr = (pi/2 + PI)/sqrt(297);     % Time from 0 to 1.0 (if zeta >= 1, 0.1 ~ 0.9)

str_rise = ['Rise tr : ', num2str(tr)];

%% Settling Time

Designated_Error = 0.02;

ts = -log(Designated_Error*sqrt(297)/19)/8;        % Timr for reaching its steady-state value include designated error

str_settle = ['Settle ts : ', num2str(ts)];

%% Percent Overshoot

osp = exp(-zeta*pi/sqrt(1-zeta^2))*100;     % (peak value - 1) / 1

str_percent_os = ['OS% : ', num2str(osp), '%'];

%% Bounded

boundary = ones(size(t));

Upper_Bound = boundary.*(1 + Designated_Error);

Lower_Bound = boundary.*(1 - Designated_Error);

%% Plot

figure
plot(t, CT(t))
xlabel('time(sec)')
ylabel('c(t)')
title('Closed-Loop Transfer Function')
grid on
hold on

plot(tp, CT(tp) , 'r*')
text(tp+0.02,CT(tp), str_peak)
hold on

plot(tr, CT(tr), 'r*')
text(tr+0.02, CT(tr), str_rise)
hold on

plot(ts, CT(ts), 'r*')
text(ts+0.02, CT(ts)-0.03, str_settle)
hold on

text(0.6, 1.2, str_percent_os, 'FontSize', 15)
hold on

plot(t, Upper_Bound, '--')
hold on

plot(t, Lower_Bound, '--')
hold on

%% Closed-Loop Transfer Function
function[y] = CT(t)

    y = 1-exp(-8.*t).*(cos(sqrt(297).*t) + 8*sin(sqrt(297).*t)/sqrt(297));

end

%% Closed-Loop Transfer Function Derivative
function[dy] = CT_D(t)

    dy = 1-exp(-8*tp).*(cos(sqrt(297)*tp) + 8*sin(sqrt(297)*tp)/sqrt(297));

end