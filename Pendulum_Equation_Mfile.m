%% Initialize
clear
close all
clc


%% Pendulum Equation
% 0.1223 x'' + 0.012 x' + 2 x = 0
% -> x'' + 0.0981 x' + 16.3532 = 0
% v = x' 

func = @(tVal, xVal, vVal) -0.0981*vVal -16.3532*xVal;


%% Main

del_t = 0.01;
time = 0:del_t:10;

x = zeros(size(time));
v = zeros(size(time));

% Initial Condition
x(1) = 1;
v(1) = 0;


for idx = 1:(length(time)-1)

   [x(idx+1), v(idx+1)] = RK4(func, time(idx), x(idx), v(idx), del_t);
    
end


%% Plot

figure
subplot(1,2,1)
plot(time, x)
title('theta')
xlabel('sec')
ylabel('position')
grid on
hold on

subplot(1,2,2)
plot(time, v)
title('thetaprime')
xlabel('sec')
ylabel('gradient')
grid on
hold on


%% Runge-Kutta 4th

function[x_next, v_next] = RK4(func, time, x, v, del_t)
    
    kx1 = v;
    kv1 = func( time, x, v );

    kx2 = v + del_t*kv1/2;
    kv2 = func( time + del_t/2, x + del_t*kx1/2, v + del_t*kv1/2 );

    kx3 = v + del_t*kv2/2;
    kv3 = func( time + del_t/2, x + del_t*kx2/2, v + del_t*kv2/2 );

    kx4 = v + del_t*kv3;
    kv4 = func( time + del_t, x + del_t*kx3, v + del_t*kv3 );

    dx = del_t*(kx1 + 2*kx2 + 2*kx3 + kx4)/6;
    dv = del_t*(kv1 + 2*kv2 + 2*kv3 + kv4)/6;

    x_next = x + dx;
    v_next = v + dv;

end