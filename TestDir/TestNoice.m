% this script is used to test noice functions

close all; clear all; clc;
pct_error=10;

t=0:1e-3:4;
signal=1*sin(pi*t);

signal_w_error=add_nor_error(signal,5);
signal_w_error2 = signal + randn(size(signal)).*0.01.*abs(signal);

plot(t,signal)
hold on
plot(t,signal_w_error)
legend('Original Signal','Signal with Error')
hold off

figure
plot(t,signal)
hold on
plot(t,signal_w_error2)
legend('Original Signal','Signal with Error')
