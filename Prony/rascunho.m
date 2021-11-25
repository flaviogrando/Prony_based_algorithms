


clc
clear all
close all

f0 = 50;
N = 25;
Fs = N*f0;
Ts = 1/Fs;

f=50;
R = 0.01;

for n=0:25
    
    arg_lin(n+1) = 2*pi*f*Ts*n/N;
    arg_exp(n+1) = 1-1/(exp(1.2*2*pi*f*Ts*n/N)); %-2*pi;
    
end

figure, hold on, grid on
plot(arg_lin)
plot(arg_exp)