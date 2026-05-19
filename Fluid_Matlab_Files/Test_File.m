% This  is a test. 
clear all
close all
clc 

% Time Discretization: 
N     = 1e4;
t_0   = 0; 
t_end = 10;
t     = linspace(t_0, t_end, N);

%  Define a function: 
f = @(t) 1/30 * exp(-t.^2) .* cos(t);

% Plot: 
figure(1)
plot(t,f(t),'k','linewidth',5)
title('My First Plot','fontsize',20)
grid on
grid minor



