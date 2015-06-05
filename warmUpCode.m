clc
close all;
clear all;
x = -2.9:0.1:2.9;
y = randn(10000,1);
subplot(1,3,1); hist(y,x)
y = randn(10000,1);
subplot(1,3,2); hist(y,x)
y = randn(10000,2);
subplot(1,3,3); hist(y,x)


