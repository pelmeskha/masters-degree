clear all; clc; close all;
A=-2; B=3;
a=-3; b=1;
space=(logspace(a,b,10)-10^a)*(B-A)/(10^b-10^a)+A;
space=(logspace(b,a,10)-10^b)*(B-A)/(10^a-10^b)+A;
plot(space,'k.')