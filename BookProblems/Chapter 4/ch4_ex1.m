% Chapter 4, exercise 1
clc; clear all; close all;
%z=     [I:1_2	I :2_1 	I:1_3 	I:3_1	I:2_3	I:3_2	V:1     V:2     V:3]^T
z=      [52,	-49,	-15,	16,     -82,	80,     11, 	6,  	14].';
V=z(7:9);
om=0.1;

H_DC=[ 1/om,   -1/om,  0;
    -1/om,  1/om,   0;
    1/om,   0,      -1/om;
    -1/om,  0,      1/om;
    0,      1/om,   -1/om;
    0,      -1/om,  1/om;
    1,      0,      0;
    0,      1,      0;
    0,      0,      1;];

fprintf('The true V values are:\n')
disp(V')
% for part a, R=identity
R=eye(length(z));
V=(transpose(H_DC)*R^-1*H_DC)^-1*transpose(H_DC)*R^-1*z;
fprintf('for part a the V values were:\n')
disp(V')
%for part b, the variance 100 times smaller than the variance values used
%in a:
R(1,1)=R(1,1)/100; R(3,3)=R(3,3)/100;
V=(transpose(H_DC)*R^-1*H_DC)^-1*transpose(H_DC)*R^-1*z;
fprintf('for part a the b values were:\n')
disp(V')
%pg 106