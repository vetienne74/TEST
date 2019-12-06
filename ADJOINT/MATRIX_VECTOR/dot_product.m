
close all; clear all;

% input parameters

N = 100 ;
SCALE = 1.e0 ;

A = sprand(N,N,1) .* SCALE ;

figure
hold on
imagesc(A)
title('Operator') ;
colorbar
axis tight

figure
hold on
imagesc(A')
title('Adjoint Operator') ;
colorbar
axis tight

X = sprand(N,1,1) ;
Y = sprand(N,1,1) ;

AX = A * X ;
A_ADJY = A' * Y ;

EQ1 = Y' * AX
EQ2 = A_ADJY' * X
RES=EQ1 - EQ2
RES/SCALE




