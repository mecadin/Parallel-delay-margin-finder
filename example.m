% Delay Margin Finder for Neutral TDS

% This script calls the parDMF function, which returns the DM of the
% system: x'(t)-Cx'(t-tau)=Ax(t)+Bx(t-tau)

% Assumption 1: rho(C)<1
% Assumption 2: (I-C)^-1*(A+B) is Hurwitz

clc
close all
profile on


% System matrices
A = [0        1.0000;
     0       -0.3444];
B = [0        0;
    -34.4407 -4.8485];
C = [ 0         0;
         0   -0.8610];

     
% Delay margin approximation
g    = 2e-3;    % coarse grid step 
a    = 2e-7;    % accuracy for the solution

tic;
DM = parDMF(g,a,A,B,C);
t=toc;

% Display the results

fprintf('###      Delay margin: %10.7s',DM);
fprintf('\n###  Computation time: %10.5f',t);
fprintf('\n###          Accuracy: %10.1s\n\n',a);




   
