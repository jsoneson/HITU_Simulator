%% Authored by Joshua Soneson 2018
function[P1,P2] = BuildPade12operators(A,kk,dz,k,JJ)

I = speye(JJ);
kkk = k*kk;
A = A/kkk/kkk;
s = i*kkk*dz;
P1 = I + 0.25*(1-s)*A;
P2 = I + 0.25*(1+s)*A;
