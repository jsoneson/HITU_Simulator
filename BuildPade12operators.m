%% Authored by Joshua Soneson 2018
function[P1,P2,P3] = BuildPade12operators(A,kk,dz,k,JJ)

I = speye(JJ);
kkk = k*kk;
A = A/kkk/kkk;
s = i*kkk*dz;
muplus = (3-2*s*s+i*sqrt((((2*s+6)*s-6)*s-18)*s-9))/12/(1+s);
muminus = (3-2*s*s-i*sqrt((((2*s+6)*s-6)*s-18)*s-9))/12/(1+s);
epsilon = ((s+3)*s+3)/6/(1+s);
P1 = I + muplus*A;
P2 = I + muminus*A;
P3 = I + epsilon*A;
