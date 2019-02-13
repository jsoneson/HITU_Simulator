%% Authored by Joshua Soneson 2007
function[A] = matrixize(v,A,J,K)
%% Converts J*Kx1 vector into JxK matrix

v=conj(v');
for j=1:J
  A(j,1:K) = v(K*(j-1)+[1:K]);
end
