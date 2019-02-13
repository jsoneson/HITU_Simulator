%% Authored by Joshua Soneson 2007
function[v] = vektorize(v,A,J,K)
%% converts JxK matrix into a column-stacked J*Kx1 vector

for j=1:J
  v(K*(j-1)+[1:K]) = conj(A(j,:)');
end
