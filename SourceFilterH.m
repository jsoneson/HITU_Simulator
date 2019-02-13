%% Authored by Joshua Soneson 2018
function[Af] = SourceFilterH(x,A,k)

% get some vital parameters:
JJ = length(x);
R = x(JJ);

% find zeros of BesselJ(0,r):
c = zeros(1,JJ+1);
for jj=1:JJ+1
  y = pi*(4*jj-1)/4;				% first guess based on asymptotic approximation
  for ii=1:5
    c(jj) = y + besselj(0,y)/besselj(1,y);      % Newton iteration
    if abs(c(jj)-y) < eps break;                % check for convergence
    else y = c(jj); end
  end
end

V = c(JJ+1)/(2*pi*R);    % Maximum spatial frequency
r = c(1:JJ)'*R/c(JJ+1);  % vector of radial nodes (nonuniform)
v = c(1:JJ)'/(2*pi*R);   % vector of spatial frequencies

[Jn,Jm] = meshgrid(c(1:JJ),c(1:JJ));
C = (2/c(JJ+1))*besselj(0,Jn.*Jm/c(JJ+1)) ...
  ./ (abs(besselj(1,Jn)).*abs(besselj(1,Jm)));

m1 = (abs(besselj(1,c(1:JJ)))/R)';   
m2 = m1*R/V;                            
clear Jn
clear Jm

% filter:
q = 40;
s = 1.15;
F = (1-tanh(q*(v/k-s/2/pi)))/2;

% Look at filter function:
%figure
%plot(2*pi*v/k,F)

% perform transforms and filtering:
Ahat = C*(A./m1);		% Apply Hankel transform
Ahatf = F.*Ahat;		% apply filter 
Af = (C*Ahatf).*m1;		% apply inverse Hankel transform

