%% Authored by Joshua Soneson 2018
function[CN1,CN2,Hvec,Grid2] = BuildBHTperipherals(Grid,Layer,Q,dt)


% heating grid:
rskip = 2;		% heat grid is twice as coarse as propagation grid in r
zskip = 3;		% three times as coarse in z
Grid2.r = Grid.r(1:rskip:length(Grid.r));
dr = Grid2.r(2);
Grid2.JJ = length(Grid2.r);
Grid2.z = Grid.z(1:zskip:length(Grid.z));
dz = Grid2.z(2)-Grid2.z(1);
Grid2.NN = length(Grid2.z);
Grid2.JN = Grid2.JJ*Grid2.NN;
Q2 = Q(1:rskip:length(Grid.r),1:zskip:length(Grid.z));
Qvec = zeros(Grid2.JN,1);
Qvec = vektorize(Qvec,Q2,Grid2.JJ,Grid2.NN);

II = length(Layer);

% reporting
fprintf('\tGrid stepsize\n')
fprintf('\t\tdz = %2.2f mm\n',10*dz)
fprintf('\t\tdr = %2.2f mm\n',10*dr)

% compute diffusivity d, the reciprocal of the perfusion time constant v,
% and  the coefficient for the power density for all the tissue layers:
for ii=1:II
  Layer(ii).d = 1e4*Layer(ii).kappa/Layer(ii).Cp/Layer(ii).rho;	% units cm^2/s
  Layer(ii).v = Layer(ii).w/Layer(ii).rho;	% units 1/s
  Layer(ii).coef = 1e2/Layer(ii).Cp/Layer(ii).rho; % units K cm s^2/kg
end


% build matrix operator's vector "bands" 
alpha_plus = ones(Grid2.JN,1)/dz/dz;
alpha_minus = ones(Grid2.JN,1)/dz/dz;
bp = zeros(Grid2.JJ,1);
bm = zeros(Grid2.JJ,1);
bp(1) = 2/dr/dr;
bp(2:Grid2.JJ) = (1/dr + 0.5./Grid2.r(2:Grid2.JJ))/dr;
bm(2:Grid2.JJ) = (1/dr - 0.5./Grid2.r(2:Grid2.JJ))/dr;
beta_plus = zeros(Grid2.JN,1);
beta_minus = zeros(Grid2.JN,1);
for jn=1:Grid2.JN
  if(mod(jn,Grid2.NN)==0) alpha_plus(jn) = 0; end
  if(mod(jn,Grid2.NN)==1) alpha_minus(jn) = 0; end
  qq = ceil(jn/Grid2.NN);
  beta_plus(jn) = bp(qq);
  beta_minus(jn) = bm(qq);
end
gamma = -2*(1/dr/dr + 1/dz/dz)*ones(Grid2.JN,1);



% build diffusivity and perfusion coefficient matrices (describe layers)
for ii=1:II
  z_t(ii) = Layer(ii).z;
end
if z_t(1) < Grid2.z(1)
  ii = 1;
else
  ii = 0;
end
d = zeros(1,Grid2.NN);
for nn=1:Grid2.NN
  if min(abs(Grid2.z(nn)-z_t)) < dz/2
    ii = ii + 1;
  end
  d(nn) = Layer(ii).d;
  v(nn) = Layer(ii).v;
  c(nn) = Layer(ii).coef;
end
D = d;
V = v;
C = c;
for jj=1:Grid2.JJ-1
  D = [D d];
  V = [V v];
  C = [C c];
end
D = D';
V = V';
C = C';
D = spdiags(D,0,Grid2.JN,Grid2.JN);
V = spdiags(V,0,Grid2.JN,Grid2.JN);



% put it all together
A = spdiags([beta_plus alpha_plus gamma alpha_minus beta_minus],[-Grid2.NN -1 0 1 Grid2.NN],Grid2.JN,Grid2.JN);
A = A';
I = speye(Grid2.JN);
%CN1 = I - dt*(D*A - V)/2;	% Crank-Nicolson operators w/o Factorize (slow)
CN1 = factorize(I - dt*(D*A - V)/2);	% Crank-Nicolson operators
CN2 = I + dt*(D*A - V)/2;

Hvec = C.*Qvec;			% rescale power density to obtain heating rate
