%% Authored by Joshua Soneson 2018, %%%%%%%%%%%%%%%%%%%
function[Grid,Layer,Q] = WAKZK_planar()
%% Implementation of wide-angle parabolic method for axisymmetric HITU
%% beams.
tic;


%% Transducer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tx.f = 1e6;%2e6;		% frequency (Hz)
Tx.a1 = 0;%1;		% inner radius (cm)
Tx.a2 = 1;%3;		% outer radius (cm)
Tx.d = Inf;%5;		% geometric focus (cm) [= Inf if planar array]
Tx.P = 10;		% total acoustic power (W)

 

%% Computational Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Tx.d==Inf
  z_start = 0;
else
  z_start = Tx.d-sqrt(Tx.d^2-Tx.a2^2);	% axial location of equivalent source (cm)
end
Grid.Z = 10;%7;		% max axial distance (cm)
Grid.KK = 1;%128;		% max number of harmonics in simulation (use power of 2)
ppw_r = 15;		% grid resolution in r-direction; points per wavelength
ppw_z = 10;		% and z-direction


%% Spatial averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hd = 0;%0.2;		% diameter of hydrophone element (mm).  [hd = 0 no averaging] 


%% Graphical output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_output = [0 6.5];%[4.65 4.95 5.25];	% locations on z-axis where plots are produced
LL = length(z_output);	% code determines number of plot locations
ll = 1;			% initialize index


%% Layered media %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II = 1;%3;			% number of layers

% material 1 parameters:
Layer(1).z = 0;		% material transition distance (cm)	
Layer(1).c = 1482;	% small-signal sound speed (m/s)
Layer(1).rho = 1000;	% mass density (kg/m^3)
Layer(1).alpha = 0.217;	% attenuation at 1MHz (dB/m)
Layer(1).fraction = 0;	% fracion of attenuation due to absorption	
Layer(1).eta = 2;	% exponent in attenuation power law
Layer(1).beta = 3.5;	% nonlinear parameter
Layer(1).Cp = 4180;	% heat capacity (J/kg/K)
Layer(1).kappa = 0.6;	% thermal conductivity (W/m/K)
Layer(1).w = 0;		% perfusion rate (kg/m^3/s)

% material 2 parameters:
Layer(2).z = 4;		% material transition distance (cm)	
Layer(2).c = 1629;
Layer(2).rho = 1000;
Layer(2).alpha = 58;
Layer(2).fraction = 0.9;
Layer(2).eta = 1;
Layer(2).beta = 4.5;
Layer(2).Cp = 4180;
Layer(2).kappa = 0.6;
Layer(2).w = 20;

% material 3 parameters:
Layer(3).z = 6;		% material transition distance (cm)	
Layer(3).c = 1482;	
Layer(3).rho = 1000;	
Layer(3).alpha = 0.217;	
Layer(3).fraction = 0;
Layer(3).eta = 2;	
Layer(3).beta = 3.5;	
Layer(3).Cp = 4180;	
Layer(3).kappa = 0.6;	
Layer(3).w = 0;	

% User may add more Layers if necessary: Layer(4), Layer(5), etc..
% Be sure to change the value of II above.
%
%
%
%
%















% calculate wavenumber in each layer:
for ii=1:II
  Layer(ii).k = 2*pi*Tx.f/(100*Layer(ii).c);	% wavenumber (cm^-1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% More computational domain stuff %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Tx.d==Inf
  a2p = Tx.a2;
else
  a2p = Tx.a2*(Tx.d-z_start)/sqrt(Tx.d^2-Tx.a2^2);	% radius of equivalent source (a_2')
end
w = 1.5*a2p;%1.05*a2p;		% width (radius) of physical domain (Tx radius + 5% pad)
lambda = 2*pi/Layer(1).k;	% wavelength (cm)
th = 20*lambda;%2*lambda;			% PML thickness
Grid.R = w + th;		% max radius of computational domain (cm)
Grid.JJ = ceil(ppw_r*Grid.KK^0.35*Grid.R/lambda);	% Gridpoints in r-dir
Grid.NN = ceil(ppw_z*(Grid.Z-z_start)/lambda); % Gridpoints in z-direction
% node vectors
Grid.r = linspace(0,Grid.R,Grid.JJ)';
Grid.z = linspace(z_start,Grid.Z,Grid.NN);
dr = Grid.r(2);
dz = Grid.z(2)-Grid.z(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equivalent source %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This source is a converging spherical wave at z=z_start bounded by a2p. 
% Note: it's dimensionless:
if Tx.d==Inf
  A = ones(Grid.JJ,1).*Grid.r<a2p;
  a1p = Tx.a1;
else
  A = Tx.d*exp(-i*Layer(1).k*sqrt(Grid.r.^2+(Tx.d-z_start)^2)).*(Grid.r<a2p) ...
     ./sqrt(Grid.r.^2+(Tx.d-z_start)^2);
  if(Tx.a1 ~= 0)		% if there's a hole in the Tx
    a1p = Tx.a1*(Tx.d-z_start)/sqrt(Tx.d^2-Tx.a1^2);	
    A = A.*(Grid.r>a1p);
  end
end



% The user could specify a custom source here [any complex function A=A(r)].
%A = ;



% Apply a low-pass filter to the source:
A = SourceFilterH(Grid.r,A,Layer(1).k);

% Next scale the source by the appropriate pressure coefficient so that it has
% the proper total acoustic power: 
integral = 2*pi*dr*trapz(abs(A).^2.*Grid.r);
integral = 1e-4*integral;		% convert units from cm^2 to m^2
p0 = sqrt(2*Layer(1).rho*Layer(1).c*Tx.P/integral);
A = p0*A;	% dimensionalize the boundary condition (units of Pa)



%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial averaging %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hr = 0.1*hd/2;		% convert to radius in cm 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate attenuation, dispersion %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = Tx.f*[1:Grid.KK]'/1e6;		% vector of frequencies
for ii=1:II
  Layer(ii).alpha = Layer(ii).alpha/8.686/100;% convert to Np/cm
  if Layer(ii).eta == 1                                 % linear media
    Layer(ii).alpha = Layer(ii).alpha*v.*(1+2*i*log(v)/pi);
  else                                                  % everything else
    Layer(ii).alpha = Layer(ii).alpha*(v.^Layer(ii).eta ...
                    - i*tan(pi*Layer(ii).eta/2)*(v.^Layer(ii).eta - v));
end; end


%% some reporting
fprintf('\n\tWavelength = %2.2f mm\n',10*lambda)
fprintf('\tNode count\n')
fprintf('\t\tAxial %d\n',Grid.NN)
fprintf('\t\tRadial %d\n',Grid.JJ)
fprintf('\tGrid stepsize\n')
fprintf('\t\tdz = %2.2f mm\n',10*dz)
fprintf('\t\tdr = %2.2f mm\n',10*dr)



%% dependent variable (pressure) matrices:
p = zeros(Grid.JJ,Grid.KK);	% new pressure (n+1 th step)
q = zeros(Grid.JJ,Grid.KK);	% old pressure (n th step)
q(:,1) = A;			% apply boundary condition (source)

dr_max = Grid.r(Grid.JJ)-Grid.r(Grid.JJ-1);         % mesh spacing near PML
JJ_ = Grid.JJ - ceil(hr/dr_max);      % Index of radial limit where spatial averaging occurs
p_r = zeros(1,Grid.NN);	% peak rarefactional pressure
p_c = zeros(1,Grid.NN);	% peak compressional pressure
[p_r(1),p_c(1),p5(:,1)] = SynthAxScan(Grid.r,q,hr,JJ_);
I = zeros(Grid.JJ,Grid.NN);	% intensity
Q = zeros(Grid.JJ,Grid.NN);	% power density




%% PML stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
phi = 2*(Grid.r-Grid.R+th)/th.*(Grid.r>Grid.R-th);
phi = phi+(1-phi).*(Grid.r>Grid.R-th/2);
u = exp(-i*pi*phi/4);
Du2 = zeros(Grid.JJ,1);
Du2 = (-i*pi/2/th).*u.*(Grid.r>Grid.R);
Du2 = Du2.*(Grid.r<Grid.R+th/2);
ur = sparse(diag(u./Grid.r));
ur(1,1) = 0;
u = sparse(diag(u));
Du = sparse(diag(Du2));

%% Build transverse Laplacian operator w/PML:
clear A;
e = ones(Grid.JJ,1);
D1 = spdiags([-e e],[-1 1],Grid.JJ,Grid.JJ)/2/dr;
D2 = spdiags([e -2*e e],-1:1,Grid.JJ,Grid.JJ)/dr/dr;
A = u*((ur+Du)*D1 + u*D2);		
A(1,2) = 2*A(1,2);			% zero flux BC at r=0;


%% peripherals for nonlinear integrator:
Ppos = zeros(Grid.JJ,Grid.NN);
Ppos(:,1) = abs(q(:,1));
Pneg = zeros(Grid.JJ,Grid.NN);
Pneg(:,1) = -abs(q(:,1));
p5 = zeros(min(Grid.KK,5),Grid.NN);
p5(1,1) = abs(q(1,1));
w = zeros(Grid.NN,2*Grid.KK);                % waveform data vectors
Y = zeros(1,2*Grid.KK);
I_td = zeros(Grid.JJ,2);			% change in intensity
dt = 1/Tx.f/(2*Grid.KK-1);			% in seconds
t = 1e6*[0:dt:1/Tx.f];				% in us

%% more reporting:
fprintf('\t\tdt = %2.2f ns\n',1e9*dt)


%% find indices of first Gridpoint in each Layer:
Layer(1).index = 1;
for ii=2:II
  Layer(ii).index = ceil((Layer(ii).z-z_start)/dz)+1;
end
Layer(II+1).index = Grid.NN;



%% integration loop:
for ii=1:II
  %build operators for Layer ii:
  for kk=1:Grid.KK		
    %[M(kk).P1,M(kk).P2,M(kk).P3] = ...
    %  BuildPade12operators(A,kk,dz,Layer(ii).k,Grid.JJ);
    [M(kk).P1,M(kk).P2] = ...
      BuildPade11operators(A,kk,dz,Layer(ii).k,Grid.JJ);
  end
  mu = (Layer(ii).beta/2/Layer(ii).rho/Layer(ii).c^3)*(0.01*dz/dt);
  cutoff = Layer(ii).alpha(1)*Layer(ii).rho*Layer(ii).c^2 ...
         / Layer(ii).beta/Layer(ii).k; % cutoff for nonlinearity
  for nn=Layer(ii).index:Layer(ii+1).index-1
    %integrate nonlinear term:	  
    [p,w(nn+1,:),Ppos(:,nn+1),Pneg(:,nn+1),I_td(:,1)] = ...
      TDNL(q,w(nn+1,:),Y,Grid.KK,Grid.JJ,mu,cutoff,Ppos(:,nn),Pneg(:,nn),I_td(:,1));	
    %attenuation/dispersion term and diffraction term:
    for kk=1:Grid.KK
      p(:,kk) = p(:,kk).*exp(-Layer(ii).alpha(kk)*dz);
      %p(:,kk) = M(kk).P1 \ (M(kk).P2 \ (M(kk).P3*p(:,kk)));	% for Pade 12
      p(:,kk) = M(kk).P1 \ (M(kk).P2*p(:,kk));			% for Pade 11
    end
    Norm = norm(p(:,1));
    if(Norm==NaN | Norm==Inf)	% stop if something goes wrong
      fprintf('\tNaN or Inf detected! Simulation stopped at z = %d cm.\n',Grid.z(nn))
    end
    q = p;
    for jj=1:Grid.JJ	% calculate intensity I and power density H
      I(jj,nn+1) = sum(abs(p(jj,:)).^2)/2/Layer(ii).rho/Layer(ii).c;
      Q(jj,nn+1) = (sum(Layer(ii).fraction*real(Layer(ii).alpha).*abs(p(jj,:)').^2) ...
                 + sum(I_td(jj,:))/dz/(2*Grid.KK-1))/Layer(ii).rho/Layer(ii).c;
    end
    %collect/process data:
    if hd==0;
      p5(:,nn+1) = p(1,1:min(Grid.KK,5))';
    else
      [p_r(nn+1),p_c(nn+1),p5(:,nn+1)] = SynthAxScan(Grid.r,p,hr,JJ_);
    end
    if ll<=LL
      if Grid.z(nn+1)>z_output(ll) & Grid.z(nn)<=z_output(ll)  % find special output locations
        if hd==0
          SpecOut(ll).pr = Pneg(:,nn+1);
          SpecOut(ll).pc = Ppos(:,nn+1);
          SpecOut(ll).p5 = abs(p(:,1:min(5,Grid.KK)));
          SpecOut(ll).w = w(nn+1,:);		
        else
          SpecOut(ll) = SynthRadScan(Grid.r,p,hr,JJ_);
        end
        SpecOut(ll).I = I(:,nn+1);
        ll = ll + 1;
      end
    end
  end
  %rescale pressure due to transmission at interface between layers ii and ii+1:
  if ii < II
    q = 2*Layer(ii+1).rho*Layer(ii+1).c*q ...
      / (Layer(ii).rho*Layer(ii).c + Layer(ii+1).rho*Layer(ii+1).c);
  end
end
fprintf('\tTook %2.1f minutes.\n\n',toc/60)

%LinearHeating(Layer(2),Grid,I);	% assumes a focused beam and that
					% the focus is in layer 2
    
Layer = Layer(1:II);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot routine %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure	% 3d plot showing axial and radial pressure amplitudes
r_ones = ones(1,length(SpecOut(1).p5));
z_zeros = zeros(1,length(Grid.z));
plot3(Grid.z,z_zeros,abs(p5)/1e6,'LineWidth',2)
hold on
for ll=1:LL
  plot3(z_output(ll)*r_ones,Grid.r(1:length(SpecOut(ll).p5)),SpecOut(ll).p5/1e6,'LineWidth',2)
end
xlabel('z (cm)')
ylabel('r (cm)')
zlabel('|p| (MPa)')
grid on

figure	% axial plots of amplitude of first 5 harmonics and intensity
subplot(211)
plot(Grid.z,abs(p5)/1e6,'LineWidth',2)
ylabel('|p| (MPa)')
grid on
subplot(212)
plot(Grid.z,I(1,:)/1e4,'LineWidth',2)
xlabel('z (cm)')
ylabel('I (W/cm^2)')
grid on

for ll=1:LL		% build plot label
  V(ll).str = strcat('z= ',num2str(z_output(ll)),' cm');
end

if(Grid.KK>1)
  figure	% temporal waveforms at specified axial locations
  hold on
  for ll=1:LL
    plot(t,SpecOut(ll).w/1e6,'LineWidth',2)
  end
  xlim([t(1) t(length(t))])
  legend(V.str)
  grid on
  xlabel('t (us)')
  ylabel('p (MPa)')
  figure	% axial plots of compressional and rarefactional pressure
  if(hd==0)
    p_c = Ppos(1,:);
    p_r = Pneg(1,:);
  end
  plot(Grid.z,p_c/1e6,Grid.z,p_r/1e6,'LineWidth',2)
  xlabel('z (cm)')
  ylabel('p (MPa)')
  grid on
end

for ll=1:LL
  figure	% radial plots of amplitude of first 5 harmonics and intensity
  subplot(211)	% at specified axial locations
  plot(Grid.r(1:length(SpecOut(ll).p5)),SpecOut(ll).p5/1e6,'LineWidth',2)
  ylabel('|p| (MPa)')
  title(V(ll).str)
  grid on
  subplot(212)
  plot(Grid.r,SpecOut(ll).I/1e4,'LineWidth',2)
  xlabel('r (cm)')
  ylabel('I (W/cm^2)')
  grid on
end

figure	% spatial distribution of field emphasizing low-amplitude variations
r = [-Grid.r(Grid.JJ:-1:2);Grid.r];
I = [I(Grid.JJ:-1:2,:);I];
imagesc(Grid.z,r,I.^0.2)
xlabel('z (cm)')
ylabel('r (cm)')
grid
