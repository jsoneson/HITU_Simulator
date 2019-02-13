%% Authored by Joshua Soneson 2018
function[] = BHT(Grid,Layer,Q)
tic;
%% set path to include Factorize for shorter run times:
path(path,'/home/josh/Downloads/Factorize/')


% set heating and cooling durations:
t_h = 2;	% heating time (s)
t_c = 3;	% cooling time (s)

% set equilibrium temperature:
Teq = 37;	% celcius degrees

% Thermal damage thresholds:
safety = 80;	% dose (CEM43) where 0% cell necrosis occurs; safety threshold
efficacy = 240;	% dose at which 100% cell nectosis occurs

% temporal grid setup:
dt = min(0.01,t_h/100);		% estimate timestep dt
MM = ceil((t_h+t_c)/dt);	% find number of timesteps
t = linspace(0,t_h+t_c,MM);	
dt = t(2);			% refined dt
T_sp = zeros(MM,1);		% spatial peak temp as a function of time

% Get grid and operators:
[CN1,CN2,Hvec,Grid2] = BuildBHTperipherals(Grid,Layer,Q,dt);

% reporting:
fprintf('\t\tdt = %2.2f s\n',dt)
fprintf('\tNode count\n')
fprintf('\t\tAxial %d\n',Grid2.NN)
fprintf('\t\tRadial %d\n',Grid2.JJ)
fprintf('\t\tTemporal %d\n',MM)

Tvec = zeros(Grid2.JN,1);	% Temp values (overwritten at each timestep
Tvec_max = zeros(Grid2.JN,1);	% Max temp at each spatial location
Tmat = zeros(Grid2.JJ,Grid2.NN);
Tmat_max = zeros(Grid2.JJ,Grid2.NN);
Dvec = zeros(Grid2.JN,1);	% Dose values
Dmat = zeros(Grid2.JJ,Grid2.NN);

% Integration loop:
for mm=1:MM-1
  if t(mm)<t_h
    Tvec = CN1 \ (CN2*Tvec + dt*Hvec);
  else
    Tvec = CN1 \ (CN2*Tvec);
  end
  %Tmat = matrixize(Tvec,Tmat,Grid2.JJ,Grid2.NN);
  %mesh(Grid2.z,Grid2.r,Tmat)
  %pause(0.01)
  T_sp(mm+1) = max(Tvec);
  s = find(Tvec>Tvec_max);
  Tvec_max(s) = Tvec(s);
  Dvec = Dvec + dt*2.^(Tvec+Teq-43)/60;	% accrue thermal dose
end

% convert vectors to matrices for presentation
Tmat_max = matrixize(Tvec_max,Tmat_max,Grid2.JJ,Grid2.NN);
Dmat = matrixize(Dvec,Dmat,Grid2.JJ,Grid2.NN);

% plot routines:
figure
plot(t,T_sp+Teq,'LineWidth',2)
xlabel('Time (s)')
ylabel('Peak Temp (^{\circ}C)')
grid on

figure
Tmat_max = [Tmat_max(Grid2.JJ:-1:2,:);Tmat_max];
r = [-Grid2.r(Grid2.JJ:-1:2);Grid2.r];
imagesc(Grid2.z,r,Tmat_max+Teq);
xlabel('z (cm)')
ylabel('r (cm)')
colorbar
title('Max temp during simulation')
grid on

figure
Dmat = [Dmat(Grid2.JJ:-1:2,:);Dmat];
[c,h] = contour(Grid2.z,r,Dmat,[safety efficacy],'LineWidth',2);
clabel(c,h)
xlabel('z (cm)')
ylabel('r (cm)')
title('Thermal dose accumulation')
grid on

fprintf('\tTook %2.1f seconds.\n',toc)




  



