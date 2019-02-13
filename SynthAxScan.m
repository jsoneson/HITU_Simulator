%% Authored by Joshua Soneson 2018
function[p_r,p_c,p5] = SynthAxScan(r,p,b,JJ_)
% r = radial node vector (cm)
% p = pressure matrix (radial x harmonic spectrum)
% b = hydrophone element radius (cm)

% returns p_r and p_c, peak averaged rarefactional and compressional pressure
% and amplitude of first (up to) 5 averaged pressure harmonics 

[JJ,KK] = size(p);		% JJ = number of radial nodes; KK = number of harmonics
dr_min = r(2);			% mesh spacing near axis
p_h = zeros(KK,1);		% vector of spatially averaged pressure values
NN = max(4,ceil(10*b/dr_min));	% number of points over which to spatially average
dr = b/(NN-1);
x = linspace(0,b,NN);
for kk=1:KK
  q(:,kk) = interp1(r,p(:,kk),x,'linear');
  p_h(kk) = dr*trapz(q(:,kk).*x');
end
p_h = 2*p_h/b/b;
p5 = p_h(1:min(5,KK));

% determine peak compressional p_c and rarefactional p_r pressure
if(KK==1)        % linear case - do nothing
  p_c = abs(p_h(1));
  p_r = -p_c;
else		% nonlinear case - transform to time domain
  U(2:KK+1) = conj(p_h);
  U(2*KK:-1:KK+2) = p_h(1:KK-1);
  U(1) = 0;
  % transform to time domain:
  U = KK*real(ifft(U));
  p_r = min(U);
  p_c = max(U);
end


