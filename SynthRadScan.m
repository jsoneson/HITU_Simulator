%% Authored by Joshua Soneson 2018
function[SpecOut] = SynthRadScan(r,p,b,JJ_)
% r = radial node vector (cm)
% p = pressure matrix (radial x harmonic spectrum)
% b = hydrophone element radius (cm)

% returns SpecOut, a structure containing p_r and p_c, vectors of averaged 
% peak rarefactional and compressional pressure and first (up to) 5 
% averaged harmonic pressure amplitudes, all as a function of radius.
% Also contains the averaged waveform on axis. 

[JJ,KK] = size(p);		% JJ = number of radial nodes; KK = number of harmonics
dr_min = r(2);			% mesh spacing near axis
p_h = zeros(JJ_,KK);		% matrix of spatially averaged pressure values
NN = max(4,ceil(10*b/dr_min));		% number of points over which to spatially average
dr = b/(NN-1);
x = linspace(0,b,NN);
for kk=1:KK
  q(:,kk) = interp1(r,p(:,kk),x,'linear');
  p_h(1,kk) = dr*trapz(q(:,kk).*x');
end

for jj=2:JJ_
  if(r(jj)<b)		% if element overlays central axis
    lowerlimit = 0;	% setup for "inner circle"
    upperlimit = b-r(jj);
    NN = ceil(10*(upperlimit-lowerlimit)/dr_min);
    dr = (upperlimit-lowerlimit)/(NN-1);
    x = linspace(lowerlimit,upperlimit,NN);
    clear q;
    for kk=1:KK
      q(:,kk) = interp1(r,p(:,kk),x,'linear');
      p_h(jj,kk) = dr*trapz(conj(q(:,kk)').*x);
    end
    lowerlimit = b-r(jj);	% setup for outer crescent
    upperlimit = r(jj)+b;
    NN = ceil(10*(upperlimit-lowerlimit)/dr_min);
    dr = (upperlimit-lowerlimit)/(NN-1);
    x = linspace(lowerlimit,upperlimit,NN);
    clear q;
    for kk=1:KK
      q(:,kk) = interp1(r,p(:,kk),x,'linear');
      p_h(jj,kk) = p_h(jj,kk) + dr*trapz(conj(q(:,kk)').*weight(x,r(jj),b).*x);
    end
  else
    lowerlimit = r(jj)-b;
    upperlimit = r(jj)+b;
    NN = ceil(10*(upperlimit-lowerlimit)/dr_min);
    dr = (upperlimit-lowerlimit)/(NN-1);
    x = linspace(lowerlimit,upperlimit,NN);
    clear q;
    for kk=1:KK
      q(:,kk) = interp1(r,p(:,kk),x,'linear');
      W = weight(x,r(jj),b);
      p_h(jj,kk) = dr*trapz(conj(q(:,kk)').*W.*x);
    end
  end
end
p_h = 2*p_h/b/b;
p5 = abs(p_h(:,1:min(5,KK)));

% determine peak compressional p_c and rarefactional p_r pressure
if(KK==1)        % linear case - do nothing
  for jj=1:JJ_
    p_c(jj) = abs(p_h(jj,1));
  end
  p_r = -p_c;
else		% nonlinear case - transform to time domain
  %for jj=1:JJ_
  for jj=JJ_:-1:1
    U(2:KK+1) = conj(p_h(jj,:));
    U(2*KK:-1:KK+2) = p_h(jj,1:KK-1);
    U(1) = 0;
    % transform to time domain:
    U = KK*real(ifft(U));
    p_r(jj) = min(U);
    p_c(jj) = max(U);
  end
  SpecOut.w = U;
end
SpecOut.pr = p_r;
SpecOut.pc = p_c;
SpecOut.p5 = p5;
SpecOut.I = p_r;	% placeholder; intensity is assigned in WAKZK.m
%figure
%plot(r(1:JJ_),SpecOut.p5,'-o')

function[x] = weight(r,r0,b)
  x = real(acos((r.*r + r0*r0 - b*b)./(2*r0*r))/pi);
  x(isnan(x)) = 0; 

