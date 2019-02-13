%% Authored by Joshua Soneson 2018
function[u,U,Ppos,Pneg,I_td] = TDNL(u,U,X,KK,JJ,c,cutoff,Ppos,Pneg,I_td)
% converts spectrum to one cycle of the time-domain waveform and 
% integrates the invicid Burger's equation using upwind/downwind 
% method with periodic boundary conditions. TDNL stands for Time
% Domain NonLinear.
% set peak and trough values to zero; in case they are not assigned later
if(KK==1)	% linear case - do nothing
  for jj=1:JJ
    Ppos(jj) = abs(u(jj,1));
  end
  Pneg = -Ppos;
else		% nonlinear case - enter loop
  for jj=JJ:-1:1
    % execute nonlinear step only if amplitude is above cutoff;
    % row jj=1 is always computed so plots look nice
    I_td(jj) = 0;
    if(abs(u(jj,1)) < cutoff/20)	% if pressure is too low, skip nonlinear step
    else
      % convert from sin & cos representation to complex exponential 	
      U(2:KK+1) = conj(u(jj,:));
      U(2*KK:-1:KK+2) = u(jj,1:KK-1);
      U(1) = 0;
      % transform to time domain:
      U = KK*real(ifft(U));
      I_td(jj) = sum(U.^2);
      % determine how many steps necessary for CFL<1 (CFL<0.9 to be safe).
      PP = ceil(2*c*max(abs(U))/0.9);
      %if(j==1) 
      %  plot(U);
      %  hold on 
      %end
      % Nonlinear integration (upwind/downwind) algorithm 
      for pp=1:PP
        for kk=1:2*KK		
          if(U(kk)<0)
            if(kk==1)
              X(kk) = U(kk) + c*(U(kk)*U(kk) - U(2*KK)*U(2*KK))/PP;
            else
              X(kk) = U(kk) + c*(U(kk)*U(kk) - U(kk-1)*U(kk-1))/PP;
            end
          else
            if(kk==2*KK)
              X(kk) = U(kk) + c*(U(1)*U(1) - U(kk)*U(kk))/PP;
            else
              X(kk) = U(kk) + c*(U(kk+1)*U(kk+1) - U(kk)*U(kk))/PP;
            end
          end
        end
        U = X;
      end
      % account for nonlinear losses:
      I_td(jj) = I_td(jj) - sum(X.^2); 
      % transform back to frequency domain:
      Ppos(jj) = max(X);
      Pneg(jj) = min(X);
      X = fft(X)/KK;
      % convert back to sin & cos representation:
      u(jj,:) = conj(X(2:KK+1));
    end
  end
end
