%% Authored by Joshua Soneson 2018
function[] = LinearHeating(Layer,Grid,I)

r = Grid.r;
z = Grid.z;
magicNNindex = find(max(I)==max(max(I)));
magicJJindex = find(max(I')==max(max(I)));
if(magicJJindex>1)
  fprintf('\tLinear heating failed--maximum intensity does not occur on central axis.\n')
end
I0 = I(1,magicNNindex);
a = interp1(I(:,magicNNindex),r,exp(-1)*I0)/100;		% effective beam radius (m)
z1 = interp1(I(1,1:magicNNindex),z(1:magicNNindex),exp(-1)*I0)/100;
z2 = interp1(I(1,magicNNindex:length(z)),z(magicNNindex:length(z)),exp(-1)*I0)/100;
b = (z2-z1)/2;
c1 = 100*real(Layer.alpha(1))*Layer.fraction*I0*a*a*b/Layer.kappa/sqrt(a*a-b*b);
c2 = b*b/(a*a-b*b);
c3 = 4*Layer.kappa/Layer.rho/Layer.Cp/(a*a-b*b);

Teq = 37;
t = linspace(0,1,30);
T = Teq + c1*(atan(sqrt(c2))-atan(sqrt(c2+c3*t)));
figure
plot(t,T,'LineWidth',2)
xlabel('t (sec)')
ylabel('\Delta T (^{\circ}C)')
title('Linear Heating Estimate')
grid on


