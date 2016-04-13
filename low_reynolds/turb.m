% example of how to read DNS data.
%
% Channel data at Re_delta=7890, Re_tau=395, Re_theta=700.
% All quantites are normalized by u_tau and nu unless stated otherwise; delta denotes
% the half-channel width.  Data compiled from an unpublished work of Kim (1989).
%
% See also Kim, Moin & Moser (1987, JFM, vol 177) and Mansour, Kim & Moin (1988, JFM, vol 194).
%
%
%close all
clear
clc
format long

%Defining simulation constants
c1 = 1.45;
c2 = 2.0;
cMu = 0.09;
sigmaK = 1.00;
sigmaEps = 1.30;
visc=1/395;
urC = 0.7;
BCU = [2 2];
BCk = [2 2];
BCeps = [2 2];
% for the new model:
% at the wall, both eps and k should be kept at zero
% at the middle, both gradients should be zero

% wall friction velocity
ustar=1;
% read u from DNS data base
load u_dns.dat

% read y from DNS data base
 load y_dns.dat
yc=y_dns;          %yc is a vector contains the faces coordinates (grid)
nj=length(yc)+1; % nj: number of node points  
y_node = zeros(nj,1);
u_dns(nj)=u_dns(nj-1);
y_node(1)=(yc(1)-yc(2))/2;
for j=2:nj-1   
    y_node(j)= (yc(j)+yc(j-1))/2;
end
y_node(nj)=yc(nj-1)+(yc(nj-1)-yc(nj-2))/2;

cMu = ones(nj,1) * cMu;

%Calculating dY and deltaY
dY(:,1) = [0 diff(y_node(2:end)') 0]';
dY(:,2) = [0 diff(y_node(1:end-1)') 0]';
deltaY = [0 diff(yc') 0]';


% init velocity, vist, k & eps
U = zeros(nj,1);
k = zeros(nj,1);
eps = zeros(nj,1);
vist = zeros(nj,1);
dudy = zeros(nj,1);
d2udy2 = zeros(nj,1);
dsqrtkdy = zeros(nj,1);
U(1)=0;
k(1)=0;
for j=2:nj-1
   U(j)=1;
   k(j)=10^(-5);
   eps(j)=10^(-5);
   vist(j)=10^(-5);
end
k(nj)=k(nj-1);
U(nj)=U(nj-1);
eps(nj)=eps(nj-1);

kappa=0.41;
error=1000;
count=0;
max_error=0.001;
urf=0.8;

disp('Go!')
UStore = [];
kStore = [];
epsStore = [];

%%
old_error = error;

%while error > max_error 
while count < 1
   
    count = count+1;
    % implementing boundary conditions
    U(end) = U(end-1);
    k(end) = k(end-1);
    k(1) = 0;
    eps(end) = eps(end-1);
    eps(1) = 0;
    
    % Compute the gradients du/dy, d^2u/dy^2 and d(sqr(k))/dy
    % by fitting a quadratic function to three points and then
    % using the derivative of the quadratic function
   for j=2:nj-1
      %dudy(j)= (U(j+1) - U(j-1)) / (dY(j,1) + dY(j,2));
      
      dN = dY(j,1);
      dS = dY(j,2);
      factor = dS^2/(2*dN*dS + dN^2);
      
      b = (U(j) - U(j-1) + ((U(j) - U(j+1))*factor))/(dS - dN*factor);
      a = -(U(j) - U(j+1) + dN*b)*factor/dS^2;

      dudy(j) = 2*a*dS + b;

      b = (dudy(j) - dudy(j-1) + ((dudy(j) - dudy(j+1))*factor))/(dS - dN*factor);
      a = -(dudy(j) - dudy(j+1) + dN*b)*factor/dS^2;
      
      d2udy2(j) = 2*a*dS + b;

      b = (sqrt(k(j)) - sqrt(k(j-1)) + ((sqrt(k(j)) - sqrt(k(j+1)))*factor))/(dS - dN*factor);
      a = -(sqrt(k(j)) - sqrt(k(j+1)) + dN*b)*factor/dS^2;
      
      dsqrtkdy(j) = 2*a*dS + b;      
   end
   
   %Calculating cMu and c2
   Rt = k.^2 ./ (visc .* eps);
   cMu = 0.09  .* exp(-2.5 ./ (1 + 0.02*Rt));
   c2 = 2.0 .* (1 - 0.3 .* exp(-Rt.^2));
    
%
%
%  Often it can be tricky to start the simulations. They often diverge.
%  It can then be useful to compute the turbulent viscosity from the 
%  mixing-length model for the first 2000 iterations.
%  In this way the K and EPS (or OMEGA) eqns. are de-coupled from U 
%  for the first 2000 iterations, i.e. the K and EPS (or OMEGA) do not
%  influence U since the viscosity is taken from the mixing-length model.

   vist_old=vist;

   
   if count < 0   % use mixing length model for turbulent viscosity if count >2000
      for j=2:nj-1
% compute turbulent viscosity
         yplus=ustar*y_node(j)/visc;
         damp=1-exp(-yplus/26);
         ell=min(damp*kappa*y_node(j),0.09);
         vist(j)=urf*abs(dudy(j))*ell^2+(1-urf)*vist_old(j);
      end
   else
       vist(2:nj-1) =  cMu(2:nj-1) .* (k(2:nj-1).^2)./eps(2:nj-1);
       vist(1) = vist(2);
       vist(end) = vist(end-1);
   end
     
   %Calculating source terms
   uSp = zeros(nj,1);
   uSu = ones(nj,1) .* deltaY;
   
   kSp = (-deltaY./k) .* (eps + dsqrtkdy.^2);
   kSu = deltaY .* vist .* dudy.^2;

   epsSp = -deltaY .* c2 .* eps ./ k;
   epsSu = deltaY .* vist .* (c1 .* eps .* dudy.^2 ./ k + 2*visc .* d2udy2.^2); 
 

   %Calculating coefficients
   UCoeff = CalcCoeffs( 1, dY, deltaY, visc, vist, uSp, nj, BCU);
   kCoeff = CalcCoeffs( sigmaK, dY, deltaY, visc, vist, kSp, nj, BCk);
   epsCoeff = CalcCoeffs( sigmaEps, dY, deltaY, visc, vist, epsSp, nj,BCeps);
   
   %Gauss-Seidel iteration
   U_new = GaussSeidel(U,uSu,UCoeff);
   k_new = GaussSeidel(k,kSu,kCoeff);
   eps_new = GaussSeidel(eps,epsSu,epsCoeff);
   
   U_new(end,:) = U_new(end-1,:);
   k_new(end,:) = k_new(end-1,:);
   eps_new(end,:) = eps_new(end-1,:);
   
   %disp([epsCoeff.point, epsCoeff.south,epsCoeff.north])
   
% after having computed ap and su, use under-relaxation (see lecture notes) 
%  Compute the velocity U
   
%  compute U
    U(2:nj-1) = U(2:nj-1) + urC.*(U_new(2:nj-1) - U(2:nj-1));
   

%  Compute the turbulent kinetic energy k
%  compute k  
     k(2:nj-1) = k(2:nj-1) + urC*(k_new(2:nj-1) - k(2:nj-1));

%  Compute the turbulent dissipation epsilon
%  compute epsi
     eps(2:nj-1) = eps(2:nj-1) + urC*(eps_new(2:nj-1) - eps(2:nj-1));
  
% Convergence criterian (Check the error)
% compute residuals R
     R = ComputeResidual(U,UCoeff,uSu,k,kCoeff,kSu,eps,epsCoeff,epsSu,nj);
% compute the flux F
     F = ComputeFlux(U,deltaY,nj);

  error = abs(R/F);
  if(mod(count,10) == 0)
      UStore = [UStore U];
      kStore = [kStore k];
      epsStore = [epsStore eps];
      
%       if(error > old_error)
%           urC = max(0.1*urC , 0.1);
%       else
%           urC = min(1.01*urC , 0.99999);
%       end
      fprintf('Iteration: %d, %e,  %f \n',count,error, urC);
      old_error = error;
  end
  
%   disp([R,F,max(diff(vist))])
%   disp([ max(UCoeff.point) max(kCoeff.point) max(epsCoeff.point)])
%   disp(max(epsSp))
%   disp(' ')
  
  
end  %while
%
% plot
% compare with DNS
%
% plot k

%%


figure(1)

hold on
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat

k_dns=0.5*(u2_dns+v2_dns+w2_dns);
plot(y_node,k,'rx')
plot(y_dns,k_dns,'bo')
xlabel('x')
ylabel('turbulent kinetic energy, k')
legend('Calc. k','DNS')
print k.ps -deps

% plot epsi
figure(2)
hold on
load dns_data.dat

% eps is normalized by ustar^4/viscos
% your computed eps is normalized with ustar^3/delta=1
%
eps_dns=dns_data(:,2)*ustar^4/visc;
plot(y_node,eps,'rx');
plot(y_dns,eps_dns,'bo')
xlabel('x')
ylabel('dissipation of k')
legend('Calc. eps','DNS')
axis([0 1 0 90])
print eps.ps -deps

% plot shear stress
figure(3)
hold on
load u_dns.dat
plot(y_node,U,'rx')
plot(y_dns,u_dns,'bo')

xlabel('x')
ylabel('U')
legend('Calc.','DNS')
print u.ps -deps

% Compare also with the different terms in the k-eq. 
% Read DNS data from file 'dns_data.dat'
%
% 6 coulumns as below:
%
%      y+         Diss        prod     vel_p_grad   Turb_diff   Visc_diff
%
% Please note that all terms are normalized by ustar^4/viscos
%


%close all


% %
% 
% figure(1)
% contourf(UStore)
% contourf(UStore(1:80,:))
% colorbar
% figure(2)
% contourf(kStore)
% contourf(kStore(1:80,:))
% colorbar
% figure(3)
% contourf(epsStore)
% contourf(epsStore(1:80,:))
% colorbar

