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
cMuTilde = 0.09;
sigmaR = 1;     % just guessing...
c1Eps = 1.44;
c2Eps = 1.92;
visc=1/395;
C1 = 0.48;
C2 = 0.08;
urC = 0.1;
BCU = [2 2];
BCR = [2 2];
BCDamping = [2 2];  % Not sure...


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
y_node(1)=yc(1);
for j=2:nj-1   
    y_node(j)= (yc(j)+yc(j-1))/2;
end
y_node(nj)=yc(nj-1);

%Calculating dY and deltaY
dY(:,1) = [0 diff(y_node(2:end)') 0]';
dY(:,2) = [0 diff(y_node(1:end-1)') 0]';
deltaY = [1 diff(yc') 1]';


% init velocity, vist, k & eps
U = zeros(nj,1);
k = zeros(nj,1);
eps = zeros(nj,1);
R = zeros(nj,1);
damping = 0.5 * ones(nj,1);
vist = zeros(nj,1);
dudy = zeros(nj,1);
dRtildedy = zeros(nj,1);
U(1)=0;
k(1)=0;
for j=2:nj-1
   U(j)=1;
   k(j)=10^(-5);
   eps(j)=10^(-5);
   vist(j)=10^(-5);
   R(j) = k(j)^2/eps(j);
end
k(nj)=k(nj-1);
U(nj)=U(nj-1);
eps(nj)=eps(nj-1);
eps(1) = eps(2);
R(nj) = R(nj-1);
damping(1) = 0;
damping(nj) = damping(nj-1);


kappa=0.41;
error=1;
count=0;
max_error=0.001;
urf=0.8;

disp('Go!')
UStore = [];
kStore = [];
epsStore = [];

%%
old_error = error;
while error > max_error 
%while count < 200
    
    count = count+1;
    
    % Compute the velocity gradient du/dy
   for j=2:nj-1
      %dudy(j)= (U(j+1) - U(j-1)) / (dY(j,1) + dY(j,2));
      
      dN = dY(j,1);
      dS = dY(j,2);
      factor = dS^2/(2*dN*dS + dN^2);
      
      b = (U(j) - U(j-1) + ((U(j) - U(j+1))*factor))/(dS - dN*factor);
      a = -(U(j) - U(j+1) + dN*b)*factor/dS^2;

      dudy(j) = 2*a*dS + b;
   end
   
   %Computing cMu
   S = dudy .* sqrt(2);
   W = 2*dudy;
   Re = sqrt(2);%abs(W./S);
   Tt = sqrt(k.^2 ./ eps.^2 + 2*visc./eps);
   psi = Tt .* S .* max(1,Re);
   cNy = 1 ./ (2*(1 + Tt .* S .* sqrt(1 + Re.^2)));
   Pib = cNy.^2 .* psi.^2;
   g = (1 + 2 * cNy .* psi.^2).^(-1);
   alpha1 = g .* (0.25 + (2/3)* Pib.^(1/2));
   alpha2 = 3 * g ./ (8*sqrt(2));
   alpha3 = 3 * alpha2/sqrt(2);
   eta = alpha2 .* Tt .* S;
   xi = alpha3 .* Tt .* W;
   cMu = 3 * (1 + eta.^2).*alpha1 ./ (3 + eta.^2 + 6*eta.^2 .* xi.^2 + 6*xi.^2);
   
   %Computing Rtilde gradient dRtildedy
   Rtilde = k .* Tt;
   for j=2:nj-1      
      dN = dY(j,1);
      dS = dY(j,2);
      factor = dS^2/(2*dN*dS + dN^2);
      
      b = (Rtilde(j) - Rtilde(j-1) + ((Rtilde(j) - Rtilde(j+1))*factor))/(dS - dN*factor);
      a = -(Rtilde(j) - Rtilde(j+1) + dN*b)*factor/dS^2;

      dRtildedy(j) = 2*a*dS + b;
   end
   
   %Computing other coefficients
   C1 = 2 * cMu .* psi .* (1 - cMu .* psi);
   C2 = min(2*cMu,cMuTilde * sqrt(1 + (C1/(6*cMuTilde))));
   sigma = cMu + damping./sqrt(2);
   Slambda = sqrt(eps .* cNy .* psi.^2 .* Rtilde);
   
%
%
%  Often it can be tricky to start the simulations. They often diverge.
%  It can then be useful to compute the turbulent viscosity from the 
%  mixing-length model for the first 2000 iterations.
%  In this way the K and EPS (or OMEGA) eqns. are de-coupled from U 
%  for the first 2000 iterations, i.e. the K and EPS (or OMEGA) do not
%  influence U since the viscosity is taken from the mixing-length model.

   vist_old=vist;

   
   if count < 2000   % use mixing length model for turbulent viscosity if count >2000
      for j=2:nj-1
% compute turbulent viscosity
         yplus=ustar*y_node(j)/visc;
         damp=1-exp(-yplus/26);
         ell=min(damp*kappa*y_node(j),0.09);
         vist(j)=urf*abs(dudy(j))*ell^2+(1-urf)*vist_old(j);
      end
   else
       vist(2:nj-1) =  cMu .* (k(2:nj-1).^2)./eps(2:nj-1);
   end
     
   %Calculating source terms
   %Pk = (vist .* (dudy).^2);

   uSp = zeros(nj,1);
   uSu = ones(nj,1) .* deltaY;
   %uSp(2) = -(cMu)^(1/4)*k(2)^(1/2);%*deltaY(2); % /U(2);
   
%   Perhaps useful reference for b.c. implementation
%    kSp = (-eps./ k) .* deltaY;
%    kSu = Pk .* deltaY;
%    kSp(2) = -cMu^(3/4)*k(2)^(1/2);%*U(2);
%    kSu(2) = U(2);
% 
%    epsSp = ((c1Eps .* Pk - c2Eps .* eps) ./ k) .* deltaY;
%    epsSu = zeros(nj,1);
%    epsSp(2) = -1e10;
%    epsSu(2) = cMu^(3/4)*k(2)^(3/2)*1e10/(kappa*deltaY(2));
%    
   RSp = -C2 .* (dRtildedy).^2 .* deltaY ./ R;
   RSu = C1 .* Slambda .* deltaY;
 
   dampingSu = 1;

   %Calculating coefficients
   UCoeff = CalcCoeffs( 1, dY, deltaY, visc, vist, uSp, nj, BCU);
   RCoeff = CalcRCoeffs(U, sigma, dY, deltaY, visc, vist,RSp,nj, BCR);
   dampingCoeff = CarcDampingCoeff(L_sq, dY, nj, BCDamping);
   % using visc and vist since rho = 1 --> kin_visc = dyn_visc
   
   %Gauss-Seidel iteration
   U_new = GaussSeidel(U,uSu,UCoeff);
   R_new = GaussSeidel(R,RSu,RCoeff);
   damping_new = GaussSeidel(damping, dampingSu, dampingCoeff);
   
   U_new(end) = U_new(end-1);
   R_new(end) = R_new(end-1);
   damping(end) = damping(end-1);
   
   %disp([epsCoeff.point, epsCoeff.south,epsCoeff.north])
   
% after having computed ap and su, use under-relaxation (see lecture notes) 
%  Compute the velocity U
   
%  compute U & R
    U(2:nj-1) = U(2:nj-1) + urC.*(U_new(2:nj-1) - U(2:nj-1));
    R(2:nj-1) = R(2:nj-1) + urC.*(R_new(2:nj-1) - R(2:nj-1));
    damping(2:nj-1) = damping(2:nj-1) + urC.*(damping_new(2:nj-1) - damping(2:nj-1));
   
% Create k and eps again
% ...
    
% Convergence criterian (Check the error)
% compute residuals R
     residual = ComputeResidual(U,UCoeff,uSu,R,RCoeff,RSu,nj);
% compute the flux F
     F = ComputeFlux(U,deltaY,nj);

  error = abs(residual/F);
  if(mod(count,200) == 0)
      UStore = [UStore U];
      kStore = [kStore k];
      epsStore = [epsStore eps];
      
      if(error > old_error)
          urC = max(0.1*urC , 0.1);
      else
          urC = min(1.01*urC , 0.99999);
      end
      fprintf('%e,  %f \n',error, urC);
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
print eps.ps -deps

% plot shear stress
figure(3)
hold on
load uv_dns.dat
plot(y_dns,-uv_dns,'bo')
plot(y_node,dudy,'rx')
xlabel('x')
ylabel('turbulent shear stress -uv')
legend('DNS','Calc.')
print uv.ps -deps

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


%%

figure(1)
contourf(UStore)
%contourf(UStore(1:80,:))
colorbar
figure(2)
contourf(kStore)
%contourf(kStore(1:80,:))
colorbar
figure(3)
contourf(epsStore)
%contourf(epsStore(1:80,:))
colorbar

