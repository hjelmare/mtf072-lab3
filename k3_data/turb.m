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
format long

%Defining simulation constants
cMu = 0.09;
sigmaK = 1.00;
sigmaEps = 1.30;
c1Eps = 1.44;
c2Eps = 1.92;
viscos=1/395;

% wall friction velocity
ustar=1;
% read u from DNS data base
load u_dns.dat

% read y from DNS data base
load y_dns.dat
yc=y_dns;          %yc is a vector contains the faces coordinates (grid)
nj=length(yc)+1; % nj: number of node points  
u_dns(nj)=u_dns(nj-1);
y_node(1)=yc(1);
for j=2:nj-1
    y_node(j)= (yc(j)+yc(j-1))/2;
end
y_node(nj)=yc(nj-1);

%Calculating dY and deltaY
dY(:,1) = [0 diff(y_node(2:end)) 0]';
dY(:,2) = [0 diff(y_node(1:end-1)) 0]';


% init velocity, vist, k & eps
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
error=1;
count=0;
max_error=0.001;
urf=0.8;
while error > max_error 
    
    count = count+1;
    
    % Compute the velocity gradient du/dy
   for j=2:nj-1
      dudy(j)= (U(j+1) - U(j-1)) / (dY(j,1) + dY(j,2));
   end
    
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
         
         yplus=ustar*y_node(j)/viscos;
         damp=1-exp(-yplus/26);
         ell=min(damp*kappa*y_node(j),0.09);
         vist(j)=urf*abs(dudy(j))*ell^2+(1-urf)*vist_old;
      end
   else
       for j = 2:nj-1
         vist(j) =  cMu * (k(j)^2)/epsi(j);
       end
   end
   
   %Calculating source terms
   Pk = (vist .* (dudy).^2);

   kSp = (-eps./ k) .* deltaY;
   kSu = Pk .* deltaY;

   uSp = zeros(length,1);
   uSu = ones(length,1) .* deltaY;

   epsSp = (-c2Eps .* eps) ./ k;
   epsSu = (eps ./ k) * c1Eps .* Pk;
%
%
% ....
% ....
% ....
% ....
%  your finite volume code
% ....
% ....
% ....
% ....
% ....
%
% after having computed ap and su, use under-relaxation (see lecture notes) 
%  Compute the velocity U
   for j=2:nj-1
%  compute U
%    U(j)=...
   end

%  Compute the turbulent kinetic energy k
   for j=2:nj-1
%  compute k  
%     k(j)=....
   end

%  Compute the turbulent dissipation epsilon
   for j=2:nj-1
%  compute epsi
%     epsi(j)=...
   end

  
% Convergence criterian (Check the error)
   for j=2:nj-1
% compute residuals R
%     R(j)=...
% compute the flux F
%     F(j)=....
   end
  error = sum(R)/sum(F);
  

end  %while
%
% plot
% compare with DNS
%
% plot k
figure(1)

hold on
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat

k_dns=0.5*(u2_dns+v2_dns+w2_dns);
plot(y_dns,k_dns,'bo')
hold
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
eps_dns=dns_data(:,2)*ustar^4/viscos;
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
