% cartersian grid
[X,Y,Z] = meshgrid(-10:.5:10,-10:.5:10,-10:.5:10);
% cylindrical coordinates
[phi,s,z] = cart2pol(X,Y,Z);
% physical parameters/constants
a = 5;                % radius of solenoid
mu_0 = 1.25663706e-6; % permeability of free space
L = 0.5;              % length of the solenoid
I = 1;                % current 1 A
gamma1 = z - L/2;
gamma2 = z + L/2;
h = sqrt(4.*a.*s./(s + a).^2);
k1 = sqrt(4.*a.*s./((s + a).^2 + gamma1.^2));
k2 = sqrt(4.*a.*s./((s + a).^2 + gamma2.^2));
% elliptic integrals
[K1,E1] = ellipke(k1.^2); % complete elliptic integral of 1st and 2nd kind of k1.^2
[K2,E2] = ellipke(k2.^2); % complete elliptic integral of 1st and 2nd kind of k2.^2
T1 = ellipticPi(h.^2,k1.^2); % complete elliptic integral of the third 
T2 = ellipticPi(h.^2,k2.^2); % complete elliptic integral of the third 

% magnitude of the radial component of the mag. field
Bs_mid = (mu_0./4.*pi).*I.*(1./L).*sqrt(a./s).*((((k1.^2 -2)./k1).*K1 + (2./k1).*E1) ... 
    - (((k2.^2 - 2)./k2).*K2 + (2./k2).*E2));
% magnitude of the z-component of the mag. field
Bz_mid = -1.*(mu_0./4.*pi).*((2.*L.*sqrt(a.*s)).^(-1)).*((gamma1.*k1.*(K1 + ((a-s)./(a+s)).*T1)) ...
    - (gamma2.*k2.*(K2 + ((a-s)./(a+s)).*T2)));
% magnitude of the phi-component of the mag. field (cylindrical symmetry)
Bphi_mid = zeros(size(Bz_mid));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BOTTOM SOLENOID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cartersian grid
[X1,Y1,Z1] = meshgrid(-10:.5:10,-10:.5:10,-10:.5:10 + 5*sqrt(3/7));
% cylindrical coordinates
[phi,s,z] = cart2pol(X1,Y1,Z1);
% physical parameters/constants
a = a*sqrt(4/7);      % radius of solenoid
I = I*(49/64);        % current 1 A
gamma1 = z - L/2;
gamma2 = z + L/2;
h = sqrt(4.*a.*s./(s + a).^2);
k1 = sqrt(4.*a.*s./((s + a).^2 + gamma1.^2));
k2 = sqrt(4.*a.*s./((s + a).^2 + gamma2.^2));
% elliptic integrals
[K1,E1] = ellipke(k1.^2); % complete elliptic integral of 1st and 2nd kind of k1.^2
[K2,E2] = ellipke(k2.^2); % complete elliptic integral of 1st and 2nd kind of k2.^2
T1 = ellipticPi(h.^2,k1.^2); % complete elliptic integral of the third 
T2 = ellipticPi(h.^2,k2.^2); % complete elliptic integral of the third 

% magnitude of the radial component of the magnetic field
Bs_bot = (mu_0./4.*pi).*I.*(1./L).*sqrt(a./s).*((((k1.^2 -2)./k1).*K1 + (2./k1).*E1) ... 
    - (((k2.^2 - 2)./k2).*K2 + (2./k2).*E2));
% magnitude of the radial component of the magnetic field
Bphi_bot = zeros(size(Bs_bot));
% magnitude of the z-component of the magnetic field
Bz_bot = -1.*(mu_0./4.*pi).*((2.*L.*sqrt(a.*s)).^(-1)).*((gamma1.*k1.*(K1 + ((a-s)./(a+s)).*T1)) ...
    - (gamma2.*k2.*(K2 + ((a-s)./(a+s)).*T2)));
% resize
Bs_bot(:,:,1:6) = [];
Bphi_bot(:,:,1:6) = [];
Bz_bot(:,:,1:6) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP SOLENOID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cartersian grid
[X2,Y2,Z2] = meshgrid(-10:.5:10,-10:.5:10,-10 - 5*sqrt(3/7):.5:10);
% cylindrical coordinates
[phi,s,z] = cart2pol(X2,Y2,Z2);
% physical parameters/constants
gamma1 = z - L/2;
gamma2 = z + L/2;
h = sqrt(4.*a.*s./(s + a).^2);
k1 = sqrt(4.*a.*s./((s + a).^2 + gamma1.^2));
k2 = sqrt(4.*a.*s./((s + a).^2 + gamma2.^2));
% elliptic integrals
[K1,E1] = ellipke(k1.^2); % complete elliptic integral of 1st and 2nd kind of k1.^2
[K2,E2] = ellipke(k2.^2); % complete elliptic integral of 1st and 2nd kind of k2.^2
T1 = ellipticPi(h.^2,k1.^2); % complete elliptic integral of the third 
T2 = ellipticPi(h.^2,k2.^2); % complete elliptic integral of the third 

% magnitude of the radial component of the magnetic field
Bs_top = (mu_0./4.*pi).*I.*(1./L).*sqrt(a./s).*((((k1.^2 -2)./k1).*K1 + (2./k1).*E1) ... 
    - (((k2.^2 - 2)./k2).*K2 + (2./k2).*E2));
% magnitude of the radial component of the magnetic field
Bphi_top = zeros(size(Bs_top));
% magnitude of the z-component of the magnetic field
Bz_top = -1.*(mu_0./4.*pi).*((2.*L.*sqrt(a.*s)).^(-1)).*((gamma1.*k1.*(K1 + ((a-s)./(a+s)).*T1)) ...
    - (gamma2.*k2.*(K2 + ((a-s)./(a+s)).*T2)));
% resize
Bs_top(:,:,42:47) = [];
Bphi_top(:,:,42:47) = [];
Bz_top(:,:,42:47) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUPERPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bs = Bs_mid + Bs_top + Bs_bot;
Bphi = Bphi_mid + Bphi_top + Bphi_bot;
Bz = Bz_mid + Bz_top + Bz_bot;
figure()
quiver3(X,Y,Z,Bs,Bphi,Bz)
field_strength = sqrt(Bs.^2 + Bphi.^2 + Bz.^2);
% figure()
% coneplot(X,Y,Z,Bs,Bphi,Bz,X,Y,Z,field_strength,1)
% view(2)
% shading interp
% colorbar
% 
% hsurfaces = slice(X,Y,Z,field_strength,[xmin,xmax],ymax,zmin);
% set(hsurfaces,'FaceColor','interp','EdgeColor','none')
% hold off