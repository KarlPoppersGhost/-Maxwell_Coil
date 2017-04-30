% cartersian grid
[X,Y,Z] = meshgrid(-10:.5:10,-10:.5:10,-10:.5:10);
% cylindrical coordinates
[phi,s,z] = cart2pol(X,Y,Z);
% physical parameters/constants
a = 5;                % radius of solenoid
mu_0 = 1.25663706e-6; % permeability of free space
L = 1;                % length of the solenoid
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

% phi-component of the magnetic vector potential
A_phi = (mu_0./4.*pi).*I.*(1./L).*sqrt(a./s).*((gamma1.*k1.*((((k1.^2 + h.^2 - (h.^2.*k1.^2))./(h.^2.*k1.^2)).*K1)...
    - ((k1.^-2).*E1) + ((h.^2 - 1)./h.^2).*T1)) - (gamma2.*k2.*((((k2.^2 + h.^2 - (h.^2.*k2.^2))./(h.^2.*k2.^2)).*K2)...
    - ((k2.^-2).*E2) + ((h.^2 - 1)./h.^2).*T2)));
% magnitude of the radial component of the mag. field
Bs_mid = (mu_0./4.*pi).*I.*(1./L).*sqrt(a./s).*((((k1.^2 -2)./k1).*K1 + (2./k1).*E1) ... 
    - (((k2.^2 - 2)./k2).*K2 + (2./k2).*E2));
% magnitude of the z-component of the mag. field
Bz_mid = -1.*(mu_0./4.*pi).*((2.*L.*sqrt(a.*s)).^(-1)).*((gamma1.*k1.*(K1 + ((a-s)./(a+s)).*T1)) ...
    - (gamma2.*k2.*(K2 + ((a-s)./(a+s)).*T2)));
% magnitude of the phi-component of the mag. field (cylindrical symmetry)
Bphi_mid = zeros(size(B_z));