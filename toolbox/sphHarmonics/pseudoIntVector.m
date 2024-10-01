function I = pseudoIntVector(P,k,Omega,R)
%I = pseudoIntVector(P,k,Omega,R) Obtain the pseudo-intensity vector for an
%Eigenmike.
%   Input:
%       - P         : sound pressure in freq domain. [Nf x M]
%       - k         : wavenumber. [Nf x 1]
%       - Omega     : microphone coordinates (omega, phi). [M x 2]
%       - R         : array radius. Scalar
%   Output:
%       - I         : pseudo-intensity vector. [Nf x 3]
%
% Author: Antonio Figueroa-Duran
% Date: June 2024

% Genreate SH
Nmax = 1;
Y = sh(Omega,Nmax);

% Gen rotation coefficients alpha
alpha_x = sh([pi/2 0],1);       alpha_x = alpha_x(2:end);
alpha_y = sh([pi/2 pi/2],1);    alpha_y = alpha_y(2:end);
alpha_z = sh([0 0],1);          alpha_z = alpha_z(2:end);

% SHT in the LS sense
j_1 = radFunction(k,R,1,"open");
b_1 = 4*pi./((k*R).^2.*j_1(:,2).');
p_SH = (4*pi/size(P,2))*pinv(Y)*P.';

% Cartesian proyection
p_a(1,:) = (1./b_1).*(alpha_x*p_SH(2:4,:));
p_a(2,:) = (1./b_1).*(alpha_y*p_SH(2:4,:));
p_a(3,:) = (1./b_1).*(alpha_z*p_SH(2:4,:));

% Pseudo-intensity vector
I = 0.5*(real(conj(p_SH(1,:)).*p_a))';

end