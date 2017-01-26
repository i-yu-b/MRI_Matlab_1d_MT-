function result = superlorentzian(delta,T2)
% superlorentzian lineshape function
% inputs: T2 - relaxation time of restiricted pool
%         delta - frequency offset, linear frequency=nu

integrand = @(u,delta,T2) sqrt(2/pi)*T2./...
    abs(3*u.^2-1).*exp(-2*(2*pi*delta*T2./(3*u.^2-1)).^2);
result = integral(@(u)integrand(u,delta,T2),0,1); 

% % gaussian
% result = T2/sqrt(2*pi)*exp(-(2*pi*delta*T2)^2/2);