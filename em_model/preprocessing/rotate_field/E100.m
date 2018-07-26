function [res] = E100(E,phi,theta,psi)
 res=(-1).*sin(theta).*E(1,1,2)+cos(theta).*(sin(psi).*E(1,2,1)+cos(psi).*E(2,1,1));
end