function [res] = dVd020(x,y,z)
 res=(1/4).*pi.^(-1).*z.*((-2).*((-1)+x).*((-1)+y).*(1+(-2).*y+y.^2+z.^2).^( ...
  -1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+(-4).*((-1)+x).*((-1)+y) ...
  .*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -1/2)+2.*(1+x).*((-1)+y).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2) ...
  .*y+y.^2+z.^2).^(-3/2)+4.*(1+x).*((-1)+y).*(1+(-2).*y+y.^2+z.^2).^(-2).* ...
  (2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+2.*((-1)+x).*(1+y).*(1+2.*y+ ...
  y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+4.*((-1)+x).*( ...
  1+y).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+ ...
  (-2).*(1+x).*(1+y).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-3/2)+(-4).*(1+x).*(1+y).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2));
end