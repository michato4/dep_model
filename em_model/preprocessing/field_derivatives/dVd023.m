function [res] = dVd023(x,y,z)
 res=(1/2).*pi.^(-1).*z.*(180.*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-7/2)+216.*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1) ...
  .*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -5/2)+96.*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2) ...
  .*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2)+(-90).* ...
  ((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2)+288.*((-1)+x).*( ...
  (-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3) ...
  .*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+128.*((-1)+x).*((-1)+y) ...
  .^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+( ...
  -2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+(-120).*((-1)+x).*((-1)+y).*z.*( ...
  1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+32.*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).* ...
  x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).* ...
  y+y.^2+z.^2).^(-3/2)+(-72).*((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2) ...
  .^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2) ...
  .^(-3/2)+576.*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+ ...
  (-2).*y+y.^2+z.^2).^(-4).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+ ...
  256.*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+ ...
  y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+(-240).*(( ...
  -1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2) ...
  .^(-3).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+64.*((-1)+x).*((-1)+ ...
  y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+( ...
  -2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+(-144).*((-1)+x).*((-1)+y).*z.*( ...
  1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+(-48).*((-1)+x).*((-1)+y).*z.*(1+(-2).* ...
  x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).* ...
  y+y.^2+z.^2).^(-1/2)+(-180).*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2) ...
  .^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -7/2)+(-216).*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).* ...
  y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2)+(-96).*(1+x) ...
  .*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1) ...
  .*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2)+90.*(1+x).*((-1)+y).*z.*(1+2.* ...
  x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-5/2)+(-288).*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^( ...
  -1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -3/2)+(-128).*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).* ...
  y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+120.*(1+x).* ...
  ((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+ ...
  2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+(-32).*(1+x).*((-1)+y).^3.*z.*(1+ ...
  2.*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).* ...
  y+y.^2+z.^2).^(-3/2)+72.*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-2).*( ...
  1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+( ...
  -576).*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-4).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+(-256).*(1+x).*(( ...
  -1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+ ...
  2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+240.*(1+x).*((-1)+y).*z.*(1+2.*x+ ...
  x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-1/2)+(-64).*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^( ...
  -3).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -1/2)+144.*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+ ...
  y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+48.*(1+x).*(( ...
  -1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.* ...
  x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+(-180).*((-1)+x).*(1+y).^3.*z.*(1+(-2) ...
  .*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-7/2)+(-216).*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2) ...
  .^(-1).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -5/2)+(-96).*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.* ...
  y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2)+90.*((-1)+x) ...
  .*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+( ...
  -2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2)+(-288).*((-1)+x).*(1+y).^3.*z.*(1+( ...
  -2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-3/2)+(-128).*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2) ...
  .^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -3/2)+120.*((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+ ...
  y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+(-32).*((-1)+ ...
  x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-1).* ...
  (2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+72.*((-1)+x).*(1+y).*z.*(1+(-2) ...
  .*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-3/2)+(-576).*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2) ...
  .^(-1).*(1+2.*y+y.^2+z.^2).^(-4).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -1/2)+(-256).*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.* ...
  y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+240.*((-1)+ ...
  x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+ ...
  (-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+(-64).*((-1)+x).*(1+y).^3.*z.*(1+( ...
  -2).*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-1/2)+144.*((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-2) ...
  .*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+48.* ...
  ((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^( ...
  -1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+180.*(1+x).*(1+y).^3.*z.*( ...
  1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-7/2)+216.*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*( ...
  1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-5/2)+96.*(1+x) ...
  .*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-5/2)+(-90).*(1+x).*(1+y).*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -5/2)+288.*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+ ...
  z.^2).^(-3).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+128.*(1+x).*(1+y).^3.* ...
  z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.* ...
  y+y.^2+z.^2).^(-3/2)+(-120).*(1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).* ...
  (1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+32.*(1+x) ...
  .*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-1).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+(-72).*(1+x).*(1+y).*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -3/2)+576.*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+ ...
  z.^2).^(-4).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+256.*(1+x).*(1+y).^3.* ...
  z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+2.* ...
  y+y.^2+z.^2).^(-1/2)+(-240).*(1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).* ...
  (1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+64.*(1+x) ...
  .*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-2).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+(-144).*(1+x).*(1+y).*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -1/2)+(-48).*(1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+ ...
  z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+(-105).*((-1)+x).*(( ...
  -1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).* ...
  (2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-9/2).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+2.*z.^2)+(-120).*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^( ...
  -1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -7/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-60).*((-1)+x).*((-1)+y) ...
  .^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+( ...
  -2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-7/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.* ...
  z.^2)+45.*((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).* ...
  y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-7/2).*(2+(-2) ...
  .*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-144).*((-1)+x).*((-1)+y).^3.*z.*(1+(-2) ...
  .*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2) ...
  .*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-72).*(( ...
  -1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+2.*z.^2)+54.*((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-5/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-24).*((-1)+x).*(( ...
  -1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-1).* ...
  (2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+2.*z.^2)+36.*((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*( ...
  1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).* ...
  (2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-192).*((-1)+x).*((-1)+y).^3.*z.* ...
  (1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-4).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+( ...
  -96).*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+ ...
  y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+(-2).* ...
  x+x.^2+(-2).*y+y.^2+2.*z.^2)+72.*((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-3/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-32).*((-1)+x).*(( ...
  -1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-2).* ...
  (2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+2.*z.^2)+48.*((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*( ...
  1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).* ...
  (2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+24.*((-1)+x).*((-1)+y).*z.*(1+(-2) ...
  .*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2) ...
  .*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-384).*(( ...
  -1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-5).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+2.*z.^2)+(-192).*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+ ...
  x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-4).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+144.*((-1)+x) ...
  .*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-4) ...
  .*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+2.*z.^2)+(-64).*((-1)+x).*((-1)+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^( ...
  -3).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -1/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+96.*((-1)+x).*((-1)+y).*z.* ...
  (1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+ ...
  48.*((-1)+x).*((-1)+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+2.*z.^2)+105.*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-9/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+120.*(1+x).*((-1)+y) ...
  .^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-7/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+60.* ...
  (1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2) ...
  .^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-7/2).*(2+2.*x+x.^2+(-2).*y+ ...
  y.^2+2.*z.^2)+(-45).*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+( ...
  -2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-7/2).*(2+2.* ...
  x+x.^2+(-2).*y+y.^2+2.*z.^2)+144.*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-5/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+72.*(1+x).*((-1)+y) ...
  .^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+( ...
  -54).*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+(-2) ...
  .*y+y.^2+2.*z.^2)+24.*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-3).*( ...
  1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+ ...
  2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-36).*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-5/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+192.*(1+x).*((-1)+y) ...
  .^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-4).*(2+2.*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+96.* ...
  (1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2) ...
  .^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+(-2).*y+ ...
  y.^2+2.*z.^2)+(-72).*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+( ...
  -2).*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+2.* ...
  x+x.^2+(-2).*y+y.^2+2.*z.^2)+32.*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-3/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-48).*(1+x).*((-1)+y) ...
  .*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+( ...
  -24).*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+(-2) ...
  .*y+y.^2+2.*z.^2)+384.*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).* ...
  (1+(-2).*y+y.^2+z.^2).^(-5).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*( ...
  2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+192.*(1+x).*((-1)+y).^3.*z.*(1+2.*x+ ...
  x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-4).*(2+2.*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-144).*(1+x).*(( ...
  -1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-4).*(2+2.* ...
  x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+ ...
  64.*(1+x).*((-1)+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+(-2) ...
  .*y+y.^2+2.*z.^2)+(-96).*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+z.^2).^(-2).*( ...
  1+(-2).*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+ ...
  2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-48).*(1+x).*((-1)+y).*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-3).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-1/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+105.*((-1)+x).*(1+y) ...
  .^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).* ...
  x+x.^2+2.*y+y.^2+z.^2).^(-9/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+120.* ...
  ((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2) ...
  .^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-7/2).*(2+(-2).*x+x.^2+2.*y+ ...
  y.^2+2.*z.^2)+60.*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*( ...
  1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-7/2).*(2+( ...
  -2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-45).*((-1)+x).*(1+y).*z.*(1+(-2).*x+ ...
  x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-7/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+144.*((-1)+x).*(1+y) ...
  .^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+(-2).* ...
  x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+72.*( ...
  (-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^( ...
  -2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+2.*y+y.^2+ ...
  2.*z.^2)+(-54).*((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.* ...
  y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+ ...
  x.^2+2.*y+y.^2+2.*z.^2)+24.*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2) ...
  .^(-3).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -5/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-36).*((-1)+x).*(1+y).*z.*(1+ ...
  (-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.* ...
  y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+192.*((-1)+x).* ...
  (1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-4).*(2+( ...
  -2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+ ...
  96.*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+ ...
  z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+ ...
  2.*y+y.^2+2.*z.^2)+(-72).*((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^( ...
  -1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2).* ...
  (2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+32.*((-1)+x).*(1+y).^3.*z.*(1+(-2).* ...
  x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-48).*((-1)+x).*(1+y) ...
  .*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-24).* ...
  ((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^( ...
  -1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+ ...
  2.*z.^2)+384.*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.* ...
  y+y.^2+z.^2).^(-5).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+ ...
  x.^2+2.*y+y.^2+2.*z.^2)+192.*((-1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-4).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2) ...
  .^(-1/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-144).*((-1)+x).*(1+y).* ...
  z.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-4).*(2+(-2).*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+64.*(( ...
  -1)+x).*(1+y).^3.*z.*(1+(-2).*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^( ...
  -3).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+2.*y+y.^2+ ...
  2.*z.^2)+(-96).*((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.* ...
  y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+ ...
  x.^2+2.*y+y.^2+2.*z.^2)+(-48).*((-1)+x).*(1+y).*z.*(1+(-2).*x+x.^2+z.^2) ...
  .^(-3).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -1/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-105).*(1+x).*(1+y).^3.*z.*( ...
  1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-9/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-120).*(1+x).*(1+y) ...
  .^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-7/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-60).*(1+ ...
  x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-7/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+45.*( ...
  1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-7/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+( ...
  -144).*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2) ...
  .^(-3).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+2.*y+y.^2+2.* ...
  z.^2)+(-72).*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+ ...
  z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+2.*y+ ...
  y.^2+2.*z.^2)+54.*(1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+ ...
  y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+2.* ...
  y+y.^2+2.*z.^2)+(-24).*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+ ...
  2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+2.*z.^2)+36.*(1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-2).*( ...
  1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+2.*z.^2)+(-192).*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2) ...
  .^(-1).*(1+2.*y+y.^2+z.^2).^(-4).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-3/2).* ...
  (2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-96).*(1+x).*(1+y).^3.*z.*(1+2.*x+x.^2+ ...
  z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+72.*(1+x).*(1+y).*z.*(1+2.*x+ ...
  x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+2.*y+y.^2+z.^2) ...
  .^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-32).*(1+x).*(1+y).^3.*z.*(1+ ...
  2.*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+48.*(1+x).*(1+y).*z.*(1+ ...
  2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+24.*(1+x).*(1+y).*z.*(1+ ...
  2.*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-384).*(1+x).*(1+y).^3.* ...
  z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-5).*(2+2.*x+x.^2+2.* ...
  y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-192).*(1+x).*(1+ ...
  y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-4).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+144.*(1+x) ...
  .*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-4).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-64).*(1+ ...
  x).*(1+y).^3.*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-3).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+96.*( ...
  1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-3).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+48.*( ...
  1+x).*(1+y).*z.*(1+2.*x+x.^2+z.^2).^(-3).*(1+2.*y+y.^2+z.^2).^(-2).*(2+ ...
  2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2))+( ...
  1/2).*pi.^(-1).*((-24).*((-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+z.^2).^( ...
  -1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -5/2)+(-32).*((-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2) ...
  .*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+(-8).*( ...
  (-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2) ...
  .^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+18.*((-1)+x).*((-1)+ ...
  y).*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).* ...
  x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+(-64).*((-1)+x).*((-1)+y).^3.*(1+(-2) ...
  .*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2) ...
  .*y+y.^2+z.^2).^(-1/2)+(-16).*((-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-1/2)+36.*((-1)+x).*((-1)+y).*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+( ...
  -2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+12.* ...
  ((-1)+x).*((-1)+y).*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2) ...
  .^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+24.*(1+x).*((-1)+y) ...
  .^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-5/2)+32.*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-3/2)+8.*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).* ...
  y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+(-18).*(1+x) ...
  .*((-1)+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+ ...
  2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2)+64.*(1+x).*((-1)+y).^3.*(1+2.*x+ ...
  x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-1/2)+16.*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+ ...
  (-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+(-36) ...
  .*(1+x).*((-1)+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2) ...
  .*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2)+(-12).*(1+x).*((-1)+y).*(1+2.* ...
  x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-1/2)+24.*((-1)+x).*(1+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).* ...
  (1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2)+32.*(( ...
  -1)+x).*(1+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-2) ...
  .*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+8.*((-1)+x).*(1+y).^3.*(1+(-2) ...
  .*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-3/2)+(-18).*((-1)+x).*(1+y).*(1+(-2).*x+x.^2+z.^2).^(-1).* ...
  (1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+64.*(( ...
  -1)+x).*(1+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3) ...
  .*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+16.*((-1)+x).*(1+y).^3.*(1+( ...
  -2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+ ...
  y.^2+z.^2).^(-1/2)+(-36).*((-1)+x).*(1+y).*(1+(-2).*x+x.^2+z.^2).^(-1).* ...
  (1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+(-12) ...
  .*((-1)+x).*(1+y).*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1) ...
  .*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+(-24).*(1+x).*(1+y).^3.*(1+2.* ...
  x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-5/2)+(-32).*(1+x).*(1+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+ ...
  y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+(-8).*(1+x).*(1+ ...
  y).^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+ ...
  2.*y+y.^2+z.^2).^(-3/2)+18.*(1+x).*(1+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+ ...
  2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-3/2)+(-64).*(1+x) ...
  .*(1+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2)+(-16).*(1+x).*(1+y).^3.*(1+2.*x+x.^2+z.^2) ...
  .^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+ ...
  36.*(1+x).*(1+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-2).*( ...
  2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+12.*(1+x).*(1+y).*(1+2.*x+x.^2+z.^2) ...
  .^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+2.*y+y.^2+z.^2).^(-1/2)+ ...
  15.*((-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+ ...
  z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-7/2).*(2+(-2).*x+ ...
  x.^2+(-2).*y+y.^2+2.*z.^2)+18.*((-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-5/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+6.*((-1)+x).*((-1)+ ...
  y).^3.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2) ...
  .*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.* ...
  z.^2)+(-9).*((-1)+x).*((-1)+y).*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+ ...
  y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+(-2).* ...
  x+x.^2+(-2).*y+y.^2+2.*z.^2)+24.*((-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-3/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+8.*((-1)+x).*((-1)+ ...
  y).^3.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+(-2) ...
  .*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.* ...
  z.^2)+(-12).*((-1)+x).*((-1)+y).*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).* ...
  y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+(-2) ...
  .*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-6).*((-1)+x).*((-1)+y).*(1+(-2).*x+ ...
  x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+48.*((-1)+x).* ...
  ((-1)+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-4).*( ...
  2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+ ...
  2.*z.^2)+16.*((-1)+x).*((-1)+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2) ...
  .*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+( ...
  -2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-24).*((-1)+x).*((-1)+y).*(1+(-2).*x+ ...
  x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-12).*((-1)+ ...
  x).*((-1)+y).*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).* ...
  (2+(-2).*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+(-2).*y+ ...
  y.^2+2.*z.^2)+(-15).*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+( ...
  -2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-7/2).*(2+2.* ...
  x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-18).*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-5/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-6).*(1+x).*((-1)+y) ...
  .^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+ ...
  x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+9.*( ...
  1+x).*((-1)+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-1).*( ...
  2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.* ...
  z.^2)+(-24).*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+ ...
  y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+ ...
  (-2).*y+y.^2+2.*z.^2)+(-8).*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+z.^2).^(-2) ...
  .*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-3/2).* ...
  (2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+12.*(1+x).*((-1)+y).*(1+2.*x+x.^2+ ...
  z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+ ...
  z.^2).^(-3/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+6.*(1+x).*((-1)+y).*( ...
  1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+y.^2+z.^2).^(-1).*(2+2.*x+x.^2+(-2) ...
  .*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-48).*(1+x) ...
  .*((-1)+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+(-2).*y+y.^2+z.^2).^(-4).*( ...
  2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.* ...
  z.^2)+(-16).*(1+x).*((-1)+y).^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+(-2).*y+ ...
  y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+ ...
  (-2).*y+y.^2+2.*z.^2)+24.*(1+x).*((-1)+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+ ...
  (-2).*y+y.^2+z.^2).^(-3).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^(-1/2).*(2+ ...
  2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+12.*(1+x).*((-1)+y).*(1+2.*x+x.^2+z.^2) ...
  .^(-2).*(1+(-2).*y+y.^2+z.^2).^(-2).*(2+2.*x+x.^2+(-2).*y+y.^2+z.^2).^( ...
  -1/2).*(2+2.*x+x.^2+(-2).*y+y.^2+2.*z.^2)+(-15).*((-1)+x).*(1+y).^3.*(1+ ...
  (-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.* ...
  y+y.^2+z.^2).^(-7/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-18).*((-1)+x) ...
  .*(1+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-2).*(2+( ...
  -2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+ ...
  (-6).*((-1)+x).*(1+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+ ...
  z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+(-2).*x+x.^2+ ...
  2.*y+y.^2+2.*z.^2)+9.*((-1)+x).*(1+y).*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+ ...
  2.*y+y.^2+z.^2).^(-1).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+(-2) ...
  .*x+x.^2+2.*y+y.^2+2.*z.^2)+(-24).*((-1)+x).*(1+y).^3.*(1+(-2).*x+x.^2+ ...
  z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2) ...
  .^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+(-8).*((-1)+x).*(1+y).^3.*( ...
  1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+ ...
  2.*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+12.*((-1)+x) ...
  .*(1+y).*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2) ...
  .*x+x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+6.* ...
  ((-1)+x).*(1+y).*(1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).* ...
  (2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.* ...
  z.^2)+(-48).*((-1)+x).*(1+y).^3.*(1+(-2).*x+x.^2+z.^2).^(-1).*(1+2.*y+ ...
  y.^2+z.^2).^(-4).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+ ...
  x.^2+2.*y+y.^2+2.*z.^2)+(-16).*((-1)+x).*(1+y).^3.*(1+(-2).*x+x.^2+z.^2) ...
  .^(-2).*(1+2.*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+y.^2+z.^2).^( ...
  -1/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+24.*((-1)+x).*(1+y).*(1+(-2).* ...
  x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+(-2).*x+x.^2+2.*y+y.^2+ ...
  z.^2).^(-1/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+12.*((-1)+x).*(1+y).*( ...
  1+(-2).*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+(-2).*x+x.^2+ ...
  2.*y+y.^2+z.^2).^(-1/2).*(2+(-2).*x+x.^2+2.*y+y.^2+2.*z.^2)+15.*(1+x).*( ...
  1+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-7/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+18.*(1+x) ...
  .*(1+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+6.*(1+x).* ...
  (1+y).^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-9).*(1+ ...
  x).*(1+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-5/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+24.*(1+x) ...
  .*(1+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+8.*(1+x).* ...
  (1+y).^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-12).*(1+ ...
  x).*(1+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-6).*(1+ ...
  x).*(1+y).*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-1).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-3/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+48.*(1+x) ...
  .*(1+y).^3.*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-4).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+16.*(1+x) ...
  .*(1+y).^3.*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-24).*(1+ ...
  x).*(1+y).*(1+2.*x+x.^2+z.^2).^(-1).*(1+2.*y+y.^2+z.^2).^(-3).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2)+(-12).*(1+ ...
  x).*(1+y).*(1+2.*x+x.^2+z.^2).^(-2).*(1+2.*y+y.^2+z.^2).^(-2).*(2+2.*x+ ...
  x.^2+2.*y+y.^2+z.^2).^(-1/2).*(2+2.*x+x.^2+2.*y+y.^2+2.*z.^2));
end