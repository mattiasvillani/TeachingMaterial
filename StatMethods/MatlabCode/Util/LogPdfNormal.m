function  LogDens = LogPdfNormal(x,mu,sigma)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % PURPOSE:     Computing the log height of the N(mu,sigma^2) density
 %
 % INPUT:       x          (scalar)        Point where the density is evaluated  
 %              mu         (scalar)        Mean
 %              sigma      (scalar)        Standard deviation
 %      
 % OUTPUT:      LogDens    (scalar)        Log height of the N(mu,sigma^2) distribution
 %
 % AUTHOR:      Mattias Villani, Research Department, Sveriges Riksbank and
 %              Department of Statistics, Stockholm University. 
 %              E-mail: mattias.villani@riksbank.se
 %
 % REVISED:     2004-03-04
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 LogDens = -log(sigma) -0.5*log(2*pi) -0.5*((x-mu)./sigma).^2;

