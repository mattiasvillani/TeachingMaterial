function [XRepStatSample,XObsStat,MuRes,SigmaRes]=PostPredAnalysisNormal(Data,PriorMean,PriorStd,NIter,NRepIter)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % PURPOSE:     Evaluating fit of the usual iid normal model. Normal prior.
 %
 % INPUT:       Data        (NObs-by-1)     Observed data
 %              PriorMean   (scalar)        Prior mean
 %              PriorStd    (scalar)        Prior standard deviation
 %              NIter  (scalar)             Number of Gibbs sampling draws from posterior
 %              NRepIter    (scalar)        Number of data replications for each posterior draw.
 %
 % AUTHOR:      Mattias Villani, Research Department, Sveriges Riksbank and
 %              Department of Statistics, Stockholm University. 
 %              E-mail: mattias.villani@riksbank.se
 %
 % REVISED:     2004-12-15
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 XBar=mean(Data);
 SigmaData=std(Data);
 n=length(Data);
 Sqrtn=sqrt(n);
 Sigma=SigmaData; % Initial value
 Sigma2Data=SigmaData^2;
 XObsStat=max(abs(Data));
 MuRes=zeros(NIter,1);
 SigmaRes=zeros(NIter,1);
 XRepStatSample=zeros(NIter*NRepIter,1);
 for i=1:NIter
     Mu=XBar+randn*(Sigma/Sqrtn);
     Sigma=1/sqrt(gamrnd((n+3)/2,2/((n-1)*Sigma2Data)));
     MuRes(i)=Mu;
     SigmaRes(i)=Sigma;
     for j=1:NRepIter
         RepSample=Mu*ones(n,1)+Sigma*randn(n,1);
         XRepStatSample(j+(i-1)*NRepIter)=max(abs(RepSample));
     end
 end
