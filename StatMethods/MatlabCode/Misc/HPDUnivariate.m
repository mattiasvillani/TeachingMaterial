function HPDRegions=HPDUnivariate(GridPoints,Density,HPDCovProb)

NGridPoints=length(GridPoints);
[Sorted,Index]=sort(Density);
NCutOff=sum(cumsum(Sorted)<(1-HPDCovProb));
HPDGridPoints=sort(GridPoints(Index(NCutOff+1:end)));
Jumps=find(diff(sort(Index(NCutOff+1:end)))>1);
Jumps(end+1)=length(HPDGridPoints);
LowerEndPoints=[];
UpperEndPoints=[];
NJumps=length(Jumps);
Start=1;

for j=1:NJumps
    LowerEndPoints=[LowerEndPoints min(HPDGridPoints(Start:Jumps(j)))];
    UpperEndPoints=[UpperEndPoints max(HPDGridPoints(Start:Jumps(j)))];
    Start=Jumps(j)+1;
end

HPDRegions=[LowerEndPoints' UpperEndPoints']; % NHPDRegions-by-2 matrix, i:th row contains the lower and upper limit of the i:th region