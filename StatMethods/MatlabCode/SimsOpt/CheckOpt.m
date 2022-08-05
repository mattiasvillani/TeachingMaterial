function []=CheckOpt(ObjFuncStr,CenterPoint,MinMaxGrid,InvHessian,NumberOfEval,SameVerticalScale,SubplotStructure,EstParTeXNames,varargin)

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % PURPOSE:     Computing and plotting perturbations of the log posterior around a given point (CenterPoint). The pertubations
 %              are one parameter at the time keeping the other parameters fixed at CenterPoint. 
 %              
 %              NaN is allowed in both CenterPoint and MinMaxGrid, but not at the same time for a given row (parameter). 
 %              If CenterPoint(i)=nan, then CenterPoint(i) becomes the midpoint of the interval specified by MinMaxGrid(i,:)
 %              If any of the bounds in MinMaxGrid(i,:) is nan, then the bounds are +- 4 posterior standard deviations, using the InvHessian 
 %              to compute approximate posterior standard deviations.
 %
 % INPUT:       ObjFuncStr      (string)            String with name of .m file containing the objective function
 %              CenterPoint     (1-by-Npara)        Point where the normal approximation of the log posterior is centered.
 %              MinMaxGrid      (Npara-by-2)        MinMaxGrid(i,:) contains the bounds for the i:th parameter. May contain nan's or MinMaxGrid=[]. See PURPOSE. 
 %              InvHessian      (Npara-by-Npara)    InvHessian of the log posterior.
 %              NumberOfEval    (scalar)            Number of point in the grid.
 %              SameVerticalScale (scalar)          If SameVerticalScale==1, then the same scale is used on the vertical axis for each parameter.
 %              EstParTeXNames  (Npara-by-Something)EstParTeXNames(i,:) contains the i:ths parameter's name in TeX code.
 %              SubplotStructure (1-by-2)           Determines the structure of the Subplots (the number of rows and cols)
 %              varargin                            List of additional arguments needed to compute the objective function in ObjFuncStr
 %
 % AUTHOR:      Mattias Villani, Research Department, Sveriges Riksbank and
 %              Department of Statistics, Stockholm University. 
 %              E-mail: mattias.villani@riksbank.se
 %
 % OLDER:       2005-11-03
 % REVISED:     2008-01-22
 %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CenterPoint=CenterPoint(:); % Make sure it is a col vector.
NumberOfParameters=length(CenterPoint);
if isempty(SubplotStructure)
    SubplotStructure=[ceil(sqrt(NumberOfParameters)) ceil(sqrt(NumberOfParameters))];
end
NplotsPerFigure=prod(SubplotStructure);
NFigures=ceil(NumberOfParameters/NplotsPerFigure);
PriorFigHandles=cell(NFigures,1);

if NFigures==1
    [Nrows,Ncols]=ConstructOptimalSubplot(NumberOfParameters);
else
    Nrows=SubplotStructure(1);
    Ncols=SubplotStructure(2);
end

% Computing approximate conditional standard deviations
ApproxStd=zeros(NumberOfParameters,1);
S=InvHessian;
FullSet=1:NumberOfParameters;
for i=1:NumberOfParameters
    iExcluded=setdiff(FullSet,i);
    ApproxStd(i)=sqrt(S(i,i)-S(i,iExcluded)*inv(S(iExcluded,iExcluded))*S(iExcluded,i));
end

if isempty(MinMaxGrid)
    LowerBoundGrid=CenterPoint-4*ApproxStd;
    UpperBoundGrid=CenterPoint+4*ApproxStd;
    StepLength=(8*ApproxStd)/(NumberOfEval-1);
else
    for i=1:NumberOfParameters
        if any(isnan(MinMaxGrid(i,:)))
            LowerBoundGrid(i)=CenterPoint(i)-4*ApproxStd(i);
            UpperBoundGrid(i)=CenterPoint(i)+4*ApproxStd(i);
            StepLength(i)=(8*ApproxStd(i))/(NumberOfEval-1);
            if isnan(CenterPoint(i))
                errordlg(['Neither CenterPoint or the bounds for parameter ',EstParTeXNames(i,:),' are supplied. Exiting'])
                return;
            end
        else
            LowerBoundGrid(i)=MinMaxGrid(i,1);
            UpperBoundGrid(i)=MinMaxGrid(i,2);
            StepLength(i)=(UpperBoundGrid(i)-LowerBoundGrid(i))/(NumberOfEval-1);
            if isnan(CenterPoint(i))
                CenterPoint(i)=(LowerBoundGrid(i)+UpperBoundGrid(i))/2;
            end
            if CenterPoint(i)<LowerBoundGrid(i) | CenterPoint(i)>UpperBoundGrid(i)
                disp(['CenterPoint of parameter ',EstParTeXNames(i,:),'is outside user specified bounds. Setting CenterPoint to midpoint of bounds.'])
                CenterPoint(i)=(LowerBoundGrid(i)+UpperBoundGrid(i))/2;
            end
        end
    end
end

Minimum=inf;
Maximum=-inf;
for j=1:NFigures
    PriorFigHandles{j}=figure;
    k=0;
    if j<NFigures
        PlotCounter=1+(j-1)*NplotsPerFigure:j*NplotsPerFigure;
    else
        PlotCounter=NplotsPerFigure*(NFigures-1)+1:NumberOfParameters;
    end
    fprintf(1,'Slicing posterior for parameter: ')
    for i=PlotCounter
        fprintf(1,[int2str(i),',']) % Print to screen
        k=k+1;
        subplot(Nrows,Ncols,k)
        % Computes the log posterior around the supposed optimum
        GridAroundOpt=LowerBoundGrid(i):StepLength(i):UpperBoundGrid(i);
        if isempty(MinMaxGrid)
            if rem(NumberOfEval,2)==0
                GridAroundOpt=[GridAroundOpt(1:NumberOfEval/2) CenterPoint(i) GridAroundOpt(NumberOfEval/2+1:NumberOfEval)];
                CenterPointLoc=NumberOfEval/2+1;
            else
                GridAroundOpt=[GridAroundOpt(1:floor(NumberOfEval/2)) CenterPoint(i) GridAroundOpt(floor(NumberOfEval/2)+1:NumberOfEval)];
                CenterPointLoc=floor(NumberOfEval/2)+1;
            end
        else
            NLowerThanCenter=sum(GridAroundOpt<=CenterPoint(i));
            NHigherThanCenter=sum(GridAroundOpt>=CenterPoint(i));
            if NLowerThanCenter+NHigherThanCenter>length(GridAroundOpt) % Center is one of the grid points
                GridAroundOpt(GridAroundOpt~=CenterPoint(i)); % Exclude centerpoint
            end
            GridAroundOpt=[GridAroundOpt(1:NLowerThanCenter) CenterPoint(i) GridAroundOpt(NLowerThanCenter+1:length(GridAroundOpt))];
            CenterPointLoc=NLowerThanCenter+1;
        end
        LogPostOnGrid=zeros(1,length(GridAroundOpt));
        TempPoint=CenterPoint;
        u=0;
         
        LogPostOnGridMode=feval(ObjFuncStr,TempPoint,varargin{:});
        for q=GridAroundOpt
            u=u+1;
            TempPoint(i)=q;
            LogPostDens=feval(ObjFuncStr,TempPoint,varargin{:});
            LogPostOnGrid(u)=LogPostDens;
        end
        
        CloseToMode=(LogPostOnGridMode-LogPostOnGrid<10); % Plotting only grid point with log posterior
        LogPostOnGrid=LogPostOnGrid(CloseToMode==1);      % not more than 10 units smaller than Log posterior mode 
        GridAroundOpt=GridAroundOpt(CloseToMode==1);
        Minimum=min(Minimum,min(LogPostOnGrid));
        Maximum=max(Maximum,max(LogPostOnGrid));
        plot(GridAroundOpt,LogPostOnGrid,'.-')
        ApproxLogPDFMode=LogPdfNormal(CenterPoint(i),CenterPoint(i),ApproxStd(i));
        hold on
        if isempty(InvHessian)==0
            ApproxLogPDF=LogPdfNormal(GridAroundOpt,CenterPoint(i),ApproxStd(i));
            plot(GridAroundOpt,ApproxLogPDF-ApproxLogPDFMode+LogPostOnGridMode,'m:')
        end
        drawnow
     end
end
ParamCount=1;
for j=1:NFigures
    figure(PriorFigHandles{j})
    set(PriorFigHandles{j},'name',strcat('Checking Curvature at Posterior Mode no. ',int2str(j)))
    k=0;
    for Row=1:Nrows
        for Col=1:Ncols
            if ParamCount<=NumberOfParameters
                k=k+1;
                subplot(Nrows,Ncols,k)
                if SameVerticalScale
                    set(gca,'ylim',[Minimum Maximum]);
                else
                    YLimits=get(gca,'ylim');
                    Minimum=YLimits(1);
                    Maximum=YLimits(2);
                end
                LineHandle=line([CenterPoint(ParamCount) CenterPoint(ParamCount)],[Minimum Maximum]);
                set(LineHandle,'color','r')
                title(EstParTeXNames{ParamCount})
                set(gca,'fontsize',10)
                box off
                axis tight
            end
            ParamCount=ParamCount+1;
            if k==1
                legend('Log Conditional PDF - Perturbed','Log Conditional PDF - Normal Approx')
            end
        end
    end
end


fprintf(1,'Done!\n')