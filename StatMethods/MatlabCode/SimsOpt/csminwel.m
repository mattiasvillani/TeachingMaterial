function [fh,xh,gh,H,itct,fcount,retcodeh] = csminwel(fcn,x0,H0,grad,crit,nit,varargin)
%[fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel(fcn,x0,H0,grad,crit,nit,varargin)
% fcn:   string naming the objective function to be minimized
% x0:    initial value of the parameter vector
% H0:    initial value for the inverse Hessian.  Must be positive definite.
% grad:  Either a string naming a function that calculates the gradient, or the null matrix.
%        If it's null, the program calculates a numerical gradient.  In this case fcn must
%        be written so that it can take a matrix argument and produce a row vector of values.
% crit:  Convergence criterion.  Iteration will cease when it proves impossible to improve the
%        function value by more than crit.
% nit:   Maximum number of iterations.
% varargin: A list of optional length of additional parameters that get handed off to fcn each
%        time it is called.
%        Note that if the program ends abnormally, it is possible to retrieve the current x,
%        f, and H from the files g1.mat and H.mat that are written at each iteration and at each
%        hessian update, respectively.  (When the routine hits certain kinds of difficulty, it
%        write g2.mat and g3.mat as well.  If all were written at about the same time, any of them
%        may be a decent starting point.  One can also start from the one with best function value.)

%% Mattias Villani 080107 - Monitor progress graphically
PlotSeq=1;
if PlotSeq
    ColorMapBasic=[1 0 0;0 1 0;0 0 1;0 0 0;0.5 0.5 0.5;1 1 0;0 1 1;1 0 1];
    LineaStyleBasic='-|:|--|.-';
    ZoomFactor=1.1;
    fCollect=nan;
    xCollect=[nan(1,length(x0))];
    IterCountMV=0;
    PlotfFigHandle=figure('name','Evolution of function evaluations');
    subplot(2,1,1)
    subplot(2,1,2)
    set(gcf,'DefaultAxesColorOrder',ColorMapBasic,'DefaultAxesLineStyleOrder',LineaStyleBasic);
    hold on
    ParamLegend=cell(1,length(x0));
    for i=1:length(x0)
        ParamLegend{i}=strcat('x',int2str(i));
    end
end
%%

[nx,no]=size(x0);
nx=max(nx,no);
Verbose=0;
NumGrad= isempty(grad);
done=0;
itct=0;
fcount=0;
snit=100;
%tailstr = ')';
%stailstr = [];
% Lines below make the number of Pi's optional.  This is inefficient, though, and precludes
% use of the matlab compiler.  Without them, we use feval and the number of Pi's must be
% changed with the editor for each application.  Places where this is required are marked
% with ARGLIST comments
%for i=nargin-6:-1:1
%   tailstr=[ ',P' num2str(i)  tailstr];
%   stailstr=[' P' num2str(i) stailstr];
%end
f0 = feval(fcn,x0,varargin{:});
%ARGLIST
%f0 = feval(fcn,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
% disp('first fcn in csminwel.m ----------------') % Jinill on 9/5/95
if f0 > 1e50, disp('Bad initial parameter.'), return, end
if NumGrad
   if length(grad)==0
      [g badg] = numgrad(fcn,x0, varargin{:});
      %ARGLIST
      %[g badg] = numgrad(fcn,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
   else
      badg=any(find(grad==0));
      g=grad;
   end
   %numgrad(fcn,x0,P1,P2,P3,P4);
else
   [g badg] = feval(grad,x0,varargin{:});
   %ARGLIST
   %[g badg] = feval(grad,x0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13);
end
retcode3=101;
x=x0;
f=f0;
H=H0;
cliff=0;
while ~done
   g1=[]; g2=[]; g3=[];
   %addition fj. 7/6/94 for control
   %disp('-----------------')
   %disp('-----------------')
   %disp('f and x at the beginning of new iteration')
   disp(sprintf('f at the beginning of new iteration, %20.10f',f))
   
   
   %% Added MV 080107 - plot of function evalutions
   if PlotSeq
       fCollect=[fCollect f];
       xCollect=[xCollect;x'];
       IterCountMV=[IterCountMV IterCountMV(end)+1];
       figure(PlotfFigHandle)
       NotExtreme=(fCollect/min(fCollect)<ZoomFactor);

       subplot(2,1,1)
       hold off
       plot(IterCountMV(NotExtreme),fCollect(NotExtreme))
       hold on
       plot(IterCountMV(end),fCollect(end),'ro')
       text(IterCountMV(end),fCollect(end),num2str(fCollect(end)),'horizontalalignment','right','verticalalignment','bottom','color','red')
       xlabel('Optimization iterations')
       ylabel('Function values')

       subplot(2,1,2)
       hold off
       plot(IterCountMV(NotExtreme),xCollect(NotExtreme,:))
       xlabel('Optimization iterations')
       ylabel('Parameter values')
       legend(ParamLegend,'location','westoutside')
   end
   %% End MV

   %-----------Comment out this line if the x vector is long----------------
      %disp([sprintf('x = ') sprintf('%15.8g %15.8g %15.8g %15.8g\n',x)]);
   %-------------------------
   itct=itct+1;
   [f1 x1 fc retcode1] = csminit(fcn,x,f,g,badg,H,varargin{:});
   %ARGLIST
   %[f1 x1 fc retcode1] = csminit(fcn,x,f,g,badg,H,P1,P2,P3,P4,P5,P6,P7,...
   %           P8,P9,P10,P11,P12,P13);
   % itct=itct+1;
   fcount = fcount+fc;
   % erased on 8/4/94
   % if (retcode == 1) | (abs(f1-f) < crit)
   %    done=1;
   % end
   % if itct > nit
   %    done = 1;
   %    retcode = -retcode;
   % end
   if retcode1 ~= 1
      if retcode1==2 | retcode1==4
         wall1=1; badg1=1;
      else
         if NumGrad
            [g1 badg1] = numgrad(fcn, x1,varargin{:});
            %ARGLIST
            %[g1 badg1] = numgrad(fcn, x1,P1,P2,P3,P4,P5,P6,P7,P8,P9,...
            %                P10,P11,P12,P13);
         else
            [g1 badg1] = feval(grad,x1,varargin{:});
            %ARGLIST
            %[g1 badg1] = feval(grad, x1,P1,P2,P3,P4,P5,P6,P7,P8,P9,...
            %                P10,P11,P12,P13);
         end
         wall1=badg1;
         % g1
         save g1 g1 x1 f1 varargin;
         %ARGLIST
         %save g1 g1 x1 f1 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
      end
      if wall1 & (length(H) > 1)% 
         % Bad gradient or back and forth on step length.  Possibly at
         % cliff edge.  Try perturbing search direction if problem not 1D
         %
         %fcliff=fh;xcliff=xh;
         Hcliff=H+diag(diag(H).*rand(nx,1));
         disp('Cliff.  Perturbing search direction.')
         [f2 x2 fc retcode2] = csminit(fcn,x,f,g,badg,Hcliff,varargin{:});
         %ARGLIST
         %[f2 x2 fc retcode2] = csminit(fcn,x,f,g,badg,Hcliff,P1,P2,P3,P4,...
         %     P5,P6,P7,P8,P9,P10,P11,P12,P13);
         fcount = fcount+fc; % put by Jinill
         if  f2 < f
            if retcode2==2 | retcode2==4
                  wall2=1; badg2=1;
            else
               if NumGrad
                  [g2 badg2] = numgrad(fcn, x2,varargin{:});
                  %ARGLIST
                  %[g2 badg2] = numgrad(fcn, x2,P1,P2,P3,P4,P5,P6,P7,P8,...
                  %      P9,P10,P11,P12,P13);
               else
                  [g2 badg2] = feval(grad,x2,varargin{:});
                  %ARGLIST
                  %[g2 badg2] = feval(grad,x2,P1,P2,P3,P4,P5,P6,P7,P8,...
                  %      P9,P10,P11,P12,P13);
               end
               wall2=badg2;
               % g2
%                badg2
               save g2 g2 x2 f2 varargin
               %ARGLIST
               %save g2 g2 x2 f2 P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 P12 P13;
            end
            if wall2
               disp('Cliff again.  Try traversing')
               if norm(x2-x1) < 1e-13
                  f3=f; x3=x; badg3=1;retcode3=101;
               else
                  gcliff=((f2-f1)/((norm(x2-x1))^2))*(x2-x1);
                  if(size(x0,2)>1), gcliff=gcliff', end
                  [f3 x3 fc retcode3] = csminit(fcn,x,f,gcliff,0,eye(nx),varargin{:});
                  %ARGLIST
                  %[f3 x3 fc retcode3] = csminit(fcn,x,f,gcliff,0,eye(nx),P1,P2,P3,...
                  %         P4,P5,P6,P7,P8,...
                  %      P9,P10,P11,P12,P13);
                  fcount = fcount+fc; % put by Jinill
                  if retcode3==2 | retcode3==4
                     wall3=1; badg3=1;
                  else
                     if NumGrad
                        [g3 badg3] = numgrad(fcn, x3,varargin{:});
                        %ARGLIST
                        %[g3 badg3] = numgrad(fcn, x3,P1,P2,P3,P4,P5,P6,P7,P8,...
                        %                        P9,P10,P11,P12,P13);
                     else
           