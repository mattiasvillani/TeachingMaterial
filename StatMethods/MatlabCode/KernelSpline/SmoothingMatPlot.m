function DF = SmoothingMatPlot(L,x)


    % Hat matrix plot
    n = size(x,1);
    [SortX,SortIndex]=sort(x);
    L=L(SortIndex,SortIndex);
    DF=trace(L);
    
    FigHandle=figure;
    set(FigHandle,'name','Smoothing matrix')
    subplot(1,2,1)
    surfHandle = surf(L,'CDataMapping','scaled');
    set(surfHandle,'linestyle','none')
    axis square
    ylabel('yHat (obs are ordered wrt x')
    ylabel('y (obs are ordered wrt x')

    subplot(1,2,2)
    image(L,'CDataMapping','scaled')
    axis square
    colorbar
    box off
    ylabel('yHat (obs are ordered wrt x')
    ylabel('y (obs are ordered wrt x')

    % Equivalent kernels
    FigHandle=figure;
    set(FigHandle,'name','Equivalent Kernels - Smoothing matrix')
    NRows=24;
    PlotRows=ceil(1:(n/NRows):n);
    if PlotRows(end)<n
        PlotRows=[PlotRows n];
    end
    Count=0;
    for ii=PlotRows
        Count=Count+1;
        subplot(5,5,Count)
        plot(SortX,L(ii,:),'r')
        hold on
        LineHandle=line([SortX(ii) SortX(ii)],get(gca,'ylim'));
        set(LineHandle,'color','k','linestyle',':')
        LineHandle=line(get(gca,'xlim'),[0 0]);
        set(LineHandle,'color','k','linestyle',':')
        title(['Obs. No. ',int2str(ii)])
        set(gca,'xlim',[min(x) max(x)])
    end
