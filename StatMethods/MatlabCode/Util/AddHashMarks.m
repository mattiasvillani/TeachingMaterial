function []=AddHashMarks(FigHandle,HashMarks,HashColor,HashWidth)

switch nargin
    case 2
        HashColor=[0.6 0 0.4];
        HashWidth=1;
    case 3
        HashWidth=1;
    otherwise
end

if ~isempty(HashMarks)
    figure(FigHandle)
    hold on
    YLimits=get(gca,'ylim');
    HashLength=(YLimits(2)-YLimits(1))*0.02;
    for i=1:length(HashMarks)
        line([HashMarks(i),HashMarks(i)],[YLimits(1) YLimits(1)+HashLength],'color',HashColor,'linewidth',HashWidth);
    end
    %set(gca,'ylim',[YLimits(1)-HashLength/2 YLimits(2)])
end