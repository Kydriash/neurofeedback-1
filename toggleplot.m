function toggleplot(src,eventdata) %#ok
global selected
%global goodplots; 

%cellnum = getappdata(src,'PlotNumber'); 
%goodplots(cellnum) = ~goodplots(cellnum); 
%if select
select = 1;
for k = 1:length(src.Parent.Children)
    if strcmp(src.Parent.Children(k).Tag,'Selection')
        delete(src.Parent.Children(k))
        parent_idx = src.Parent.Title.String;
        selected = selected(~strcmp(selected,parent_idx));
        select = 0;
        break;
    end
end
if select
    selected{end+1} = src.Parent.Title.String;
    ang=0:0.01:2*pi;
    x = (src.Parent.XLim(2) + src.Parent.XLim(1))/2;
    y = (src.Parent.YLim(2) + src.Parent.YLim(1))/2;
    xp=src.Parent.XLim(2)*cos(ang);
    yp=src.Parent.YLim(2)*sin(ang);
    plot(src.Parent,x+xp,y+yp,'Color',[0.5 0 1],'LineWidth',3,'Tag','Selection');
end


%else
    
%end
%set(src,'FaceAlpha',0.25 - get(src,'FaceAlpha')); 
% set(src,'EdgeAlpha',0.25 - get(src,'EdgeAlpha'));
end