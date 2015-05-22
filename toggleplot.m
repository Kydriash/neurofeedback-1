function toggleplot(src,eventdata) %#ok

%global goodplots; 

%cellnum = getappdata(src,'PlotNumber'); 
%goodplots(cellnum) = ~goodplots(cellnum); 
src.Color = src.Color*0.75;
%set(src,'FaceAlpha',0.25 - get(src,'FaceAlpha')); 
% set(src,'EdgeAlpha',0.25 - get(src,'EdgeAlpha'));
