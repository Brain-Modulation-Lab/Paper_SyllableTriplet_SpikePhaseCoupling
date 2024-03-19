function [] = export_graph(fig_name, fig_ext)
% EXPORT_GRAPH - saves figures to the hard-drive, if you want to adjust
% export settings (e.g. use PDF instead of PNG) change this here
% INPUT
%  fig_h       - figure handle
%  fig_name    - name of the figure, including full path

set(fig_h, 'PaperOrientation', 'portrait');
print('-r300', '-dpng', [fig_name, fig_ext]);
