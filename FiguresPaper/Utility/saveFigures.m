function saveFigures(fh, figname, ext)

if ~exist('ext', 'var')
    saveas(fh,[figname, '.fig'])
    saveas(fh,[figname, '.png'])
    saveas(fh,[figname, '.pdf'])
    saveas(fh,[figname, '.svg'])
else
    for ext_i = 1 : numel(ext)
        saveas(fh,[figname, ext{ext_i}])
    end
end
