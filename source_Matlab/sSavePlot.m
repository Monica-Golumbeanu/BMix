% This scripts saves the figure, Fig, in the address, OutFig, in .eps, and
% .fig formats and closes it.
%
% Written by Pejman Mohammadi
% pejman.mohammadi@bsse.ethz.ch
%

function sSavePlot(Fig, OutFig, Boxoff)
if nargin <3; Boxoff = false;end

saveas(Fig, OutFig)
if Boxoff;
    h = get(Fig,'Children');
    for i = 1:length(h)
        try
            subplot(h(i));
            legend boxoff;
        catch ignore
            
        end
    end
end
try
set(Fig, 'PaperPositionMode', 'auto');
catch err
end
print(Fig, '-depsc2', OutFig);
close(Fig)
disp(['Figure saved: ' OutFig] )
end