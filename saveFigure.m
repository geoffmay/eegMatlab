function [ output_args ] = saveFigure( handle, fileName )
% saveFigure 
%   Saves figure specified by `handle` as `fileName` in fullscreen
%   as to get around the stupid behavior.

    screen_size = get(0, 'ScreenSize');
    origSize = get(handle, 'Position'); % grab original on screen size
    set(handle, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to scren size
    set(handle,'PaperPositionMode','auto') %set paper pos for printing
    saveas(handle, fileName) % save figure
    set(handle,'Position', origSize) %set back to original dimensions

end