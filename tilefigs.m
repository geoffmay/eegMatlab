function tilefigs(dims)
%TILEFIGS Tile all open figure windows around on the screen.
% TILEFIGS places all open figure windows around on the screen with no
% overlap. TILEFIGS(FIGS) can be used to specify which figures that
% should be tiled. Figures are not sorted when specified.
%
% 1998-01-22 19:40:40 Peter J. Acklam <jacklam@math.uio.no>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the handles to the figures to process.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figs = findobj( 'Type', 'figure' ); % ...find all figures.
numbers = [figs.Number];
titles = cell(length(figs), 1);
for i = 1:length(figs)
    a = figs(i).Name;
    chd = figs(i).Children;
    for j = 1:length(chd)
        if(isa(chd(j), 'matlab.graphics.axis.Axes'))
            a = get(chd(j), 'title');
            titles{i} = a.String;
        end
    end
end
[~, indT] = sort(titles);
[~, ind] = sort(numbers);
figs = figs(ind);

if isempty( figs )
   disp( 'No open figures.' );
   return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set miscellaneous parameter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfigs = length(figs); % Number of figures.

if(~exist('dims', 'var'))
nh = ceil( sqrt(nfigs) ); % Number of figures horisontally.
nv = ceil( nfigs/nh ); % Number of figures vertically.

nh = max( nh, 2 );
nv = max( nv, 2 );
else
    if(~isnumeric(dims))
        error('dimensions must be numeric: [horizontal, vertical]');
    end
    nh = dims(1);
    nv = dims(2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the screen size.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set( 0, 'Units', 'pixels' ); % Set root units.
scrdim = get( 0, 'ScreenSize' ); % Get screen size.
scrwid = scrdim(3); % Screen width.
scrhgt = scrdim(4); % Screen height.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The elements in the vector specifying the position.
% 1 - Window left position
% 2 - Window bottom position
% 3 - Window horizontal size
% 4 - Window vertical size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hspc = 50; % Horisontal space.
topspc = 80; % Space above top figure.
medspc = 80; % Space between figures.
botspc = 40; % Space below bottom figure.

figwid = ( scrwid - (nh+1)*hspc )/nh;
fighgt = ( scrhgt - (topspc+botspc) - (nv-1)*medspc )/nv;
% figwid = ( scrwid - (nh)*hspc )/nh;
% fighgt = ( scrhgt - (topspc+botspc) - (nv)*medspc )/nv;

figwid = round(figwid);
fighgt = round(fighgt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Put the figures where they belong.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for row = 1:nv
   for col = 1:nh
      idx = (row-1)*nh + col;
      if ( idx <= nfigs )
         figlft = col*hspc + (col-1)*figwid;
         figbot = scrhgt - topspc - row*fighgt - (row-1)*medspc;
         figpos = [ figlft figbot figwid fighgt ]; % Figure position.
         fighnd = figs(idx); % Figure handle.
         units = get( fighnd, 'Units' ); % Get units.
         set( fighnd, 'Units', 'pixels' ); % Set new units.
         set( fighnd, 'Position', figpos ); % Set position.
         figure( fighnd ); % Raise figure.
      end
   end
end