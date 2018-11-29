function fig = plotChannelPairs(labelPairs, pairWeights, colorLimits, drawColorbar, fig)
%PLOTCHANNELPAIRS Plots pairs of channels using the 'line' function.
%   The function attempts to minimize collinearity by using "jitter", which
%   skews the x values away from the center based on how far the y values
%   are from the center, and vice-versa.

%setting this to zero executes new code that adapts to larger electrode
%sets.
fixedElectrodeSet = 0; 

if(fixedElectrodeSet)
    baseLabels = antChannelLocs;
    mastoids = [find(strcmp(baseLabels, 'M1')), find(strcmp(baseLabels, 'M2'))];
    baseLabels(mastoids) = [];
    baseLabels(31:end) = [];
    x1s = [4 5 6 2 3 5 7 8 2 4 6 8 1 3 5 7 9 2 4 6 8 2 3 5 7 8 5 4 5 6];
    y1s = [1 1 1 2 2 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 6 6 6 6 6 7 8 8 8];
    y1s = 9 - y1s;
else
    %this file contains electrode locations in the 10-5 system, in the
    %'chanlocs' variable.
    load('chanlocs.mat');
    
    %find unique channel labels, assign to "baseLabels"
    counter = 1;    
    for i = 1:length(labelPairs)
        pair = labelPairs{i};
        %remove frequency info (delimited by a space)
        space = strfind(pair, ' ');
        pair(space(2):end) = [];
        pair(1:space(1)) = [];
        %get individual channels (separated by a '-')
        channels = strsplit(pair, '-');
        for j = 1:length(channels)
            allChannels{counter} = channels{j};
            counter = counter + 1;
        end
    end
    baseLabels = unique(allChannels);
    
    %x1s and y1s are the coordinates without any offset; they have a one 
    %to one pairing with baseLabels    
    for i = 1:length(baseLabels)
        index = find(strcmp({chanlocs.labels}, baseLabels{i}));
        x1s(i) = chanlocs(index).X;
        y1s(i) = chanlocs(index).Y;
    end
end

for i = 1:9
    jitter(i) = 16.25 + sqrt(16.25*16.25 - ((i-5)*(i-5))) - 32;
end
if(true)
    if(fixedElectrodeSet)
        offset = 2;
        xMin = 1;
        xMax = 9;
        xCenter = 5;
    else
        xMin = min(x1s);
        xMax = max(x1s);
        yMin = min(y1s);
        yMax = max(y1s);
        xCenter = (xMin + xMax) / 2;
        yCenter = (yMin + yMax) / 2;
        offset = (xMax - xMin) / 10;
    end
    if(fixedElectrodeSet)
        [centre radius] = calc_circle([xMin offset], [xMax offset], [xCenter yCenter]);
        y0 = radius;
        x0 = xCenter;
        for i = 1:9
            jitter(i) = sqrt(y0*y0 - ((i-x0)*(i-x0))) - y0;
        end
        for i = 1:length(x1s)
            x2s(i) = x1s(i) + jitter(y1s(i)) * .2*(x1s(i)-4.5);
            y2s(i) = y1s(i) + jitter(x1s(i)) * .2*(y1s(i)-4.5);
        end
    else
        [centre radius] = calc_circle([xMin offset], [xMax offset], [xCenter yCenter]);
        for i = 1:length(x1s)
            %'jitter' is additional deviation from the center; a circle is
            %drawn with its center at the outermost edge of the plot
            circleOffset = centre(2);
            delX = x1s(i) - xCenter;
            delY = y1s(i) - yCenter;
            xJitter = sqrt(radius * radius - delY * delY) - circleOffset;
            yJitter = sqrt(radius * radius - delY * delY) - circleOffset;
            %scale numbers span from 1 to -1. jitter closer to the center
            %is dampened.
            xScale = -1 + (1 - -1)/(xMax - xMin)*(x1s(i) - xMin);
            yScale = -1 + (1 - -1)/(yMax - yMin)*(y1s(i) - yMin);
            x2s(i) = x1s(i) + xJitter * xScale;
            y2s(i) = y1s(i) + yJitter * yScale;
        end
        %debug
%         tab = table(baseLabels', x1s', y1s', x2s', y2s');
        %end debug
        %flip.
        temp = x2s;
        x2s = -y2s;
        y2s = temp;
        %debug
        x2s = -y1s;
        y2s = x1s;
        %end debug
    end
end

if(~exist('fig','var'))
fig = figure;
end
axis off;


%make colormap
finished = false;
if(exist('pairWeights', 'var'))
    colorMap = jet(256);
    maxWeight = max(pairWeights);
    minWeight = min(pairWeights);
    if(~exist('colorLimits', 'var'))
        colorLimits = [minWeight maxWeight];
    end
    rangeWeight = colorLimits(2) - colorLimits(1);
    colors = NaN(length(pairWeights),3);
    for i = 1:length(pairWeights)
        mapIndex = floor((pairWeights(i) - colorLimits(1)) / rangeWeight * 255) + 1;
        if(mapIndex > size(colorMap,1))
            mapIndex = size(colorMap,1);
        elseif(mapIndex < 1)
            mapIndex = 1;
        end
        colors(i,:)  = colorMap(mapIndex, :);
    end
    
    caxis(colorLimits);
    %     colorMap = colormap('jet');
    if(false)
        colorbar('Location', 'north');
        print('colorbar', '-dpng');
    end
    if(~exist('drawColorbar','var'))
        drawColorbar = true;
    end
    if(drawColorbar)
      colormap('jet');
        colorbar();        
    end
end


for i = 1:length(x2s)
    text(x2s(i), y2s(i), baseLabels{i});
end
test = labelPairs{1};
if(ischar(test))
    toPlot = labelPairs;
    color = [0 1 1];
else
    toPlotIndex = 1;
    toPlot = labelPairs{toPlotIndex};
    colorMap = jet(length(labelPairs));
    color = colorMap(toPlotIndex,:);
end
% toPlot = strrep(toPlot, 'T3', 'T7');
% toPlot = strrep(toPlot, 'T4', 'T8');
% toPlot = strrep(toPlot, 'T5', 'P7');
% toPlot = strrep(toPlot, 'T6', 'P8');


%drop excess from labels
for i = 1:length(toPlot)
    label = toPlot{i};
    space = strfind(label, ' ');
    if(length(space) > 1)
        label(space(2):end) = [];
        label(1:space(1)) = [];
        toPlot{i} = label;
    end
end

%draw actual plot
while(~finished)
    for i = 1:length(toPlot)
        if(exist('pairWeights', 'var'))
            color = colors(i,:);
        end
        if(length(strfind(toPlot{i}, '-')) > 0)
          labels = strsplit(toPlot{i}, '-');
          li1 = find(strcmp(baseLabels, labels(1)));
          li2 = find(strcmp(baseLabels, labels(2)));
          line([x2s(li1), x2s(li2)], [y2s(li1), y2s(li2)], 'Color', color);
        end
    end
    
    if(ischar(test))
        finished = true;
    else
        toPlotIndex = toPlotIndex + 1;
        if(toPlotIndex > length(labelPairs))
            finished = true;
        else
            toPlot = labelPairs{toPlotIndex};
            color = colors(toPlotIndex,:);
        end
    end
end

end


function [centre radius] = calc_circle(pt1, pt2, pt3)

delta_a = pt2 - pt1;
delta_b = pt3 - pt2;

ax_is_0 = abs(delta_a(1)) <= 0.000000001;
bx_is_0 = abs(delta_b(1)) <= 0.000000001;

% check whether both lines are vertical - collinear
if (ax_is_0 && bx_is_0)
    centre = [0 0];
    radius = -1;
    return
end

% make sure delta gradients are not vertical
% re-arrange points to change deltas
if (ax_is_0)
    [centre radius] = calc_circle(pt1, pt3, pt2);
    return
end
if (bx_is_0)
    [centre radius] = calc_circle(pt2, pt1, pt3);
    return
end

grad_a = delta_a(2) / delta_a(1);
grad_b = delta_b(2) / delta_b(1);

% check whether the given points are collinear
if (abs(grad_a-grad_b) <= 0.000000001)
    centre = [0 0];
    radius = -1;
    return
end

% swap grads and points if grad_a is 0
if abs(grad_a) <= 0.000000001
    tmp = grad_a;
    grad_a = grad_b;
    grad_b = tmp;
    tmp = pt1;
    pt1 = pt3;
    pt3 = tmp;
end

% calculate centre - where the lines perpendicular to the centre of
% segments a and b intersect.
centre(1) = ( grad_a*grad_b*(pt1(2)-pt3(2)) + ...
grad_b*(pt1(1)+pt2(1)) - grad_a*(pt2(1)+pt3(1)) ) /...
(2*(grad_b-grad_a));
centre(2) = ((pt1(1)+pt2(1))/2 - centre(1)) / grad_a +...
(pt1(2)+pt2(2))/2;

% calculate radius
radius = norm(centre - pt1); 
end