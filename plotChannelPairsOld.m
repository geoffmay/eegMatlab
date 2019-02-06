function plotChannelPairs(labelPairs, pairWeights, colorLimits, drawColorbar)

baseLabels = antChannelLocs;
mastoids = [find(strcmp(baseLabels, 'M1')), find(strcmp(baseLabels, 'M2'))];
baseLabels(mastoids) = [];
x1s = [4 5 6 2 3 5 7 8 2 4 6 8 1 3 5 7 9 2 4 6 8 2 3 5 7 8 5 4 5 6];
y1s = [1 1 1 2 2 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5 6 6 6 6 6 7 8 8 8];
y1s = 9 - y1s;

for i = 1:9
    jitter(i) = 16.25 + sqrt(16.25*16.25 - ((i-5)*(i-5))) - 32;
end
if(true)
    offset = 2;
    [centre radius] = calc_circle([1 offset], [9 offset], [5 0]);
    y0 = radius;
    x0 = 5;
    for i = 1:9
        jitter(i) = sqrt(y0*y0 - ((i-x0)*(i-x0))) - y0;
    end    
end
for i = 1:length(x1s)
        x2s(i) = x1s(i) - jitter(y1s(i)) * .2*(x1s(i)-4.5);
        y2s(i) = y1s(i) - jitter(x1s(i)) * .2*(y1s(i)-4.5);
end

fig = figure;
axis off;

if(exist('pairWeights', 'var'))
    caxis(colorLimits);
    colors = colormap('jet');
    if(false)
        colorbar('Location', 'north');
        print('colorbar', '-dpng');
    end
    if(~exist('drawColorbar','var'))
        drawColorbar = true;
    end
    if(drawColorbar)
        colorbar;        
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
    colors = jet(length(labelPairs));
    color = colors(toPlotIndex,:);
end
toPlot = strrep(toPlot, 'T3', 'T7');
toPlot = strrep(toPlot, 'T4', 'T8');
toPlot = strrep(toPlot, 'T5', 'P7');
toPlot = strrep(toPlot, 'T6', 'P8');


finished = false;
if(exist('pairWeights', 'var'))
    map = jet(256);
    maxWeight = max(pairWeights);
    minWeight = min(pairWeights);
    if(~exist('colorLimits', 'var'))
        colorLimits = [minWeight maxWeight];
    end
    rangeWeight = colorLimits(2) - colorLimits(1);
    colors = NaN(length(pairWeights),3);
    for i = 1:length(pairWeights)
        mapIndex = floor((pairWeights(i) - colorLimits(1)) / rangeWeight * 255) + 1;
        colors(i,:)  = map(mapIndex, :);
    end
end

while(~finished)
    for i = 1:length(toPlot)
        if(exist('pairWeights', 'var'))
            color = colors(i,:);
        end
        labels = strsplit(toPlot{i}, '-');
        li1 = find(strcmp(baseLabels, labels(1)));
        li2 = find(strcmp(baseLabels, labels(2)));
        line([x2s(li1), x2s(li2)], [y2s(li1), y2s(li2)], 'Color', color);
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