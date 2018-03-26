function plotCoherenceCluster(clustMat, clustLabels, basefilename)

fig = figure;
handles.baseMatrix = clustMat;
handles.labels = clustLabels;
handles.basefilename = basefilename;
handles.viewMinX = 1;
handles.viewMaxX = length(clustMat);

btn = uicontrol('String', 'zm in');
set(btn, 'Callback', @zoomIn);

btn2 = uicontrol('String', 'zm out');
btn2.Position(1) = btn.Position(1) + btn.Position(3);
set(btn2, 'Callback', @zoomOut);

btn3 = uicontrol('String', 'plot');
btn3.Position(1) = btn2.Position(1) + btn2.Position(3);
set(btn3, 'Callback', @plotConnections);

im = imagesc(clustMat);
set(gca, 'visible', 'off');
colorbar;
set(im, 'ButtonDownFcn', @handleImgClick);

guidata(fig, handles);

end

function plotConnections(src, evt)
handles = guidata(src.Parent);
minX = min(handles.x1, handles.x2);
maxX = max(handles.x1, handles.x2);
x1 = ceil(minX);
x2 = floor(maxX);

outsideSum = 0;
outsideCount = 0;
insideSum = 0;
insideCount = 0;
overlapSum = 0;
overlapCount = 0;

minX = min(x1, x2);
maxX = max(x1, x2);
labelPairs = handles.labels(minX:maxX);
matLength = length(handles.baseMatrix);

for i = 1:matLength
    iHit = (i >= minX && i <= maxX);
    for j = (i+1):matLength
        jHit = (j >= minX && j <= maxX);
        if(iHit && jHit)
            insideSum = insideSum + handles.baseMatrix(i,j);
            insideCount = insideCount + 1;
        elseif(xor(iHit, jHit))
            overlapSum = overlapSum + handles.baseMatrix(i,j);
            overlapCount = overlapCount + 1;
        else
            outsideSum = outsideSum + handles.baseMatrix(i,j);
            outsideCount = outsideCount + 1;
        end
    end
end

clustInfo.insideAvg = insideSum / insideCount;
clustInfo.outsideAvg = outsideSum / outsideCount;
clustInfo.overlapAvg = overlapSum / overlapCount;
clustInfo.segregation = clustInfo.insideAvg / clustInfo.overlapAvg;
clustInfo.count = insideCount;
clustInfo.channelPairLabels = labelPairs;
clustInfo.baseFilename = handles.basefilename;

outputFolder = '/home/data/EEG/processed/Robi/coherenceAnticorrelatedClusters/';
filename = handles.basefilename;
filename = filename(1:(strfind(filename, '_63')-1));
path = sprintf('%s%s(%d-%d).mat', outputFolder, filename, minX, maxX);
save(path, 'clustInfo');

plotChannelPairs(labelPairs);

guidata(src.Parent, handles);
clustInfo
end

function zoomIn(src, evt)
handles = guidata(src.Parent);
minX = min(handles.x1, handles.x2);
maxX = max(handles.x1, handles.x2);
xRange = maxX - minX;
viewMinX = round(minX - xRange * .3);
viewMaxX = round(maxX + xRange * .3);
viewMinX = max(1, viewMinX);
viewMaxX = min(viewMaxX, length(handles.baseMatrix));
handles.viewMinX = viewMinX;
handles.viewMaxX = viewMaxX;
winSize = viewMaxX - viewMinX;
line1 = minX - viewMinX + 1;
line2 = maxX - viewMinX + 1;
ind = viewMinX:viewMaxX;
im = imagesc(handles.baseMatrix(ind, ind));
set(im, 'ButtonDownFcn', @handleImgClick);
set(gca, 'visible', 'off');
handles.vLin1 = line([.5, winSize+5], [line1 line1], 'Color', 'r');
handles.hLin1 = line([line1 line1], [.5, winSize+5], 'Color', 'r');
handles.vLin2 = line([.5, winSize+5], [line2 line2], 'Color', 'r');
handles.hLin2 = line([line2 line2], [.5, winSize+5], 'Color', 'r');

guidata(src.Parent, handles);
dummy = 1;
end

function zoomOut(src, evt)
handles = guidata(src.Parent);
minX = min(handles.x1, handles.x2);
maxX = max(handles.x1, handles.x2);
handles.viewMinX = 1;
handles.viewMaxX = length(handles.baseMatrix);
winSize = handles.viewMaxX - handles.viewMinX;
line1 = minX - handles.viewMinX + 1;
line2 = maxX - handles.viewMinX + 1;
ind = handles.viewMinX:handles.viewMaxX;
im = imagesc(handles.baseMatrix(ind, ind));
set(im, 'ButtonDownFcn', @handleImgClick);
set(gca, 'visible', 'off');
handles.vLin1 = line([.5, winSize+5], [line1 line1], 'Color', 'r');
handles.hLin1 = line([line1 line1], [.5, winSize+5], 'Color', 'r');
handles.vLin2 = line([.5, winSize+5], [line2 line2], 'Color', 'r');
handles.hLin2 = line([line2 line2], [.5, winSize+5], 'Color', 'r');

guidata(src.Parent, handles);
dummy = 1;
end

function handleImgClick(src, evt)
x = floor(evt.IntersectionPoint(1)) + .5;
handles = guidata(src);
if(isfield(handles, 'x1'))
    if(handles.nextX == 2)
        handles.nextX = 1;
        handles.x2 = x + handles.viewMinX - 1;
        if(isfield(handles, 'vLin2'))
            delete(handles.vLin2);
            delete(handles.hLin2);
        end
    else
        handles.nextX = 2;
        handles.x1 = x + handles.viewMinX - 1;
        if(isfield(handles, 'vLin1'))
            delete(handles.vLin1);
            delete(handles.hLin1);
        end
    end
else
    handles.x1 = x + handles.viewMinX - 1;
    handles.nextX = 2;
end
vLin = line([x x], [.5, length(src.CData) + .5]);
set(vLin, 'Color', 'r');
hLin = line([.5, length(src.CData) + .5], [x x]);
set(hLin, 'Color', 'r');

if(handles.nextX == 1)
    handles.vLin2 = vLin;
    handles.hLin2 = hLin;
else
    handles.vLin1 = vLin;
    handles.hLin1 = hLin;
end

guidata(src, handles);
dummy = 1;
end

function testSearch(src)

handles = guidata(src);


matLength = length(handles.baseMatrix);
bestSeg.minX = -1;
bestSeg.maxX = -1;
bestSeg.segregation = realmin;
segs = NaN(matLength);
for minX = 1:matLength
    for maxX = (minX + 1):matLength
        fprintf('\n%d %d', minX, maxX);
        labelPairs = handles.labels(minX:maxX);
        outsideSum = 0;
        outsideCount = 0;
        insideSum = 0;
        insideCount = 0;
        overlapSum = 0;
        overlapCount = 0;
                
        for i = 1:matLength
            iHit = (i >= minX && i <= maxX);
            for j = (i+1):matLength
                jHit = (j >= minX && j <= maxX);
                if(iHit && jHit)
                    insideSum = insideSum + handles.baseMatrix(i,j);
                    insideCount = insideCount + 1;
                elseif(xor(iHit, jHit))
                    overlapSum = overlapSum + handles.baseMatrix(i,j);
                    overlapCount = overlapCount + 1;
                else
                    outsideSum = outsideSum + handles.baseMatrix(i,j);
                    outsideCount = outsideCount + 1;
                end
            end
        end
        
        clustInfo.insideAvg = insideSum / insideCount;
        clustInfo.outsideAvg = outsideSum / outsideCount;
        clustInfo.overlapAvg = overlapSum / overlapCount;
        clustInfo.segregation = clustInfo.insideAvg / clustInfo.overlapAvg;
        clustInfo.count = insideCount;
        clustInfo.channelPairLabels = labelPairs;
        clustInfo.baseFilename = handles.basefilename;
        segs(minX, maxX) = clustInfo.segregation;
        if(clustInfo.segregation > bestSeg.segregation)
            bestSeg.minX = minX;
            bestSeg.maxX = maxX;
            bestSeg.segregation = clustInfo.segregation;
            bestSeg.clustInfo = clustInfo;
        end
    end
end
dummy = 1;
scaledSegs = segs;
for i = 1:length(scaledSegs)
    for j = (i+1):length(scaledSegs)
        scaledSegs(i,j) = scaledSegs(i,j) * (j - i);
    end
end
logScaledSegs = log10(scaledSegs);
figure;
imagesc(logScaledSegs);

end
