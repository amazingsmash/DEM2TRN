%TRN adjustment to DEM data
%Code by Jose Miguel Santana Nunez

function [mx, my, mz] = TRN2DEM(x, y, z, ... %DEM
                                TRNresX, TRNresY, ... %TRN size in mesh vertices
                                tileSize, ... %Length of tile side in mesh vertices (including frame)
                                f) % Frame width in mesh vertices
                            
showAnimation = 0; %Showing intermediate results

%Defining TRN
xs = linspace(min(x(:)), max(x(:)), TRNresX);
ys = linspace(min(y(:)), max(y(:)), TRNresY);
[mx, my] = meshgrid(xs, ys);
mz = NaN(TRNresY, TRNresX);

%We store here added values for accurate averaging
summz = zeros(TRNresY, TRNresX);
countmz = summz;

if (showAnimation)
    animFigure = figure;
    myaxis = [min(x(:)) max(x(:)) min(y(:)) max(y(:))];
    axis(myaxis);
end

innerWindowSize = tileSize-2*f;

for i = 1:innerWindowSize-1:TRNresX
    
    tileXs = max(i - f, 1): min(i + innerWindowSize - 1 + f, TRNresX);
    centerTileXs = find(tileXs == i) : find(tileXs == i + innerWindowSize - 1);
    if isempty(centerTileXs)
        centerTileXs = find(tileXs == i) : length(tileXs);
    end
    
    for j = 1:innerWindowSize-1:TRNresY
        
        tileYs = max(j - f, 1): min(j + innerWindowSize - 1 + f, TRNresY);
        centerTileYs = find(tileYs == j) : find(tileYs == j + innerWindowSize - 1);
        if isempty(centerTileYs)
            centerTileYs = find(tileYs == j) : length(tileYs);
        end
        
        %Mesh tile
        pmx = mx(tileXs, tileYs);
        pmy = my(tileXs, tileYs);
        
        %DEM tile
        select = x >= min(pmx(:)) & x <= max(pmx(:)) & ...
              y >= min(pmy(:)) & y <= max(pmy(:));
        demSizeX = max(sum(select));
        demSizeY = max(sum(select, 2));
        
        DEMtileX = reshape(x(select), demSizeX, demSizeY);
        DEMtileY = reshape(y(select), demSizeX, demSizeY);
        DEMtileZ = reshape(z(select), demSizeX, demSizeY);

        %Adjusting tile
        pmz = TRN2DEMTileAdjustment(DEMtileX, DEMtileY, DEMtileZ, pmx, pmy);

        %Extracing inner window
        pmx = pmx(centerTileXs,centerTileYs);
        pmy = pmy(centerTileXs,centerTileYs);
        pmz = pmz(centerTileXs,centerTileYs);

        %Finding window in supermatrix
        pos = ismember(mx, pmx) & ismember(my, pmy);
        pos = find(pos);
        
        %Display
        if (showAnimation)
            mz(pos) = pmz; %Placing values for display
            clf
            figure(animFigure);
            surf(mx, my, mz);
            axis(myaxis);
            pause(0.3); %Frame pause
        end
        
        %Data aggregation
        summz(pos) = summz(pos) + pmz(:);
        countmz(pos) = countmz(pos) + 1;
        
    end
end

%Calculating final result
mz = summz ./ countmz; %Averaging
end

%Tile adjustment algorithm
%TRN heights as output
function mz = TRN2DEMTileAdjustment(x, y, z, mx, my)

if std(z) == 0
    mz = repmat(z(1,1), size(mx));
    return;
end

[A, B] = createLinearSystemForDEMAdjustment(x, y, z, mx, my);

%First approximation X0 = BILINEAR
x0 = interp2(x,y,z, mx, my, 'linear');
x0 = x0(:);
x0(isnan(x0)) = 0;

%Approximating with LSQR
X = lsqr(A, B,1e-10, 100, [], [], x0); 
mz = reshape(X, size(mx));

%Iterative refinement of the error
tol = 1e-10;

iniError = mean(abs(A * X - B));
fprintf('Initial Mean Error: %f\n', iniError);

anyModified = 1;
while anyModified
    
    anyModified = 0;
    errors = getPointError(A, X, B);
    [~,is] = sort(errors);
    is = is(end:-1:1)';
    is = is(abs(errors(is)) > tol)';
    fprintf('-----\n');
    
    for i = is
        [X, meanErrorGainI] = updateVertex(A, B, X, errors, i);
        if meanErrorGainI > tol
            anyModified = 1;
        end
    end
    currentError = mean(abs(A*X-B));
    fprintf('Mean Error After Point By Point Displacement: %f\n', currentError);
end

mz = reshape(X, [size(mz,1), size(mz,2)]);

finalError = mean(abs(A * X - B));
fprintf('Final Mean Error: %f, (Gain: %E) \n', mean(abs(A * X - B)), iniError -finalError);
end

%Update vertex height considering previous error
%New mesh and new mean error as result
function [X, meanErrorGain] = updateVertex(A, B, X, errors, i)
lastError = mean(abs(A * X - B));

prevValueI = X(i);

X(i) = prevValueI - errors(i);
minError = mean(abs(A * X - B));

if minError < lastError
    meanErrorGain = lastError - minError;
else
    X(i) = prevValueI;
    meanErrorGain = 0;
end

end

function errors = getPointError(A, X, B)

e = (A * X - B);
errors = zeros(1,size(A, 2));

for i = 1:size(A, 2) %All points in MZ
    %Error in DEM point i multiplied by the contribution of that point to mesh vertex
    errors(i) = sum(A(:,i) .* e) / sum(A(:,i));
end

end

%Creating linear system A x X = B describing adjustment errors
function [A, B] = createLinearSystemForDEMAdjustment(x, y, z, mx, my)

if sum(size(z) ~= size(x)) || sum(size(z)~= size(y)) || sum(size(mx) ~= size(my))
    warning('Dimension mistmach');
end

%Quadrants and displacement in quadrant
%Only if DEM is regular
% [qx, qy] = meshgrid(linspace(1,size(mx,1), size(z,1)),... 
%                    linspace(1,size(mx,2), size(z,2))); 

qx = (((x - min(min(x))) / (max(max(x)) - min(min(x)))) * (size(mx,2)-1))+1;
qy = (((y - min(min(y))) / (max(max(y)) - min(min(y)))) * (size(my,1)-1))+1;
                
dqx = qx;
dqy = qy;
qx = floor(qx);
qy = floor(qy);
dqx = dqx - qx;
dqy = dqy - qy;
dqx(:,end) = 1.0;
dqy(end,:) = 1.0;
qx(:,end) = qx(:,end-1);
qy(end,:) = qy(end-1,:);

%One quadrant triangulation
P = [0, 0;
    1, 0;
    0, 1;
    1, 1];
T = [1 2 4;
     1 3 4];
TR = triangulation(T,P);

%Upper or lower triangle in quadrant
ut = dqx < dqy;
ut = ut + 1; %1 = upper triangle, 2 = lower triangle

%Barycentric coordinates of all DEM points
ps = [dqx(:), dqy(:)];
bary = cartesianToBarycentric(TR,ut(:),ps);

%All barycentrinc coordinates must be >=0

%Composing A
axs = 1:length(bary); %Vertex index
axs = [axs' axs' axs'];

ays = zeros(length(bary), 3); %Triangle vertices indexes
mrl = size(mx,1); %Mesh col length
vind = qy + (qx-1) * mrl;
for i = 1:length(bary)
    vi = vind(i); %Upper left mesh vertex of quadrant
    if ut(i) == 1 %Upper triangle
        ays(i, :) = [vi, vi + mrl, vi + 1 + mrl]; %1 2 4
    else
        ays(i, :) = [vi, vi + 1, vi + 1 + mrl];% 1 3 4
    end
end

nrow = numel(z); %one row per DEM point
ncol = numel(mx); %one col per mesh vertex
gb = bary ~= 0;
axs = axs(gb);
ays = ays(gb);
bary = bary(gb);
A = sparse(axs, ays, bary, nrow, ncol); %Using matlabs sparse matrix

%A*X=B
B = sparse(1:numel(z), 1, z(:));

end