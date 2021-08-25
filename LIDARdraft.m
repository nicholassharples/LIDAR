Files = dir("50cmDTMWindsor/*.asc");

%MDX Uni College Building Front Door OS grid reference = TQ 22828 89360
% TQ  is 51 so in numeric grid reference this is 522828 189360
% This is in tile tq2289.


%%
cellsize = 0.5;

X = zeros(0,1000/cellsize);
Y = zeros(0,1000/cellsize);
data = zeros(0,1000/cellsize);


for k = 1:length(Files)
    FileName = fullfile(Files(k).folder, Files(k).name) %% ToDo: construct the path rather than filename.
    
    opts = detectImportOptions(FileName, "FileType", "delimitedtext");
    opts.Delimiter = " ";
    opts.VariableTypes = "double";

    opts.DataLines = [1 6];
    metadata = readmatrix(FileName, opts)

    opts.DataLines = 7;
    datatemp = readmatrix(FileName, opts); 

    datatemp(datatemp==-9999) = NaN;

    x = metadata(3,5) + (0:cellsize:(1000-cellsize));
    y = metadata(4,5) + (0:cellsize:(1000-cellsize));

    [Xtemp,Ytemp] = meshgrid(x,y);    

    X = [X Xtemp];
    Y = [Y Ytemp];
    data = [data datatemp];
end  

%%
file = "50cm/tq2289_DSM_50cm.asc";


opts = detectImportOptions(file, "FileType", "delimitedtext");
opts.Delimiter = " ";
opts.VariableTypes = "double";

opts.DataLines = [1 6];
metadata = readmatrix(file, opts);

opts.DataLines = 7;
data = readmatrix(file, opts); 

data(data==-9999) = NaN;

%% Get points for plotting
x = metadata(3,5) + (0:0.5:999.5);
y = metadata(4,5) + (0:0.5:999.5);

[X,Y] = meshgrid(x,y);

%% Plot
s = surf(X,Y, data);
s.EdgeColor = 'none';

%% Define the kernel for smoothing
kernel = ones(10,10);

%% Plot smoothed
s = surf(X,Y,conv2(data,kernel,'same'));
s.EdgeColor = 'none';


%%
for k = 1:length(Files)
    FileName = fullfile(Files(k).folder, Files(k).name)
    tmp = readlines(FileName);
    tmp(1:6)
end

%% Now GeoTIFFs
file = "1M2020DTM/TQ28nw_DTM_1m.tif";
[A,R,cmap] = readgeoraster(file, "OutputType","Double");

%[B,RB] = mapcrop(A,R,[522828-300 522828+100], [189360-200 189360+200]);
%[B,RB] = mapcrop(A,R,[522588 522736], [189000 189500]); % Close up of campus

[B,RB] = mapcrop(A,R,[522000 523200], [188400 189600]); % Hill up Hendon way


lon = RB.XWorldLimits(1):(RB.XWorldLimits(2)-1);
lat = RB.YWorldLimits(1):(RB.YWorldLimits(2)-1);
lat = flip(lat); %% It looks like the raster import array B doesn't index in the "obvious" way.
[X, Y] = meshgrid(lon,lat);

%% ShapeFile from OpenStreetMapping
% Need to grab shapes inside our ROI. But: GeoTIFF is in OS Geoid
% Coordinates (meters) and shapefile is in lat/lon.
% Solution: convert OSM geoid values to lat/lon
[latmin,lonmin] = projinv(R.ProjectedCRS, 522000, 188500);
[latmax,lonmax] = projinv(R.ProjectedCRS, 523000, 189500);


file = "greater-london-latest-free.shp/gis_osm_buildings_a_free_1.shp"; %% All of London - massive!
BuildingShapes = shaperead(file, "BoundingBox", [lonmin, latmin; lonmax, latmax], "UseGeoCoords", true);

file = "greater-london-latest-free.shp/gis_osm_roads_free_1.shp"; %% All of London - massive!
RoadShapes = shaperead(file, "BoundingBox", [lonmin, latmin; lonmax, latmax], "UseGeoCoords", true);



% There's a crazy discrepency here - the shapefile Lat/Lon is accurate, but
% if we convert it to x,y using the LIDAR/OS coordinate system R then we
% get a significant drift compared to the LIDAR data. Is the LIDAR data
% inaccurate?





%%
%mapshow(S)
%%
% shapelon = [S.Lon];
% shapelat = [S.Lat];
% shapetype = {S(:).type};
% 
% [shapex,shapey] = projfwd(R.ProjectedCRS, shapelat, shapelon);


%% Project ShapeFile Lat/Lon onto LIDAR/OS grid references
for index = 1:length(BuildingShapes)
    [BuildingShapes(index).X, BuildingShapes(index).Y] = projfwd(RB.ProjectedCRS, BuildingShapes(index).Lat, BuildingShapes(index).Lon);
end

%% Correct drift
% College building front entrance is sent to 522720 189418
% It should be 522831 189361
for index = 1:length(BuildingShapes)
    BuildingShapes(index).X = BuildingShapes(index).X + 111;
    BuildingShapes(index).Y = BuildingShapes(index).Y - 57;
end

%% Interpolate a Z-coordinates for the shapefile
%shapez = interp2(X,Y,B,shapex,shapey);
for index = 1:length(BuildingShapes)
    BuildingShapes(index).Z = interp2(X,Y, B, BuildingShapes(index).X, BuildingShapes(index).Y);
end

%% Project ShapeFile Lat/Lon onto LIDAR/OS grid references
for index = 1:length(RoadShapes)
    [RoadShapes(index).X, RoadShapes(index).Y] = projfwd(RB.ProjectedCRS, RoadShapes(index).Lat, RoadShapes(index).Lon);
end

%% Correct drift
% College building front entrance is sent to 522720 189418
% It should be 522831 189361
for index = 1:length(RoadShapes)
    RoadShapes(index).X = RoadShapes(index).X + 111;
    RoadShapes(index).Y = RoadShapes(index).Y - 57;
end

%% Interpolate a Z-coordinates for the shapefile
%shapez = interp2(X,Y,B,shapex,shapey);
for index = 1:length(RoadShapes)
    RoadShapes(index).Z = interp2(X,Y, B, RoadShapes(index).X, RoadShapes(index).Y);
end

%% Prepare MDX and Hendon Central coordinates
MDX.X = 522831;
MDX.Y = 189361;
MDX.Z = interp2(X,Y, B, MDX.X, MDX.Y);

HenC.X = 522980;
HenC.Y = 188658;
HenC.Z = interp2(X,Y, B, HenC.X, HenC.Y);

%% Plot
h = figure;
hold on;
mapshow(B,RB,"DisplayType","surface");
ax = gca;
demcmap(B);
BuildingMap = mapshow(BuildingShapes);
RoadMap = mapshow(RoadShapes);
%M.ZData = shapez ;
plot3(MDX.X, MDX.Y, MDX.Z + 20 ,"o", "Color", "#e30613","MarkerFaceColor", "#e30613", "MarkerSize",10);
plot3([MDX.X MDX.X], [MDX.Y MDX.Y], [MDX.Z MDX.Z+20], "Color", "#e30613");
plot3(HenC.X, HenC.Y, HenC.Z + 20 ,"o", "Color", "b", "MarkerFaceColor", "b", "MarkerSize",10);
plot3([HenC.X HenC.X], [HenC.Y HenC.Y], [HenC.Z HenC.Z+20],"Color", "b");
axis normal
%mapshow(double(B),RB,"DisplayType","contour","LineColor","k","ShowText","on")
view(3)
camorbit(ax,0,15)
hold off;

%% Add Zdata
for index = 1:length(BuildingShapes)
    %M.Children(index).ZData = repelem(max(S(index).Z), length(S(index).Z)) + 2;
    BuildingMap.Children(index).ZData = BuildingShapes(index).Z + 2;
end

for index = 1:length(RoadShapes)
    %M.Children(index).ZData = repelem(max(S(index).Z), length(S(index).Z)) + 2;
    RoadMap.Children(index).ZData = RoadShapes(index).Z + 2;
end

%% Animate
videofile = "Hendon.gif";
for framenumber = 1:360
    camorbit(ax,1,0)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256); 
    if framenumber == 1 
          imwrite(imind,cm,videofile,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,videofile,'gif','WriteMode','append'); 
    end    
end 

 for framenumber = 1:30
    camorbit(ax,0,1)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);  
    imwrite(imind,cm,videofile,'gif','WriteMode','append'); 
 end

 for framenumber = 1:30
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);  
    imwrite(imind,cm,videofile,'gif','WriteMode','append'); 
 end
 
%% Get 2D images  
view(2);
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,"Hendon.png","png");
delete(BuildingMap);
delete(RoadMap);
hold on;
contour = mapshow(B,RB,"DisplayType","contour","LineColor","k","ShowText","on");
hold off;

%% Get 2D images
h = figure;
hold on;
mapshow(B,RB,"DisplayType","texturemap");
ax = gca;
demcmap(B);
BuildingMap = mapshow(BuildingShapes);
RoadMap = mapshow(RoadShapes);
plot(MDX.X, MDX.Y, "o", "Color", "#e30613","MarkerFaceColor", "#e30613", "MarkerSize",10);
plot(HenC.X, HenC.Y, "o", "Color", "b", "MarkerFaceColor", "b", "MarkerSize",10);
axis normal
%mapshow(double(B),RB,"DisplayType","contour","LineColor","k","ShowText","on")
view(2)
xlim([522000,523200]);
ylim([188400 189600]);

frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,"Hendon.png","png");
delete(BuildingMap);
delete(RoadMap);
contour = mapshow(B,RB,"DisplayType","contour","LineColor","k","ShowText","on");
hold off;
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,"HendonContour.png","png");

 %%
mapshow(double(B),RB,"DisplayType","texturemap");

%%
hold on
mapshow(double(B),RB,"DisplayType","surface"); %% Correct orientation
view(3)
axis normal
surf(X,Y,double(B)) %% Correct orientation
hold off


%% Fourier fitting

F = fft2(double(B));

%% Take first n modes by magnitude
n = 1000000;
sortedabsF = sort(reshape(F,[],1),"descend");

G = ifft2((abs(F)>=sortedabsF(n)).*F);

%% Plot
hold on
terrain = surf(double(B),"EdgeColor","none","FaceAlpha",0.5,"FaceColor","g");
model = surf(G,"EdgeColor","none","FaceAlpha", 0.5, "FaceColor", "r");
view(3);
hold off

%% Polynomial fitting
f = fit([reshape(X,[],1),reshape(Y,[],1)], reshape(B,[],1), "poly44");
g = fit([reshape(X,[],1)-MDX.X, reshape(Y,[],1)-MDX.Y], reshape(B,[],1), "poly44");

%%
%centredg = 

%% Plot
h = figure;
hold on
terrain = mapshow(B,RB,"DisplayType","surface");
ax= gca;
demcmap(B);
%terrain = surf(X,Y, B,"EdgeColor","none","FaceAlpha",0.5,"FaceColor","g");
%model = plot(f);
model2 = plot3(reshape(X,[],1), reshape(Y,[],1), g([reshape(X,[],1)-MDX.X, reshape(Y,[],1)-MDX.Y]) );
view(3);
terrain.FaceAlpha = 0.5;
%model.FaceAlpha = 0.5;
%model.FaceColor = "#2f2552";
model2.Color = "#2f2552"
xlim([522000,523200]);
ylim([188400 189600]);
camorbit(ax,0,15)
hold off

%% Animate
videofile = "HendonModel.gif";
for framenumber = 1:360
    camorbit(ax,1,0)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256); 
    if framenumber == 1 
          imwrite(imind,cm,videofile,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,videofile,'gif','WriteMode','append'); 
    end    
end 

 for framenumber = 1:30
    camorbit(ax,0,1)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);  
    imwrite(imind,cm,videofile,'gif','WriteMode','append'); 
 end

 for framenumber = 1:60
    camorbit(ax,0,-1)
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);  
    imwrite(imind,cm,videofile,'gif','WriteMode','append'); 
 end

 %%
 f
 
 