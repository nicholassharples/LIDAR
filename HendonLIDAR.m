%% Uses Matlab's Mapping Toolbox

%% Load LIDAR data
% First we get the LIDAR data. This is in the form of a GeoTIFF which can be found at
% https://environment.data.gov.uk/DefraDataDownload/?Mode=survey
% Draw the region of interest, select "LIDAR Composite DTM" and choose the
% most recent data. Data after 2019 is in the required GeoTIFF format.
% Older data uses ASCII.

file = "1M2020DTM/TQ28nw_DTM_1m.tif";
[A,R,cmap] = readgeoraster(file, "OutputType","Double");

%This is a large 5km*5km region, so we crop it down to our area of
%interest. All coordinates are in Ordnance Survery Geoid coordinates
%https://getoutside.ordnancesurvey.co.uk/guides/beginners-guide-to-grid-references/


[B,RB] = mapcrop(A,R,[522000 523200], [188400 189600]); % Hill up Hendon way


Xpoints = RB.XWorldLimits(1):(RB.XWorldLimits(2)-1);
Ypoints = RB.YWorldLimits(1):(RB.YWorldLimits(2)-1);
% Accessing the array directly the top-left element corresponds with the
% bottom-left geographic point. So we flip the array.
Ypoints = flip(Ypoints); %%
[X, Y] = meshgrid(Xpoints,Ypoints);

%% Load building/road data from OpenStreetMapping
% This data is in the form of a shapefile, which can be found at, for
% example,
% http://download.geofabrik.de/europe/great-britain/england/greater-london.html

% This file has OSM data for all of London, which is huge!
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
% get a significant drift compared to the LIDAR data.
% Is the projection information inside the LIDAR GeoTIFF inaccurate?




%% We now prepare the ShapeFile data for plotting ontop of the LIDAR data
% First the Building Shapes
% Project ShapeFile Lat/Lon onto LIDAR/OS grid references
for index = 1:length(BuildingShapes)
    [BuildingShapes(index).X, BuildingShapes(index).Y] = projfwd(RB.ProjectedCRS, BuildingShapes(index).Lat, BuildingShapes(index).Lon);
end

% Correct drift
% College building front entrance is sent to 522720 189418
% It should be 522831 189361
for index = 1:length(BuildingShapes)
    BuildingShapes(index).X = BuildingShapes(index).X + 111;
    BuildingShapes(index).Y = BuildingShapes(index).Y - 57;
end

% Interpolate a Z-coordinates for the shapefile
for index = 1:length(BuildingShapes)
    BuildingShapes(index).Z = interp2(X,Y, B, BuildingShapes(index).X, BuildingShapes(index).Y);
end


% Now the RoadShapes
% Project ShapeFile Lat/Lon onto LIDAR/OS grid references
for index = 1:length(RoadShapes)
    [RoadShapes(index).X, RoadShapes(index).Y] = projfwd(RB.ProjectedCRS, RoadShapes(index).Lat, RoadShapes(index).Lon);
end

% Correct drift
% College building front entrance is sent to 522720 189418
% It should be 522831 189361
for index = 1:length(RoadShapes)
    RoadShapes(index).X = RoadShapes(index).X + 111;
    RoadShapes(index).Y = RoadShapes(index).Y - 57;
end

% Interpolate a Z-coordinates for the shapefile
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
plot3(MDX.X, MDX.Y, MDX.Z + 20 ,"o", "Color", "#e30613","MarkerFaceColor", "#e30613", "MarkerSize",10); % Middlesex red
plot3([MDX.X MDX.X], [MDX.Y MDX.Y], [MDX.Z MDX.Z+20], "Color", "#e30613"); % Middlesex red
plot3(HenC.X, HenC.Y, HenC.Z + 20 ,"o", "Color", "b", "MarkerFaceColor", "b", "MarkerSize",10);
plot3([HenC.X HenC.X], [HenC.Y HenC.Y], [HenC.Z HenC.Z+20],"Color", "b");
axis normal
view(3)
camorbit(ax,0,15)
hold off;

% Add Zdata
for index = 1:length(BuildingShapes)
    %M.Children(index).ZData = repelem(max(S(index).Z), length(S(index).Z)) + 2;
    BuildingMap.Children(index).ZData = BuildingShapes(index).Z + 2;
end

for index = 1:length(RoadShapes)
    %M.Children(index).ZData = repelem(max(S(index).Z), length(S(index).Z)) + 2;
    RoadMap.Children(index).ZData = RoadShapes(index).Z + 2;
end

%% Animate
% This will take a while...
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

% Get 2D images
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



%% Polynomial fitting
% We want to describe the altitudes using a function, so that we can do
% some calculus on it. To do this we will fit a degree 4 multivariable
% polynomial to this data.
g = fit([reshape(X,[],1)-MDX.X, reshape(Y,[],1)-MDX.Y], reshape(B,[],1), "poly44");

% Note that we're recentred the data so that the origin is at the College
% Building main entrance. This means that the x,y data are in the range
% [-1000,1000] which is easier to manage.

%Get details of the model.
g

%% Plot
h = figure;
hold on
terrain = mapshow(B,RB,"DisplayType","surface");
ax= gca;
demcmap(B);
model2 = plot3(reshape(X,[],1), reshape(Y,[],1), g([reshape(X,[],1)-MDX.X, reshape(Y,[],1)-MDX.Y]) );
view(3);
terrain.FaceAlpha = 0.5;
model2.Color = "#2f2552" % Middlesex Indigo
xlim([522000,523200]);
ylim([188400 189600]);
camorbit(ax,0,15)
hold off

%% Animate
% This will take a while... 
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
