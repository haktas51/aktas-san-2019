%%Automatic Generating Seed Cell Polygons Algorithm for 2LRS method
%%Code by: Hakan Aktaþ and Bekir Taner San
%%email: hakanaktas5151@gmail.com
%%The algorithm has 3 main parts. First part is the preprocessing part. 
%%This part ends in line 156. The second part of the code is the main part 
%%where the seed ceel extracting algorithm works. This part finish in line 472. 
%%Final part of the code is the 2 Level Random Sampling Method. 
%%In the last part the final aim is to generate the array_learner, 
%%arrayImage_test_ls and arrayImage_test_non_ls data sets. array_learner is the 
%%data set for training process. By opening matlab classification learner toolbox
%%user can use this train data set during training operations.
%%Many figure codes are commented, the user uncomment the figure code lines
%%and display the results.
clc
close all
clear all
%Reading Multispectral Image and Generating DEM file
ms  = imread('0765_aster_sub.tif'); % the folder in which ur files exists
DEM = ms(:,:,15);
% Plesae enter the percantage of Train and Test points 
Train_points_percentage = 50;
Test_poits_percentage = 100 - Train_points_percentage;
point_percentage = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading all landslides from a current folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirName = 'C:\Users\VA\Desktop\Generating_Seed_Cell\shp_files';  % the folder in which ur files exists
srcFiles= dir( fullfile(dirName,'*.shp') ); 
for i = 1 : length(srcFiles)
    filename = strcat('C:\Users\VA\Desktop\Generating_Seed_Cell\shp_files\',srcFiles(i).name);
    Slides(i) = shaperead (filename); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading all landslides buffers from a current folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirName_2 = 'C:\Users\VA\Desktop\Generating_Seed_Cell\shp_files_buf';  % the folder in which ur files exists
srcFiles_2 = dir( fullfile(dirName_2,'*.shp')); 
for i = 1 : length(srcFiles_2)
    filename = strcat('C:\Users\VA\Desktop\Generating_Seed_Cell\shp_files_buf\',srcFiles_2(i).name);
    Buffers(i) = shaperead (filename); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-generating buffer files. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : length(Buffers)
     for j = 1 : length(Buffers(i).X)-1
         if isnan(Buffers(i).X(j))
             break
         end
         j;
         Only_Buffers(i).X(j) = Buffers(i).X(j);
         Only_Buffers(i).Y(j) = Buffers(i).Y(j);
     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lat_begin = 254260.798; 
lon_begin = 4083499.701;
size_of_y = 1927;
size_of_x = 778;
pixel_length = 15;
number_of_polygons = length(Buffers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating working_area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_working_area = shaperead('working_area_borders.shp');
lat_working_area = [S_working_area.X];  
lon_working_area = [S_working_area.Y]; 

format longG
working_area = zeros(size_of_y,size_of_x);
i_limit = length(lat_working_area)-1 ;
for i = 1 : i_limit-1
    lat_diff = (lat_working_area(i)-lat_begin) / pixel_length;
    if (lat_diff < 1)
        lat_diff = 1;
    end
    lat_diff = round(lat_diff);
    lon_diff = ( lon_begin - lon_working_area(i)) / pixel_length;
    if (lon_diff < 1)
        lon_diff = 1;
    end
    lon_diff = round(lon_diff);
    working_area (lon_diff,lat_diff) = 1;
end
% Display and saving working area, uncomment the lines if needed. 
% imshow(working_area)
% title('Working Area')
% print(gcf,'figure1.tif','-dtiffn','-r600');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating zeros of slide polygon, buffer_polygon and only_buf_structs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slide_lat_lon = zeros(size_of_y,size_of_x);
for j=1:length(Buffers)
    shape_struct(j).slides = slide_lat_lon;
end
for j=1:length(Buffers)
    only_buf_shape_struct(j).slides = slide_lat_lon;
end
for j=1:length(Buffers)
    interested_buf_struct(j).final = slide_lat_lon;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating slide polygon points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:number_of_polygons
	S_land_slide = Slides(j); 
	l = length(S_land_slide.X);
	x = (S_land_slide.X);
	y = (S_land_slide.Y);

	for i = 1:l-1
		slide_lat_diff = ( x(i) - lat_begin ) / pixel_length;
		if (slide_lat_diff < 1)
			slide_lat_diff = 1;
		end
		slide_lat_diff = round(slide_lat_diff);  
		slide_lon_diff = ( lon_begin - y(i)) / pixel_length;
		if (slide_lon_diff < 1)
			slide_lon_diff = 1;
		end
		slide_lon_diff = round(slide_lon_diff);       
		shape_struct(j).slides (slide_lon_diff,slide_lat_diff) = 1;
		shape_struct(j).dem(i) = DEM(slide_lon_diff-1,slide_lat_diff);
		shape_struct(j).y(i) = slide_lon_diff;
		shape_struct(j).x(i) = slide_lat_diff;
	end	 
end
% figure
% imshow(shape_struct(j).slides + working_area)
% title([' ' num2str(j) '. landslide' ])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating landslide-buffer points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:number_of_polygons
	Buf_land_slide = Only_Buffers(j);
    l = length(Buf_land_slide.X)-1;
    x = (Buf_land_slide.X);
    y = (Buf_land_slide.Y);
	for i = 1:l
		slide_lat_diff = ( x(i) - lat_begin ) / pixel_length;
        if (slide_lat_diff < 1)
			slide_lat_diff = 1;
        end
        slide_lat_diff = round(slide_lat_diff);  
        slide_lon_diff = ( lon_begin - y(i)) / pixel_length;
        if (slide_lon_diff < 1)
			slide_lon_diff = 1;
		end
        slide_lon_diff = round(slide_lon_diff);       
        only_buf_shape_struct(j).slides (slide_lon_diff,slide_lat_diff) = 1;
        only_buf_shape_struct(j).dem(i) = DEM(slide_lon_diff,slide_lat_diff);
        only_buf_shape_struct(j).y(i) = slide_lon_diff;
        only_buf_shape_struct(j).x(i) = slide_lat_diff;
    end
end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Here is the main Process of the algorithm. The main proess is done in only
%%one for loop which is mentioned below. With this main process Generating Seed
%%Cell Polygons is done for all Landslide polygons. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
mask_interested_bufs = zeros(size_of_y,size_of_x);
landslide_landslidebuf_workingarea = zeros(size_of_y,size_of_x);

for main_tour = 1: number_of_polygons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Finding max min points of slides 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j= main_tour;
index_1 = 1;
index_2 = 1;
find_max = 0;
find_min = 100000;
for l = 4:length(shape_struct(j).dem)-3  	     
    if shape_struct(j).dem(l-3)+shape_struct(j).dem(l-2)+shape_struct(j).dem(l-1)+ shape_struct(j).dem(l)+shape_struct(j).dem(l+1)+shape_struct(j).dem(l+2)+shape_struct(j).dem(l+3) > find_max
        find_max = shape_struct(j).dem(l-3)+shape_struct(j).dem(l-2)+shape_struct(j).dem(l-1)+ shape_struct(j).dem(l)+shape_struct(j).dem(l+1)+shape_struct(j).dem(l+2)+shape_struct(j).dem(l+3);
        index_1 = l;
    else 
        find_max = find_max;
        index_1 = index_1;
    end 
    if  shape_struct(j).dem(l-3)+shape_struct(j).dem(l-2)+shape_struct(j).dem(l-1)+ shape_struct(j).dem(l)+shape_struct(j).dem(l+1)+shape_struct(j).dem(l+2)+shape_struct(j).dem(l+3) < find_min
        find_min = shape_struct(j).dem(l-3)+shape_struct(j).dem(l-2)+shape_struct(j).dem(l-1)+ shape_struct(j).dem(l)+shape_struct(j).dem(l+1)+shape_struct(j).dem(l+2)+shape_struct(j).dem(l+3);
        index_2 = l;
    else 
        find_min = find_min;
        index_2 = index_2;
    end  
	end
max_min_struct(j).max.max = shape_struct(j).dem(index_1);
max_min_struct(j).max.lat = shape_struct(j).x(index_1);
max_min_struct(j).max.lon = shape_struct(j).y(index_1);
%
max_min_struct(j).min.min = shape_struct(j).dem(index_2);
max_min_struct(j).min.lat = shape_struct(j).x(index_2);
max_min_struct(j).min.lon = shape_struct(j).y(index_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drawing a line from max to min 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_to_min = zeros(size_of_y,size_of_x); % (y,x) = (lon, lat)
theta = atan2(max_min_struct(j).max.lon - max_min_struct(j).min.lon , max_min_struct(j).max.lat - max_min_struct(j).min.lat);
r = sqrt(( max_min_struct(j).max.lon - max_min_struct(j).min.lon)^2 + (max_min_struct(j).max.lat - max_min_struct(j).min.lat)^2 );
for line = 0:0.1: r;
x = max_min_struct(j).min.lat + line*cos(theta);
x  = round(x);
y = max_min_struct(j).min.lon + line*sin(theta);
y = round(y);
max_to_min(y,x)=1;
end
% figure
% imshow(max_to_min)
% title('Verticle Line')
slide_plus_max_to_min = shape_struct(main_tour).slides + max_to_min;
% figure
% imshow(slide_plus_max_to_min)
% title([' ' num2str(main_tour) '.landslide + max to min line' ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finding line formula from max to min point... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j= main_tour;
draw_struct(j).adding_lat= abs(round(2*(max_min_struct(j).max.lat - max_min_struct(j).min.lat)/3));
draw_struct(j).adding_lon= abs(round(2*(max_min_struct(j).max.lon - max_min_struct(j).min.lon)/3));
draw_struct(j).m1 = (max_min_struct(j).max.lon - max_min_struct(j).min.lon)/ (max_min_struct(j).max.lat- max_min_struct(j).min.lat);
if draw_struct(j).m1 == 0
    draw_struct(j).m1 = 0.1;
end
if max_min_struct(j).max.lat > max_min_struct(j).min.lat
    draw_struct(j).expected_lat = max_min_struct(j).max.lat - draw_struct(j).adding_lat;
end
if max_min_struct(j).max.lat < max_min_struct(j).min.lat
        draw_struct(j).expected_lat = max_min_struct(j).max.lat + draw_struct(j).adding_lat;
end
if max_min_struct(j).max.lat == max_min_struct(j).min.lat
        draw_struct(j).expected_lat = max_min_struct(j).max.lat + draw_struct(j).adding_lat;
end
if max_min_struct(j).max.lon > max_min_struct(j).min.lon
    draw_struct(j).expected_lon = max_min_struct(j).max.lon - draw_struct(j).adding_lon;
end
if max_min_struct(j).max.lon < max_min_struct(j).min.lon
        draw_struct(j).expected_lon = max_min_struct(j).max.lon + draw_struct(j).adding_lon;
end
if max_min_struct(j).max.lon == max_min_struct(j).min.lon
        draw_struct(j).expected_lon = max_min_struct(j).max.lon + draw_struct(j).adding_lon;
end    
draw_struct(j).m2 = -1/draw_struct(j).m1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw a circle around mid point 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
centx = draw_struct(j).expected_lat;
centy = draw_struct(j).expected_lon;
r = 80; % user can define the radius
% figure
% imshow(slide_plus_max_to_min)
% title('landslide + max to min line + circle')
% hold on;
theta = 0 : (2 * pi / 1000) : (2 * pi);
circleline_x = r * cos(theta) + centx;
circleline_y = r * sin(theta) + centy;
% plot(ceil(circleline_x),ceil(circleline_y), 'r-', 'LineWidth', 2);
% plot(centx,centy,'+','MarkerSize',8);
% hold off;
% end
% line formula can be shown as in below
% y-y1 = m2(x-x1)
% y-(draw_struct(j).expected_lon) = m2(x- draw_struct(j).expected_lat )
loop = length(circleline_y);
a=1;
for i = 1:loop-1
	if ( (( circleline_y(i)- draw_struct(j).expected_lon) - (draw_struct(j).m2*(circleline_x(i) - draw_struct(j).expected_lat)) ) < 6.9 &...
    (( circleline_y(i) - draw_struct(j).expected_lon) -  (draw_struct(j).m2*(circleline_x(i) - draw_struct(j).expected_lat)) ) > -6.9)  
    % ideally 6.9 should be one where exact 2 points fit the line formula
    % however using one can cause some problems. Because of different
    % shape formats of landslides when one is used sometimes no points are
    % found that fit this line formula. User can increase this number
    % reasonably then, more than two points will fit to this line formula.
    % After that these points will be checked one by one and two best
    % fitting points will be found by the code automatically in next steps by using max_p1_p2_distance parameter.
	a1_struct(main_tour).a1(a) = circleline_x(i);
	a2_struct(main_tour).a2(a) = circleline_y(i);
	a = a+1;
 end
end
max_p1_p2_distance = 0;
for p1_p2_tour = 2:length( a1_struct(main_tour).a1)
	distance = sqrt(( a1_struct(main_tour).a1(1) - a1_struct(main_tour).a1(p1_p2_tour))^2 + (a2_struct(main_tour).a2(1) - a2_struct(main_tour).a2(p1_p2_tour))^2) ;    
	if distance > max_p1_p2_distance
	max_p1_p2_distance = distance;
	index_p1_p2 = p1_p2_tour;
	end
end   
p1 = [a2_struct(main_tour).a2(1),a1_struct(main_tour).a1(1)] ; % 
p2 = [a2_struct(main_tour).a2(index_p1_p2),a1_struct(main_tour).a1(index_p1_p2)] ; % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drawing a line from p1 to p2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawline = zeros(size_of_y,size_of_x);
theta = atan2( p2(1) - p1(1), p2(2) - p1(2));
r2 = sqrt( (p2(1) - p1(1))^2 + (p2(2) - p1(2))^2);
for line = 0:0.1: r2;
	x = p1(2) + line*cos(theta);
	x  = round(x);
	y = p1(1) + line*sin(theta);
	y = round(y);
	drawline(y,x)=1;
end
% figure
% imshow(drawline)
% title('Verticle Line')
buffer_plus_draw = only_buf_shape_struct(main_tour).slides + drawline + slide_plus_max_to_min;
% figure
% imshow(slide_plus_max_to_min + working_area ) %buffer_plus_draw   %slide_plus_max_to_min + working_area
% title([' ' num2str(main_tour) '.landslide + circle + verle line'])
% hold on;
% theta = 0 : (2 * pi / 1000) : (2 * pi);
% circleline_x = r * cos(theta) + centx;
% circleline_y = r * sin(theta) + centy;
% plot(ceil(circleline_x),ceil(circleline_y), 'r-', 'LineWidth', 2);
% plot(centx,centy,'+','MarkerSize',8);
% hold off;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acording to p1 and p2 poitns Finding interested points on land slide buffers... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interested_points_binary_struct(j).points = drawline & only_buf_shape_struct(j).slides; 
% after &(logical and) operation more than 2 points can be found. To find the best
% points some resonable operations are done in code below. 
% figure
% imshow(interested_points_binary_struct(j).points)
% title('p1 and p2 interested points')
interested_points_location_struct(j).locations = find(interested_points_binary_struct(j).points) ;
for i=1:length(interested_points_location_struct(j).locations)
    lon = mod(interested_points_location_struct(j).locations(i),size_of_y);
    lat = round(interested_points_location_struct(j).locations(i)/size_of_y)+1;
	interested_points_lat_struct(j).lat(i) = lat;
	interested_points_lon_strcut(j).lon(i) = lon;
end       
first_point_lat_struct(j).first_point_lat = interested_points_lat_struct(j).lat(1) ;
second_point_lat_strcut(j).second_point_lat = interested_points_lat_struct(j).lat(length(interested_points_lat_struct(j).lat));
first_point_lon_strcut(j).first_point_lon = interested_points_lon_strcut(j).lon(1) ;
second_point_lon_struct(j).second_point_lon = interested_points_lon_strcut(j).lon(length(interested_points_lon_strcut(j).lon));
% after &(logical and) operation more than 2 points can be found. To find the best
% points some resonable operations are done in code below. First case is
% the points should have the same lat or lon values. 
first_point_index =1;
control = 0; %control statement is used during algorithm tests 
for search_first_point =1:length(only_buf_shape_struct(main_tour).x)
    if (only_buf_shape_struct(main_tour).x(search_first_point)== first_point_lat_struct(j).first_point_lat) & ...
       (only_buf_shape_struct(main_tour).y(search_first_point)== first_point_lon_strcut(j).first_point_lon)
        finding_first_point_begin_struct(j).finding_first_point_begin = first_point_index;
        control = control+1;
    else
        first_point_index = first_point_index+1;
    end
end
% after &(logical and) operation more than 2 points can be found. To find the best
% points some resonable operations are done in code below. Second case is
% the difference of two points should be one. This shows the closest two
% points. Then one of this point can be choosen. 
first_point_index =1;
if control == 0 
for search_first_point =1:length(only_buf_shape_struct(main_tour).x)
    if abs(only_buf_shape_struct(main_tour).x(search_first_point) - first_point_lat_struct(j).first_point_lat) == 1 & ...
       abs(only_buf_shape_struct(main_tour).y(search_first_point)- first_point_lon_strcut(j).first_point_lon) == 1
        finding_first_point_begin_struct(j).finding_first_point_begin = first_point_index;
        control = control+1;
    else
        first_point_index = first_point_index+1;
    end
end
end
second_point_index =1;
control_2 = 0; %control statement is used during algorithm tests
for search_second_point =1:length(only_buf_shape_struct(main_tour).x)
    if (only_buf_shape_struct(main_tour).x(search_second_point)== second_point_lat_strcut(j).second_point_lat) & ...
       (only_buf_shape_struct(main_tour).y(search_second_point)== second_point_lon_struct(j).second_point_lon)
        finding_second_point_begin_struct(j).finding_second_point_begin = second_point_index;
        control_2 = control_2+1;
    else
        second_point_index = second_point_index+1;
    end
end
second_point_index =1;
if control_2 == 0
for search_second_point =1:length(only_buf_shape_struct(main_tour).x)
    if abs(only_buf_shape_struct(main_tour).x(search_second_point)- second_point_lat_strcut(j).second_point_lat) == 1 & ...
       abs(only_buf_shape_struct(main_tour).y(search_second_point) - second_point_lon_struct(j).second_point_lon)== 1
        finding_second_point_begin_struct(j).finding_second_point_begin = second_point_index;
        control_2 = control_2+1;
    else
        second_point_index = second_point_index+1;
    end
end
end
% After finding these two points, using matlab features the interested
% points can be found as shown in below.
interested_buf = zeros(size_of_y,size_of_x);
if  (finding_first_point_begin_struct(j).finding_first_point_begin <  finding_second_point_begin_struct(j).finding_second_point_begin)
    finding_first_point = finding_first_point_begin_struct(j).finding_first_point_begin;
    finding_second_point = finding_second_point_begin_struct(j).finding_second_point_begin;
end
if  (finding_first_point_begin_struct(j).finding_first_point_begin >  finding_second_point_begin_struct(j).finding_second_point_begin)
    finding_first_point = finding_second_point_begin_struct(j).finding_second_point_begin;
    finding_second_point = finding_first_point_begin_struct(j).finding_first_point_begin;
end
if (finding_second_point - finding_first_point) > round(length(only_buf_shape_struct(main_tour).x)/2)
    for i = finding_first_point : finding_second_point;
        interested_buf(only_buf_shape_struct(main_tour).y(i),only_buf_shape_struct(main_tour).x(i))=1;
    end
end
if (finding_second_point - finding_first_point) < round(length(only_buf_shape_struct(main_tour).x)/2)
    for i = 1: finding_first_point
        interested_buf(only_buf_shape_struct(main_tour).y(i),only_buf_shape_struct(main_tour).x(i))=1;
    end
    for i = finding_second_point:length(only_buf_shape_struct(main_tour).x)
         interested_buf(only_buf_shape_struct(main_tour).y(i),only_buf_shape_struct(main_tour).x(i))=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final arithmetic operations 
% To generate exact seed cell polygons some final arithmetic operations are
% needed. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shape_struct_filled = imfill(shape_struct(main_tour).slides,'holes');
% figure 
% imshow(shape_struct_filled)
% title('shape struct filled')
only_buf_filled = only_buf_shape_struct(main_tour).slides;
only_buf_filled_struct(main_tour).filled = imfill(only_buf_filled,'holes');
% figure
% imshow(only_buf_filled_struct(main_tour).filled)
% title('only buf filled')
% figure
% imshow(drawline + only_buf_shape_struct(main_tour).slides + max_to_min)
% title('Expected Buffer Points Location')
interested_buf_with_drawline = interested_buf+drawline;
% figure
% imshow(interested_buf_with_drawline)
% title('interested buf with drawline')
interested_buf_filled = imfill(interested_buf_with_drawline,'holes');
% figure
% imshow(interested_buf_filled)
% title('interested buf filled')
buf_minus_slide_filled = only_buf_filled_struct(main_tour).filled - shape_struct_filled;
% figure
% imshow (buf_minus_slide_filled)
% title('buf filled minus lanslide filled')
interested_buf_final = (buf_minus_slide_filled & interested_buf_filled);
% figure
% imshow(interested_buf_final)
% title(['Index: ' num2str(j) '. mask' ])
interested_buf_struct(main_tour).final = interested_buf_final;
mask_interested_bufs = mask_interested_bufs + interested_buf_struct(main_tour).final;
landslide_landslidebuf_workingarea = landslide_landslidebuf_workingarea + shape_struct(main_tour).slides + only_buf_shape_struct(main_tour).slides + working_area;
% figure 
% imshow(landslide_landslidebuf_workingarea);
% title('landslide + landslide buf + working area')
% figure
% imshow(mask_interested_bufs + working_area)
% title('All Buffer Masks in Working Area')
% figure
% imshow(landslide_landslidebuf_workingarea)
% title('landslide + landslide buf + working area')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rest of the code from here can be custom code. The rest of the code can
% be implemented by any user acording to their neeeds. Here is just an
% example how to generate train and test points by using 2LRS algorithm. 
% Generating test and training buffers with random numbers. This is first
% random sampling of the 2LRS algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Training_buffers = zeros(size_of_y,size_of_x);
Test_buffers = zeros(size_of_y,size_of_x);
mask_train_buffers_struct = zeros(size_of_y,size_of_x);
mask_train_buffers_struct_total = zeros(size_of_y,size_of_x);
mask_train_buffer_points_struct = zeros(size_of_y,size_of_x);
mask_train_buffer_points_struct_total= zeros(size_of_y,size_of_x);
%
mask_test_buffers_struct = zeros(size_of_y,size_of_x);
mask_test_buffers_struct_total = zeros(size_of_y,size_of_x);
mask_test_buffer_points_struct = zeros(size_of_y,size_of_x);
mask_test_buffer_points_struct_total= zeros(size_of_y,size_of_x);
%
N = main_tour;    % Numbers from 1 to N will be permuted
n = main_tour;    % Numbers to be extracted
x = randperm(N);    % Permute numbers between 1 and N
x = x(1:n);    % Retain first n
%
train_number = round(main_tour * (Train_points_percentage/100));
test_number = main_tour-train_number;
for i = 1: train_number
Training_buffers_struct(i).buffers = interested_buf_struct(x(i)).final;
Training_buffers = Training_buffers + interested_buf_struct(x(i)).final;
end
figure
imshow(Training_buffers)
title('Training Buffers')
for j = train_number + 1: main_tour
Test_buffers_struct(j-train_number).buffers = interested_buf_struct(x(j)).final;
Test_buffers = Test_buffers + interested_buf_struct(x(j)).final;
end
figure
imshow(Test_buffers)
title('Test Buffers')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating random test and train buffer pixels from random buffers. This is second random 
% sampling of the 2LRS algorithm.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1: train_number
train_buffers_struct(k).points = find (Training_buffers_struct(k).buffers);
for i= 1:length(train_buffers_struct(k).points)
    r = rem(train_buffers_struct(k).points(i),size_of_y);
    y_coordinate = r; %mod(train_buffer_points_struct(k).points(i),size_of_y);
    x_coordinate = (train_buffers_struct(k).points(i)-r)/size_of_y+1;
    train_buffers_struct(k).lat(i) = y_coordinate;
    train_buffers_struct(k).lon(i) = x_coordinate;
    mask_train_buffers_struct(y_coordinate,x_coordinate) = 1;
end
mask_train_buffers_struct_total = mask_train_buffers_struct_total + mask_train_buffers_struct;
end
% figure
% imshow(mask_train_buffers_struct_total)
% figure
% imshow(Training_buffers - mask_train_buffers_struct_total)
mask_train_buffers_struct_total2 = mask_train_buffers_struct_total(1:1906,1:size_of_x);
% figure
% imshow(mask_train_buffers_struct_total2)
total_train_points = 0;
%
for k = 1:train_number
N = length(train_buffers_struct(k).points);    % Numbers from 1 to N will be permuted
n = length(train_buffers_struct(k).points);    % Numbers to be extracted
x2 = randperm(N);    % Permute numbers between 1 and N
x2 = x2(1:n);    % Retain first n
train_points = round( length(train_buffers_struct(k).points) * (point_percentage/100));

for i=1:train_points
    r = rem(train_buffers_struct(k).points(x2(i)),size_of_y);
    y_coordinate = r %mod(train_buffer_points_struct(k).points(i),size_of_y);
    x_coordinate = (train_buffers_struct(k).points(x2(i))-r)/size_of_y+1;
    train_buffers_struct(k).lat(i) = y_coordinate;
    train_buffers_struct(k).lon(i) = x_coordinate;
    mask_train_buffer_points_struct(y_coordinate,x_coordinate) = 1;
end
mask_train_buffer_points_struct_total = mask_train_buffer_points_struct_total + mask_train_buffer_points_struct;
total_train_points = total_train_points + train_points;
end
mask_train_buffer_points_struct_total = mask_train_buffer_points_struct_total(1:1906,1:size_of_x);
figure
imshow(mask_train_buffer_points_struct_total)
title('mask train buffer points struct total')
% imtool(mask_train_buffer_points_struct_total)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
for k = 1:test_number
test_buffers_struct(k).points = find (Test_buffers_struct(k).buffers);
for i= 1:length(test_buffers_struct(k).points)
    r = rem(test_buffers_struct(k).points(i),size_of_y);
    y_coordinate = r; %mod(train_buffer_points_struct(k).points(i),size_of_y);
    x_coordinate = (test_buffers_struct(k).points(i)-r)/size_of_y+1;
    test_buffers_struct(k).lat(i) = y_coordinate;
    test_buffers_struct(k).lon(i) = x_coordinate;
    mask_test_buffers_struct(y_coordinate,x_coordinate) = 1;
end
mask_test_buffers_struct_total = mask_test_buffers_struct_total + mask_test_buffers_struct;
end
%  figure
%  imshow(mask_test_buffers_struct_total)
%  figure
%  imshow(Test_buffers - mask_test_buffers_struct_total)
total_test_points = 0;
for k = 1:test_number
N = length(test_buffers_struct(k).points);    % Numbers from 1 to N will be permuted
n = length(test_buffers_struct(k).points);    % Numbers to be extracted
x2 = randperm(N);    % Permute numbers between 1 and N
x2 = x2(1:n);    % Retain first n
test_points = round( length(test_buffers_struct(k).points) * (point_percentage/100));
for i=1:test_points
    r = rem(test_buffers_struct(k).points(x2(i)),size_of_y);
    y_coordinate = r; %mod(train_buffer_points_struct(k).points(i),size_of_y);
    x_coordinate = (test_buffers_struct(k).points(x2(i))-r)/size_of_y+1;
    test_buffers_struct(k).lat(i) = y_coordinate;
    test_buffers_struct(k).lon(i) = x_coordinate;
    mask_test_buffer_points_struct(y_coordinate,x_coordinate) = 1;
end
mask_test_buffer_points_struct_total = mask_test_buffer_points_struct_total + mask_test_buffer_points_struct;
total_test_points = total_test_points + test_points;
end
mask_test_buffer_points_struct_total = mask_test_buffer_points_struct_total(1:1906,1:size_of_x);
figure
imshow(mask_test_buffer_points_struct_total)
title('Mask Test Buffer Points Struct Total')
total_test_points = length(find(mask_test_buffer_points_struct_total));
%  figure
%  imshow(Test_buffers - mask_test_buffer_points_struct_total)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NON  LS points and its train masks...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
working_area_filled = imfill(working_area,'holes');
non_ls_area = zeros(size_of_y,size_of_x);
mask_non_ls_area_train_points_struct = zeros(size_of_y,size_of_x);
mask_non_ls_area_train_points_struct_total = zeros(size_of_y,size_of_x);
mask_non_ls_area_test_points_struct = zeros(size_of_y,size_of_x);
mask_non_ls_area_test_points_struct_total = zeros(size_of_y,size_of_x);
%
for i=1:99
    non_ls_area = non_ls_area + only_buf_filled_struct(i).filled; 
end
non_ls_area = (~non_ls_area & working_area_filled);
figure
imshow(non_ls_area)
title('Non LS Area in White')
%
non_ls_area_points = find (non_ls_area);
N = length(non_ls_area_points);    % Numbers from 1 to N will be permuted
n = length(non_ls_area_points);    % Numbers to be extracted
x3 = randperm(N);    % Permute numbers between 1 and N
x3 = x3(1:n);    % Retain first n
%
for i=1:total_train_points
    r = rem(non_ls_area_points(x3(i)),size_of_y);
    y_coordinate = r; %mod(train_buffer_points_struct(k).points(i),size_of_y);
    x_coordinate = (non_ls_area_points(x3(i))-r)/size_of_y+1;
    non_ls_area_train_points_struct(1).lat(i) = y_coordinate;
    non_ls_area_train_points_struct(1).lon(i) = x_coordinate;
    mask_non_ls_area_train_points_struct(y_coordinate,x_coordinate) = 1;
end
mask_non_ls_area_train_points_struct_total = mask_non_ls_area_train_points_struct_total + mask_non_ls_area_train_points_struct;
mask_non_ls_area_train_points_struct_total = mask_non_ls_area_train_points_struct_total(1:1906,1:size_of_x);
figure 
imshow(mask_non_ls_area_train_points_struct_total)
title('Mask Non Ls Area Train Points Struct Total')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NON  LS points and its test masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(non_ls_area_points);    % Numbers from 1 to N will be permuted
n = length(non_ls_area_points);    % Numbers to be extracted
x4 = randperm(N);    % Permute numbers between 1 and N
x4 = x4(1:n);    % Retain first n
%
for i=1:total_test_points
    r = rem(non_ls_area_points(x4(i)),size_of_y);
    y_coordinate = r; %mod(train_buffer_points_struct(k).points(i),size_of_y);
    x_coordinate = (non_ls_area_points(x4(i))-r)/size_of_y+1;
    non_ls_area_test_points_struct(1).lat(i) = y_coordinate;
    non_ls_area_test_points_struct(1).lon(i) = x_coordinate;
    mask_non_ls_area_test_points_struct(y_coordinate,x_coordinate) = 1;
end
mask_non_ls_area_test_points_struct_total = mask_non_ls_area_test_points_struct_total + mask_non_ls_area_test_points_struct;
mask_non_ls_area_test_points_struct_total = mask_non_ls_area_test_points_struct_total(1:1906,1:size_of_x);
% figure 
% imshow(mask_non_ls_area_test_points_struct_total)
% title('Mask Non Ls Area Test Points Struct Total')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading multispectral image
ms1 = imread('stretch_bands_A_band1-14_.tif');
ms2 = imread('stretch_bands_B_band15-28_.tif');
% seperating the MS image to its all bands 
DEM                 = ms1(:,:,1);
slope               = ms1(:,:,2);
aspect              = ms1(:,:,3);
plan_curvature      = ms1(:,:,4);
profil_curvature    = ms1(:,:,5);
convergence_index   = ms1(:,:,6);
TWI                 = ms1(:,:,7);
LS_factor           = ms1(:,:,8);
NDVI                = ms1(:,:,9);
Kaolinite           = ms1(:,:,10);
Calcite             = ms1(:,:,11);
OHI                 = ms1(:,:,12);
distance_to_fault   = ms1(:,:,13);
distance_to_channel = ms1(:,:,14);
DECORR_1            = ms2(:,:,1);
DECORR_2            = ms2(:,:,2);
DECORR_3            = ms2(:,:,3);
DECORR_4            = ms2(:,:,4);
DECORR_5            = ms2(:,:,5);
DECORR_6            = ms2(:,:,6);
DECORR_7            = ms2(:,:,7);
DECORR_8            = ms2(:,:,8);
DECORR_9            = ms2(:,:,9);
DECORR_10           = ms2(:,:,10);
DECORR_11           = ms2(:,:,11);
DECORR_12           = ms2(:,:,12);
DECORR_13           = ms2(:,:,13);
DECORR_14           = ms2(:,:,14);
% 
ms = zeros(1906,778,28);
ms(:,:,1:14)=ms1;
ms(:,:,15:28)=ms2;

% As the mask will be used for logical AND operation all values of the mask
% must be zero or one. 
for i = 1:1906
    for j =1:size_of_x
    if mask_train_buffer_points_struct_total(i,j) >= 1 
        mask_train_buffer_points_struct_total(i,j) = 1;
    end
    end
end
% As the mask will be used for logical AND operation all values of the mask
% must be zero or one. 
for i = 1:1906
    for j =1:size_of_x
    if mask_test_buffer_points_struct_total(i,j) >= 1 
        mask_test_buffer_points_struct_total(i,j) = 1;
    end
    end
end
% If the user want to use whole image as a test data set, user can use the notation shown in below.    
% test_image_full = [DEM(:),slope(:),aspect(:), plan_curvature(:), profil_curvature(:), convergence_index(:),TWI(:),LS_factor(:),NDVI(:),Kaolinite(:),Calcite(:),...
% OHI(:),distance_to_fault(:),distance_to_channel(:),DECORR_1(:),DECORR_2(:),DECORR_3(:),DECORR_4(:),DECORR_5(:),DECORR_6(:),DECORR_7(:),DECORR_8(:),...
%     DECORR_9(:),DECORR_10(:),DECORR_11(:),DECORR_12(:),DECORR_13(:),DECORR_14(:)];
%    
% test_image_full = [DEM(:),slope(:),aspect(:), plan_curvature(:)];
%
% To generate train points from all of the bands,
% mask_train_buffer_points_struct_total is multiplied with all bands.
how_many_bands = 28;
for j = 1:how_many_bands
train_points_layer_struct(j).layers = mask_train_buffer_points_struct_total .*ms(:,:,j)
end
% figure
% imshow(train_points_layer_struct(7).layers)
%
% To generate non_ls train points from all of the bands,
% mask_non_ls_area_train_points_struct_total is multiplied with all bands.
for j = 1:how_many_bands
train_points_layer_non_ls_struct(j).layers = mask_non_ls_area_train_points_struct_total .*ms(:,:,j)
end
% figure
% imshow(train_points_layer_struct(1).layers)
total_train_points = length(find(train_points_layer_non_ls_struct(14).layers));
arrayImage_train = zeros(total_train_points*2-200, how_many_bands); 
% Generating a train data set
for bands =1:how_many_bands
a = 1;
arrayImage_first = zeros (total_train_points-100,1);
for i =1: 1906
    for j =1:size_of_x
        if (train_points_layer_struct(bands).layers(i,j) > 0 & a<(total_train_points-99))
        arrayImage_first(a) = train_points_layer_struct(bands).layers(i,j);
        a = a+1;
        end
    end
end
arrayImage_train(1:(total_train_points-100),bands) = arrayImage_first;
end
% Generating a train data set
for bands =1:how_many_bands
a = 1;
arrayImage_first = zeros (total_train_points-100,1);
for i =1: 1906
    for j =1:size_of_x
        if (train_points_layer_non_ls_struct(bands).layers(i,j) > 0 & a<(total_train_points-99))
        arrayImage_first(a) = train_points_layer_non_ls_struct(bands).layers(i,j);
        a = a+1;
        end
    end
end
arrayImage_train(total_train_points-100+1:total_train_points*2-200,bands) = arrayImage_first;
end
% Generating a train data set
class = zeros((total_train_points*2-200),1);
for i = 1:total_train_points*2-200
    if i<total_train_points-99
        class(i)=1;
    else
        class(i)=0;
    end
end
% Generating a train data set
array_learner = [arrayImage_train, class];% array_learner is the final train data set that will be used in matlab classification learner toolbox.
% To generate test_points_layer_struct from all of the bands,
% mask_test_buffer_points_struct_total is multiplied with all bands
for j = 1:how_many_bands
test_points_layer_struct(j).layers = mask_test_buffer_points_struct_total .*ms(:,:,j)
end
% To generate test_points_layer_non_ls_struct from all of the bands,
% mask_non_ls_area_test_points_struct_total is multiplied with all bands
for j = 1:how_many_bands
test_points_layer_non_ls_struct(j).layers = mask_non_ls_area_test_points_struct_total .*ms(:,:,j)
end
% Generating a test data set from ls points. 
arrayImage_test_ls = zeros(total_test_points-100, how_many_bands);
for bands =1:how_many_bands
a = 1;
arrayImage_first = zeros (total_test_points-100,1);
for i =1: 1906
    for j =1:size_of_x
        if (test_points_layer_struct(bands).layers(i,j) > 0 & a<total_test_points-100+1)
        arrayImage_first(a) = test_points_layer_struct(bands).layers(i,j);
        a = a+1;
        end
    end
end
arrayImage_test_ls(1:total_test_points-100,bands) = arrayImage_first;
end
%Generating a test data set from non ls points. 
arrayImage_test_non_ls = zeros(total_test_points-100, how_many_bands);
for bands =1:how_many_bands
a = 1;
arrayImage_first = zeros (total_test_points-100,1);
for i =1: 1906
    for j =1:size_of_x
        if (test_points_layer_non_ls_struct(bands).layers(i,j) > 0 & a<total_test_points-100+1)
        arrayImage_first(a) = test_points_layer_non_ls_struct(bands).layers(i,j);
        a = a+1;
        end
    end
end
arrayImage_test_non_ls(1:total_test_points-100,bands) = arrayImage_first;
end
