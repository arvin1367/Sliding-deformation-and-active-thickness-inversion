%%%%%% Extract deformation time series in the sliding direction
%%%%%% The code is written by Huiyuan Luo.
%% Add this code package to the path first. Contains a parent folder 
%  and three sub folders, as well as a time table of image pairs for 
%  ascending and descending.
% example:
%   ascending_los      (dir) 
%   descending_los     (dir)
%   supporting_documents   (dir)
%   date_ascending.xls     (table)
%   date_descending.xls    (table)
% Notice: Please ensure that the folder and file names are consistent, 
% except for the sbas-pot processing results, which are the files in the 
% first four folders.
% supporting_documents include the smoothed dem_aspect.tif (local aspect angle) 
% and dem_slope.tif (local solpe angle) obtained.
clear
clc
%% Begin
% Folder Name
dirc='xiongba';
obj_arvin=myPreprocess_sliding(dirc,36.6714,-12.7536,42.2848,192.7446); 
% Input 2: the incident angle of ascending
% Input 3: the azimuth angle of ascending
% Input 4: the incident angle of descending
% Input 5: the azimuth angle of descending
save obj_arvin
%Calculate and obtain results
obj_data=myProcess_sliding(obj_arvin,0.1);
% Input 2: regularization parameter value lamda
%%  

%%%%%% To invert the thickness of a landslide using the three-dimensional 
%%%%%% flow velocity of the landslide.
% Directly call the function.
thickness_cal('ve.tif','vn.tif','vu.tif',30,0.5,0.1,150);
% For specific introduction, please refer to the function
%Input 1: East-west deformation rate .tif
%Input 2: North-south deformation rate .tif 
%Input 3: Vertical deformation rate .tif
%Input 4: Grid resolution 
%Input 5: rheological parameter
%Input 6: Regularization coefficient     recommend 0.1
%Input 7: Maximum constraint.






