 function shrink_movie_b(day, BIN_SIZE)  
%Bins the data so you can run ClusterVermis and Cluster_roi_lever to get
%TCs
 
 if(~exist('BIN_SIZE', 'var'))
     BIN_SIZE =20;
 end
 DATA_DIR =  'S:\Data\';

ROI_x =[];
ROI_y = [];

mati_code_cd = 'C:\Users\minhal\Documents\Repositories\Imaging\Jake\LeverAnalysis';
image_dest  = [DATA_DIR day '\'];
old_cd = cd(image_dest);
% get all files

%------ assume the file has a a unique string [_MMStack_[???].ome] and ?? is the number of the file -1
img_str_indicator = 'MMStack';
all_files = dir(['*' img_str_indicator '*'] );
% ---- get order of imaging
file_order = [];
for i=1:length(all_files)
    beg_inx = findstr(all_files(i).name, img_str_indicator) + length(img_str_indicator);
    end_no_inx = strfind(all_files(i).name , '.ome');
    if(end_no_inx == beg_inx)
        file_order(i) = 1;
    else
        file_order(i) = str2num(all_files(i).name(beg_inx+1:end_no_inx-1))+1;
    end
end
all_files(file_order) = all_files;

for i=1:length(all_files)
    info{i} = imfinfo(all_files(i).name);
    
    if(i==1)
        all_info= info{i};
    else
        all_info = cat(1, all_info,info{i});
    end
end
frame_times = get_frame_time_by_movie_info(all_info);

dest =  [image_dest  day '_ROI'];
shrink_img = [];
sz = {};
for i=1:length(all_files)
    cd(mati_code_cd);
    [img, sz] = get_movie_by_ROI(info{i}(1).Filename, info{i}, ROI_x, ROI_y, BIN_SIZE, 1, length(info{i}));
    shrink_img = [shrink_img, img]; 
end

% sterch values --- keep values betwee 0-255 
% min_v = min(shrink_img(:));
% max_v = max(shrink_img(:));
% shrink_img = 255*(shrink_img - min_v)/(max_v-min_v); 

    % ----- write tif file with ROI only
for i=1:size(shrink_img,2)    
    c_img = shrink_img(:,i);
    s_img = uint8(reshape(c_img,sz));
    mode= 'append';
    if(i==1)
        mode = 'overwrite';
    end
    imwrite(s_img, [dest 'shrink.tif'], 'WriteMode', mode);    
end



cd(old_cd);