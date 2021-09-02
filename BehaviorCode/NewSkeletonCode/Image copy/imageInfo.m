clear;
imagesInfo = struct;
for i = 1 : 23
    filename = strcat(int2str(i), '.jpg');
    I = imread(filename);
    imagesInfo(i).max = max(I(:));
    imagesInfo(i).min = min(I(:));
    imagesInfo(i).all = I;
end