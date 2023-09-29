%% Version III -- XZ and YZ through local max , XY through masking

clc ; clear ; %close all

tic

tmin = 1 ;
tmax = 30 ; 

cd '/Users/Sharr/Dropbox/PostDoc/Projects/Morphogenesis/Codes and ImageProcessing/Stress_Flow_Analysis/DynamicalAnalysis/Codes'
orgdir = pwd;

cd '/Users/Sharr/Dropbox/PostDoc/Projects/Morphogenesis/Codes and ImageProcessing/Stress_Flow_Analysis/DynamicalAnalysis/Set3Data/SegtFilter_z15_31'
savedir = pwd;

% cd '/Users/Sharr/Dropbox/PostDoc/Projects/Morphogenesis/Codes and ImageProcessing/Stress_Flow_Analysis/Set3_Myosin/Myo_Stacks_Symmetrized'
cd '/Users/Sharr/Dropbox/PostDoc/Projects/Morphogenesis/Codes and ImageProcessing/Stress_Flow_Analysis/Set5_Membrane/Membrane30/Mem30'
imgdir = pwd ;

fsample = sprintf('Set5_mem_%02d.tif' , tmin);
Inf = imfinfo(fsample);
Isamp = imread(fsample , 1);

[Ny , Nx] = size(Isamp);
Nz = size(Inf , 1);
xmin = 1;%601;
xmax = Nx;%1650;
N_x = xmax - xmin + 1;

ymin = 1;%151;
ymax = Ny;%Ny - 150;
N_y = ymax - ymin + 1;

zmin = 2;%15;%4;%15;7       % this is the minimum z-level above which the 
zmax = 9;%31;%20;
N_z = zmax - zmin + 1;

bg_rad = 30;

hisTiles = [64 , 64];
nBin = 128;
clipLim = 0.1;
gaussFiltRad = 4;

z_WL = 10 ;                      % Distance between peaks in z-dir;
y_WL = 20 ;                      % Distance between peaks in y-dir;
x_WL = 15 ;                      % Distance between peaks in x-dir;

cy = 30 ;                      % Distance between peaks in y-dir;
cx = 30 ;                      % Distance between peaks in x-dir;

sy = 15;
sx = 15;

for tt = 20 : 20
    
    disp('time = ')
    disp(tt);
    
    tic
    
    cd(imgdir);

%     stName = sprintf('S2_Myo_%02d.tif' , tt);
    stName = sprintf('Set5_mem_%02d.tif' , tt);
    IB = zeros(N_z , N_y , N_x);      % Image Block
    Iraw = zeros(N_y , N_x , N_z);      % Image Block
    
    for zz = zmin:zmax
        
        Iz = (imread(stName , zz));
        Iraw(:,:,zz) = Iz;
        
    end
    
    for zz = 1
        Izz = adapthisteq(Iz);
        Izz0 = imreducehaze(Izz);
        Izz1 = imflatfield(Iz , 5);
        figure ; imshow(Izz,[]);
        figure ; imshow(Izz0,[]);
        figure ; imshow(Izz1,[]);

        I2 = adapthisteq(Iz , 'NumTiles' , hisTiles , 'NBins' , nBin, 'ClipLimit' , clipLim);
        I3g = imgaussfilt(I2 , 2);
        figure ; imshow(I3g,[]);

        bg_t = imopen(I3g , strel('disk' ,  bg_rad) );
        IzzBG = imsubtract(I3g , bg_t);
        figure ; imshow(adapthisteq(IzzBG),[]);
        IB(zz - zmin + 1 , 1:N_y , 1:N_x) = IzzBG(ymin:ymax , xmin:xmax);
        
        II = IzzBG;
        I0 = imgaussfilt(II , 4);
        Iws = ~logical(watershed(I0));
        IwsOver = labeloverlay(I3g , Iws);
        figure ; imshow(IwsOver)
        
    end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    Ixy = zeros(N_z , N_y , N_x);       % Reconstructed Image along x;
    Iyz = zeros(N_z , N_y , N_x);       % Reconstructed Image along y;
    Izx = zeros(N_z , N_y , N_x);       % Reconstructed Image along y;
    IbwZ = zeros(N_z , N_y , N_x);       % Reconstructed Image along y;

    for yy = 1:N_y
        Irx = squeeze(IB(: , yy , :));
        Iy = imgaussfilt(Irx , 1);
        for zz = 1:N_z
            Irx = Iy(zz , :);
            Irx = smooth(Irx , 'moving');
            [~ , xMax] = findpeaks( Irx , 'MinPeakDistance' , x_WL );
            Iyz(zz , yy , xMax) = 1;
        end
    end
    
    disp('1');
    Ixx = Iyz;
    
    for xx = 1:N_x
        Iry = squeeze(IB(: , : , xx));
        Ix = imgaussfilt(Iry , 1);
        for zz = 1:N_z    
            Iry = Ix(zz , :);
            Iry = smooth(Iry , 'moving');
            [~ , yMax] = findpeaks( Iry , 'MinPeakDistance' , y_WL );
            Izx(zz , yMax , xx) = 1;
        end
    end
    
    disp('2');    
    Iyy = Izx;
    
    for zz = 1:N_z
        J = imflatfield(squeeze(IB(zz , : , :)) , 10);
        se2 = strel('disk' , 5);
        J_bhf = imbothat(J,se2);
        Jdiff0 = imsubtract(imadd(J , imtophat(J,se2)) , J_bhf);
        Jdiff = imgaussfilt(Jdiff0 , 4);
        Jws = ~logical(watershed((Jdiff)));
        Irz = squeeze(IB(zz , : , :));
        II = Irz;
        I0 = imgaussfilt(II , gaussFiltRad);
        Iws = watershed(I0);
        IbwZ(zz , : , :) = Jws;%~logical(Iws);
    end

    disp('3');

    cd(savedir); 

    saveNameTot = sprintf('SmoothWS_TotzRec_Myo_%02d.tif' , tt) ; 

    I_xx = permute(Ixx , [2,3,1]);
    I_yy = permute(Iyy , [2,3,1]);
    I_zz = permute(IbwZ , [2,3,1]);
    I_tot = I_xx + I_yy + I_zz ; 

    
    for zz = 1:N_z
        II = squeeze(I_tot(:,:,zz));
        II = bwmorph(II , 'diag');
        II = bwmorph(II , 'bridge');
        II = bwmorph(II , 'clean');
        II = bwmorph(II , 'spur' , 10);
        If = bwmorph(II , 'fill' , 10);
        If = bwareaopen(If , 50);
        If0 = If;
        cc = bwconncomp(If);
        for reg = 1:cc.NumObjects
            if length(cc.PixelIdxList{reg}) < 50
                If0(cc.PixelIdxList{reg}) = 0;
            end
        end

        IfFat = bwmorph(If0 , 'fatten');
        IfFat = bwmorph(IfFat , 'majority');

        In1 = ~IfFat;
        In0 = bwmorph(In1 , 'break');
        ccn = bwconncomp(In0);
        for reg = 1:ccn.NumObjects
            if length(ccn.PixelIdxList{reg}) < 20
                In0(ccn.PixelIdxList{reg}) = 0;
            end
        end
        Iff = ~In0;
        Iff = bwmorph(Iff , 'fill');
        Iff = bwmorph(Iff , 'thin' , 10);
        Iff = bwmorph(Iff , 'spur' , 5);
        I_tot(:,:,zz) = Iff;
    end

    for zz = 1 : N_z
        imwrite(I_tot(:,:,zz) , saveNameTot , 'WriteMode', 'append', 'Compression', 'none');                
    end
    
   toc
   
end

cd(orgdir); 

%%

xrange = 701:1000 ; yrange = 1601:1850 ; zrange = 2:7;%5:17 ;
nx = length(xrange); ny = length(yrange); nz = length(zrange);

IT = I_tot(yrange , xrange , zrange);
IT2 = zeros(ny , nx , nz);
volz = zeros(ny , nx , 1);    
cover = ones(ny , nx);
for zz = 1:nz
    II = squeeze(IT(:,:,zz));
    II = bwareaopen(II , 50);
    II = bwmorph(II , 'fatten' , 2);
    II = ~bwareaopen(~II , 50);
    II = bwmorph(II , 'diag');
    II = bwmorph(II , 'bridge');
    II = ~bwareaopen(~II , 50);
    II = bwmorph(II , 'diag');
    II = bwmorph(II , 'shrink', 3);
    II = bwmorph(II , 'diag');

    volz(:,:,zz) = II;
    vol_z = volz(: , : , 1:zz);
    vol_z = bwmorph3(vol_z , 'majority');
    II = bwmorph(logical(II) , 'thin' , 5);
    II = imfill(II , 'holes');
    figure ; imshow(II)
    IT2(:,:,1:zz) = vol_z;
    
end

Vol = (IT2);%(yrange , xrange , zrange);

figure ; volview = volshow(imcomplement(Vol));
volview.ScaleFactors = [1,1,1];

Vol(:,:,zz + 1) = cover;

Vol = bwmorph3(Vol , 'majority');
Vol = bwmorph3(Vol , 'majority');
Vol = bwmorph3(Vol , 'majority');
Vol = bwmorph3(Vol , 'majority');

%%

II1 = squeeze(IT2(:,:,end));
figure ; imshow(II1);
II1 = bwmorph(II1 , 'shrink' , 10);
figure ; imshow(II1);

%%

xrange = 651:900 ; yrange = 1701:1800 ; zrange = 3:10 ;
for zz = 1:8

    volz = squeeze(I_tot(yrange,xrange,zz));
    figure ; imshow(volz , []);
    imtool(volz)
end

%%

figure ; volview = volshow(imcomplement(Vol));
volview.ScaleFactors = [1,1,5];

figure ; volview = volshow(imcomplement(IT2));
volview.ScaleFactors = [1,1,5];


%%
volumeViewer(Vol)
volumeViewer(imcomplement(Vol))

volview = volshow(Vol);
volview.ScaleFactors = [1,1,5];

volview = volshow(imcomplement(Vol));
volview.ScaleFactors = [1,1,5];

IB = permute(IB , [2,3,1]);

Imb = Iraw(yrange , xrange , 2:9);

figure ; volview = volshow(imcomplement(Ifilt2));
volview.ScaleFactors = [1,1,5];

figure ; volview = volshow(rescale(imgaussfilt(Imb , 4)));
volview.ScaleFactors = [1,1,5];

%%
Imb = IB(yrange , xrange , :);
Ifilt = imboxfilt3(Imb);

nx = length(xrange);
ny = length(yrange);
nz = 10;

Ifilt2 = zeros(ny , nx , nz);
IT = zeros(ny , nx , nz);
for zz = 1:nz
    II = squeeze(Imb(:,:,zz));
    Ifilt2(:,:,zz) = II;
    II = bwareaopen(II , 50);
    II = bwmorph(II , 'fatten' , 2);
    II = ~bwareaopen(~II , 50);
    IT(:,:,zz) = II;
end

volview = volshow(Vol);
volview.ScaleFactors = [1,1,5];

figure ; volview = volshow(imcomplement(Ifilt2));
volview.ScaleFactors = [1,1,5];

figure ; volview = volshow(imcomplement(Imb));
volview.ScaleFactors = [1,1,5];

% %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Defining the ROI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ilast = imread(SkelName , Time);

if length(size(Ilast)) > 2
    segbw = rgb2gray(Ilast) ;
else
    segbw = Ilast ;
end

seg0 = bwmorph(segbw , 'shrink') ;
reg0 = ~seg0;
cc0 = bwconncomp(reg0 , 4);

for id = 1:cc0.NumObjects
    if length(cc0.PixelIdxList{id}) < 100
        seg0(cc0.PixelIdxList{id}) = 1;
    end
end

% This is only added to have the vertices be extracted from a seg0
% version without any white blob
seg0 = bwmorph(seg0 , 'shrink', Inf) ;  
seg1 = bwmorph(seg0 , 'bridge') ;
seg2 = bwmorph(seg1 , 'fill') ;
seg3 = bwmorph(seg2 , 'diag') ;
seg4 = bwmorph(seg3 , 'thick') ;
seg5 = bwmorph(seg4 , 'majority') ;
seg6 = bwmorph(seg5 , 'thin'  , Inf) ;    
seg7 = bwmorph(seg6 , 'shrink' , Inf) ;
seg8 = bwmorph(seg7 , 'clean' , Inf) ;
segt_last = seg8;

% Coordinates of the final boundaries of ROI. H: Head, and, G: Gut.
Xh_f = [60 , 470 , 576 , 470 , 60];
Yh_f = [130 , 130 , 1025 , 1880 , 1880];

Xg_f = [610 , 1510 , 1510 , 1120 , 1120 , 1380 , 1140 , 580 , 740];
Yg_f = [130 , 130 , 200 , 570 , 1330 , 1660 , 1950 , 1950 , 1025];

Xg_r = [725 , 1270 , 1045 , 770 , 645];
Yg_r = [144 , 144 , 485 , 485 , 340];

Xg_l = [725 , 790 , 1080 , 1325 , 1165];
Yg_l = [1815 , 1320 , 1320 , 1600 , 1815];

% Coordinates of the boundaries of ROI in time. H: Head, and, G: Gut.
Xh = zeros(length(Xh_f) , Time);
Yh = zeros(length(Yh_f) , Time);

Xg = zeros(length(Xg_f) , Time);
Yg = zeros(length(Yg_f) , Time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xh(: , end) = Xh_f ; Yh(: , end) = Yh_f ; %
Xg(: , end) = Xg_f ; Yg(: , end) = Yg_f ; %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maskT = zeros(Ny , Nx , Time);
maskedSeg = zeros(Ny , Nx , Time);

maskT(: , : , Time) = roipoly(maskT(: , : , Time) , Xh_f , Yh_f);
maskT(: , : , Time) = maskT(: , : , Time) + roipoly(maskT(: , : , Time) , Xg_r , Yg_r);
maskT(: , : , Time) = maskT(: , : , Time) + roipoly(maskT(: , : , Time) , Xg_l , Yg_l);
maskT(801:1300 , : , Time) = 0;

maskLast = squeeze(maskT(: , : , Time));

 
ver0 = bwmorph(segt , 'branchpoints');
ver0 = bwmorph(ver0 , 'thicken' , 1);
segt1 = segt - (segt & ver0);

segMasked0 = bwmorph(segMasked0 , 'shrink', Inf);
segMasked = bwmorph(segMasked0 , 'clean' );
% segMasked = maskLast .* double(segt_last);

fsave = sprintf('Sources/NewROI/Laterals/TotalMask_t%02d.tif' , tmax);
fsegmask = sprintf('Sources/NewROI/Laterals/SegtMask_t%02d.tif' , tmax);
imwrite(double(maskLast) , fsave);
imwrite(double(segMasked) , fsegmask);

for t = 2:Time
   
    im1  = double(imread(MemName , Time - t + 2)) ;
    im2  = double(imread(MemName , Time - t + 1)) ;

    Iseg = imread(SkelName , Time - t + 1);

    if length(size(Iseg)) > 2
        segbw = rgb2gray(Iseg) ;
    else
        segbw = Iseg ;
    end

    seg0 = bwmorph(segbw , 'shrink') ;
    reg0 = ~seg0;
    cc0 = bwconncomp(reg0 , 4);
    
    for id = 1:cc0.NumObjects
        if length(cc0.PixelIdxList{id}) < 100
            seg0(cc0.PixelIdxList{id}) = 1;
        end
    end
    
    % This is only added to have the vertices be extracted from a seg0
    % version without any white blob
    seg0 = bwmorph(seg0 , 'shrink', Inf) ;      
    seg1 = bwmorph(seg0 , 'bridge') ;
    seg2 = bwmorph(seg1 , 'fill') ;
    seg3 = bwmorph(seg2 , 'diag') ;
    seg4 = bwmorph(seg3 , 'thick') ;
    seg5 = bwmorph(seg4 , 'majority') ;
    seg6 = bwmorph(seg5 , 'thin'  , Inf) ;    
    seg7 = bwmorph(seg6 , 'shrink' , Inf) ;
    seg8 = bwmorph(seg7 , 'clean' , Inf) ;
    Segm = seg8;

    mask_old = maskT(:,:, Time - t + 2);
    
    tform = imregdemons(im1 , im2 , [300 , 100 , 50]);
    mask_new = imwarp(mask_old , tform , 'cubic');
    
    mask_new = imfill(mask_new);
    maskT(:,:, Time - t + 1) = mask_new;

    mask_new = squeeze(maskT(: , : , Time - t + 1));

    segMasked0 = bwmorph(mask_new .* double(Segm) , 'shrink', Inf);
    segMasked = mask_new .* double(Segm);
    
    fsave = sprintf('Sources/NewROI/Laterals/TotalMask_t%02d.tif' , tmax - t + 1);
    imwrite(double(mask_new) , fsave);
     
    fsegmask = sprintf('Sources/NewROI/Laterals/SegtMask_t%02d.tif' , tmax - t + 1);
    imwrite(double(segMasked) , fsegmask);
    
end



%%

clear
clc

fsample = sprintf('S2_Myo_%02d.tif' , 30);
Inf = imfinfo(fsample);
Isamp = imread(fsample , 1);
I = imread(fsample , 6);
J = adapthisteq(I);
figure ; imshow(J)

Iw0 = imgaussfilt(I , 4);
Jw0 = imgaussfilt(J , 4);

figure ; imshow(rescale(J))

Iw = ~logical(watershed(Iw0));
Jw = ~logical(watershed(Jw0));

figure ; imshow(Iw);
figure ; imshow(Jw);

Jf = bwmorph(Jw , 'fatten');
If = bwmorph(Iw , 'fatten');

figure ; imshowpair(Jf , If);




