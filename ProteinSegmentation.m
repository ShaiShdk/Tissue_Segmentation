%% Time Evolution of segmented pictures and correcting the bad regions/edges.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The velocity field is first found using the PIV of the membrane data.
%   Then the velocity field is used to evolve the Segmented image at time
%   to next timepoint t+1.
%   
%   At each iteration two rounds of warp between evolved version and actual
%   segmented at t+1 is performed for further accuracy. 
%
%   The second section consists of a for loop that changes the initial time
%   point. For each time point, the initial frame is taken and moved along
%   with the velocity field that is obtained using the phase correlation of
%   the evolved image and the main image at time t+1. 
%
%   A map is defined between each two successive time points that connects  
%   the label of the matrices. First the centroids of each time point is
%   transferred according to the flow field obtained from raw images. The
%   new position of the centroid is now assumed to be close to that of the
%   corresponding region (cell) at the next time point. Using the
%   properties of the second image we compare the area (total number of pixels)
%   between the two putative corresponding cells. We then decide whether
%   the map is correct or not. A level of confidence is needed. The map is
%   a time-global map, in the sense that in the end all the images should
%   agree with each other. However in reallity there might exist
%   topological changes that forbids the exact correspondence. We aim at
%   identifying these too, using the adjacency matrix. 
% 
%                 Shahriar Shadkhoo , Dec. 2018, KITP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc ; clear ; close all

clear CellEd ; clear EdCell ; clear VerCell ; clear CellVer ;

cd 'images_directory'
orgdir = pwd;

MyoName = sprintf('Sources/MyoNormal_Stacks4_8.tif');
MemName = sprintf('Sources/MemBGSubt.tif');
MaskedSkelName = sprintf('Sources/C1-Edited2.tif');
metr = load('Sources/Metric_4_Sam.mat');

im1 = imread(MemName , 1);
Ny = size(im1 , 1);    Nx = size(im1 , 2);

metric0 = metr.g;  % metric from 3D (real body) to 2D (image plane)
metr1 = zeros(Ny , Nx , 2 , 2);
metr1(:,:,1,1) = metric0{1,1};
metr1(:,:,1,2) = metric0{1,2};
metr1(:,:,2,1) = metric0{2,1};
metr1(:,:,2,2) = metric0{2,2};
metr2 = permute(metr1 , [3,4,1,2]);

cd(orgdir);

invMetr = multinv(metr2);
Metric = permute(invMetr , [3,4,1,2]);  % metric from 2D (image plane) to 3D (real space)

cd ../

% 'APpart' defines the part of the embryo under analysis, 'head','tail','full'
APpart = 'tail';
Ceph_x = 570;

% 'LatPart' defines the parts: Left, Right, Lateral = Left + Right, and Dorsal 
LatPart = 'Lateral';
Dors_y = 1025;

RightLat = 1301:1800;  LeftLat = 301:800;

% lower and upper limits of number of pixels in an edge
EdPxl_LL = 5 ; EdPxl_UL = 30 ; 

tminSeries = 5 ; tmaxSeries = 14;

tmin = 5 ;    tmax = 14 ;   
step = 1 ;    % step in timeframes.

global minTP
minTP = 10; % minimum timepoints for analysis

global Time
Time = tmax - tmin + 1;     % length of time frames

IterNum = 1 ;        % number of iterations to match the evolved image with the actual image
EdgeLength  = 40;   % Length of box edges in pixels; 

marX = 20 ; marY = 20 ;

isf         = 1;    % 0.4;   % image scaling factor. 
smooth      = 1;    % set to 1 if gaussian smoothing is desired
sigma       = 1.7;  % standard deviation of gaussian kernel

im1 = imresize(im1,isf,'bicubic');

% define the grid on which to compute the flow field
[X1,Y1] = meshgrid(EdgeLength/2 : EdgeLength : size(im1,2) - EdgeLength/2 , EdgeLength/2 : EdgeLength : size(im1,1) - EdgeLength/2); 

%%  

Tissue = struct;

for t = 1 : Time
    
    clear cellc; clear edc; clear verc; clear EdAng; clear EdAng0;
    clear CellEd; clear EdCell; clear VerCell; clear CellVer; clear EdLen; 
    
    seg = logical(imread(MaskedSkelName , t));

    if length(size(seg)) > 2
        segbw = rgb2gray(seg) ;
    else
        segbw = seg ;
    end
    
    if strcmp(APpart , 'head') == 1
        segbw(: , Ceph_x:end) = false;
    elseif strcmp(APpart , 'tail') == 1
        segbw(: , 1 : Ceph_x) = false;
    end

    if strcmp(LatPart , 'Right') == 1
        segbw(801:end , :) = false;
    elseif strcmp(LatPart , 'Left') == 1
        segbw(1:1300 , :) = false;
    elseif strcmp(LatPart , 'Lateral') == 1
        segbw(801:1300 , :) = false;        
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
        
    cc1 = bwconncomp( ~ Segm , 4);

    px1 = cat(1 , cc1.PixelIdxList);
    cell_area = zeros(cc1.NumObjects,1);

    for ce = 1:cc1.NumObjects
        cell_area(ce) = length(cc1.PixelIdxList{ce});
    end
    
    [MaxArea , locMax] = max(cell_area);    
    
    % removing the background region
    px1(locMax) = [];
    cc1.PixelIdxList = px1;
    cc1.NumObjects = cc1.NumObjects - 1 ;

    cen1 = regionprops(cc1 , 'centroid');
    cent1 = round(cat(1 , cen1.Centroid));
    per = [2 , 1];
    cent1 = cent1(: , per);    
    L1 = labelmatrix(cc1);

%%%%%%%%%%%%%%%%%%% Modifying Vertices and New ConnComp %%%%%%%%%%%%%%%%%%%    

    ver0 = bwmorph(Segm , 'branchpoints');
    cv0 = bwconncomp(ver0);
    VerCent = regionprops(cv0 , 'centroid');
    V_loc = round(cat(1 , VerCent.Centroid));

    per = [2 , 1];
    V_loc = V_loc(: , per);

    ver0 = bwmorph(ver0, 'diag');
    ver = bwmorph(ver0 , 'thicken' , 1);

    verc = bwconncomp(ver);
    vPxl = verc.PixelIdxList;

%%%%%%%%%%%%%% Array of Verts and their corresponding Cells %%%%%%%%%%%%%%%

    % This is the cell array of Cells with their corresponding vertices.
    % Initially creating an array with maximum number of Cells so that an
    % array of that size is created and we can address all the Cells in the
    % below vertex-assignment.

    VerCell = cell(verc.NumObjects , 1); 
    CellVer = cell(cc1.NumObjects , 1);
    
    for nv = 1:verc.NumObjects

        [ymin , xmin] = ind2sub([Ny , Nx] , verc.PixelIdxList{nv}(1));
        [ymax , xmax] = ind2sub([Ny , Nx] , verc.PixelIdxList{nv}(end));

        yrange_v = ymin : ymax ;
        xrange_v = xmin : xmax ;

        Cell_Ver_overlap = zeros(1 , 10) ;
        N_Cells = length( unique( L1(yrange_v , xrange_v) ) ) ;

        Cell_Ver_overlap(1 : N_Cells) = sort(unique( L1(yrange_v , xrange_v) ) , 'descend');
        CellInds = find( Cell_Ver_overlap ~= 0 ) ; % Cell_Ver_overlap ; %
        VerCell{nv} = Cell_Ver_overlap(CellInds)' ; %Cell_Ver_overlap ; %
        
        for ncv = 1:length(VerCell{nv})
            CellVer{VerCell{nv}(ncv)}(length(CellVer{VerCell{nv}(ncv)}) + 1) = nv ;
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    ed = logical(Segm) - (logical(Segm) & logical(ver)) ;
    edc = bwconncomp(ed);
    
    %%%%%% Here we should throw away all the edges that are in the marX
    %%%%%% and marY, namely they have even one pixel in these regions.
    
    edPxl = edc.PixelIdxList;
    endpts = bwmorph(ed , 'endpoints');
    endpt_ind = find(endpts == 1);

%%%%%%%% Array of Edges and their corresponding Cells and Vertices %%%%%%%%

    CellEd = cell(cc1.NumObjects + 1 , 1) ;
    EdCell = cell(length(edPxl) , 1) ;
    
    EdLen = zeros(length(edPxl) , 1) ;

    for ned = 1:length(edPxl)
        if length(edPxl{ned}) > EdPxl_LL && length(intersect(edPxl{ned} , endpt_ind)) == 2 && length(edPxl{ned}) < EdPxl_UL

            intsct = intersect(edPxl{ned} , endpt_ind);
            [y1 , x1] = ind2sub([Ny , Nx] , intsct(1));
            [y2 , x2] = ind2sub([Ny , Nx] , intsct(2));
            
            midpt_ind = edPxl{ned}(ceil(length(edPxl{ned})/2));
            %%%%% midpt_r is midpt_row and midpt_c is midpt_column
            [midpt_r , midpt_c] = ind2sub([Ny , Nx] , midpt_ind);
            
            g11 = Metric(midpt_r , midpt_c , 1,1); g12 = Metric(midpt_r , midpt_c , 1,1); 
            g21 = Metric(midpt_r , midpt_c , 2,1); g22 = Metric(midpt_r , midpt_c , 2,2);
            
            EdLen(ned , 1) = sqrt( g11*(y1 - y2)^2 + (g12+g21)*(y1 - y2)*(x1 - x2) ...
                + g22*(x1 - x2)^2 ); 

            FatMidpt_y = midpt_r-1 : midpt_r+1 ; 
            FatMidpt_x = midpt_c-1 : midpt_c+1 ;

            Cell_Ed_overlap = unique(L1(FatMidpt_y , FatMidpt_x)) ;
            CellInds = find(Cell_Ed_overlap ~= 0);
            EdCell{ned} = Cell_Ed_overlap(CellInds)' ;
                
            for ce = 1:length(EdCell{ned})
                CellEd{EdCell{ned}(ce)}(length(CellEd{EdCell{ned}(ce)}) + 1) = ned;
            end
        else
             EdLen(ned , 1) = 0  ;
             EdCell{ned} = [] ;
        end
    end
    
    VerEd = cell(verc.NumObjects , 1);
    EdVer = cell(edc.NumObjects , 1);
    
    for nv = 1:verc.NumObjects
        
        nCv = VerCell{nv};
        if length(nCv) > 1
            nCv(length(nCv) + 1) = nCv(1);
            intsVC = zeros(length(nCv) - 1 , 1);

            for ii = 1:length(nCv) - 1
                if isempty(intersect( CellEd{nCv(ii)} , CellEd{nCv(ii + 1)} ))
                    intsVC(ii) = 0;
                else
                    intsVC(ii) = intersect( CellEd{nCv(ii)} , CellEd{nCv(ii + 1)} );
                end
            end

            VerEd{nv} = unique(intsVC);
        end
        
        for edv = 1:length(VerEd{nv})
            if VerEd{nv}(edv) ~= 0
                EdVer{VerEd{nv}(edv)}(length(EdVer{VerEd{nv}(edv)})+1) = nv;
            end
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%% Array of Angles of Edges %%%%%%%%%%%%%%%%%%%%%%%%%
    
    edge_t = edc;
    angArray = regionprops(edge_t , 'Orientation') ;
    EdgeAng0 = cat(1 , angArray.Orientation) ;
    EdgeAng = EdgeAng0 ; 
    
    Tissue(t).EdAng = cell(edge_t.NumObjects , 1) ;

    for edd = 1:edge_t.NumObjects
        
        [yed , xed] = ind2sub([Ny , Nx] , edge_t.PixelIdxList{edd});
        if yed < Dors_y
            EdgeAng(edd) = - EdgeAng0(edd);
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%% Array of Angles of Cells %%%%%%%%%%%%%%%%%%%%%%%%%
    
    cell_t = cc1;
    angArray = regionprops(cell_t , 'Orientation') ;
    CellAng0 = cat(1 , angArray.Orientation) ;
    CellAng = CellAng0 ; 
    
    Tissue(t).CCAng = cell(cell_t.NumObjects , 1) ;

    for ccl = 1:cell_t.NumObjects
        
        [ycc , xcc] = ind2sub([Ny , Nx] , cell_t.PixelIdxList{ccl});
        if ycc < Dors_y
            CellAng(ccl) = - CellAng0(ccl);
        end
        
    end
    
%%%%%%%%%%% Array of Eccentricity, Major and Minor Axes of Cells %%%%%%%%%%
    
    CexArray = regionprops(cell_t , 'Eccentricity') ;
    Cex0 = cat(1 , CexArray.Eccentricity) ;
    
    CMajAx = regionprops(cell_t , 'MajorAxisLength') ;
    Cmj0 = cat(1 , CMajAx.MajorAxisLength) ;
    
    CMinAx = regionprops(cell_t , 'MinorAxisLength') ;
    Cmi0 = cat(1 , CMinAx.MinorAxisLength) ;

%%%%%%%%%%% Making the Structure of Tissue based on Properties %%%%%%%%%%%%

    Tissue(t).Cell = cc1        ;
    Tissue(t).Edge = edc        ;
    Tissue(t).Vert = verc       ;
    Tissue(t).CV = CellVer      ;
    Tissue(t).VC = VerCell      ;
    Tissue(t).EC = EdCell       ;
    Tissue(t).CE = CellEd       ;
    Tissue(t).VE = VerEd        ;
    Tissue(t).EV = EdVer        ;    
    Tissue(t).EL = EdLen        ;
    Tissue(t).EdAng = EdgeAng   ;
    Tissue(t).CCAng = CellAng   ;
    Tissue(t).Cex = Cex0        ;
    Tissue(t).Cmia = Cmj0       ;
    Tissue(t).Cmja = Cmi0       ;

end

%% Tracking the cells and edges based on the topology and connectivity of the network

cd(orgdir);

for t = 1 : Time - step

    clear CCMap_F ; clear CCMap_B ; clear EEMap_F ; clear EEMap_B ;
    
    cc1 = Tissue(t).Cell ;
    cc2 = Tissue(t + 1).Cell ;
    ec1 = Tissue(t).Edge ;
    ec2 = Tissue(t + 1).Edge ;
    
%%%%%%%%%%%%%%%%%%%%%% Flow Field From Raw Images %%%%%%%%%%%%%%%%%%%%%%%%%
    
    im1 = double(imread(MemName , t)) ;
    im2 = double(imread(MemName , t+1)) ;

    im1     = imresize(im1 , isf , 'bicubic')       ;
    im2     = imresize(im2 , isf , 'bicubic')       ;
    
    cd(orgdir)
    
    for iter = 1 : IterNum

        [VY1 , VX1] = GetPIV(im1,im2,Y1,X1,EdgeLength); 
        [VY2 , VX2] = GetPIV(im2,im1,Y1,X1,EdgeLength); 

        % smooth if desired
        if smooth == 1
            VX1  = imgaussfilt(VX1 , sigma);
            VY1  = imgaussfilt(VY1 , sigma);
            VX2  = imgaussfilt(VX2 , sigma);
            VY2  = imgaussfilt(VY2 , sigma);
        end

        [mesh_y , mesh_x] = size(VX1);

        Vx01 = zeros(size(im1 , 1) , size(im1 , 2))  ;
        Vy01 = zeros(size(im1 , 1) , size(im1 , 2))  ;
        Vx02 = zeros(size(im1 , 1) , size(im1 , 2))  ;
        Vy02 = zeros(size(im1 , 1) , size(im1 , 2))  ;

        for xx = 1:mesh_x
            for yy = 1:mesh_y
                Vx01( (yy-1)*EdgeLength+1 : yy*EdgeLength , (xx-1)*EdgeLength+1 : xx*EdgeLength ) = VX1(yy , xx) ;       
                Vy01( (yy-1)*EdgeLength+1 : yy*EdgeLength , (xx-1)*EdgeLength+1 : xx*EdgeLength ) = VY1(yy , xx) ;       
                Vx02( (yy-1)*EdgeLength+1 : yy*EdgeLength , (xx-1)*EdgeLength+1 : xx*EdgeLength ) = VX2(yy , xx) ;       
                Vy02( (yy-1)*EdgeLength+1 : yy*EdgeLength , (xx-1)*EdgeLength+1 : xx*EdgeLength ) = VY2(yy , xx) ;       
            end
        end
        
        if smooth == 1
            Vx01  = imgaussfilt(Vx01 , EdgeLength);
            Vy01  = imgaussfilt(Vy01 , EdgeLength);
            Vx02  = imgaussfilt(Vx02 , EdgeLength);
            Vy02  = imgaussfilt(Vy02 , EdgeLength);
        end
    
    end

    cd ../
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cen1 = regionprops(cc1 , 'centroid');
    cent1 = round(cat(1 , cen1.Centroid));
    per = [2 , 1];
    cent1 = cent1(:,per);
    L1 = labelmatrix(cc1);
    reg1 = logical(L1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cen2 = regionprops(cc2 , 'centroid');
    cent2 = round(cat(1 , cen2.Centroid));  
    per = [2 , 1];
    cent2 = cent2(:,per);        
    L2 = labelmatrix(cc2);
    reg2 = logical(L2);
    
%%%%%%%%% Evolving the centroids according to the velocity field %%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%% This following method screws up the matching for some reason! %%%%%%
%     VyId1 = zeros(cc1.NumObjects) ; VxId1 = zeros(cc1.NumObjects);
%     for c = 1:cc1.NumObjects
%         VyId1(c) = Vy01(cent1(c,1) , cent1(c,2));
%         VxId1(c) = Vx01(cent1(c,1) , cent1(c,2));
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    cenid1 = sub2ind([Ny , Nx] , cent1(:,1) , cent1(:,2));
    VyId1 = Vy01(cenid1);        VxId1 = Vx01(cenid1);
    V_vec12 = cat(2 , VyId1 , VxId1);

    cent1wrpd = cent1 + V_vec12 ;

    cenid2 = sub2ind([Ny , Nx] , cent2(:,1) , cent2(:,2));        
    VyId2 = Vy02(cenid2);        VxId2 = Vx02(cenid2);
    V_vec21 = cat(2 , VyId2 , VxId2);

    cent2wrpd = cent2 + V_vec21 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    CCMap_F = cell(cc1.NumObjects + 1 , 1) ;
    CCMap_B = cell(cc2.NumObjects + 1 , 1) ;

    for nccf = 1 : cc1.NumObjects
        CCMap_F{nccf} = L2( round(cent1wrpd(nccf,1)) , round(cent1wrpd(nccf,2)) );
    end

    for nccb = 1 : cc2.NumObjects    
        CCMap_B{nccb} = L1( round(cent2wrpd(nccb , 1)) , round(cent2wrpd(nccb , 2)) );
    end


    EEMap_F = cell(ec1.NumObjects + 1 , 1) ;
    EEMap_B = cell(ec2.NumObjects + 1 , 1) ;

    % cell1_(1 or 2) are the two cells sharing the edge 'neef' at
    % time 't'. cell2_(1 or 2) are the two cells sharing the 
    % corresponding edge at time 't+1'. 
    for neef12 = 1 : ec1.NumObjects
        if length(Tissue(t).EC{neef12}) == 2
    
            cell1_1 = Tissue(t).EC{neef12}(1) ; cell1_2 = Tissue(t).EC{neef12}(2) ;
            cell2_1 = CCMap_F{cell1_1} ; cell2_2 = CCMap_F{cell1_2} ;
  
            if cell2_1 ~= 0 && cell2_2 ~= 0 &&  cell2_1 ~= cell2_2
                EEMap_F{neef12} = intersect( Tissue(t+1).CE{cell2_1} , Tissue(t+1).CE{cell2_2} ) ;
            else 
                EEMap_F{neef12} = {};
            end
            
        end
    end

    for neef21 = 1 : ec2.NumObjects
        if length(Tissue(t+1).EC{neef21}) == 2
    
            cell2_1 = Tissue(t+1).EC{neef21}(1) ; cell2_2 = Tissue(t+1).EC{neef21}(2) ;
            cell1_1 = CCMap_B{cell2_1} ; cell1_2 = CCMap_B{cell2_2} ;
            % the following if can be removed if EC is constructed more
            % carefully to avoid double assignment to one edge or assiging
            % zero to them. 
            if cell1_1 ~= 0 && cell1_2 ~= 0 && cell1_1 ~= cell1_2
                EEMap_B{neef21} = intersect( Tissue(t).CE{cell1_1} , Tissue(t).CE{cell1_2} ) ;
            else
                EEMap_B{neef21} = {};
            end
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % These velocities are not to acted on by the metric. Only defined in
    % the 2D plane for tracking cells. 
    Tissue(t).VxFrwrd = Vx01;
    Tissue(t).VyFrwrd = Vy01;
    Tissue(t+1).VxBck = Vx02;
    Tissue(t+1).VyBck = Vy02;

    Tissue(t).CCmapF = CCMap_F ;
    Tissue(t+1).CCmapB = CCMap_B ;
    Tissue(t).EEmapF = EEMap_F ;
    Tissue(t+1).EEmapB = EEMap_B ;

end


% We need to first find the total number of paths to be able to assign to
% the edges. 
% to each edge (all of the edges of all time), we associate a lable which
% determines to what path it belongs. 
for edid = 1:Tissue(1).Edge.NumObjects
    Tissue(1).EdgeLabel0{edid} = edid ;
end    
% initial label value after labeling all the edges in t = 1;
Ed_label = Tissue(1).Edge.NumObjects + 1;

% The same strategy for cells is also used:
for cid = 1:Tissue(1).Cell.NumObjects
    Tissue(1).CellLabel0{cid} = cid ;
end
% initial label value after labeling all the edges in t = 1;
CC_label = Tissue(1).Cell.NumObjects + 1;

% After the initial labeling in the above line, we ask each edge at later
% times. If it is preceded by another edge that is already labeled, the new
% edge too, adopts the same label. If there exist no corresponding edge we
% assign a new label to it and raise the label by one. 
for t = 1 : tmax - tmin + 1 - 1
    
    clear Elabel_t
    clear Clabel_t

    Elabel_t = cell(Tissue(t+1).Edge.NumObjects , 1);
    
    for edid = 1:Tissue(t+1).Edge.NumObjects
        % Elable_t and Clabel_t are local temporary variables to store the
        % corresponding information at time 't', and transfer to an
        % actuall array of Tissue(t+1).
        if ~isempty(Tissue(t+1).EEmapB{edid})       % this 0 or 1 loop of logical propositions!
            Elabel_t{edid} = Tissue(t).EdgeLabel0{Tissue(t+1).EEmapB{edid}} ;    
        else
            Elabel_t{edid} = Ed_label;
            Ed_label = Ed_label + 1 ;
        end
        
    end

    Tissue(t+1).EdgeLabel0 = Elabel_t ;
    
    Clabel_t = cell(Tissue(t+1).Cell.NumObjects , 1);

    for cid = 1:Tissue(t+1).Cell.NumObjects
        
        if Tissue(t+1).CCmapB{cid} ~= 0       % this 0 or 1 loop of logical propositions!
            Clabel_t{cid} = Tissue(t).CellLabel0{Tissue(t+1).CCmapB{cid}} ;    
        else
            Clabel_t{cid} = CC_label;
            CC_label = CC_label + 1 ;
        end
        
    end

    Tissue(t+1).CellLabel0 = Clabel_t ;

end

Npath0 = Ed_label;
Cpath0 = CC_label;

% creating the Ed_Paths and CC_Paths cell-array which stores the (local) label of edges
% that belong to a certain path. Local here means the label of an edge in
% its own timepoint. 
% {N_EdPaths} = {};
% for lb = 1:N_EdPaths
%     Ed_Paths{lb} = zeros(1 , Time);
% end

Ed_Paths0 = zeros(Npath0 , Time) ;
CC_Paths0 = zeros(Cpath0 , Time) ;

for t = 1:Time
    for ed_t = 1:Tissue(t).Edge.NumObjects
        Ed_Paths0(Tissue(t).EdgeLabel0{ed_t} , t) = ed_t;        
    end    
    for cc = 1:Tissue(t).Cell.NumObjects
        CC_Paths0(Tissue(t).CellLabel0{cc} , t) = cc;        
    end
end

% Ed_Paths{ed} and {cell} are each a vector of length Time, so that we can
% know the interval within which each edge and cell appears. 

%%%%%%%%%%%%%%%% condition on min{timepoints} of tracking %%%%%%%%%%%%%%%%%

EdPaths = Ed_Paths0;

ee = 0;
for pt = 1:Npath0
    edpt = EdPaths(pt-ee , :);
    if length(edpt(edpt ~= 0)) < minTP
        EdPaths(pt-ee , :) = [];
        ee = ee + 1;        
    end
end

Npath = size(EdPaths , 1);

for tt = 1 : tmax - tmin + 1
    
    NedT = Tissue(tt).Edge.NumObjects;
    edT = EdPaths(: , tt) ;
    for ee = 1:NedT
        if ~isempty(find( edT == ee , 1 ))
            Tissue(tt).EdgeLabel{ee} = find( edT == ee );
        else
            Tissue(tt).EdgeLabel{ee} = 0;
        end
    end
    
end

disp(Npath)

global NpathF
NpathF = Npath;

%%%%%%%%%%%%%%%%%%%%% Minimum Time of Cell Tracking %%%%%%%%%%%%%%%%%%%%%%%

CCPaths = CC_Paths0;

cc = 0;
for cpt = 1:Cpath0
    edpt = CCPaths(cpt-cc , :);
    if length(edpt(edpt ~= 0)) < minTP
        CCPaths(cpt-cc , :) = [];
        cc = cc + 1;        
    end
end

Cpath = size(CCPaths , 1);

for tt = 1 : tmax - tmin + 1
    
    NcT = Tissue(tt).Cell.NumObjects;
    cT = CCPaths(: , tt) ;
    for cc = 1:NcT
        if ~isempty(find( cT == cc , 1 ))
            Tissue(tt).CellLabel{cc} = find( cT == cc );
        else
            Tissue(tt).CellLabel{cc} = 0;
        end
    end
    
end

disp(Cpath)

global CpathF
CpathF = Cpath;

%% This script is supposed to associate the Myosin cables to edges and find
% the projection of Myosin intensity along each edge. There is a level of
% uncertainty as the cables do NOT match perfectly with the edges, often
% thicker and/or with an offset w.r.t. the edges. 

meanMyo = zeros(Time , 1);
minMyo = zeros(Time , 1);
maxMyo = zeros(Time , 1);

myMax = 0;

bg_rad = 20;
bg0 = zeros(Time , 1);
for t = 1:Time
    
    MaskName = sprintf('Sources/NewROI/Laterals/LateralMask/TotalMask_t%02d.tif' , t+4);
    mask = double(imread(MaskName))./255;
    mask(: , 1:600) = 0;
    NZind = find(mask ~= 0) ;
    
    myo = imread(MyoName , t) ;
    
    bg_t = imopen(myo , strel('disk' ,  bg_rad) );
    bg0(t) = mean(bg_t(NZind));
   
end

bgMax = max(bg0);

minMy = 1000;

for t = 1:Time
    
    clear EdMyo_arr;
    EdMyoDensity = zeros(Ny , Nx);

    myo = double(imread(MyoName , t)) ;
    myo = myo .* double(bgMax ./ bg0(t));
    bg = imopen(myo , strel('disk' ,  bg_rad) );
    myo = double(imsubtract(myo , bg));
    
    meanMyo(t) = mean(mean(myo));
    minMyo(t) = min(min(myo));
    maxMyo(t) = max(max(myo));
    
    ed0 = Tissue(t).Edge ;
    Led0 = labelmatrix(ed0) ; 
    Ed0Seg = logical( labelmatrix(ed0) ) ;
    EdDiag = bwmorph(Ed0Seg , 'diag') ;

    myowid = 2;
    EdFat = bwmorph(EdDiag , 'thick' , myowid) ;
    eF = bwconncomp(EdFat , 8) ;
    LF = labelmatrix(eF) ;
    
    thin_thick_ratio = sum(logical(Led0(:)))/sum(logical(LF(:)));
    
    EdMyo0 = (double(EdFat)) .* double(myo);
    EdMyo_arr = zeros(ed0.NumObjects , 1);

    for e = 1:ed0.NumObjects
        if Tissue(t).EL(e) ~= 0
            edPxl = eF.PixelIdxList{ LF(ed0.PixelIdxList{e}(ceil(end/2))) }; 
            MyoInt = sum( EdMyo0( edPxl ) ) .* thin_thick_ratio ./ (Tissue(t).EL(e))  ;% * (2*myowid + 1)) ;%(length(edF)/(2*wid + 1))  ;
            EdMyo_arr(e , 1) = MyoInt;
            EdMyoDensity(ed0.PixelIdxList{e}) = MyoInt ;
        else 
            EdMyo_arr(e , 1) = 0;
            EdMyoDensity(ed0.PixelIdxList{e}) = 0 ;
        end 
    end
        
    Tissue(t).EdFat = eF ;
    Tissue(t).EdMyo = EdMyo_arr ;
    Tissue(t).EdMyoImg = EdMyoDensity;

    myMax = max(myMax , max(EdMyoDensity(:)));

    EdMinMy = min(EdMyo_arr(EdMyo_arr ~= 0));
    minMy = min(minMy , EdMinMy);
    
end


%%

cd(orgdir)

global wid
wid = 5;

%%%%%%%%%%%%%%%%%%% All the data including conditions %%%%%%%%%%%%%%%%%%%%%

[Ed_H , Ed_V , Ed_P , Ed_N , Lpt , Mpt , Apt] = M_R_condition(Tissue, EdPaths);

%%%%%%%%%%%%%%%%%%%%% Intersection of the conditions %%%%%%%%%%%%%%%%%%%%%%

Lpt = Lpt./mean(Lpt(:));
Mpt = Mpt./mean(Mpt(:));

dMpt = diff(log(Mpt) , 1 , 2);
dLpt = diff(log(Lpt) , 1 , 2);

PH = intersect(Ed_P , Ed_H);
PV = intersect(Ed_P , Ed_V);
NH = intersect(Ed_N , Ed_H);
NV = intersect(Ed_N , Ed_V);

PP = union(PH , PV);
NN = union(NH , NV);
VV = union(NV , PV);
HH = union(NH , PH);

AllCond = HH;
AllCond = union(AllCond , VV);
AllCond = union(AllCond , PP);
AllCond = union(AllCond , NN);

%%

EpsXX = zeros(Ny , Nx , Time - 1);
EpsXY = zeros(Ny , Nx , Time - 1);
EpsYY = zeros(Ny , Nx , Time - 1);

for tt = 1:Time - 1

    Uy = Tissue(tt).VyFrwrd;
    Ux = Tissue(tt).VxFrwrd;

    [Exx,Exy] = gradient(Ux);
    [Eyx,Eyy] = gradient(Uy);

    EpsXX(: , : , tt) = Exx;
    EpsXY(: , : , tt) = (Exy + Eyx)/2;
    EpsYY(: , : , tt) = Eyy;

end

%%      Coarse-grained continuum analysis of strain-stress

xl_lat = 601:1300;
yl_lat_l = 1301:1850;
yl_lat_r = 201:800;

PsizeY = 100 ; PsizeX = 100 ;
StepY  = 100 ; StepX  = 100 ; 

Xs = xl_lat(1) + floor(PsizeX/2) : StepX : xl_lat(end) ; 
Ysl = yl_lat_l(1) + floor(PsizeY/2) : StepY : yl_lat_l(end);
Ysr = yl_lat_r(1) + floor(PsizeY/2) : StepY : yl_lat_r(end);
Ys = union(Ysl , Ysr);

nX = length(Xs) ; nY = length(Ys) ;
[Xpts , Ypts] = meshgrid(Xs , Ys) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% t = Time %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tissue(Time).PatchBkw.AngLen{nY}{nX} = [];

for yy = 1:nY
    for xx = 1:nX
        xblock = Xpts(yy , xx) - floor(PsizeX/2) + 1 : Xpts(yy , xx) + floor(PsizeX/2);
        yblock = Ypts(yy , xx) - floor(PsizeY/2) + 1 : Ypts(yy , xx) + floor(PsizeY/2);
        [xb , yb] = meshgrid(xblock , yblock);
        ptsXY = cat(2 , yb(:) , xb(:));
        pts = sub2ind([Ny , Nx] , ptsXY(:,1), ptsXY(:,2) );
        Tissue(Time).PatchBkw.Pts{yy}{xx} = pts;        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% t = 1 : Time - 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tt = 1:Time - 1
    for yy = 1:nY
        for xx = 1:nX
            pts_t1 = Tissue(Time - tt + 1).PatchBkw.Pts{yy}{xx};
            [ptsXY_t1y , ptsXY_t1x] = ind2sub([Ny , Nx] , pts_t1);
            ptsXY_t1 = cat(2 , ptsXY_t1y , ptsXY_t1x);

            VyId = Tissue(Time - tt + 1).VyBck ;   VyB = VyId(pts_t1) ;
            VxId = Tissue(Time - tt + 1).VxBck ;   VxB = VxId(pts_t1) ;
            V_vec21 = cat(2 , VyB , VxB) ;

            ptsXY_t2 = round(ptsXY_t1 + V_vec21) ;
            pts_t2 = sub2ind( [Ny , Nx] , ptsXY_t2(:,1) , ptsXY_t2(:,2) );
            pts_t2 = round(pts_t2);             % index of pixels in the box
            % shifting +/- 1, i.e. up-down , and +/- Ny, i.e. right and
            % left to remove the gaps
            pts_t2 = union(pts_t2 , pts_t2 + 1);
            pts_t2 = union(pts_t2 , pts_t2 - 1);
            pts_t2 = union(pts_t2 , pts_t2 + Ny);
            pts_t2 = union(pts_t2 , pts_t2 - Ny);

            Tissue(Time - tt).PatchBkw.Pts{yy}{xx} = pts_t2;
        end
    end
end

%  Checking that forward moving of patches moves them back to last time point

ed1 = Tissue(1).Edge;
EdL0 = labelmatrix(ed1) ;
nED0 = zeros(nY , nX , 1);

for yy = 1:nY
    for xx = 1:nX
        pts = Tissue(1).PatchBkw.Pts{yy}{xx};
        unqs = unique(EdL0(pts)); 
        unqEds = unqs(unqs ~= 0);

        Tissue(1).PatchFrw.Pts{yy}{xx} = Tissue(1).PatchBkw.Pts{yy}{xx};        
        Tissue(1).PatchFrw.Eds{yy}{xx} = unqEds;
        nED0(yy , xx , 1) = length(unqEds);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% t = 1 : Time - 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for tt = 2:Time

    edt = Tissue(tt).Edge;
    EdL0 = labelmatrix(edt);

    for yy = 1:nY
        for xx = 1:nX
            pts_t1 = Tissue(tt-1).PatchFrw.Pts{yy}{xx};
            [ptsXY_t1y , ptsXY_t1x] = ind2sub([Ny , Nx] , pts_t1);
            ptsXY_t1 = cat(2 , ptsXY_t1y , ptsXY_t1x);

            VyId = Tissue(tt-1).VyFrwrd ;   VyB = VyId(pts_t1) ;
            VxId = Tissue(tt-1).VxFrwrd ;   VxB = VxId(pts_t1) ;
            V_vec21 = cat(2 , VyB , VxB) ;

            ptsXY_t2 = round(ptsXY_t1 + V_vec21) ;
            pts_t2 = sub2ind( [Ny , Nx] , ptsXY_t2(:,1) , ptsXY_t2(:,2) );
            pts_t2 = round(pts_t2);             % index of pixels in the box
            % shifting +/- 1, i.e. up-down , and +/- Ny, i.e. right and
            % left to remove the gaps
            pts_t2 = union(pts_t2 , pts_t2 + 1);
            pts_t2 = union(pts_t2 , pts_t2 - 1);
            pts_t2 = union(pts_t2 , pts_t2 + Ny);
            pts_t2 = union(pts_t2 , pts_t2 - Ny);
            unqs_t2 = unique(EdL0((pts_t2))) ; 
            unqEds_t2 = unqs_t2(unqs_t2 ~= 0) ; % label of edges in the box

            Tissue(tt).PatchFrw.Eds{yy}{xx} = unqEds_t2;
            Tissue(tt).PatchFrw.Pts{yy}{xx} = pts_t2;
            nED0(yy , xx , tt) = length(unqEds_t2);
        end
    end

end

NEDGES = nED0(nED0(:) ~= 0);
nED = nED0;
avNed = mean(NEDGES(:));

%%

PtcMyo = zeros(nY , nX , Time);
yyStress = zeros(nY , nX , Time);
xyStress = zeros(nY , nX , Time);
xxStress = zeros(nY , nX , Time);

Myy = zeros(1 , Time);
Mxx = zeros(1 , Time);
Mxy = zeros(1 , Time);

Eyy = zeros(1 , Time - 1);
Exx = zeros(1 , Time - 1);
Exy = zeros(1 , Time - 1);

TrMM = zeros(nY , nX , Time);

for tt = 1:Time
    MM = zeros(nY , nX , 2 , 2);
    for yy = 1:nY
        for xx = 1:nX
            if mean(nED0(yy , xx , :)) > avNed
                edsptc = Tissue(tt).PatchFrw.Eds{yy}{xx};
                nCond = 0;
                Len_ed = 0;
                for eds = 1:length(edsptc)
                    % ALL angles must be included, otherwise the patch's
                    % distortion in x direction won't include the
                    % corresponding stress!!!
                    if abs(Tissue(tt).EdAng(edsptc(eds))) >= 0 && Tissue(tt).EL(edsptc(eds)) ~= 0 && Tissue(tt).EdMyo(edsptc(eds)) ~= 0
                        nCond = nCond + 1;
                        Nem11 = cosd(Tissue(tt).EdAng(edsptc(eds))).^2;
                        Nem12 = sind(Tissue(tt).EdAng(edsptc(eds))) .* cosd(Tissue(tt).EdAng(edsptc(eds))); 
                        Nem21 = Nem12 ; 
                        Nem22 = sind(Tissue(tt).EdAng(edsptc(eds))).^2;

                        Nem = 2*[Nem11 , Nem12 ; Nem21 , Nem22]; 
                        KK = zeros(1,1,2,2);
                        KK(1,1,:,:) = Nem;
                        
                        Len_ed = Len_ed + Tissue(tt).EL(edsptc(eds));
                        PtcMyo(yy , xx , tt) = PtcMyo(yy , xx , tt) + Tissue(tt).EdMyo(edsptc(eds));
                        MM(yy , xx , : , :) = MM(yy , xx , : , :) + (Tissue(tt).EdMyo(edsptc(eds)) * Tissue(tt).EL(edsptc(eds)) .* KK(1,1,:,:));
                    end
                end
                if nCond ~= 0
                    MM(yy , xx , : , :) = (MM(yy , xx , : , :) / Len_ed);
                    PtcMyo(yy , xx , tt) = PtcMyo(yy , xx , tt) / nCond;
                end
                TrMM(yy , xx , tt) = trace(squeeze(squeeze(MM(yy , xx , : , :))));
            end
        end
    end

    xxStress(: , : , tt) = MM(:, : , 1 , 1);
    xyStress(: , : , tt) = MM(:, : , 1 , 2);
    yyStress(: , : , tt) = MM(:, : , 2 , 2);
 
    ind = 0;
    for ii = 1:nX
        for jj = 1:nY 
            if yyStress(jj , ii , tt) ~= 0
                ind = ind + 1;
                Myy(ind , tt) = yyStress(jj , ii , tt);
                Mxx(ind , tt) = xxStress(jj , ii , tt);
                Mxy(ind , tt) = xyStress(jj , ii , tt);
            end
        end
    end

end

xxStressArr = (Mxx(:,1:end-1) + Mxx(:,2:end))/2;
xyStressArr = (Mxy(:,1:end-1) + Mxy(:,2:end))/2;
yyStressArr = (Myy(:,1:end-1) + Myy(:,2:end))/2;

Myy = zeros(1 , Time - 1);
Mxx = zeros(1 , Time - 1);
Mxy = zeros(1 , Time - 1);

NzPatches = zeros(1 , 2);

majEps = zeros(1 , Time - 1);
majSig = zeros(1 , Time - 1);
minEps = zeros(1 , Time - 1);
minSig = zeros(1 , Time - 1);

Projection = zeros(1 , Time - 1);

for tt = 1:Time - 1
    
    EpsYY_t = squeeze(EpsYY(:,:,tt));
    EpsXY_t = squeeze(EpsXY(:,:,tt));
    EpsXX_t = squeeze(EpsXX(:,:,tt));
    
    ind = 0;
    for ii = 1:nX
        for jj = 1:nY 
            if yyStress(jj , ii , tt) ~= 0                
                ind = ind + 1;
                NzPatches(ind , 1) = jj;
                NzPatches(ind , 2) = ii;
                
                Myy(ind , tt) = yyStress(jj , ii , tt);
                Mxx(ind , tt) = xxStress(jj , ii , tt);
                Mxy(ind , tt) = xyStress(jj , ii , tt);
                
                Eyy(ind , tt) = mean(EpsYY_t(Tissue(tt).PatchFrw.Pts{jj}{ii}));
                Exy(ind , tt) = mean(EpsXY_t(Tissue(tt).PatchFrw.Pts{jj}{ii}));
                Exx(ind , tt) = mean(EpsXX_t(Tissue(tt).PatchFrw.Pts{jj}{ii}));
                
                Sigma(1,1 , ind , tt) = Mxx(ind , tt);
                Sigma(1,2 , ind , tt) = Mxy(ind , tt);
                Sigma(2,1 , ind , tt) = Mxy(ind , tt);
                Sigma(2,2 , ind , tt) = Myy(ind , tt);

                Epsilon(1,1 , ind , tt) = Exx(ind , tt);
                Epsilon(1,2 , ind , tt) = Exy(ind , tt);
                Epsilon(2,1 , ind , tt) = Exy(ind , tt);
                Epsilon(2,2 , ind , tt) = Eyy(ind , tt);
                
                Eps_t = [squeeze(Exx(ind , tt)) , squeeze(Exy(ind , tt)) ; squeeze(Exy(ind , tt)) , squeeze(Eyy(ind , tt)) ];
                Sig_t = [squeeze(Mxx(ind , tt)) , squeeze(Mxy(ind , tt)) ; squeeze(Mxy(ind , tt)) , squeeze(Myy(ind , tt)) ];

                [Veps , Deps] = eig(Eps_t);
                [Vsig , Dsig] = eig(Sig_t);
                                
                majEps(ind , tt) = Deps(2,2);%max(eig(Eps_t));
                minEps(ind , tt) = Deps(1,1);

                majSig(ind , tt) = Veps(:,2)' * Sig_t * Veps(:,2);%Dsig(2,2);%
                minSig(ind , tt) = Veps(:,1)' * Sig_t * Veps(:,1);%Dsig(1,1);%

                Projection(ind , tt) = abs(Veps(:,2)' * Vsig(:,2));
                
            end
        end
    end
end

reg = 10;

figure ; plot(1:Time - 1 , majEps(reg , :)./mean(majEps(reg , :)) , 'r-o' , 1:Time - 1 , majSig(reg , :)./mean(majSig(reg , :)) , 'b-o');% , 1:Time - 1 , Exy(reg , :) , 'g-*')
title('Major Axes');
figure ; plot(1:Time - 1 , minEps(reg , :)./mean(minEps(reg , :)) , 'r-o' , 1:Time - 1 , minSig(reg , :)./mean(minSig(reg , :)) , 'b-o');% , 1:Time - 1 , Mxy(reg , :) , 'g-*')
title('Minor Axes');

t1 = 7;

figure ; plot(1:ind , majEps(: , t1)./max(majEps(:)) , 'b-o' , 1:ind , majSig(: , t1)./max(majSig(:)) , 'r-o');% , 1:Time - 1 , Exy(reg , :) , 'g-*')
title('Major Axes');
figure ; plot(1:ind , minEps(: , t1)./max(minEps(:)) , 'b-o' , 1:ind , minSig(: , t1)./max(minSig(:)) , 'r-o');% , 1:Time - 1 , Mxy(reg , :) , 'g-*')
title('Minor Axes');

TmajCorr = zeros(ind , 1);
TminCorr = zeros(ind , 1);
for reg = 1:ind
    
    majX = corrcoef(majEps(reg,:) , majSig(reg,:));
    TmajCorr(reg , 1) = majX(1,2);
    minX = corrcoef(minEps(reg,:) , minSig(reg,:));
    TminCorr(reg , 1) = minX(1,2);
    
end
    
figure ; plot(1:ind , TmajCorr , 'b-o' , 1:ind , TminCorr , 'r-*');

majCorr = zeros(Time - 1 , 1);
minCorr = zeros(Time - 1 , 1);
for tt = 1:Time-1
    
    majX = corrcoef(majEps(:,tt) , majSig(:,tt));
    majCorr(tt , 1) = majX(1,2);
    minX = corrcoef(minEps(:,tt) , minSig(:,tt));
    minCorr(tt , 1) = minX(1,2);
    
end
    
figure ; plot(1:Time - 1 , majCorr , 'b-o' , 1:Time - 1 , minCorr , 'r-*');

maxEpsyy = max(Eyy(:));
maxEpsxx = max(Exx(:));
maxEps = max(maxEpsxx , maxEpsyy);
for tt = 1:Time - 1
    AA = zeros(Ny , Nx);
    BB = zeros(Ny , Nx);
    CC = zeros(Ny , Nx);
    for reg = 1:ind    
        AA(Tissue(tt).PatchFrw.Pts{NzPatches(reg , 1)}{NzPatches(reg , 2)}) = Eyy(reg , tt)/maxEps;
        BB(Tissue(tt).PatchFrw.Pts{NzPatches(reg , 1)}{NzPatches(reg , 2)}) = Myy(reg , tt)/maxEps;
          CC(Tissue(tt).PatchFrw.Pts{NzPatches(reg , 1)}{NzPatches(reg , 2)}) = Projection(reg , tt);%/maxEps;
    end
    figure ; imshow(CC);%imshowpair(rescale(AA) , rescale(BB));
end

%%

xxStressRateArr = diff(log(Mxx) , 1 , 2);
xyStressRateArr = diff(log(Mxy) , 1 , 2);
yyStressRateArr = diff(log(Myy) , 1 , 2);

totStress = yyStressArr + xxStressArr;
totStressRate = yyStressRateArr + xxStressRateArr;

totptcs = size(yyStressRateArr , 1);
shuffptcs1 = randperm(totptcs , totptcs);
shuffptcs2 = randperm(totptcs , totptcs);

yyStressRateShuffled = yyStressRateArr(shuffptcs1 , :);
xxStressRateShuffled = xxStressRateArr(shuffptcs1 , :);
totStressRateShuffled = totStressRate(shuffptcs1 , :);

yyStressShuffled = yyStressArr(shuffptcs1 , :);
xxStressShuffled = xxStressArr(shuffptcs1 , :);
totStressShuffled = totStress(shuffptcs1 , :);

figure ; plot(1:Time-1 , mean(yyStressRateArr) , 'r-o');
title('d\langle \sigma_{yy} \rangle /dt')
figure ; plot(1:Time-1 , mean(yyStrainRateArr) , 'b-o');
title('d\langle \epsilon_{yy} \rangle /dt')
figure ; plot(1:Time-1 , mean(xxStressRateArr) , 'r-*');
title('d\langle \sigma_{xx} \rangle /dt')
figure ; plot(1:Time-1 , mean(xxStrainRateArr) , 'b-*');
title('d\langle \epsilon_{xx} \rangle /dt')

figure ; plot(1:Time-1 , mean(yyStressArr) , 'r-o');
title('\langle \sigma_{yy} \rangle')
figure ; plot(1:Time , mean(yyStrainArr) , 'b-o');
title('\langle \epsilon_{yy} \rangle')
figure ; plot(1:Time-1 , mean(xxStressArr) , 'r-*');
title('\langle \sigma_{xx} \rangle')
figure ; plot(1:Time , mean(xxStrainArr) , 'b-*');
title('\langle \epsilon_{xx} \rangle')

figure ; plot(1:Time-1 , mean(yyStressArr + xxStressArr) , 'r-o');
title('\langle \sigma_{tot} \rangle')
figure ; plot(1:Time , mean(yyStrainArr + xxStrainArr) , 'b-o');
title('\langle \epsilon_{tot} \rangle')

%%

CArea = zeros(Cpath , Time);
CMyo = zeros(Cpath , Time);
CPer = zeros(Cpath , Time);
Cexx = zeros(Cpath , Time);

for tt = 1:Time
    for cc = 1:Cpath

        cell_t = CCPaths(cc , tt);
        NumEd = length(Tissue(tt).CE{cell_t});
        cellmyo = 0;
        cellper = 0;
        
        for edcell = 1:NumEd
            cellper = cellper + Tissue(tt).EL(Tissue(tt).CE{cell_t}(edcell)) ;
        end    
        for edcell = 1:NumEd
            cellmyo = cellmyo + Tissue(tt).EdMyo(Tissue(tt).CE{cell_t}(edcell))...
                      .* Tissue(tt).EL(Tissue(tt).CE{cell_t}(edcell))./cellper ;
        end
        
        CArea(cc , tt) = length(Tissue(tt).Cell.PixelIdxList{CCPaths(cc , tt)});
        CMyo(cc , tt) = cellmyo/2; % divide by 2 because edges are shared
        CPer(cc , tt) = cellper;
        Cexx(cc , tt) = Tissue(tt).Cex(CCPaths(cc , tt));
    end
end

CArea = smoothdata(CArea , 2 , 'Gaussian' , wid);
CMyo = smoothdata(CMyo , 2 , 'Gaussian' , wid);
CPer = smoothdata(CPer , 2 , 'Gaussian' , wid);

CsizeMean = mean(CArea(:));
CMyoMean = mean(CMyo(:));

figure ; plot(1:Time , mean(CMyo)/mean(CMyo(:)));
figure ; plot(1:Time , (mean(CPer)/mean(CPer(:))));
figure ; plot(1:Time , mean(CArea)/CsizeMean , 1:Time , mean(CMyo)/CMyoMean);
figure ; plot(1:Time , mean(Cexx));

%%

dCellM = diff(log(CMyo) , 1 , 2);
dCellP = diff(log(CPer) , 1 , 2);
dCellS = diff(log(CArea) , 1 , 2);
dCellX = diff(Cexx , 1 , 2);

nMp = 0;
nMn = 0;
nSp = 0;
nSn = 0;
eMp = zeros(1,1);
eMn = zeros(1,1);
eSp = zeros(1,1);
eSn = zeros(1,1);

for cc = 1:Cpath    
    if CMyo(cc,end) - CMyo(cc,1) > 0.2 * mean(CMyo(cc,:)) && all(abs(dCellM(cc,:)) < 0.5) 
        nMp = nMp + 1;
        eMp(nMp) = cc;
    elseif CMyo(cc,end) - CMyo(cc,1) < - 0.2 * mean(CMyo(cc,:)) && all(abs(dCellM(cc,:)) < 0.5)
        nMn = nMn + 1;
        eMn(nMn) = cc;
    end
    if CArea(cc,end) - CArea(cc,1) > 0.2 * mean(CArea(cc,:)) && all(abs(dCellS(cc,:)) < 0.5)
        nSp = nSp + 1;
        eSp(nSp) = cc;
    elseif CArea(cc,end) - CArea(cc,1) < - 0.2 * mean(CArea(cc,:)) && all(abs(dCellS(cc,:)) < 0.5)
        nSn = nSn + 1;
        eSn(nSn) = cc;
    end
end

CondCell = eMp;
CondCell = union(CondCell , eMn);
CondCell = union(CondCell , eSp);
CondCell = union(CondCell , eSn);


%% Coarse Graining

xl_lat = 601:1400;
yl_lat_r = 301:800;
yl_lat_l = 1301:1800;

PsizeX = 100; 
PsizeY = 100;
StepX = 100; 
StepY = 100;

Xs = xl_lat(1) + floor(PsizeX/2) : StepX : xl_lat(end);
Ysl = yl_lat_l(1) + floor(PsizeY/2) : StepY : yl_lat_l(end);
Ysr = yl_lat_r(1) + floor(PsizeY/2) : StepY : yl_lat_r(end);
Ys = union(Ysl , Ysr);
nX = length(Xs); nY = length(Ys);

cc0 = Tissue(1).Cell;
Lc0 = labelmatrix(cc0);

ec0 = Tissue(1).Edge;
Le0 = labelmatrix(ec0);

Cbox = zeros(nY , nX);
CPatchArray = cell(nY,nX);
EPatchArray = cell(nY,nX);

for yy = 1:nY
    for xx = 1:nX

        Ybox = Ys(yy) - round(PsizeY/2) : Ys(yy) + round(PsizeY/2);
        Xbox = Xs(xx) - round(PsizeX/2) : Xs(xx) + round(PsizeX/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cell Patches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

        PtchCell = unique(Lc0( Ybox , Xbox ));
        PtchCell(PtchCell==0) = [];
        PtchCellCopy = PtchCell;
        [Ymesh , Xmesh] = meshgrid(Ybox , Xbox);
        Inds = sub2ind([Ny , Nx] , Ymesh , Xmesh);

        for c1 = 1:length(PtchCell)
            cPart = Tissue(1).Cell.PixelIdxList{PtchCell(c1)};
            portC1 = intersect(cPart , Inds(:));
            if length(portC1) < round(length(cPart))/2 || ismember(PtchCell(c1) , CCPaths(:,1)) == 0
                PtchCellCopy(PtchCellCopy == PtchCell(c1)) = [];
            end
        end

        CPatchArray{yy , xx} = PtchCellCopy;
        clear PtchCell;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Edge Patches %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PtchEd = unique(Le0( Ybox , Xbox ));
        PtchEd(PtchEd==0) = [];
        PtchEdCopy = PtchEd;
        [Ymesh , Xmesh] = meshgrid(Ybox , Xbox);
        Inds = sub2ind([Ny , Nx] , Ymesh , Xmesh);

        for c1 = 1:length(PtchEd)
            ePart = Tissue(1).Edge.PixelIdxList{PtchEd(c1)};
            portE1 = intersect(ePart , Inds(:));
            if length(portE1) < round(length(ePart))/2 || ismember(PtchEd(c1) , EdPaths(:,1)) == 0
                PtchEdCopy(PtchEdCopy == PtchEd(c1)) = [];
            end
        end

        EPatchArray{yy , xx} = PtchEdCopy;
        clear PtchEd;

    end
end

CCPATCH = cell(nY , nX , Time);
CMptch_t = cell(nY , nX , Time);
CSptch_t = cell(nY , nX , Time);
CPptch_t = cell(nY , nX , Time);
CAptch_t = cell(nY , nX , Time);
CMiptch_t = cell(nY , nX , Time);
CMjptch_t = cell(nY , nX , Time);
CExptch_t = cell(nY , nX , Time);

EDPATCH = cell(nY , nX , Time);
EMptch_t = cell(nY , nX , Time);
ELptch_t = cell(nY , nX , Time);
EAptch_t = cell(nY , nX , Time);

CMptch = zeros(nY , nX , Time);
CPptch = zeros(nY , nX , Time);
CSptch = zeros(nY , nX , Time);
CAptch = zeros(nY , nX , Time);
CMiptch = zeros(nY , nX , Time);
CMjptch = zeros(nY , nX , Time);
CExptch = zeros(nY , nX , Time);

EMptch = zeros(nY , nX , Time);
ELptch = zeros(nY , nX , Time);
EAptch = zeros(nY , nX , Time);

nE_XY = zeros(nY , nX);
nC_XY = zeros(nY , nX);

for tt = 1:Time
    
    for yy = 1:nY 
        for xx = 1:nX
            if isempty(CPatchArray{yy , xx}) == 0
                cptch = CPatchArray{yy , xx};
                clabel_t = CCPaths(cell2mat(Tissue(1).CellLabel(cptch)) , tt);
                CSptch_t{yy , xx , tt} = CArea(cell2mat(Tissue(1).CellLabel(cptch)) , tt);
                CPptch_t{yy , xx , tt} = CPer(cell2mat(Tissue(1).CellLabel(cptch)) , tt);
                CMptch_t{yy , xx , tt} = CMyo(cell2mat(Tissue(1).CellLabel(cptch)) , tt);
                CAptch_t{yy , xx , tt} = Tissue(tt).CCAng(clabel_t);
                CExptch_t{yy , xx , tt} = Tissue(tt).Cex(clabel_t);
                CMiptch_t{yy , xx , tt} = Tissue(tt).Cmia(clabel_t);
                CMjptch_t{yy , xx , tt} = Tissue(tt).Cmja(clabel_t);
            end    
            CSptch(yy , xx , tt) = mean(CSptch_t{yy , xx , tt});
            CMptch(yy , xx , tt) = mean(CMptch_t{yy , xx , tt});
            CPptch(yy , xx , tt) = mean(CPptch_t{yy , xx , tt});
            CAptch(yy , xx , tt) = mean(CAptch_t{yy , xx , tt});
            CExptch(yy , xx , tt) = mean(CExptch_t{yy , xx , tt});
            CMiptch(yy , xx , tt) = mean(CMiptch_t{yy , xx , tt});
            CMjptch(yy , xx , tt) = mean(CMjptch_t{yy , xx , tt});
            ncXY = size(CSptch_t{yy , xx , tt});
            nC_XY(yy , xx , tt) = ncXY(1);
        end
    end
    
    for yy = 1:nY 
        for xx = 1:nX
            if isempty(EPatchArray{yy , xx}) == 0
                eptch = EPatchArray{yy , xx};
                elabel_t = EdPaths(cell2mat(Tissue(1).EdgeLabel(eptch)) , tt);
                ELptch_t{yy , xx , tt} = Lpt(cell2mat(Tissue(1).EdgeLabel(eptch)) , tt);
                EMptch_t{yy , xx , tt} = Mpt(cell2mat(Tissue(1).EdgeLabel(eptch)) , tt);
                EAptch_t{yy , xx , tt} = Tissue(tt).EdAng(elabel_t);
            end   
            ELptch(yy , xx , tt) = mean(ELptch_t{yy , xx , tt});
            EMptch(yy , xx , tt) = mean(EMptch_t{yy , xx , tt});
            EAptch(yy , xx , tt) = mean(EAptch_t{yy , xx , tt});
            neXY = size(ELptch_t{yy , xx , tt});
            nE_XY(yy , xx , tt) = neXY(1);
        end
    end
    
end

nC_XY0 = squeeze(nC_XY(:,:,1));
nE_XY0 = squeeze(nE_XY(:,:,1));

Ncell = sum(nC_XY0(:));
Nptch = length(find(nC_XY0 ~= 0));
meanNcell = round(Ncell/Nptch);

%%

ELarray = reshape(ELptch , [nY*nX , Time]);
EMarray = reshape(EMptch , [nY*nX , Time]);

dELarray = diff(log(ELarray) , 1 , 2);
dEMarray = diff(log(EMarray) , 1 , 2);

CSarray = reshape(CSptch , [nY*nX , Time]);
CMarray = reshape(CMptch , [nY*nX , Time]);

dCSarray = diff(log(CSarray) , 1 , 2);
dCMarray = diff(log(CMarray) , 1 , 2);

index = 0;
for pt = 1:size(dEMarray , 1)
    if isnan(mean(dELarray(pt - index , :) , 2))
        ELarray(pt - index , :) = [];
        EMarray(pt - index , :) = [];
        dELarray(pt - index , :) = [];
        dEMarray(pt - index , :) = [];
        index = index + 1;
    end
end

figure ; plot(1:Time - 1 , mean(dELarray , 1)./mean(dELarray(:)) , 'b-o' , 1:Time - 1 , mean(dEMarray , 1)./mean(dEMarray(:)) , 'r-*')
title('R^{-1}dR_{ed}/dt (blue) , M^{-1}dM_{ed}/dt (red) vs. Time')
figure ; plot(1:Time , mean(ELarray , 1)./mean(ELarray(:)) , 'b-o' , 1:Time , mean(EMarray , 1)./mean(EMarray(:)) , 'r-*')
title('R_{ed} (blue) , M_{ed} (red) vs. Time')

[cX , lgX] = xcorr(mean(dELarray , 1)./abs(mean(dELarray(:))) , mean(dELarray , 1)./abs(mean(dELarray(:))) , 3 , 'coef');
cX = (Time - 1) .* cX ./ (Time - 1 - abs(lgX));
figure ; plot(lgX , cX , 'b-*')

index = 0;
for pt = 1:size(dCSarray , 1)
    if isnan(mean(dCSarray(pt - index , :) , 2))
        CSarray(pt - index , :) = [];
        CMarray(pt - index , :) = [];        
        dCSarray(pt - index , :) = [];
        dCMarray(pt - index , :) = [];
        index = index + 1;
    end
end

figure ; plot(1:Time - 1 , mean(dCSarray , 1)./abs(mean(dCSarray(:))) , 'b-o' , 1:Time - 1 , mean(dCMarray , 1)./abs(mean(dCMarray(:))) , 'r-*')
title('A^{-1}dA_{cell}/dt (blue) , M^{-1}dM_{cell}/dt (red) vs. Time')
figure ; plot(1:Time , mean(CSarray , 1)./abs(mean(CSarray(:))) , 'b-o' , 1:Time , mean(CMarray , 1)./abs(mean(CMarray(:))) , 'r-*')
title('A_{cell} (blue) , M_{cell} (red) vs. Time')

%% Average over time analysis

DeltaM = mean(dMpt , 2);
DeltaL = mean(dLpt , 2);

AllEdges = 1:Npath;
EdSet = AllEdges;

DeltaM = DeltaM(EdSet,:);
DeltaL = DeltaL(EdSet,:);

rndind = randperm(size(DeltaM , 1) , size(DeltaM , 1));
DeltaLr = DeltaL(rndind,:);

TimeAV_ML_X = corrcoef(DeltaM , DeltaL);    % Correlation of overall changes of different edges 
disp(TimeAV_ML_X)

TimeAV_ML_Xrand = corrcoef(DeltaM , DeltaLr);    % Correlation of overall changes of different edges 
disp(TimeAV_ML_Xrand)

figure ; histogram(DeltaL) ; hold on ; histogram(DeltaM)
figure ; histogram(DeltaLr) ; hold on ; histogram(DeltaM)

figure ; scatter(DeltaL(:) , DeltaM(:)) ; xlabel('d(log(r))/dt', 'FontSize' , 15) ; ylabel('d(log(m))/dt', 'FontSize' , 15)

%% constant M analysis

MctsInd = find(abs(DeltaM) < 0.05);

Mconst = Mpt(MctsInd , :);
dLconst = dLpt(MctsInd , :);

Mcst_dLX = corrcoef(mean(Mconst , 2) , mean(dLconst , 2));

disp(Mcst_dLX)

%% correlations over time
MLXt = zeros(size(dMpt , 2) , 1);

for tt = 1:size(dMpt , 2)
    
    mlt = corrcoef(dMpt(HH,tt) , dLpt(HH,tt));
    
    MLXt(tt , 1) = mlt(1,2);
    
end

figure ; plot(MLXt)

%%

alphaLMrnd = (DeltaM(:)./DeltaLr(:)).^(1);
alphaLM = (DeltaM(:)./DeltaL(:)).^(1);

alphaLMrndAccept = alphaLMrnd(alphaLMrnd < 10 & alphaLMrnd > -10);
alphaLMAccept = alphaLM(alphaLM < 10 & alphaLM > -10);

figure ; histogram(DeltaL./DeltaM(rndind)); title('dL/dM shuffled')
figure ; histogram(DeltaL./DeltaM); title('dL/dM')
figure ; histogram(alphaLMrndAccept); title('dM/dL shuffled')
figure ; histogram(alphaLMAccept); title('dM/dL')

[h , p] = kstest2(alphaLMAccept , alphaLMrndAccept);
disp('p-val');
disp(p)
disp('decision');
disp(h)


%% Frequency analysis of Linearly DETRENDED functions. 

Fs = 2; % sampling frequency (1/min)

MM = dMpt(AllCond,:);   
DetOrderM = 0; % order of detrending (1 = linear)
LL = dLpt(AllCond,:);   
DetOrderL = 0; % order of detrending (1 = linear)

DeltaT = size(MM , 2);
NumEds = size(MM , 1);

% detrended L , M to quadratic order
Mdt = zeros(NumEds , DeltaT);
Ldt = zeros(NumEds , DeltaT);

nfft = 2^(nextpow2(DeltaT) + 0);

for pt = 1:NumEds

    [pm , Sm , mum] = polyfit(1:DeltaT, MM(pt , :) , DetOrderM);
    [Mtrend  , ~] = polyval(pm , 1:DeltaT , Sm , mum);
    Mdt(pt , :) = MM(pt , :) - Mtrend;

    [pl , Sl , mul] = polyfit(1:DeltaT, LL(pt , :) , DetOrderL);
    [Ltrend  , ~] = polyval(pl , 1:DeltaT , Sl , mul);
    Ldt(pt , :) = LL(pt , :) - Ltrend;
       
end

Mft = fft(Mdt - mean(Mdt , 2) , nfft , 2); % Fast Fourier Transform
Mft = abs(Mft/(nfft)); % raw power spectrum density
Mft = 2*Mft(: , 1:1+nfft/2); % half-spectrum
f_scaleM = (0:nfft/2)* Fs/nfft; % frequency scale
figure ; plot(f_scaleM , mean(Mft)); title('Frequency distribution of Myosin')

Lft = fft(Ldt - mean(Ldt , 2) , nfft , 2); % Fast Fourier Transform
Lft = abs(Lft/(nfft)); % raw power spectrum density
Lft = 2*Lft(: , 1:1+nfft/2); % half-spectrum
f_scaleL = (0:nfft/2)* Fs/nfft; % frequency scale
figure ; plot(f_scaleL , mean(Lft)); title('Frequency distribution of Length')

figure ; errorbar(f_scaleL , mean(Lft)  , std(Lft), 'b-o'); hold on
errorbar(f_scaleM , mean(Mft) , std(Mft) , 'r-*' ); hold off
title('Frequency distribution of Length (blue) , and Myosin (red)');
Xlbl = sprintf('Freq (%2d/min)' , (Fs));
xlabel(Xlbl)

%%%%%%%%%%%% Spectrum of average length and average myosin %%%%%%%%%%%%%%%%
MftAV = fft(mean(Mdt , 1) - mean(Mdt(:)) , nfft , 2); % Fast Fourier Transform
MftAV = abs(MftAV/(nfft)); % raw power spectrum density
MftAV = 2*MftAV(1:1+nfft/2); % half-spectrum
[vm,km] = max(Mft); % find maximum
f_scaleM = (0:nfft/2)* Fs/nfft; % frequency scale
df_M = f_scaleM(km); % dominant frequency estimate
figure ; plot(f_scaleM , mean(Mft)); title('Frequency distribution of Myosin')

LftAV = fft(mean(Ldt , 1) - mean(Ldt(:)) , nfft , 2); % Fast Fourier Transform
LftAV = abs(LftAV/(nfft)); % raw power spectrum density
LftAV = 2*LftAV(1:1+nfft/2); % half-spectrum
[vl,kl] = max(Lft); % find maximum
f_scaleL = (0:nfft/2)* Fs/nfft; % frequency scale
df_L = f_scaleL(kl); % dominant frequency estimate
figure ; plot(f_scaleL , mean(Lft)); title('Frequency distribution of Length')

figure ; plot(f_scaleL , (LftAV) , 'b-o' , f_scaleM , (MftAV) , 'r-*' ); 
title('Frequency distribution of Average Length (blue) , and Average Myosin (red)');
Xlbl = sprintf('Freq (%2d/min)' , (Fs));
xlabel(Xlbl)

figure ; plot(1:DeltaT , Ldt(5,:),'r' , 1:DeltaT , Lpt(5,:) - mean(Lpt(5,:)) , 'b')

%%

Lrndind = randperm(NumEds , NumEds);
Mrndind = randperm(NumEds , NumEds);
Ldt_rnd = Ldt(Lrndind , :);
Mdt_rnd = Mdt(Mrndind , :);

MLX = zeros(NumEds , 1);
MLXrnd = zeros(NumEds , 1);
MLX_time = zeros(NumEds , 1);

for pt = 1:NumEds

    Ldtpt = Ldt_rnd(pt , :);
    Mdtpt = Mdt_rnd(pt , :);
    mlxrnd = corrcoef(Ldtpt , Mdtpt);
    MLXrnd(pt , 1) = mlxrnd(1,2);
    
    Ldtpt = Ldt(pt , :);
    Mdtpt = Mdt(pt , :);
    mlx = corrcoef(Ldtpt , Mdtpt);
    MLX(pt , 1) = mlx(1,2);
    
    mlxtime = corrcoef(Ldtpt , Mdtpt(randperm(DeltaT , DeltaT)));
    MLX_time(pt , 1) = mlxtime(1,2);
            
end

figure ; histogram(MLX , 21); title('PDF of ML Correlation')
figure ; histogram(MLXrnd , 21); title('PDF of shuffled ML Correlation')
figure ; histogram(MLX_time , 21); title('PDF of time-shuffled ML Correlation')
% 
figure ; histogram(MLX , 21 , 'Normalization' , 'cdf');
hold on
histogram(MLXrnd , 'Normalization' , 'cdf');
histogram(MLX_time , 'Normalization' , 'cdf');
title('Cumulative Histogram C_0 (blue) & C_{edge} (red) & C_{time} (yello)' , 'FontSize' , 10)
hold off

[hed , ped] = kstest2(MLX , MLXrnd);
fprintf('KS-test decision = %d \n' , hed);
fprintf('KS-test p-value = %d \n' , ped);

[htime , ptime] = kstest2(MLX , MLX_time);
fprintf('KS-test decision = %d \n' , htime);
fprintf('KS-test p-value = %d \n' , ptime);


%%  Distribution of \alpha = M_i(t)/R_i(t) for EDGE (i) and 
%   comparison with some normally distributed \alpha with the same std.

EdgeSet = AllCond;

Lfunc = Lpt(EdgeSet , :) ;
Mfunc = Mpt(EdgeSet , :) ;

TotInd = size(Lfunc , 1);

Lrand = 1:TotInd;
Mrand = 1:TotInd;

alpha0 = zeros(1,1); 
R20 = zeros(1,1);
alphaR20 = zeros(1,1);
Lff = zeros(1,1);
Mff = zeros(1,1);
naccept = 0;

RestLength0 = zeros(1,1);

for nrep = 1:1

    for Ind = 1:TotInd

        lind = Lrand(Ind);    
        mind = Mrand(Ind);

        LLgr = ((Lfunc(lind , 1:end)));
        MMgr = ((Mfunc(mind , 1:end)));

        Lii = squeeze(LLgr);
        Mii = squeeze(MMgr);

        Mii = Mii(randperm(Time , Time));
        Mii = (Mii(1:end-1) + Mii(2:end))/2;

        Lii = (Lii - mean(Lii))./ mean(Lii);
        Mii = (Mii - mean(Mii))./ mean(Mii);

        [po , S] = polyfit(Lii , Mii , 1);
        [Mfit , ~] = polyval(po , Lii , S);
        M_L0 = po(2);

        S_reg = sum((Mii - Mfit).^2);
        S_tot = sum((Mii - mean(Mfit)).^2);

        if 1 - (S_reg/S_tot) > 0 && -(M_L0) / po(1) > 0 %&& -(M_L0) / po(1) < min(Lii)  % & po(1) > 0
            naccept = naccept + 1;
            R20(naccept,1) = 1 - (S_reg/S_tot);
            alpha0(naccept,1) = po(1);
            alphaR20(naccept , 1) = po(1) .* (1 - (S_reg/S_tot));
            Lff((naccept-1)*length(Lii) + 1 : (naccept)*length(Lii)) = Lii(:);
            Mff((naccept-1)*length(Lii) + 1 : (naccept)*length(Lii)) = Mii(:);
            plot(Lii(:) , Mii(:) , 'o-'); hold on
            RestLength0(naccept , 1) = - (M_L0) / po(1);
        end
        
    end
    
end

hold off

[pf , Sf] = polyfit(Lff , Mff , 1);
[Mffit , ~] = polyval(pf , Lff , Sf);

Sf_reg = sum((Mff - Mffit).^2);
Sf_tot = sum((Mff - mean(Mffit)).^2);

figure ; plot(Lff , Mffit , 'bo' , Lff , Mff , 'r*') ; %xlim([5 , 20])

disp('overall R2');
disp(1 - (Sf_reg/Sf_tot));

disp('ratio accepted');
disp(naccept/(TotInd));

YY = normrnd(0 , std(alpha0) , [length(alpha0) , 1]);

[hchi2 , pchi2] = kstest2(alpha0 , YY);
disp('1sample KS-test decision')
disp(hchi2);
disp('1sample KS p-value')
disp(pchi2);

figure ; histogram(R20(:) , min(floor(TotInd/3) , 20));% , CtrsR2);
title('hist R^2_0' , 'FontSize' , 20)

figure ; scatter(alpha0 , R20); xlabel('\alpha_0' , 'FontSize' , 20); ylabel('R^2_0' , 'FontSize' , 20);

CtrsAlpha = min(alpha0) : (max(alpha0) - min(alpha0))/20 : max(alpha0);
figure ; histogram(alpha0(:) , min(floor(TotInd/3) , 20));% , CtrsAlpha);
title('Hist \alpha_0' , 'FontSize' , 20)

figure ; histogram(RestLength0 , min(floor(TotInd/3) , 20));
title('Hist of L^0_{rest}' , 'FontSize' , 20)

figure ; histogram(YY , min(floor(TotInd/3) , 100) , 'Normalization' , 'cdf');
hold on
histogram(alpha0 , min(floor(TotInd/3) , 100) , 'Normalization' , 'cdf');
title('Histogram of Null (blue) & \alpha_0 (red)' , 'FontSize' , 20)

%%  Distribution of \alpha = M_i(t)/R_i(t) for EDGES and
%   comparison with randomly chosen pairs (i,j) for M_i , R_j.

MaxNrep = 10;
Nperm = 1;

Lfunc = Lpt(EdgeSet , :);
Mfunc = Mpt(EdgeSet , :);

TotInd = size(Lfunc , 1);

alpha = zeros(1,1); 
R2 = zeros(1,1);
alphaR2 = zeros(1,1);
Lff = zeros(1,1);
Mff = zeros(1,1);
RestLength = zeros(1 , 1);

naccept = 0;

for nrep = 1:MaxNrep
    
    Lrand = randperm(TotInd , TotInd);
    Mrand = randperm(TotInd , TotInd);
    
    for nn = 1:Nperm
        Lrand = Lrand(randperm(TotInd , TotInd));
        Mrand = Mrand(randperm(TotInd , TotInd));
    end
    
    figure ; 
    for Ind = 1:TotInd

        lind = Lrand(Ind);    
        mind = Mrand(Ind);

        LLgr = ((Lfunc(lind , 1:end)));
        MMgr = ((Mfunc(mind , 1:end)));

        Lii = squeeze(LLgr);
        Mii = squeeze(MMgr);

        Mii = Mii(randperm(Time , Time));
        Mii = (Mii(1:end-1) + Mii(2:end))/2;

        Lii = (Lii - mean(Lii))./ mean(Lii);
        Mii = (Mii - mean(Mii))./ mean(Mii);

        [po , S] = polyfit(Lii , Mii , 1);
        [Mfit , ~] = polyval(po , Lii , S);
        M_L = po(2);

        S_reg = sum((Mii - Mfit).^2);
        S_tot = sum((Mii - mean(Mfit)).^2);

        if 1 - (S_reg/S_tot) > 0 && -(M_L) / po(1) > 0 %&& -(M_L) / po(1) < min(Lii)% & po(1) > 0
            naccept = naccept + 1;
            R2(naccept,1) = 1 - (S_reg/S_tot);
            alpha(naccept,1) = po(1);
            alphaR2(naccept , 1) = po(1) .* (1 - (S_reg/S_tot));
            Lff((naccept-1)*length(Lii) + 1 : (naccept)*length(Lii)) = Lii(:);
            Mff((naccept-1)*length(Lii) + 1 : (naccept)*length(Lii)) = Mii(:);
            plot(Lii(:) , Mii(:) , '*-'); hold on
            RestLength(naccept , 1) = -(M_L) / po(1);
        end
    end
end

disp('ratio accepted');
disp(naccept/(TotInd*MaxNrep))

[pf , Sf] = polyfit(Lff , Mff , 1);
[Mffit , ~] = polyval(pf , Lff , Sf);

Sf_reg = sum((Mff - Mffit).^2);
Sf_tot = sum((Mff - mean(Mffit)).^2);
disp(1 - (Sf_reg/Sf_tot))

HistCompare = cat(2 , alpha0 , alpha);
figure ; hist(HistCompare , 10); title('Histograms of \alpha_0 (blue) & \alpha_r (yellow)')

[hchi2 , pchi2] = kstest2(alpha , alpha0);
disp('2sample KS test decision')
disp(hchi2);
disp('2sample KS p-value')
disp(pchi2);
disp(mean(nE_XY(nE_XY ~= 0)));

figure ; histogram(R2(:) , min(floor(TotInd/3) , 20));% , CtrsR2);
title('hist R^2_r (randomized)' , 'FontSize' , 20)
CtrsAlpha = min(alphaR2) : (max(alphaR2) - min(alphaR2))/20 : max(alphaR2);
figure ; histogram(alpha(:) , min(floor(TotInd/3) , 20));% , CtrsAlpha);
title('hist \alpha_r (randomized)' , 'FontSize' , 20)

figure ; histogram(alpha0 , min(floor(TotInd/3) , 100) , 'Normalization' , 'cdf');
hold on
histogram(alpha , min(floor(TotInd/3) , 100) , 'Normalization' , 'cdf');
title('Cumulative Histogram \alpha_0 (blue) & \alpha_r (red)' , 'FontSize' , 20)
hold off

figure ; scatter(alpha , R2)
figure; scatter(Lff , Mff); xlabel('Length' , 'FontSize' , 20) ; ylabel('Myosin Concentration' , 'FontSize' , 20);
title('Scatter plot of Myosin vs. Length' , 'FontSize' , 20);

X = cat(2 , alpha , R2);
figure ; 
hist3(X,'CDataMode','auto','FaceColor','interp')
xlabel('\alpha')
ylabel('R^2')

R2mapG = imgaussfilt(R2map , 40);
alphaInv = 1 - rescale(abs(alphamap - mean(alpha0))./max(abs(alphamap - mean(alpha0))));
alphamapG = imgaussfilt(alphamap , 40);

figure ; imshow(R2mapG)
figure ; imshowpair(rescale(myo) , R2map)
figure ; imshow(rescale(alphamapG))
figure ; imshowpair(rescale(alphamapG) , R2mapG)

figure ; histogram(RestLength , min(floor(TotInd/3) , 100));
title('Hist of L^r_{rest}' , 'FontSize' , 20)


%%
for tt = 1:Time
    AA = zeros(Ny , Nx);
    for yy = 1:nY
        for xx = 1:nX
            AA(Tissue(tt).PatchBkw.Pts{yy}{xx}) = 1;    
        end
    end
    Asave = sprintf('patches40/patch_%02d.tif' , tt);    
    imwrite(double(AA) , Asave);
end

%%

for tt = 1:Time
    AA = zeros(Ny , Nx);
    for ed = 1:Npath
        if ismember(ed , AllCond)
            AA(Tissue(tt).Edge.PixelIdxList{EdPaths(ed , tt)}) = 1;
        end
    end
    Asave = sprintf('EdCond/Edges_%02d.tif' , tt);    
    imwrite(double(AA) , Asave);
end

%%
for tt = 1:Time
    AA = zeros(Ny , Nx);
    for yy = 1:nY
        for xx = 1:nX
            AA(Tissue(tt).PatchFrw.Pts{yy}{xx}) = 1;    
        end
    end
    Asave = sprintf('NewBoundFrwrd40/patch_%02d.tif' , tt);    
    imwrite(double(AA) , Asave);
end


for tt = 1:Time

    Sigma = zeros(2 , 2 , Ny , Nx);
    for ee = 1:Tissue(tt).Edge.NumObjects        
        eeAng = Tissue(tt).EdAng(ee);
        eeInd = Tissue(tt).Edge.PixelIdxList{ee}(floor(end/2));
        [eeY , eeX] = ind2sub([Ny , Nx] , eeInd);
        eeMyo = Tissue(tt).EdMyo(ee);
        Sigma(1 , 1 , eeY , eeX) = eeMyo * (cosd(eeAng))^2 ;
        Sigma(1 , 2 , eeY , eeX) = eeMyo * (cosd(eeAng)) * (sind(eeAng)) ;
        Sigma(2 , 2 , eeY , eeX) = eeMyo * (sind(eeAng))^2 ;
    end    
    Sigma(2,1,:,:) = Sigma(1,2,:,:);
    
end


maxStressyy = max(yyStress(:));
maxStressxx = max(xxStress(:));
maxStress = max(maxStressxx , maxStressyy);
for tt = 1:Time
    AA = zeros(Ny , Nx);
    for reg = 1:ind    
        AA(Tissue(tt).PatchFrw.Pts{NzPatches(reg , 1)}{NzPatches(reg , 2)}) = yyStress(NzPatches(reg , 1),NzPatches(reg , 2) , tt)/maxStress;
    end
    figure ; imshow(AA)
end



Epsilon = zeros(2 , 2 , Ny , Nx , Time);
Sigma = zeros(2 , 2 , Ny , Nx , Time);

for tt = 1:Time - 1

    Uy = Tissue(tt).VyFrwrd;    Uy(Dors_y+1:end , :) = - Uy(Dors_y+1:end , :);    
    Ux = Tissue(tt).VxFrwrd;

    [Exx,Exy] = gradient(Ux);
    [Eyx,Eyy] = gradient(Uy);

    Epsilon(1 , 1 , : , : , tt) = reshape(Exx , [1 , 1 , Ny , Nx , 1]);
    Epsilon(1 , 2 , : , : , tt) = reshape((Exy + Eyx)/2 , [1 , 1 , Ny , Nx , 1]);
    Epsilon(2 , 1 , : , : , tt) = reshape((Exy + Eyx)/2 , [1 , 1 , Ny , Nx , 1]);
    Epsilon(2 , 2 , : , : , tt) = reshape(Eyy , [1 , 1 , Ny , Nx , 1]);

end


EpsPtch = zeros(nY * nX , Time , 2 , 2);
indpt = 0;
for jj = 1:nY
    for ii = 1:nX
        indpt = indpt + 1;
        EpsPtch(indpt , : , 1 , 1) = mean(Epsilon)



    end
end

%%

BoxSize = 100;
[X2,Y2] = meshgrid(BoxSize/2 : BoxSize : size(im1,2) - BoxSize/2 , BoxSize/2 : BoxSize : size(im1,1) - BoxSize/2); 

for tt = 1:Time - 1

    Sigma1 = zeros(2 , 2 , Ny , Nx);
    Sigma2 = zeros(2 , 2 , Ny , Nx);

    Epsilon = zeros(2 , 2 , Ny , Nx);

    Uy = Tissue(tt).VyFrwrd;    Uy(Dors_y+1:end , :) = - Uy(Dors_y+1:end , :);    
    Ux = Tissue(tt).VxFrwrd;

    [Exx,Exy] = gradient(Ux);
    [Eyx,Eyy] = gradient(Uy);

    Epsilon(1 , 1 , : , :) = Exx;
    Epsilon(1 , 2 , : , :) = (Exy + Eyx)/2;
    Epsilon(2 , 1 , : , :) = (Exy + Eyx)/2;
    Epsilon(2 , 2 , : , :) = Eyy;

    for ee = 1:Tissue(tt).Edge.NumObjects        
        eeAng = Tissue(tt).EdAng(ee);
        eeInd = Tissue(tt).Edge.PixelIdxList{ee}(floor(end/2));
        [eeY , eeX] = ind2sub([Ny , Nx] , eeInd);
        eeMyo = Tissue(tt).EdMyo(ee);
        Sigma1(1 , 1 , eeY , eeX) = eeMyo * (cosd(eeAng))^2 ;
        Sigma1(1 , 2 , eeY , eeX) = eeMyo * (cosd(eeAng)) * (sind(eeAng)) ;
        Sigma1(2 , 2 , eeY , eeX) = eeMyo * (sind(eeAng))^2 ;
    end    
    Sigma1(2,1,:,:) = Sigma1(1,2,:,:);

    for ee = 1:Tissue(tt+1).Edge.NumObjects        
        eeAng = Tissue(tt+1).EdAng(ee);
        eeInd = Tissue(tt+1).Edge.PixelIdxList{ee}(floor(end/2));
        [eeY , eeX] = ind2sub([Ny , Nx] , eeInd);
        eeMyo = Tissue(tt+1).EdMyo(ee);
        Sigma2(1 , 1 , eeY , eeX) = eeMyo * (cosd(eeAng))^2 ;
        Sigma2(1 , 2 , eeY , eeX) = eeMyo * (cosd(eeAng)) * (sind(eeAng)) ;
        Sigma2(2 , 2 , eeY , eeX) = eeMyo * (sind(eeAng))^2 ;
    end    
    Sigma2(2,1,:,:) = Sigma2(1,2,:,:);

    Sigma_t = Sigma1 + Sigma2;

end

[UUy , UUx] = griddata(Uy , A);
