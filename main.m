clc; clear; close all;

DELAY = 0;
%------------------------------------ change these attributes only!
puzzle_num = 1;
pieces_count = 160;
%------------------------------------------------------------------

if(puzzle_num == 1 || puzzle_num == 2)
    
    path = ['Puzzle_' num2str(puzzle_num) '_1200_1920\Unrotated' '_' num2str(pieces_count)];
    n = uint16((0.5)^(log4(pieces_count/10) - 1) * 240);
    P = pieces_count;
    size1 = uint16(2^(log4(pieces_count/10) - 1) * 5);
    size2 = uint16(2^(log4(pieces_count/10) - 1) * 8);
    border = pieces_count;
    path_C_a1 = '\Corner_1_1.tif';
    path_C_a2 = ['\Corner_1_' num2str(size2) '.tif'];
    path_C_a3 = ['\Corner_' num2str(size1) '_1.tif'];
    path_C_a4 = ['\Corner_' num2str(size1) '_' num2str(size2) '.tif'];
else
    path = ['Puzzle_' num2str(puzzle_num) '_600_960\Unrotated' '_' num2str(pieces_count)];
    n = uint16((0.5)^(log4(pieces_count/10) - 1) * 120);
    P = pieces_count;
    size1 = uint16(2^(log4(pieces_count/10) - 1) * 5);
    size2 = uint16(2^(log4(pieces_count/10) - 1) * 8);
    border = pieces_count;
    path_C_a1 = '\Corner_1_1.tif';
    path_C_a2 = ['\Corner_1_' num2str(size2) '.tif'];
    path_C_a3 = ['\Corner_' num2str(size1) '_1.tif'];
    path_C_a4 = ['\Corner_' num2str(size1) '_' num2str(size2) '.tif'];
end
% path = 'Puzzle_1_1200_1920\Unrotated_640';

% n = 60;
% P = 640;
% size1 = 20;
% size2 = 32;
% border = 160;
% marg = 240;

directo = dir(path); 
patches_arr = zeros(n,n,3,P - 4);
big_pic = im2double(imread([path '\Output.tif']));
corner_arr = zeros(n,n,3,4);

corner_arr(:,:,:,1) = im2double(imread([path path_C_a1]));
corner_arr(:,:,:,2) = im2double(imread([path path_C_a2]));
corner_arr(:,:,:,3) = im2double(imread([path path_C_a3]));
corner_arr(:,:,:,4) = im2double(imread([path path_C_a4]));

for i = 9:4 + P
    patches_arr(:,:,:,i - 8) = im2double(imread([path '\' directo(i).name]));
end
HOGs = zeros(P - 4,size(extractHOGFeatures(patches_arr(:,:,:,1)),2));
LBPs = zeros(P - 4,size(extractLBPFeatures(rgb2gray(patches_arr(:,:,:,1))),2));
Hists = zeros(P - 4,size(imhist(patches_arr(:,:,:,1)),1));

% UHOGs = zeros(P - 4,size(extractHOGFeatures(patches_arr(1:marg,:,:,1)),2));
% ULBPs = zeros(P - 4,size(extractLBPFeatures(rgb2gray(patches_arr(1:marg,:,:,1))),2));
% UHists = zeros(P - 4,size(imhist(patches_arr(1:marg,:,:,1)),1));
% 
% RHOGs = zeros(P - 4,size(extractHOGFeatures(patches_arr(:,end - marg + 1:end,:,1)),2));
% RLBPs = zeros(P - 4,size(extractLBPFeatures(rgb2gray(patches_arr(:,end - marg + 1:end,:,1))),2));
% RHists = zeros(P - 4,size(imhist(patches_arr(:,end - marg + 1:end,:,1)),1));
% 
% DHOGs = zeros(P - 4,size(extractHOGFeatures(patches_arr(end - marg + 1:end,:,:,1)),2));
% DLBPs = zeros(P - 4,size(extractLBPFeatures(rgb2gray(patches_arr(end - marg + 1:end,:,:,1))),2));
% DHists = zeros(P - 4,size(imhist(patches_arr(end - marg + 1:end,:,:,1)),1));
% 
% LHOGs = zeros(P - 4,size(extractHOGFeatures(patches_arr(1:marg,:,:,1)),2));
% LLBPs = zeros(P - 4,size(extractLBPFeatures(rgb2gray(patches_arr(1:marg,:,:,1))),2));
% LHists = zeros(P - 4,size(imhist(patches_arr(1:marg,:,:,1)),1));
for i = 1:P - 4
    HOGs(i,:) = extractHOGFeatures(patches_arr(:,:,:,i));
    LBPs(i,:) = extractLBPFeatures(rgb2gray(patches_arr(:,:,:,i)));
    Hists(i,:) = imhist(patches_arr(:,:,:,i));

%     UHOGs(i,:) = extractHOGFeatures(patches_arr(1:marg,:,:,i));
%     ULBPs(i,:) = extractLBPFeatures(rgb2gray(patches_arr(1:marg,:,:,i)));
%     UHists(i,:) = imhist(patches_arr(1:marg,:,:,i));
% 
%     RHOGs(i,:) = extractHOGFeatures(patches_arr(:,end - marg + 1:end,:,i));
%     RLBPs(i,:) = extractLBPFeatures(rgb2gray(patches_arr(:,end - marg + 1:end,:,i)));
%     RHists(i,:) = imhist(patches_arr(:,end - marg + 1:end,:,i));
% 
%     DHOGs(i,:) = extractHOGFeatures(patches_arr(end - marg + 1:end,:,:,i));
%     DLBPs(i,:) = extractLBPFeatures(rgb2gray(patches_arr(end - marg + 1:end,:,:,i)));
%     DHists(i,:) = imhist(patches_arr(end - marg + 1:end,:,:,i));
% 
%     LHOGs(i,:) = extractHOGFeatures(patches_arr(1:marg,:,:,i));
%     LLBPs(i,:) = extractLBPFeatures(rgb2gray(patches_arr(1:marg,:,:,i)));
%     LHists(i,:) = imhist(patches_arr(1:marg,:,:,i));

end

% ((i - 1)*n + 1: i * n, (j - 1)*n + 1: j * n)
%first row filling

for i = 2:size2 - 1
    i_prev = i - 1;
    target_im = big_pic(1:n,(i_prev - 1)*n + 1: i_prev * n,:);
%     target_im = target_im(:,end - marg + 1:end,:);
    index = min_norm(target_im,HOGs,LBPs,Hists,border);
    target_rc = target_im(:,end,:);
    min = inf;

    for k = 1:size(index,1)
        if(mse(patches_arr(:,1,:,index(k)),target_rc) ...
                < min)
            min = mse(patches_arr(:,1,:,index(k)),target_rc);
            final_choice = patches_arr(:,:,:,index(k));
            final_k = index(k);
        end
    end
    big_pic(1:n , ((i - 1)*n + 1): (i * n) ,:) = final_choice;
    imshow(big_pic,[])
    pause(DELAY);
    patches_arr(:,:,:,final_k) = [];
    HOGs(final_k,:) = [];
    LBPs(final_k,:) = [];
    Hists(final_k,:) = [];

%     UHOGs(index,:) = [];
%     ULBPs(index,:) = [];
%     UHists(index,:) = [];
% 
%     DHOGs(index,:) = [];
%     DLBPs(index,:) = [];
%     DHists(index,:) = [];
% 
%     RHOGs(index,:) = [];
%     RLBPs(index,:) = [];
%     RHists(index,:) = [];
% 
%     LHOGs(index,:) = [];
%     LLBPs(index,:) = [];
%     LHists(index,:) = [];
end
%imshow(big_pic,[])
for i = 2:size2 - 1
    i_prev = i - 1;
    target_im = big_pic((size1 - 1)*n + 1: size1 * n,(i_prev - 1)*n + 1: i_prev * n,:);
%     target_im = target_im(:,end - marg + 1:end,:);
    index = min_norm(target_im,HOGs,LBPs,Hists,border);
    target_rc = target_im(:,end,:);
    min = inf;
    for k = 1:size(index,1)
        if(mse(patches_arr(:,1,:,index(k)),target_rc) < min)
            min = mse(patches_arr(:,1,:,index(k)),target_rc);
            final_choice = patches_arr(:,:,:,index(k));
            final_k = index(k);
        end
    end

    big_pic((size1 - 1)*n + 1:(size1) * n,(i - 1)*n + 1: i * n,:) = patches_arr(:,:,:,final_k);
    imshow(big_pic,[])
    pause(DELAY);
    patches_arr(:,:,:,final_k) = [];
    HOGs(final_k,:) = [];
    LBPs(final_k,:) = [];
    Hists(final_k,:) = [];

%     UHOGs(index,:) = [];
%     ULBPs(index,:) = [];
%     UHists(index,:) = [];
% 
%     DHOGs(index,:) = [];
%     DLBPs(index,:) = [];
%     DHists(index,:) = [];
% 
%     RHOGs(index,:) = [];
%     RLBPs(index,:) = [];
%     RHists(index,:) = [];
% 
%     LHOGs(index,:) = [];
%     LLBPs(index,:) = [];
%     LHists(index,:) = [];
end

for i = 2:size1 - 1

    i_prev = i - 1;
    target_im = big_pic((i_prev - 1)*n + 1: i_prev * n,1:n,:);
%     target_im = target_im(end - marg + 1:end,:,:);
    index = min_norm(target_im,HOGs,LBPs,Hists,border);

    target_rc = target_im(end,:,:);
    min = inf;
    for k = 1:size(index,1)
        if(mse(patches_arr(1,:,:,index(k)),target_rc) < min)
            min = mse(patches_arr(1,:,:,index(k)),target_rc);
            final_choice = patches_arr(:,:,:,index(k));
            final_k = index(k);
        end
    end
    

    big_pic((i - 1)*n + 1: i * n,1:n,:) = patches_arr(:,:,:,final_k);
    imshow(big_pic,[])
    pause(DELAY);
    patches_arr(:,:,:,final_k) = [];
    HOGs(final_k,:) = [];
    LBPs(final_k,:) = [];
    Hists(final_k,:) = [];
%     UHOGs(index,:) = [];
%     ULBPs(index,:) = [];
%     UHists(index,:) = [];
% 
%     DHOGs(index,:) = [];
%     DLBPs(index,:) = [];
%     DHists(index,:) = [];
% 
%     RHOGs(index,:) = [];
%     RLBPs(index,:) = [];
%     RHists(index,:) = [];
% 
%     LHOGs(index,:) = [];
%     LLBPs(index,:) = [];
%     LHists(index,:) = [];
end

for i = 2:size1 - 1

    i_prev = i - 1;
    target_im = big_pic((i_prev - 1)*n + 1: i_prev * n,(size2 - 1)*n + 1:size2 *n,:);
%     target_im = target_im(end - marg + 1:end,:,:);
    index = min_norm(target_im,HOGs,LBPs,Hists,border);

    target_rc = target_im(end,:,:);
    min = inf;
    for k = 1:size(index,1)
        if(mse(patches_arr(1,:,:,index(k)),target_rc) < min)
            min = mse(patches_arr(1,:,:,index(k)),target_rc);
            final_choice = patches_arr(:,:,:,index(k));
            final_k = index(k);
        end
    end

    big_pic((i - 1)*n + 1: i * n,(size2 - 1)*n + 1:size2 *n,:) = patches_arr(:,:,:,final_k);
    imshow(big_pic,[])
    pause(DELAY);
    patches_arr(:,:,:,final_k) = [];
    HOGs(final_k,:) = [];
    LBPs(final_k,:) = [];
    Hists(final_k,:) = [];
%     UHOGs(index,:) = [];
%     ULBPs(index,:) = [];
%     UHists(index,:) = [];
% 
%     DHOGs(index,:) = [];
%     DLBPs(index,:) = [];
%     DHists(index,:) = [];
% 
%     RHOGs(index,:) = [];
%     RLBPs(index,:) = [];
%     RHists(index,:) = [];
% 
%     LHOGs(index,:) = [];
%     LLBPs(index,:) = [];
%     LHists(index,:) = [];
end
% ((i - 1)*n + 1: i * n, (j - 1)*n + 1: j * n)

guid_table = zeros(size1,size2,'uint8');
guid_table(:,1) = 2;
guid_table(1,:) = 2;
guid_table(:,end) = 2;
guid_table(end,:) = 2;
[guid_table, list_of_ones] = set_ones(guid_table);
while (size(list_of_ones,1) ~= 0)
    for l = 1:size(list_of_ones,1)
        i = list_of_ones(l,1);
        j = list_of_ones(l,2);
        
        i_prev = i - 1;
        target_im = big_pic((i_prev - 1)*n + 1: i_prev * n, (j - 1)*n + 1: j * n,:);
    %         target_im = target_im(end - marg + 1:end,:,:);
        index = min_norm(target_im,HOGs,LBPs,Hists,border);
        target_rc = target_im(end,:,:);
        target_rc2 = big_pic((i - 1)*n + 1: i * n,(j - 1) * n,:);
        min = inf;
    
        for k = 1:size(index,1)
            if(mse(patches_arr(1,:,:,index(k)),target_rc) ...
                    + mse(patches_arr(:,1,:,index(k)),target_rc2)... 
                    < min)
                min = mse(patches_arr(1,:,:,index(k)),target_rc);
                final_choice = patches_arr(:,:,:,index(k));
                final_k = index(k);
            end
        end
        
        big_pic((i - 1)*n + 1: i * n, (j - 1)*n + 1: j * n,:) = patches_arr(:,:,:,final_k);
        imshow(big_pic,[])
        pause(DELAY);
        patches_arr(:,:,:,final_k) = [];
        HOGs(final_k,:) = [];
        LBPs(final_k,:) = [];
        Hists(final_k,:) = [];
    
        guid_table(i,j) = 2;
        
    end
    [guid_table, list_of_ones] = set_ones(guid_table);
end










% for i = 2 : size1 - 1
%     for j = 2: size2 - 1
%         i_prev = i - 1;
%         target_im = big_pic((i_prev - 1)*n + 1: i_prev * n, (j - 1)*n + 1: j * n,:);
% %         target_im = target_im(end - marg + 1:end,:,:);
%         index = min_norm(target_im,HOGs,LBPs,Hists,border);
%         target_rc = target_im(end,:,:);
%         target_rc2 = big_pic((i - 1)*n + 1: i * n,(j - 1) * n,:);
%         min = inf;
% 
%         for k = 1:size(index,1)
%             if(mse(patches_arr(1,:,:,index(k)),target_rc) ...
%                     + mse(patches_arr(:,1,:,index(k)),target_rc2)... 
%                     < min)
%                 min = mse(patches_arr(1,:,:,index(k)),target_rc);
%                 final_choice = patches_arr(:,:,:,index(k));
%                 final_k = index(k);
%             end
%         end
%         
%         big_pic((i - 1)*n + 1: i * n, (j - 1)*n + 1: j * n,:) = patches_arr(:,:,:,final_k);
%         imshow(big_pic,[])
%         pause(DELAY);
%         patches_arr(:,:,:,final_k) = [];
%         HOGs(final_k,:) = [];
%         LBPs(final_k,:) = [];
%         Hists(final_k,:) = [];
% %         UHOGs(index,:) = [];
% %         ULBPs(index,:) = [];
% %         UHists(index,:) = [];
% %     
% %         DHOGs(index,:) = [];
% %         DLBPs(index,:) = [];
% %         DHists(index,:) = [];
% %     
% %         RHOGs(index,:) = [];
% %         RLBPs(index,:) = [];
% %         RHists(index,:) = [];
% %     
% %         LHOGs(index,:) = [];
% %         LLBPs(index,:) = [];
% %         LHists(index,:) = [];
% 
%     end
% end

imshow(big_pic, [])