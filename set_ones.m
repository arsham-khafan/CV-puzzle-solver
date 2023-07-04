function [GT, list_of_ones] = set_ones(gt)
list_of_ones = [];
for i = 2:size(gt,1) - 1
    for j = 2:size(gt,2) - 1
        if(gt(i - 1,j) == 2 && gt(i ,j - 1) == 2 && gt(i,j) == 0)
            gt(i,j) = 1;
            list_of_ones = [list_of_ones ; [i j]]; 
        end
    end
end
GT = gt;
end