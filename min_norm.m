function [index] = min_norm(target_im, HOGs, LBPs, Hists, border)
   the_HOG = extractHOGFeatures(target_im);
   the_LBP = extractLBPFeatures(rgb2gray(target_im));
   the_Hist = (imhist(target_im)).';

   HOGs_dif = zeros(size(HOGs,1),1);
   LBPs_dif = zeros(size(LBPs,1),1);
   Hists_dif = zeros(size(Hists,1),1);
   for i = 1:size(HOGs,1)
      HOGs_dif(i) = sum(abs(HOGs(i,:) - the_HOG));
      LBPs_dif(i) = sum(abs(LBPs(i,:) - the_LBP));
      Hists_dif(i) = sum(abs(Hists(i,:) - the_Hist));
   end
   max_HOG = max(max(HOGs_dif));
   max_LBP = max(max(LBPs_dif));
   max_Hist = max(max(Hists_dif));
   ultra_max = max(max([max_Hist max_LBP max_HOG]));
   scale1 = ultra_max/max_HOG;
   scale2 = ultra_max/max_LBP;
   scale3 = ultra_max/max_Hist;
   all_of_em = zeros(size(HOGs,1),2);
%(scale1 * HOGs_dif(i) + scale2 * LBPs_dif(i) + Hists_dif(i) * scale3 )
   for i = 1:size(HOGs,1)
      all_of_em(i,1) = (scale1 * HOGs_dif(i) + scale2 * LBPs_dif(i) + Hists_dif(i) * scale3 );
      all_of_em(i,2) = i;
   end
   all_of_em = sortrows(all_of_em);

   if(size(all_of_em,1) < border)
       index = all_of_em(1:end,2);
   else
        index = all_of_em(1:border,2);
   end
end