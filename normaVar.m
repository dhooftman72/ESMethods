function TarVar = normaVar(TarVar)
x_range_perc_low = prctile(TarVar(:,1),2.5);
x_range_norm = TarVar(:,1) - x_range_perc_low;
x_range_norm(x_range_norm<0) = 0;
x_range_perc = prctile(x_range_norm,97.5);
x_range = (x_range_norm./x_range_perc);
x_range(x_range>1) = 1;
TarVar = x_range;
clear x_range*
end