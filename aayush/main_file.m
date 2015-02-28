
%% Gaussian-filter
% myfilter = fspecial('gaussian',[5,1], 1);
% filter_data(:,1) = data_75(:,1);
% filter_data(:,2) = imfilter(data_75(:,2), myfilter);
% filter_data(:,3) = imfilter(data_75(:,3), myfilter);

gaussFilter = gausswin(50);
gaussFilter = gaussFilter/sum(gaussFilter);
filter_data(:,1) = data_75(:,1);
filter_data(:,2) = conv(data_75(:,2),gaussFilter,'same');
filter_data(:,3) = conv(data_75(:,3),gaussFilter,'same');

%% Poly-fit
%  for i = 40:59960
%  fit(i,:) = polyfit(filter_data(i-39:i+40,1),filter_data(i-39:i+40,2),1);
%  filter_data(i,4) = fit(i,1);
%  fit(i,:) = polyfit(filter_data(i-39:i+40,1),filter_data(i-39:i+40,3),1);
%  filter_data(i,5) = fit(i,1);
%  end

%% main-code
for i=1:59
% i=59;
    [p_min,ind_min] = min(filter_data((i-1)*1000+1:i*1000,2));
    [start,fin] = check(filter_data((i-1)*1000+1:i*1000,4), ind_min);
    results(i,1) = start;
    results(i,2) = fin;
    [p_min,ind_min] = min(filter_data((i-1)*1000+1:i*1000,3));
    [start,fin] = check(filter_data((i-1)*1000+1:i*1000,5), ind_min);
    results(i,3) = start;
    results(i,4) = fin;
end