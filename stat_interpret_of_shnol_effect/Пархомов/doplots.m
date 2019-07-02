function doplots(intervals, nbins)

wbin = (max(intervals) - min (intervals)) / nbins;

subplot(4,1,1)
plot(intervals);title('интервалы между импульсами');

subplot(4,1,2)
hist(intervals, nbins);title(sprintf('гистограмма интервалов между импульсами. „исло €чеек разбиени€ %d ширина €чейки %f сек', nbins, wbin));
[h,xv] = hist(intervals, nbins);

subplot(4,1,3)
plot(xv, log(h));title('гистограмма интервалов между импульсами в полулогарифмическом масштабе');

subplot(4,1,4)
hist(log(intervals), nbins); title('гистограмма логарифмов интервалов между импульсами');
