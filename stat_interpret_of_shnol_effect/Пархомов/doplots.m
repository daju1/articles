function doplots(intervals, nbins)

wbin = (max(intervals) - min (intervals)) / nbins;

subplot(4,1,1)
plot(intervals);title('��������� ����� ����������');

subplot(4,1,2)
hist(intervals, nbins);title(sprintf('����������� ���������� ����� ����������. ����� ����� ��������� %d ������ ������ %f ���', nbins, wbin));
[h,xv] = hist(intervals, nbins);

subplot(4,1,3)
plot(xv, log(h));title('����������� ���������� ����� ���������� � ������������������� ��������');

subplot(4,1,4)
hist(log(intervals), nbins); title('����������� ���������� ���������� ����� ����������');
