filename_xls = 'D:\�����\��������\interval.xls:����1';
channel = ddeinit('Excel', filename_xls);
potash = ddereq(channel, 'r11c2:r6421c2');
uran = ddereq(channel, 'r11c6:r5063c6');


nbins = 100;% ���������� �����, �� ������� ����������� �������� ����������������� �������� (� ������ ������ ���������� ����� ����������) ��� ���������� �����������

figure(1); %title ('�����');
doplots(potash, nbins);

figure(2); % title ('����');
doplots(uran, nbins);