filename_xls = 'D:\шноль\Пархомов\interval.xls:Лист1';
channel = ddeinit('Excel', filename_xls);
potash = ddereq(channel, 'r11c2:r6421c2');
uran = ddereq(channel, 'r11c6:r5063c6');


nbins = 100;% количество ячеек, на которое разбивается интервал экспериментальных значений (в данном случае интервалов между импульсами) при построении гистограммы

figure(1); %title ('Поташ');
doplots(potash, nbins);

figure(2); % title ('Уран');
doplots(uran, nbins);