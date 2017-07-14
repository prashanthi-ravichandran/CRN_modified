%% CRN_modified trace visualization
% S. Tsutsui, July 2017
% Iyer Research Group, CUMC

vm = loaddbl('CRN_modified.out.vm');
aux = loaddbl('CRN_modified.out.aux');
%vm_baseline = loaddbl('CRN_modified_baseline.out.vm');
%aux_baseline = loaddbl('CRN_modified_baseline.out.aux');

Ina = aux(1:26:end);
Ik1 = aux(2:26:end);
Ito = aux(3:26:end);
Ikur = aux(4:26:end);
Ikr = aux(5:26:end);
Iks = aux(6:26:end);
Ical = aux(7:26:end);
Inak = aux(8:26:end);
Inaca = aux(9:26:end);
Ibca = aux(10:26:end);
Ibna = aux(11:26:end);
Ipca = aux(12:26:end);
Jrel = aux(13:26:end);
Jtr = aux(14:26:end);
Jup = aux(15:26:end);
Jxfer = aux(16:26:end);
Jupleak = aux(17:26:end);
Iion = aux(18:26:end);
Cai = aux(19:26:end);
CaNSR = aux(20:26:end);
CaSS = aux(21:26:end);
CaJSR = aux(22:26:end);
Cai_imw = aux(23:26:end);
CaNSR_imw = aux(24:26:end);
CaSS_imw = aux(25:26:end);
CaJSR_imw = aux(26:26:end);

figure;
subplot(4,1,1);
plot(vm,'LineWidth',1.4);
title('V_m');
subplot(4,1,2);
plot(Ical,'LineWidth',1.4);
title('I_{Ca_L}');
subplot(4,1,3);
plot(Jrel,'LineWidth',1.4);
title('J_{rel}');
subplot(4,1,4);
plot(Jup,'LineWidth',1.4);
title('J_{up}');

figure;
subplot(4,1,1);
plot(Cai,'LineWidth',1.4);
title('Ca_i');
subplot(4,1,2);
plot(CaSS,'LineWidth',1.4);
title('Ca_{SS}');
subplot(4,1,3);
plot(CaNSR,'LineWidth',1.4);
title('Ca_{NSR}');
subplot(4,1,4);
plot(CaJSR,'LineWidth',1.4);
title('Ca_{JSR}');

figure;
subplot(4,1,1);
plot(Cai_imw,'LineWidth',1.4);
title('Ca_i IMW');
subplot(4,1,2);
plot(CaSS_imw,'LineWidth',1.4);
title('Ca_{SS} IMW');
subplot(4,1,3);
plot(CaNSR_imw,'LineWidth',1.4);
title('Ca_{NSR} IMW');
subplot(4,1,4);
plot(CaJSR_imw,'LineWidth',1.4);
title('Ca_{JSR} IMW');