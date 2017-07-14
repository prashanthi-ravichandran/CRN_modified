%% CRN_modified trace visualization
% S. Tsutsui, July 2017
% Iyer Research Group, CUMC

vm = loaddbl('CRN_modified.out.vm');
aux = loaddbl('CRN_modified.out.aux');
%vm_baseline = loaddbl('CRN_modified_baseline.out.vm');
%aux_baseline = loaddbl('CRN_modified_baseline.out.aux');

Ina = aux(1:27:end);
Ik1 = aux(2:27:end);
Ito = aux(3:27:end);
Ikur = aux(4:27:end);
Ikr = aux(5:27:end);
Iks = aux(6:27:end);
Ical = aux(7:27:end);
Inak = aux(8:27:end);
Inaca = aux(9:27:end);
Ibca = aux(10:27:end);
Ibna = aux(11:27:end);
Ipca = aux(12:27:end);
Jrel = aux(13:27:end);
Jtr = aux(14:27:end);
Jup = aux(15:27:end);
Jxfer = aux(16:27:end);
Jupleak = aux(17:27:end);
Iion = aux(18:27:end);
Cai = aux(19:27:end);
CaNSR = aux(20:27:end);
CaSS = aux(21:27:end);
CaJSR = aux(22:27:end);
Cai_imw = aux(23:27:end);
CaNSR_imw = aux(24:27:end);
CaSS_imw = aux(25:27:end);
CaJSR_imw = aux(26:27:end);
Open = aux(27:27:end);

figure;
subplot(4,1,1);
plot(vm,'LineWidth',1.4);
title('V_m');
subplot(4,1,2);
plot(Ical,'LineWidth',1.4);
title('I_{Ca_L}');
subplot(4,1,3);
plot(Open,'LineWidth',1.4);
title('Open');
subplot(4,1,4);
plot(Jrel,'LineWidth',1.4);
title('J_{rel}');

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