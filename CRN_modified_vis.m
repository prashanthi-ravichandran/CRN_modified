%% CRN_modified trace visualization
% S. Tsutsui, July 2017
% Iyer Research Group, CUMC

vm = loaddbl('CRN_modified.out.vm');
aux = loaddbl('CRN_modified.out.aux');
%vm_baseline = loaddbl('CRN_modified_baseline.out.vm');
%aux_baseline = loaddbl('CRN_modified_baseline.out.aux');

Ina = aux(1:30:end);
Ik1 = aux(2:30:end);
Ito = aux(3:30:end);
Ikur = aux(4:30:end);
Ikr = aux(5:30:end);
Iks = aux(6:30:end);
Ical = aux(7:30:end);
Inak = aux(8:30:end);
Inaca = aux(9:30:end);
Ibca = aux(10:30:end);
Ibna = aux(11:30:end);
Ipca = aux(12:30:end);
Jrel = aux(13:30:end);
Jtr = aux(14:30:end);
Jup = aux(15:30:end);
Jxfer = aux(16:30:end);
Jupleak = aux(17:30:end);
Iion = aux(18:30:end);
Cai = aux(19:30:end);
CaNSR = aux(20:30:end);
CaSS = aux(21:30:end);
CaJSR = aux(22:30:end);
Cai_imw = aux(23:30:end);
CaNSR_imw = aux(24:30:end);
CaSS_imw = aux(25:30:end);
CaJSR_imw = aux(26:30:end);
Open = aux(27:30:end);
dC_tot = aux(28:30:end);
dCC_tot = aux(29:30:end);
dTOT = aux(30:30:end);

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
plot(Open,'LineWidth',1.4);
title('Open');
subplot(4,1,2);
plot(dC_tot,'LineWidth',1.4);
title('dC Total');
subplot(4,1,3);
plot(dCC_tot,'LineWidth',1.4);
title('dCCa Total');
subplot(4,1,4);
plot(dTOT,'LineWidth',1.4);
title('dC + dCCa Total');

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