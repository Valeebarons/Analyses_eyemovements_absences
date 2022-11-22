   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Imports single patient TRC file
  % filter data and extract epochs with absence seizures (previously detected visually)
  % 
  % Estimates:
  % - averaged maximum amplitude for each seizure for all, anterior and
  % posterior channels. 
  % - averaged maximum amplitude of all seizures together for all, anterior
  % and posterior channels
  % - peak frequency for each seizure and development of peak frequency
  % through time during absences
  % - slope of frequency through time for each absences for each patient
  % and averaged slope for all seizures of each patient
  % - centre of gravity of two groups
  %
  %
  %  
  % needs readalltrcdata.m
  % 
  % VBarone Nov, 2022
  % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% read TRC file
clear all; close all

datapath = 'C:\Users\vb\Documents\MATLAB\AbsenceChildren_paper4\data\Kemp\EEG\subj15\';
filename = 'subject15_bikemonkeymonkey.TRC';

[header, data, trigger,electrode] = readalltrcdata(filename, datapath); 
fs = data.SampleRate(1,1);
t = (0:length(data.Micromed)-1)/ fs; %in seconds
data.Micromed = data.Micromed.*1000000; %in microVolts
tminutes = t./60;
%plot(tminutes, data.Micromed(15,:))

%% filter and keep only relevant channels

EEG_data = data.Micromed';
[B,A] = butter(3,[1 30]./(fs/2), 'bandpass');         
EEG_filtered = filtfilt(B,A,EEG_data);


%remove superfluos channels 
EEG_filtered(:,[2, 20,22:36]) = []; 
%left chan
signal_names_KH =  {'Fp1','Fp2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2'};
%T7 T8 = T3 T4; P8 P7 = T6 T5

%plot one absence (pat 3)
figure(5)
plot(t(1, 1357*fs :1365*fs), EEG_filtered((   1357*fs :1365*fs),3))
xlabel('Time (s)'), ylabel('Amplitude (\muV)')
title('F7')

%% add seizure limits defined visually (-1 s from start and +1 s to end).

pat3 = [507*fs 516.5*fs; 841*fs 850*fs];
pat5 = [350.5*fs 369*fs; 901*fs 922*fs; 1326*fs 1362*fs; 1895*fs 1918*fs];
pat6 = [152*fs 160*fs; 1996*fs 2005*fs];
pat7 = [273*fs 284*fs; 572.5*fs 580*fs; 826.5 *fs 836*fs; 995*fs 1006*fs; 1082*fs 1093*fs;1356*fs 1361*fs; 1475*fs 1484*fs; 1726.5*fs 1738*fs; 1860*fs 1870*fs; 2170.5*fs 2178*fs; 2335*fs 2347*fs]; 
pat8 = [ 486*fs 500*fs;589.5*fs 607*fs; 975.5*fs 991*fs; 1118*fs 1126*fs;1268.5*fs 1292*fs;1511*fs 1530*fs;1790*fs 1827*fs; 1917*fs 1945*fs];
pat10 = [111*fs 131*fs; 189*fs 199*fs;1256*fs 1262*fs; 1292*fs 1305*fs; 1389*fs 1396*fs; 1489*fs 1499*fs;1842*fs 1862*fs; 1865*fs 1872*fs; 1926*fs 1932*fs; 2062*fs 2068*fs; 2133*fs 2139*fs];
pat11 = [67*fs 73*fs; 80*fs 86*fs; 87*fs 99*fs; 108*fs 113*fs;156*fs 165*fs; 224*fs 229*fs; 389*fs 394*fs; 651*fs 657*fs; 773*fs 779*fs; 913*fs 924*fs; 996*fs 1001*fs; 1045*fs 1050*fs; 1081*fs 1087*fs; 1848*fs 1853*fs; 1939*fs 1944*fs; 1945*fs 1950*fs; 1955*fs 1960*fs; 2016*fs 2020*fs; 2067*fs 2071*fs; 2073*fs 2077*fs; 2257*fs 2261*fs; 2755*fs 2759*fs; 2904*fs 2908*fs];
pat12 = [38*fs 46*fs; 144*fs 158*fs; 267*fs 279*fs; 445*fs 456*fs; 524.5*fs 535*fs; 855*fs 866*fs; 1249*fs 1261*fs; 1646*fs 1655*fs; 1777*fs 1788*fs; 1831*fs 1840*fs; 1975*fs 1981*fs; 2142*fs 2155*fs; 2416*fs 2427*fs; 2499*fs 2514*fs];
pat13= [178*fs 185*fs; 1574*fs 1579*fs; 1648*fs 1658*fs];
pat14 = [120*fs 128*fs];
pat15= [575*fs 580*fs; 654*fs 659*fs;705*fs 711.5*fs;1161*fs 1166*fs; 1225.6*fs 1231.6*fs; 1297*fs 1302*fs; 1357*fs 1365*fs];
    
duration = mean(pat3(:,2)-pat3(:,1))/fs;


%% amplitude 
max_amp1 = [];
for ii = 1:size(pat15,1)
    for i=1:19
        max_amp1(i,ii) = max(EEG_filtered((pat15(ii,1):pat15(ii,2)),i)); %max amplitude per channel per seizures
    end
   
end
max_amp_perseiz = mean(max_amp1,1);
mapx_Aamp_perzeis = mean(max_amp1(3:7,:),1);
mapx_Pamp_perzeis = mean(max_amp1([13,17:19],:),1);


max_amp = mean(max_amp1,2);
anterior_amp = mean(max_amp(3:7,:));
posterior_amp =  mean(max_amp([13,17:19],:));

indamp = (anterior_amp-posterior_amp)/(anterior_amp+posterior_amp);

%% frequency estimate 
Nfft=4*fs;    % Nfft, Window en Noverlap zijn instellingen van de functie pwelch. Deze kun je aanpassen om het spectrum 'gladder' of minder 'glad' te maken.
Window=3*fs; %3cycles/4Hz(lower freq i am interested in)= 0.75. 
Noverlap=fs;

%get deviation from the mean, as a measure for amplitude? 
O1=EEG_filtered(:,18); O1=O1-mean(O1);
O2=EEG_filtered(:,19); O2=O2-mean(O2);
Cz=EEG_filtered(:,10); Cz=Cz-mean(Cz);
C3=EEG_filtered(:,9); C3=C3-mean(C3);
C4=EEG_filtered(:,11); C4=C4-mean(C4);
Pz=EEG_filtered(:,15); Pz=Pz-mean(Pz);
P3=EEG_filtered(:,14); P3=P3-mean(P3); 
P4=EEG_filtered(:,16); P4=P4-mean(P4); 
Fz=EEG_filtered(:,5); Fz=Fz-mean(Fz);
F3=EEG_filtered(:,4); F3=F3-mean(F3);
F4=EEG_filtered(:,6); F4=F4-mean(F4);
F8=EEG_filtered(:,7); F8=F8-mean(F8);
F7=EEG_filtered(:,3); F7=F7-mean(F7);
Fp1=EEG_filtered(:,1); Fp1=Fp1-mean(Fp1);
Fp2=EEG_filtered(:,2); Fp2=Fp2-mean(Fp2);
T3=EEG_filtered(:,8); T3=T3-mean(T3);
T4=EEG_filtered(:,12); T4=T4-mean(T4);
T5=EEG_filtered(:,13); T5=T5-mean(T5);
T6=EEG_filtered(:,17); T6=T6-mean(T6);
EEG_power = [Fp2-F8, F8-T4, T4-T6, T6-O2, Fp1-F7, F7-T3, T3-T5, T5-O1,Fp2-F4, F4-C4, C4-P4, P4-O2, Fp1-F3,F3-C3, C3-P3, P3-O1, Fz-Cz, Cz-Pz];

[Pxx,F]=pwelch(EEG_power,Window,Noverlap,Nfft,fs);
 

%% peak frequency during absences 

SZ=struct;
for  k=1:size(pat15,1) 
    seizure1 = EEG_filtered((pat15(k,1):pat15(k,2)),:);
    %seizure1 = seizure - mean(seizure); %remove mean to make it less noisy 
    [SZ(k).powerabsences, F] = pwelch(seizure1,Window,Noverlap,Nfft,fs);
    [~,loc] = max(SZ(k).powerabsences);
    SZ(k).freqestimate= mean(F(loc));
end 


%average power for all seizures
%if more than 2 absences    
power_absences= mean(cat(size(SZ, 2),SZ.powerabsences),size(SZ, 2));

%if absences are 2 
for i=1:size(SZ,2)
    if i ==1 
        power_absences= SZ(i).powerabsences;
    elseif i>1
        power_absences= (power_absences + SZ(i).powerabsences)/2;
    end
end

frex= mean([SZ.freqestimate]);
% Plot frequency spectrum
figure(5)
plot(F, power_absences(:,5), 'b-', 'LineWidth', 2);
ylabel('PSD');
xlabel('Frequency (Hz)');
grid on;
xlim([0 30])
% Get frequency estimate (spectral peak)
title(['Frequency estimate = ',sprintf('%#.1f',round(frex,1)),' Hz']);

%% peak freq throughout absence 
Nfft=4*fs;    % Nfft, Window en Noverlap zijn instellingen van de functie pwelch. Deze kun je aanpassen om het spectrum 'gladder' of minder 'glad' te maken.
Window=1.5*fs; 
Noverlap=fs;

peakfreq_all = []; 
time_plot = []; 
for kk=1:size(pat3,1)
   seizure = EEG_filtered((pat3(kk,1):pat3(kk,2)),:);
   e = 0;
   for iii = 1:512:length(seizure) %every 2 s  
      
       if iii+512<length(seizure)
           freqabs_seiz = pwelch(seizure(iii:iii+512,:),Window,Noverlap,Nfft,fs);
           e = e+1;  
           [~,loc1] = max(freqabs_seiz);
           peakfreq_all(kk,e)= mean(F(loc1));
           time_plot(kk,e) = iii/fs;
       else
           e= 0; 
       end
       
   end
end



figure(10) 
scatter(nonzeros(time_plot(1,:)), nonzeros(peakfreq_all(1,:)), 'filled','MarkerFaceColor',	[0, 0.4470, 0.7410],'MarkerEdgeColor',[0,0,0])
hold on
Fit1 = polyfit(nonzeros(time_plot(1,:)),nonzeros(peakfreq_all(1,:)),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit1,nonzeros(time_plot(1,:)));
y1 = plot(nonzeros(time_plot(1,:)),y_est,'color',	[0, 0.4470, 0.7410],'LineWidth',2);
hold on

scatter(nonzeros(time_plot(2,:)), nonzeros(peakfreq_all(2,:)),'filled') 
hold on
Fit2 = polyfit(nonzeros(time_plot(2,:)),nonzeros(peakfreq_all(2,:)),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit2,nonzeros(time_plot(2,:)));
y1 = plot(nonzeros(time_plot(2,:)),y_est);
hold on 

scatter(nonzeros(time_plot(3,:)), nonzeros(peakfreq_all(3,:)),'filled','MarkerFaceColor',	[0.6350, 0.0780, 0.1840],'MarkerEdgeColor',[0,0,0]) 
hold on
Fit2 = polyfit(nonzeros(time_plot(3,:)),nonzeros(peakfreq_all(3,:)),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit2,nonzeros(time_plot(3,:)));
y1 = plot(nonzeros(time_plot(3,:)),y_est,'color',	[0.6350, 0.0780, 0.1840],'LineWidth',2);
hold on 
% % % 
scatter(nonzeros(time_plot(4,:)), nonzeros(peakfreq_all(4,:)),'filled') 
hold on
Fit2 = polyfit(nonzeros(time_plot(4,:)),nonzeros(peakfreq_all(4,:)),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit2,nonzeros(time_plot(4,:)));
y1 = plot(nonzeros(time_plot(4,:)),y_est);
hold on 
xlim([-1 max(time_plot(2,:))+1])
ylim([1 4])
xlabel('Time (s)')
ylabel('Peak Frequency (Hz)')
subject_nr = 5;
title(strcat('Peak frequency development during seizures, subj ', string(subject_nr)))


%velocity of peak freq change 
time1 = nonzeros(time_plot(1,1));
timelast_1 = time_plot(1,length(nonzeros(time_plot(1,:))));
timelast_2 = time_plot(2,length(nonzeros(time_plot(2,:))));
timelast_3 = time_plot(3,length(nonzeros(time_plot(3,:))));
timelast_4 = time_plot(4,length(nonzeros(time_plot(4,:))));
timelast_5 = time_plot(5,length(nonzeros(time_plot(5,:))));
timelast_6 = time_plot(6,length(nonzeros(time_plot(6,:))));
timelast_7 = time_plot(7,length(nonzeros(time_plot(7,:))));
timelast_8 = time_plot(8,length(nonzeros(time_plot(8,:))));
timelast_9 = time_plot(9,length(nonzeros(time_plot(9,:))));
timelast_10 = time_plot(10,length(nonzeros(time_plot(10,:))));
timelast_11 = time_plot(11,length(nonzeros(time_plot(11,:))));
timelast_12 = time_plot(12,length(nonzeros(time_plot(12,:))));
timelast_13 = time_plot(13,length(nonzeros(time_plot(13,:))));
timelast_14 = time_plot(14,length(nonzeros(time_plot(14,:))));
firstpeak1 = peakfreq_all(1,1);
firstpeak2 = peakfreq_all(2,1);
firstpeak3 = peakfreq_all(3,1);
firstpeak4 = peakfreq_all(4,1);
firstpeak5 = peakfreq_all(5,1);
firstpeak6 = peakfreq_all(6,1);
firstpeak7 = peakfreq_all(7,1);
firstpeak8 = peakfreq_all(8,1);
firstpeak9 = peakfreq_all(9,1);
firstpeak10 = peakfreq_all(10,1);
firstpeak11 = peakfreq_all(11,1);
firstpeak12 = peakfreq_all(12,1);
firstpeak13 = peakfreq_all(13,1);
firstpeak14 = peakfreq_all(14,1);
peaklast_1 = peakfreq_all(1,length(nonzeros(peakfreq_all(1,:))));
peaklast_2 = peakfreq_all(2,length(nonzeros(peakfreq_all(2,:))));
peaklast_3 = peakfreq_all(3,length(nonzeros(peakfreq_all(3,:))));
peaklast_4 = peakfreq_all(4,length(nonzeros(peakfreq_all(4,:))));
peaklast_5 = peakfreq_all(5,length(nonzeros(peakfreq_all(5,:))));
peaklast_6 = peakfreq_all(6,length(nonzeros(peakfreq_all(6,:))));
peaklast_7 = peakfreq_all(7,length(nonzeros(peakfreq_all(7,:))));
peaklast_8 = peakfreq_all(8,length(nonzeros(peakfreq_all(8,:))));
peaklast_9 = peakfreq_all(9,length(nonzeros(peakfreq_all(9,:))));
peaklast_10 = peakfreq_all(10,length(nonzeros(peakfreq_all(10,:))));
peaklast_11 = peakfreq_all(11,length(nonzeros(peakfreq_all(11,:))));
peaklast_12 = peakfreq_all(12,length(nonzeros(peakfreq_all(12,:))));
peaklast_13 = peakfreq_all(13,length(nonzeros(peakfreq_all(13,:))));
peaklast_14 = peakfreq_all(14,length(nonzeros(peakfreq_all(14,:))));



slope(1,:) = (peaklast_1-firstpeak1)/(timelast_1-time1);
slope(2,:) = (peaklast_2-firstpeak2)/(timelast_2-time1);
slope(3,:) = (peaklast_3-firstpeak3)/(timelast_3-time1);
slope(4,:) = (peaklast_4-firstpeak4)/(timelast_4-time1);
slope(5,:) = (peaklast_5-firstpeak5)/(timelast_5-time1);
slope(6,:) = (peaklast_6-firstpeak6)/(timelast_6-time1);
slope(7,:) = (peaklast_7-firstpeak7)/(timelast_7-time1);
slope(8,:) = (peaklast_8-firstpeak8)/(timelast_8-time1);
slope(9,:) = (peaklast_9-firstpeak9)/(timelast_9-time1);
slope(10,:) = (peaklast_10-firstpeak10)/(timelast_10-time1);
slope(11,:) = (peaklast_11-firstpeak11)/(timelast_11-time1);
slope(12,:) = (peaklast_12-firstpeak12)/(timelast_12-time1);
slope(13,:) = (peaklast_13-firstpeak13)/(timelast_13-time1);
slope(14,:) = (peaklast_14-firstpeak14)/(timelast_14-time1);

%% slope for patient (both groups)


%determine average freq change, to plot later all together 
mean_perpat = mean(peakfreq_all(:,1:4));
mean(nonzeros(peakfreq_all(:,3)));

pat3_meanpeakfreq = [2.8618,2.9276,2.7105,2.4934];
pat5_meanpeakfreq = [3.4901,3.3191,3.2730,3.0526,2.9572,3.1941,2.2993,2.7434];
pat6_meanpeakfreq = [3.3618,3.1513,3.1250];
pat7_meanpeakfreq = [2.4833,2.6017,2.5066, 2.3470, 2.1868];
pat8_meanpeakfreq = [3.1069,3.3553,3.5181,3.3515, 3.4549, 3.3985, 3.1689,3.0974,3.2336, 3.0395];
pat10_meanpeakfreq = [3.1902,3.4438,3.1992,3.0263,2.9079,2.6886];
pat12_meanpeakfreq = [2.4164,3.1466, 2.9130, 2.8257, 2.5461];
pat13_meanpeakfreq = [ 3.1754,2.9123, 2.8092];
pat14_meanpeakfreq = [3.4473,3.6578,3.6578];
pat15_meanpeakfreq = [ 3.2669,2.0639, 1.6447];

time1= [0.0039,2.0039,4.0039,6.0039,8.0039,10.0039,12.0039,14.0039,16.0039,18.0039,20.0039]; %time of subj 8 with the highest num of averaged peakfreq 
figure(10) 
scatter(time1(1,1:length(pat3_meanpeakfreq)), pat3_meanpeakfreq , 'filled','MarkerFaceColor',	[0,0,0],'MarkerEdgeColor',[0,0,0],'HandleVisibility','off')
hold on
Fit1 = polyfit(time1(1,1:length(pat3_meanpeakfreq)),pat3_meanpeakfreq,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit1,time1(1,1:length(pat3_meanpeakfreq)));
y1 = plot(time1(1,1:length(pat3_meanpeakfreq)),y_est,'-k','LineWidth',2);
hold on

scatter(time1(1,1:length(pat5_meanpeakfreq)), pat5_meanpeakfreq , 'filled','MarkerFaceColor',		[1 0 0],'MarkerEdgeColor',[0,0,0],'HandleVisibility','off')
hold on
Fit1 = polyfit(time1(1,1:length(pat5_meanpeakfreq)),pat5_meanpeakfreq,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit1,time1(1,1:length(pat5_meanpeakfreq)));
y1 = plot(time1(1,1:length(pat5_meanpeakfreq)),y_est,	'-r','LineWidth',2);
hold on


% scatter(time1(1,1:length(pat6_meanpeakfreq)), pat6_meanpeakfreq , 'filled','MarkerFaceColor',	[0,0,0],'MarkerEdgeColor',[0,0,0])
% hold on
% Fit1 = polyfit(time1(1,1:length(pat6_meanpeakfreq)),pat6_meanpeakfreq,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
% y_est = polyval(Fit1,time1(1,1:length(pat6_meanpeakfreq)));
% y1 = plot(time1(1,1:length(pat6_meanpeakfreq)),y_est,'--k','LineWidth',2);
% hold on


scatter(time1(1,1:length(pat7_meanpeakfreq)), pat7_meanpeakfreq , 'filled','MarkerFaceColor',		[1 0 0],'MarkerEdgeColor',[0,0,0])
hold on
Fit1 = polyfit(time1(1,1:length(pat7_meanpeakfreq)),pat7_meanpeakfreq,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit1,time1(1,1:length(pat7_meanpeakfreq)));
y1 = plot(time1(1,1:length(pat7_meanpeakfreq)),y_est,'--r','LineWidth',2);
hold on


scatter(time1(1,1:length(pat8_meanpeakfreq)), pat8_meanpeakfreq , 'filled','MarkerFaceColor',		[1 0 0],'MarkerEdgeColor',[0,0,0])
hold on
Fit1 = polyfit(time1(1,1:length(pat8_meanpeakfreq)),pat8_meanpeakfreq,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit1,time1(1,1:length(pat8_meanpeakfreq)));
y1 = plot(time1(1,1:length(pat8_meanpeakfreq)),y_est,'-r','LineWidth',1);
hold on


xlim([-2 max(time1)+1])
ylim([1 4])
set(gca,'fontname','times','FontSize',15)
xlabel('Time (s)')
ylabel('Peak Frequency (Hz)')
title('Peak frequency development during seizures')
legend({'black lines = Preserved','red lines = Unpreserved'},'fontname','times','FontSize',12) 

% plot boxplots of slopes
tabb = xlsread('C:\Users\vb\Desktop\V\MOTION\PhD\paper4-absencecharacteristics\absences_features\absences_characteristics.xlsx','typatyp');


slope_p = [];
slope_u = []; 
meanpeakfreq_p= [];
meanpeakfreq_u= [];
e= 1;
e1 = 1;
for jj = 1:length(tabb(:,1))
    if tabb(jj,1) == 3 || tabb(jj,1) == 6 || tabb(jj,1) == 10 || tabb(jj,1) == 12 || tabb(jj,1) == 14
       slope_p(1,e) = abs(tabb(jj,8));  
       meanpeakfreq_p(1,e) = tabb(jj,7);
       e=e+1;
       
    elseif tabb(jj,1) == 5 || tabb(jj,1) == 7 ||  tabb(jj,1) == 8 ||  tabb(jj,1) == 13 ||  tabb(jj,1) == 15 
       slope_u(1,e1) = abs(tabb(jj,8));  
       meanpeakfreq_u(1,e1) = tabb(jj,7);
       e1 = e1+1;
    end
end
C1 = [slope_p slope_u];
grp = [zeros(1,length(slope_p)),ones(1,length(slope_u))];
figure(11)
boxplot(C1,grp,'Colors','k','Symbol','*k','Labels',{'Preserved','Unpreserved'})
hold on
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
ylim([-1 2])
set(gca,'fontname','times','FontSize',15)
ylabel('Slope (Hz/s)')
title('Slope of peak frequency development during seizure')

%statistics 
%[p,h] = ranksum(slope_p,slope_u) %0.03!
%try MANOVA too instead of many multiple corrections

C2 = [meanpeakfreq_u meanpeakfreq_p];
grp = [zeros(1,length(meanpeakfreq_u)),ones(1,length(meanpeakfreq_p))];
figure(11)
boxplot(C2,grp,'Colors','k','Symbol','*k','Labels',{'Unpreserved','Preserved'})
hold on
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'r');
ylim([1 4.5])
set(gca,'fontname','times','FontSize',14)
ylabel('Peak frequency (Hz)')
title('\rm Peak frequency')

%statistics 
%[p,h] = ranksum(meanpeakfreq_p,meanpeakfreq_u) %0.02!
%bonferroni 0.05/2 = 0.025 --> questo sarebbe ok e l altro no 

%% centre of gravity 
%keep monopolar for COG

A=[-1 1 -2 -1 0 1 2 -2 -1 0 1 2 -2 -1 0 1 2 -1 1];  % x-richting 
B=[2 2 1 1 1 1 1 0 0 0 0 0 -1 -1 -1 -1 -1 -2 -2];  % y-richting


%cog between 2.5 and 4.5 Hz
COG = struct;
for j=1:size(SZ,2)
    COG(j).x = mean((SZ(j).powerabsences(5:13,:)*A')./sum(SZ(j).powerabsences(5:13,:)')');
    COG(j).y = mean((SZ(j).powerabsences(5:13,:)*B')./sum(SZ(j).powerabsences(5:13,:)')');
end

XX= mean([COG.x]);
YY = mean([COG.y]);

%preserved vs unpreserved
XX =-0.0018;
YY =0.2732;
XX1= 0.0084;
YY1 = 0.4635;
%plot

figure(6)
viscircles([0 0],2, 'Color', [0.3010 0.7450 0.9330])
hold on
%viscircles([-1 1.7],0.2, 'Color', 'black')
%hold on
%viscircles([1 1.7],0.2, 'Color', 'black')
p = nsidedpoly(3, 'Center', [0 ,2.1], 'SideLength', 0.3);
plot(p);
hold on
labels = {'preserved'};
labels1 = {'unpreserved'};
plot(XX,YY, 'o','linewidth',5, 'Color', [0.6350 0.0780 0.1840])
plot(XX1,YY1, 'o','linewidth',5, 'Color', [0 0.4470 0.7410])
text(XX,YY,-30e2, labels,'VerticalAlignment','top','HorizontalAlignment','right', 'Color', [0.6350 0.0780 0.1840], 'FontSize', 14)
text(XX1,YY1,labels1,'VerticalAlignment','bottom','HorizontalAlignment','left','Color', [0 0.4470 0.7410], 'FontSize', 14 )
hold on
plot(-3:3,linspace(0,0,7), '--', 'Color', [17 17 17]/255)
plot(linspace(0,0,7), -3:3, '--', 'Color', [17 17 17]/255)
title('Center of Gravity during absences', 'FontSize', 18)


%% time-freq plot
%logarithmic normalization is better for lower freq, so we use it.
%extract power for 2 s epochs
frex = F(7:29,:); %1.5 - 7 Hz
timevec = 0:2:9; %in s, last value depends on lenght of seizure
m = 507; %start of seizure in s
Nfft=4*fs;    % Nfft, Window en Noverlap zijn instellingen van de functie pwelch. Deze kun je aanpassen om het spectrum 'gladder' of minder 'glad' te maken.
Window=2*fs; %3cycles/4Hz(lower freq i am interested in)= 0.75. 
Noverlap=fs;
Y = [];
for i = 1:length(timevec)
    
    freqabs = pwelch(EEG_filtered((m*fs:(m+2)*fs),:),Window,Noverlap,Nfft,fs);
    average_power = mean(freqabs')';
    average_power = 10*log10(average_power(7:29,:)+1); % just between 1.5-7 Hz, logarithmic scaling
    m = m +2;
    Y(:,i) = average_power;
end


figure(3)
contourf(timevec,frex,Y,11,'linecolor','none')
%set(gca,'clim',[0 1]*10000,'xlim',[0 13])
xlabel('Time (s)'), ylabel('Frequency (Hz)')
title('Logarithmic frequency scaling between 1-6 Hz')
%caxis([-1 1]);
set(gca,'fontsize',9);
%set(gca,'linewidth',2);
axis on;  
box off;
colorbar    
colormap('jet')   

% to normalize use a ratio, like M says, because dB conversion is basically
% ratio activity/baseline. 

tt = linspace(0,1,size(Y,1));
Y1 = Y';
scatter(tt, Y1(1,:))
