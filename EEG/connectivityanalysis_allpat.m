 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Imports EEG set file of each patient. 
  % determines bivariate connectivity during each absence of each patient 
  % using Phase lag index (PLI) at 3 Hz
  % PLI is averaged for each patient and averaged PLI for two groups
  % (unpreserved vs preserved) is calculated
  % graph analyses as hubness, cluster coefficient and path length is
  % determined 
  %  
  % needs fieldtrip
  % needs standard_bem.mat
  % needs standard_1005.elc
  % needs pathlength.m 
  % needs plot_cluster coef 
  
  % VBarone Nov, 2022
  % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear

%% PLI for all seizure all patients 

for j = 1:10 %num of included patients 
if j ==1 
   subject_nr = 3;
elseif j==2   
   subject_nr = 5; 
elseif j==3  
   subject_nr = 6;
elseif j==4   
   subject_nr = 7;  
elseif j==5   
   subject_nr = 8;    
elseif j==6   
   subject_nr = 10; 
elseif j==7   
   subject_nr = 12;   
elseif j==8   
   subject_nr = 13; 
elseif j==9   
   subject_nr = 14;      
elseif j==10   
   subject_nr = 15;    
end   

%% load file
EEGfilename = strcat('C:\Users\vb\Documents\MATLAB\AbsenceChildren_paper4\data\Kemp\EEG\subj',string(subject_nr),'\subject',string(subject_nr),'_ICA_KH.set');
EEGfilename = convertStringsToChars(EEGfilename);
fs = 256;

cfg = []; 
cfg.dataset = EEGfilename; 
eeg_raw = ft_preprocessing(cfg);

%% define trials 

if subject_nr == 3
     begin_seizure = [round(508.45*fs),round(841.79*fs)]; 
     end_seizure = [round(515.6*fs),round(847.52*fs)]; 
elseif subject_nr == 5  
    begin_seizure = [round(350.5*fs),round(901.1*fs),round(1326.55*fs),round(1895.15*fs)]; 
    end_seizure = [round(367.5*fs),round(920.1*fs),round(1358*fs),round(1918*fs)]; 
  
elseif subject_nr == 6 
    begin_seizure = [round(152.5*fs),round(1996.7*fs)]; 
    end_seizure = [round(158.5*fs),round(2003.9*fs)]; 
elseif subject_nr == 7     
    begin_seizure = [round(274.4*fs),round(573*fs),round(827.5*fs),round(996.1*fs),round(1083.2*fs),round(1476*fs),round(1727.8*fs),round(1860*fs),round(2171.8*fs),round(2336*fs)]; 
    end_seizure = [round(277.4*fs),round(579.5*fs),round(834.5*fs),round(1002*fs),round(1090.2*fs),round(1484*fs),round(1735.5*fs),round(1869*fs),round(2177.4*fs),round(2346*fs)]; 
  
elseif subject_nr == 8  
    begin_seizure = [round(487.8*fs),round(590.75*fs),round(976.8*fs),round(1269.5*fs),round(1512.5*fs),round(1794*fs),round(1925*fs)]; 
    end_seizure = [round(498.8*fs),round(605*fs),round(991.4*fs),round(1292*fs),round(1527.8*fs),round(1826*fs),round(1945*fs)]; 
   
elseif subject_nr == 10
     begin_seizure = [round(113*fs),round(1257.3*fs),round(1293.2*fs),round(1843.7*fs),round(2063.5*fs)]; 
     end_seizure = [round(128*fs),round(1261.2*fs),round(1308.2*fs),round(1860.7*fs),round(2066.5*fs)]; 
  
elseif subject_nr == 12
    begin_seizure = [round(39.1*fs),round(145.2*fs),round(269.1*fs),round(447.3*fs),round(525.9*fs),round(856.3*fs),round(1250.55*fs),round(1647.65*fs),round(1779*fs),round(1832*fs),round(1976.4*fs),round(2144.7*fs),round(2418.1*fs),round(2501*fs)]; 
    end_seizure = [round(44.7*fs),round(157.2*fs),round(278.2*fs),round(455*fs),round(533.5*fs),round(864.5*fs),round(1257.6*fs),round(1653.8*fs),round(1785.5*fs),round(1839*fs),round(1980.4*fs),round(2153.8*fs),round(2426.1*fs),round(2511*fs)]; 
  
elseif subject_nr == 13
    begin_seizure = [round(178.8*fs),round(1647.65*fs)]; 
    end_seizure = [round(183*fs),round(1652.7*fs)]; 
  
elseif subject_nr == 14
    begin_seizure = [round(120.4*fs)]; 
    end_seizure = [round(128.4*fs)]; 
   
elseif subject_nr == 15     
    begin_seizure = [round(575*fs),round(654*fs),round(705*fs),round(1161*fs),round(1225.6*fs),round(1297*fs),round(1357*fs)];
    end_seizure = [round(580*fs),round(659*fs),round(711*fs),round(1166*fs),round(1231.6*fs),round(1302*fs),round(1365*fs)];
end     




%% create trial
trl_startstop = [];
for yy = 1:length(begin_seizure)
    trl_startstop(yy,1) =  begin_seizure(yy);
    trl_startstop(yy,2) =  end_seizure(yy);
end

data_trl = {};
for kk=1:size(trl_startstop,1)
cfg = [];
cfg.trl = [trl_startstop(kk,:) ones(size(trl_startstop(kk,:),1),1)];
cfg.dataset= EEGfilename;
data_trl{1,kk}= {ft_preprocessing(cfg)};
end



% specify frequencies
frex = 3; %frex = [2 3 4];
nCycles = 6; %nCycles = [ 5 6 7];

% parameters for complex Morlet wavelets -- we use wavelet and not welch
% method coz wavelet has better freq resolution (welch s masks weak
% component --> spectral leakage due to windowing)
wavtime  = -1:1/fs:1-1/fs;
half_wav = (length(wavtime)-1)/2;



%% define EEG signal
for tt = 1:size(begin_seizure,2)
     EEG_signal = data_trl{1,tt}{1,1}.trial{1,1}; %data{1,num of seiz}{1,1}.trial

% useful variables for later...
     [nchans, npnts, ntrials] = size(EEG_signal);
     
     % FFT parameters
     nWave = length(wavtime);
     nData = npnts*ntrials;
     nConv = nWave+nData-1;

% and create wavelets (sine wave taper by a gaussian)
     cmwX = zeros(length(frex),nConv);
    
     s= nCycles / (2*pi*frex);
     cmw= exp(1i*2*pi*frex.*wavtime) .* exp( (-wavtime.^2) ./ (2*s^2) );
     tempX= fft(cmw,nConv);
     cmwX(1,:) = tempX ./ max(tempX);

     
     allphases = zeros(nchans,length(frex),npnts,ntrials);

% spectrum of all channels using the fft matrix input
     dataX = fft( reshape(EEG_signal,nchans,npnts*ntrials) ,nConv,2);
    
    % run convolution (= A LOT OF TIME-SERIES DOT PRODUCTS) --> we do
    % convolution in freq domain coz it is way faster than in time domain
     as = ifft( bsxfun(@times,dataX,cmwX(1,:)) ,nConv,2 ); %then back to time domain with ifft
     as = as(:,half_wav+1:end-half_wav);   
    % phase values from all trials
     allphases(:,1,:) = angle( as );

     
     % PLI connectivity
     connmat_pli_all = zeros(nchans,nchans,length(frex));
 
     for chani=1:nchans
         for chanj=1:nchans
            connmat_pli_all(chani,chanj,:) = abs(mean(sign(imag(exp(1i* squeeze(allphases(chani,:,:)-allphases(chanj,:,:)) ))),1)) ;
         end
     end
     
     %connectivities of all seizures for each subj
     if subject_nr == 3
        connmat_pli_subj3{tt,:} = connmat_pli_all;
     elseif subject_nr == 5 
        connmat_pli_subj5{tt,:} = connmat_pli_all; 
     elseif subject_nr == 6 
        connmat_pli_subj6{tt,:} = connmat_pli_all; 
     elseif subject_nr == 7 
        connmat_pli_subj7{tt,:} = connmat_pli_all;
     elseif subject_nr == 8 
        connmat_pli_subj8{tt,:} = connmat_pli_all; 
     elseif subject_nr == 10 
        connmat_pli_subj10{tt,:} = connmat_pli_all; 
     elseif subject_nr == 12
        connmat_pli_subj12{tt,:} = connmat_pli_all;  
     elseif subject_nr == 13 
        connmat_pli_subj13{tt,:} = connmat_pli_all;  
     elseif subject_nr == 14 
        connmat_pli_subj14{tt,:} = connmat_pli_all; 
     elseif subject_nr == 15
        connmat_pli_subj15{tt,:} = connmat_pli_all; 
     end
end
end

%% plot single seizure

figure(1)
clim  = [0 1];
climD = [-.4 .4];
imagesc(squeeze(connmat_pli_subj15{5,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj15, PLI ' num2str(frex(1)) ' Hz, seizure 5' ])

figure(2)
clim  = [0 1];
climD = [-.4 .4];
imagesc(squeeze(connmat_pli_subj15{6,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj15, PLI ' num2str(frex(1)) ' Hz, seizure 6' ])


%% plot first seizure of all pat 
% define color limits
clim  = [0 1];
climD = [-.4 .4];

figure(3), set(gcf,'Name','Movie magic minimizes the magic.', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);



%sgtitle('Subplot Grid Title') 
subplot(2,5,1)
imagesc(squeeze(connmat_pli_subj3{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj3, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,2)
imagesc(squeeze(connmat_pli_subj6{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj6, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,3)
imagesc(squeeze(connmat_pli_subj10{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj10, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,4)
imagesc(squeeze(connmat_pli_subj12{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj12, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,5)
imagesc(squeeze(connmat_pli_subj14{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj14, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,6)
imagesc(squeeze(connmat_pli_subj5{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj5, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,7)
imagesc(squeeze(connmat_pli_subj7{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj7, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,8)
imagesc(squeeze(connmat_pli_subj8{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj8, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,9)
imagesc(squeeze(connmat_pli_subj13{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj13, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,10)
imagesc(squeeze(connmat_pli_subj15{1,1}))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj15, PLI' num2str(frex(1)) ' Hz' ])

ha = axes('Position',[0, 0, 1, 0.5],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Unpreserved (first seizure)')
hold on
ha = axes('Position',[0, 0, 1, 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Preserved (first seizure)')

colorbar

% PLI seizures as trials cannot do coz seizures have different lengths 

%% averaged PLI per patient 

subj3_PLI = (connmat_pli_subj3{1,1} + connmat_pli_subj3{2,1})/2;
subj5_PLI = (connmat_pli_subj5{1,1} + connmat_pli_subj5{2,1}+connmat_pli_subj5{3,1}+connmat_pli_subj5{4,1})/4;
subj6_PLI = (connmat_pli_subj6{1,1} + connmat_pli_subj6{2,1})/2;
subj7_PLI = (connmat_pli_subj7{1,1} + connmat_pli_subj7{2,1}+connmat_pli_subj7{3,1} +connmat_pli_subj7{4,1} + connmat_pli_subj7{5,1} +connmat_pli_subj7{6,1} +connmat_pli_subj7{7,1} +connmat_pli_subj7{8,1} + connmat_pli_subj7{9,1} +connmat_pli_subj7{10,1} )/10;
subj8_PLI = (connmat_pli_subj8{1,1} + connmat_pli_subj8{2,1}+connmat_pli_subj8{3,1}+connmat_pli_subj8{4,1}+connmat_pli_subj8{5,1}+connmat_pli_subj8{6,1}+connmat_pli_subj8{7,1})/7;
subj10_PLI = (connmat_pli_subj10{1,1} + connmat_pli_subj10{2,1}+connmat_pli_subj10{3,1}+connmat_pli_subj10{4,1}+connmat_pli_subj10{5,1})/5;
subj12_PLI = (connmat_pli_subj12{1,1} + connmat_pli_subj12{2,1}+connmat_pli_subj12{3,1}+connmat_pli_subj12{4,1}+connmat_pli_subj12{5,1}+connmat_pli_subj12{6,1}+connmat_pli_subj12{7,1}+connmat_pli_subj12{8,1}+connmat_pli_subj12{9,1}+connmat_pli_subj12{10,1}+connmat_pli_subj12{11,1}+connmat_pli_subj12{12,1}+connmat_pli_subj12{13,1}+connmat_pli_subj12{14,1})/14;
subj13_PLI = (connmat_pli_subj13{1,1} + connmat_pli_subj13{2,1})/2;
subj14_PLI = connmat_pli_subj13{1,1};
subj15_PLI = (connmat_pli_subj15{1,1} + connmat_pli_subj15{2,1}+connmat_pli_subj15{3,1}+connmat_pli_subj15{4,1}+connmat_pli_subj15{5,1}+connmat_pli_subj15{6,1}+connmat_pli_subj15{7,1})/7;

%% plot average of seizures of all pat 
% define color limits
clim  = [0 1];
climD = [-.4 .4];

figure(4), set(gcf,'Name','Movie magic minimizes the magic.', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);



%sgtitle('Subplot Grid Title') 
subplot(2,5,1)
imagesc(squeeze(subj3_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj3, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,2)
imagesc(squeeze(subj6_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj6, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,3)
imagesc(squeeze(subj10_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj10, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,4)
imagesc(squeeze(subj12_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj12, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,5)
imagesc(squeeze(subj14_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj14, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,6)
imagesc(squeeze(subj5_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj5, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,7)
imagesc(squeeze(subj7_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj7, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,8)
imagesc(squeeze(subj8_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj8, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,9)
imagesc(squeeze(subj13_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj13, PLI' num2str(frex(1)) ' Hz' ])


subplot(2,5,10)
imagesc(squeeze(subj15_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'subj15, PLI' num2str(frex(1)) ' Hz' ])

ha = axes('Position',[0, 0, 1, 0.5],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Unpreserved (average)')
hold on
ha = axes('Position',[0, 0, 1, 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Preserved (average)')

colorbar

%% average of two groups 

preserved_PLI = (subj3_PLI + subj6_PLI+ subj10_PLI + subj12_PLI + subj14_PLI)/5;
unpreserved_PLI = (subj5_PLI + subj7_PLI+ subj8_PLI + subj13_PLI + subj15_PLI)/5;

clim  = [0 1];
climD = [-.4 .4];

figure(6), 
colorbar
hold on
%sgtitle('Subplot Grid Title') 
subplot(1,2,1)
imagesc(squeeze(preserved_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Preserved, PLI ' num2str(frex(1)) ' Hz' ])


subplot(1,2,2)
imagesc(squeeze(unpreserved_PLI))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Unpreserved, PLI' num2str(frex(1)) ' Hz' ])

%hist(unpreserved_PLI)
[h,p,ci,stats] = ttest2(preserved_PLI,unpreserved_PLI) %data is normally distributed


%% graph analysis 

% define threshold
threshold_unp = median(median(unpreserved_PLI)) +  std(std(unpreserved_PLI));
threshold_p =  median(median(preserved_PLI)) +  std(std(preserved_PLI));

% binarize the matrix

binary_unp = zeros(length(unpreserved_PLI),length(unpreserved_PLI));
%unpr
for chani=1:length(unpreserved_PLI)
    for chanj=1:length(unpreserved_PLI)
         if unpreserved_PLI(chani,chanj) > threshold_unp
            binary_unp(chani,chanj) = 1;
         elseif unpreserved_PLI(chani,chanj) < threshold_unp
            binary_unp(chani,chanj) = 0; 
         end
    end
end

binary_p = zeros(length(preserved_PLI),length(preserved_PLI));
%unpr
for chani=1:length(preserved_PLI)
    for chanj=1:length(preserved_PLI)
         if preserved_PLI(chani,chanj) > threshold_unp
            binary_p(chani,chanj) = 1;
         elseif preserved_PLI(chani,chanj) < threshold_unp
            binary_p(chani,chanj) = 0; 
         end
    end
end

clim  = [0 1];
figure(7), 

%sgtitle('Subplot Grid Title') 
subplot(1,2,1)
imagesc(squeeze(binary_unp))
colormap(gray)
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI Unpreserved, PLI ' num2str(frex(1)) ' Hz' ])


subplot(1,2,2)
imagesc(squeeze(binary_p))
colormap(gray)
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI Preserved, PLI ' num2str(frex(1)) ' Hz' ])
colorbar('northoutside')
%colorbar(axes,'Position',[0.347142857142856 0.871349207302884 0.335714285714285 0.0447619047619048]);


%sum suprathreshold connections for each channel and divide by number of
%theoretically possible connections (in my case is 19-1 , coz 1 channel is
%the same channel (Fz-Fz), so we dont wanna consider that)

%proportion of suprathreshold connections = hubness
chan_prop_p = sum(binary_p,2)/18;
chan_prop_unp = sum(binary_unp,2)/18;

ch_list = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8','T7','C3','Cz','C4', 'T8','P7','P3','Pz','P4','P8', 'O1','O2'};

addpath C:\Users\vb\Documents\MATLAB\TBIstudy_paper3\scripts\plot_topography\

figure(10)
clim = [0 1];

plot_topography(ch_list, chan_prop_p)
colormap(jet)
set(gca,'clim',clim)
title('Preserved','fontname','times', 'FontSize', 15)
c = colorbar;
c.Limits = clim;

figure(11)

plot_topography(ch_list, chan_prop_unp)
colormap(jet)
set(gca,'clim',clim)
title('Unpreserved','fontname','times', 'FontSize', 15)
c = colorbar;
c.Limits = clim;

[p,h] = ranksum(chan_prop_p,chan_prop_unp) %>0.05
[p,h] = ranksum(mean(preserved_PLI),mean(unpreserved_PLI)) %>0.05

%same for average of single patients 
threshold_subj3 = median(median(subj3_PLI)) +  std(std(subj3_PLI));
threshold_subj5 = median(median(subj5_PLI)) +  std(std(subj5_PLI));
threshold_subj6 = median(median(subj6_PLI)) +  std(std(subj6_PLI));
threshold_subj7 = median(median(subj7_PLI)) +  std(std(subj7_PLI));
threshold_subj8 = median(median(subj8_PLI)) +  std(std(subj8_PLI));
threshold_subj10 = median(median(subj10_PLI)) +  std(std(subj10_PLI));
threshold_subj12 = median(median(subj12_PLI)) +  std(std(subj12_PLI));
threshold_subj13 = median(median(subj13_PLI)) +  std(std(subj13_PLI));
threshold_subj14 = median(median(subj14_PLI)) +  std(std(subj14_PLI));
threshold_subj15 = median(median(subj15_PLI)) +  std(std(subj15_PLI));


% binarize the matrix

binary_subj3 = zeros(length(subj3_PLI),length(subj3_PLI));
%unpr
for chani=1:length(subj3_PLI)
    for chanj=1:length(subj3_PLI)
         if subj3_PLI(chani,chanj) > threshold_subj3
            binary_subj3(chani,chanj) = 1;
         elseif subj3_PLI(chani,chanj) < threshold_subj3
            binary_subj3(chani,chanj) = 0; 
         end
    end
end
binary_subj5 = zeros(length(subj5_PLI),length(subj5_PLI));
%unpr
for chani=1:length(subj5_PLI)
    for chanj=1:length(subj5_PLI)
         if subj5_PLI(chani,chanj) > threshold_subj5
            binary_subj5(chani,chanj) = 1;
         elseif subj5_PLI(chani,chanj) < threshold_subj5
            binary_subj5(chani,chanj) = 0; 
         end
    end
end
binary_subj6 = zeros(length(subj6_PLI),length(subj6_PLI));
%unpr
for chani=1:length(subj6_PLI)
    for chanj=1:length(subj6_PLI)
         if subj6_PLI(chani,chanj) > threshold_subj6
            binary_subj6(chani,chanj) = 1;
         elseif subj6_PLI(chani,chanj) < threshold_subj6
            binary_subj6(chani,chanj) = 0; 
         end
    end
end
binary_subj7 = zeros(length(subj7_PLI),length(subj7_PLI));
%unpr
for chani=1:length(subj7_PLI)
    for chanj=1:length(subj7_PLI)
         if subj7_PLI(chani,chanj) > threshold_subj7
            binary_subj7(chani,chanj) = 1;
         elseif subj7_PLI(chani,chanj) < threshold_subj7
            binary_subj7(chani,chanj) = 0; 
         end
    end
end
binary_subj8 = zeros(length(subj8_PLI),length(subj8_PLI));
%unpr
for chani=1:length(subj8_PLI)
    for chanj=1:length(subj8_PLI)
         if subj8_PLI(chani,chanj) > threshold_subj8
            binary_subj8(chani,chanj) = 1;
         elseif subj8_PLI(chani,chanj) < threshold_subj8
            binary_subj8(chani,chanj) = 0; 
         end
    end
end
binary_subj10 = zeros(length(subj10_PLI),length(subj10_PLI));
%unpr
for chani=1:length(subj10_PLI)
    for chanj=1:length(subj10_PLI)
         if subj10_PLI(chani,chanj) > threshold_subj10
            binary_subj10(chani,chanj) = 1;
         elseif subj10_PLI(chani,chanj) < threshold_subj10
            binary_subj10(chani,chanj) = 0; 
         end
    end
end
binary_subj12 = zeros(length(subj12_PLI),length(subj12_PLI));
%unpr
for chani=1:length(subj12_PLI)
    for chanj=1:length(subj12_PLI)
         if subj12_PLI(chani,chanj) > threshold_subj12
            binary_subj12(chani,chanj) = 1;
         elseif subj12_PLI(chani,chanj) < threshold_subj12
            binary_subj12(chani,chanj) = 0; 
         end
    end
end
binary_subj13 = zeros(length(subj13_PLI),length(subj13_PLI));
%unpr
for chani=1:length(subj13_PLI)
    for chanj=1:length(subj13_PLI)
         if subj13_PLI(chani,chanj) > threshold_subj13
            binary_subj13(chani,chanj) = 1;
         elseif subj13_PLI(chani,chanj) < threshold_subj13
            binary_subj13(chani,chanj) = 0; 
         end
    end
end
binary_subj14 = zeros(length(subj14_PLI),length(subj14_PLI));
%unpr
for chani=1:length(subj14_PLI)
    for chanj=1:length(subj14_PLI)
         if subj14_PLI(chani,chanj) > threshold_subj14
            binary_subj14(chani,chanj) = 1;
         elseif subj14_PLI(chani,chanj) < threshold_subj14
            binary_subj14(chani,chanj) = 0; 
         end
    end
end
binary_subj15 = zeros(length(subj15_PLI),length(subj15_PLI));
%unpr
for chani=1:length(subj15_PLI)
    for chanj=1:length(subj15_PLI)
         if subj15_PLI(chani,chanj) > threshold_subj15
            binary_subj15(chani,chanj) = 1;
         elseif subj15_PLI(chani,chanj) < threshold_subj15
            binary_subj15(chani,chanj) = 0; 
         end
    end
end

%plot binary matrices
clim  = [0 1];

figure(12), set(gcf,'Name','Movie magic minimizes the magic.', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colormap(gray)


%sgtitle('Subplot Grid Title') 
subplot(2,5,1)
imagesc(squeeze(binary_subj3))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj3, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,2)
imagesc(squeeze(binary_subj6))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj6, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,3)
imagesc(squeeze(binary_subj10))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj10, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,4)
imagesc(squeeze(binary_subj12))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj12, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,5)
imagesc(squeeze(binary_subj14))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj14, PLI ' num2str(frex(1)) ' Hz'])


subplot(2,5,6)
imagesc(squeeze(binary_subj5))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj5, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,7)
imagesc(squeeze(binary_subj7))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj7, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,8)
imagesc(squeeze(binary_subj8))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj8, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,9)
imagesc(squeeze(binary_subj13))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj13, PLI ' num2str(frex(1)) ' Hz' ])


subplot(2,5,10)
imagesc(squeeze(binary_subj15))
axis square
set(gca,'clim',clim,'xtick',1:nchans,'ytick',1:nchans)
title([ 'Binarized PLI subj15, PLI ' num2str(frex(1)) ' Hz' ])

ha = axes('Position',[0, 0, 1, 0.5],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Unpreserved (average)')
hold on
ha = axes('Position',[0, 0, 1, 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 0.98,'Preserved (average)')

colorbar


% topoplots

%proportion of suprathreshold connections = hubness
chan_prop_subj3 = sum(binary_subj3,2)/18;
chan_prop_subj5 = sum(binary_subj5,2)/18;
chan_prop_subj6 = sum(binary_subj6,2)/18;
chan_prop_subj7 = sum(binary_subj7,2)/18;
chan_prop_subj8 = sum(binary_subj8,2)/18;
chan_prop_subj10 = sum(binary_subj10,2)/18;
chan_prop_subj12 = sum(binary_subj12,2)/18;
chan_prop_subj13 = sum(binary_subj13,2)/18;
chan_prop_subj14 = sum(binary_subj14,2)/18;
chan_prop_subj15 = sum(binary_subj15,2)/18;

ch_list = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8','T7','C3','Cz','C4', 'T8','P7','P3','Pz','P4','P8', 'O1','O2'};

addpath C:\Users\vb\Documents\MATLAB\TBIstudy_paper3\scripts\plot_topography\

figure(13)
clim = [0 1];
plot_topography(ch_list, chan_prop_subj15)
title('subj15 Hubness')
colormap(jet)
c = colorbar;
c.Limits = clim;
set(gca,'clim',clim)

%% topographical patterns (which connections are bigger) 
%print -depsc graphunp.eps

hub_unp = find(chan_prop_unp == max(chan_prop_unp)); 
hub_pres = find(chan_prop_p == max(chan_prop_p)); 

connect_1_un= find(binary_unp(hub_unp(1),:)==1)' ;
connect_2_un= find(binary_unp(hub_unp(2),:)==1)' ;
connect_p= find(binary_p(hub_pres(1),:)==1)';


pattern_conn1_un = [];
pattern_conn2_un = [];
pattern_conn_p = [];
for i = 1:19   
  for j = 1:length(connect_1_un)  
    if ismember(i,connect_1_un,'rows')
       pattern_conn1_un(i,1) = 1;
    else
       pattern_conn1_un(i,1) = 0;
    end
  end
end
    
    
for i = 1:19   
  for j = 1:length(connect_2_un)  
    if ismember(i,connect_2_un,'rows')
       pattern_conn2_un(i,1) = 1;
    else
       pattern_conn2_un(i,1) = 0;
    end
  end
end

for i = 1:19   
  for j = 1:length(connect_p)  
    if ismember(i,connect_p,'rows')
       pattern_conn_p(i,1) = 1;
    else
       pattern_conn_p(i,1) = 0;
    end
  end
end

ch_list = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8','T7','C3','Cz','C4', 'T8','P7','P3','Pz','P4','P8', 'O1','O2'};

addpath C:\Users\vb\Documents\MATLAB\TBIstudy_paper3\scripts\plot_topography\

figure(13)
clim = [0 1];
plot_topography(ch_list, pattern_conn2_un)
title('Unreserved topographical pattern of connectivity hub2')
colormap(jet)
c = colorbar;
c.Limits = clim;
set(gca,'clim',clim)

%% clustering coefficient 

%create local network = 

connectivity_degree = []; 
for i = 1:19
    connectivity_degree = find(binary_unp(i,:)==1);
    n = length(connectivity_degree);
    if n>3
                % "local" network of neighbors
                localnetwork = binary_unp(connectivity_degree,connectivity_degree);
                % localnetwork is symmetric; remove redundant values by replacing with NaN
                localnetwork = localnetwork + tril(nan(n));
                
                % compute cluster coefficient (neighbor connectivity scaled)
                clustcoef_unp(i) = 2*nansum(localnetwork(:)) / ((n-1)*n);
    end
end



connectivity_degree = []; 
for i = 1:19
    connectivity_degree = find(binary_p(i,:)==1);
    n = length(connectivity_degree);
    if n>3
                % "local" network of neighbors
                localnetwork = binary_p(connectivity_degree,connectivity_degree);
                % localnetwork is symmetric; remove redundant values by replacing with NaN
                localnetwork = localnetwork + tril(nan(n));
                
                % compute cluster coefficient (neighbor connectivity scaled)
                clustcoef_p(i) = 2*nansum(localnetwork(:)) / ((n-1)*n);
    end
end

ch_list = {'Fp1', 'Fp2', 'F7', 'F3', 'Fz', 'F4', 'F8','T7','C3','Cz','C4', 'T8','P7','P3','Pz','P4','P8', 'O1','O2'};

addpath C:\Users\vb\Documents\MATLAB\TBIstudy_paper3\scripts\plot_topography\

figure(14)
clim = [0 1];
plot_topography(ch_list, clustcoef_unp')
title('Unpreserved cluster coefficient')
colormap(jet)
c = colorbar;
c.Limits = clim;
set(gca,'clim',clim)

figure(15)
clim = [0 1];
plot_topography(ch_list, clustcoef_p')
title('Preserved cluster coefficient')
colormap(jet)
c = colorbar;
c.Limits = clim;
set(gca,'clim',clim)

%measure of network interaction

%connectivity degree (number of 1 in binary matrix) reflects long-range
%connectivity (both close and far networks), while clustering coefficient
%reflects integrity and interconnectedness of a smaller network
toplot_pres = zeros(19);
connectivity_degree = []; 
for i = 1:19
    connectivity_degree = find(binary_p(i,:)==1);
    n = length(connectivity_degree);
    if n>3
    for j = 1:19
        if j ~= i &&  ismember(j,connectivity_degree)
           con2 = find(binary_p(j,:)==1);
           [val,pos]=intersect(connectivity_degree,con2);
           if length(val)>3 && length(val)<8
               toplot_pres(i,j)=1; 
           end
        end     
    end           
    end
end

toplot_unpres = zeros(19);
connectivity_degree = []; 
for i = 1:19
    connectivity_degree = find(binary_unp(i,:)==1);
    n = length(connectivity_degree);
    if n>3
    for j = 1:19
        if j ~= i &&  ismember(j,connectivity_degree)
           con2 = find(binary_p(j,:)==1);
           [val,pos]=intersect(connectivity_degree,con2);
           if length(val)>3 && length(val)<8
               toplot_unpres(i,j)=1; 
           end
        end     
    end           
    end
end

figure(1)
plot_clustercoef(ch_list, toplot_unpres,'unp')
figure(2)
plot_clustercoef(ch_list, toplot_pres,'pres')

%% path length
% probabilities of rewiring
pathlength(binary_p)

probs = logspace(log10(0.0001),log10(1),20);

cp = zeros(size(probs));
lp = zeros(size(probs));

for neti=1:10
    for probi=1:length(probs)
        % rewire
        % find which edges to rewire
        real_edges = find(tril(binary_p));
        real_edges = real_edges(randperm(length(real_edges)));
        edges2rewire = real_edges(1:round(probs(probi)*length(real_edges)));
        
        % rewired connectivity matrix
        connmat_rewired = binary_p;
        
        % loop through edges and change target
        for ei=1:length(edges2rewire)
            
            % find XY coordinates
            [x,y] = ind2sub(size(binary_p),edges2rewire(ei));
            
            % find possible edges to change (cannot already be an edge)
            edges2change = find(~connmat_rewired(x,:));
            
            % rewire
            y2rewire = randsample(edges2change,1);
            connmat_rewired(x,y2rewire) = 1;
            connmat_rewired(y2rewire,x) = 1;
            
            % set original to zero
            connmat_rewired(x,y) = 0;
            connmat_rewired(y,x) = 0;
        end
        
        % mirror matrix
        connmat_rewired = logical(tril(connmat_rewired) + tril(connmat_rewired)');
        
        
        clustcoef_rewired  = zeros(1,19);
        pathlength_rewired = zeros(1,19);
        for chani=1:19
            
            % cluster coefficient
            neighbors = find(connmat_rewired(chani,:));
            n = length(neighbors);
            
            if n>1
                % "local" network of neighbors
                localnetwork = connmat_rewired(neighbors,neighbors);
                % localnetwork is symmetric; remove redundant values by replacing with NaN
                localnetwork = localnetwork + tril(nan(n));
                % compute cluster coefficient (neighbor connectivity scaled)
                clustcoef_rewired(chani) = 2*nansum(localnetwork(:)) / ((n-1)*n);
            end
            
        end
        % average clustering coefficient over channels
        cp(probi) = cp(probi) + nanmean(clustcoef_rewired);
        
        % average path length (remove zeros and Inf's)
        temppathlengths = nonzeros(pathlength(double(connmat_rewired)));
        lp(probi) = lp(probi) + mean(temppathlengths(isfinite(temppathlengths)));
        
        
        % save example networks from select probabilities
        if probi==1
            network1 = connmat_rewired;
        elseif probi==10
            network10 = connmat_rewired;
        elseif probi==length(probs)
            networkend = connmat_rewired;
        end
    end
end


cp = cp./neti;
lp = lp./neti;

figure
plot(probs,cp./cp(1),'-o')
hold on
plot(probs,lp./lp(1),'r-o')
xlabel('Probability of rewiring')
legend({'C_r/C';'L_r/L'})
set(gca,'xlim',[probs(1)/1.5 5],'ylim',[0 5],'xscale','lo')
plot([probs(1) probs(1)],[0 1],'k:')
plot([probs(10) probs(10)],[0 1],'k:')
plot([probs(end) probs(end)],[0 1],'k:')



xylims = [1 19];


figure
subplot(131)
imagesc(network1)
set(gca,'xlim',xylims,'ylim',xylims), axis square, colormap(1-gray)

subplot(132)
imagesc(network10)
set(gca,'xlim',xylims,'ylim',xylims), axis square

subplot(133)
imagesc(networkend)
set(gca,'xlim',xylims,'ylim',xylims), axis square

%% 

%for all
frontocentral_x =x([180, 140, 100, 60,20, 720,680,640,600,540,520,480 ]);
frontocentral_y=y([180, 140, 100, 60,20, 720,680,640,600,540,520,480 ]);
parietooccipital_x = x([  440, 400, 360, 340, 300, 260, 220]);
parietooccipital_y = y([  440, 400, 360, 340, 300, 260, 220]);
% Plot a circle. PRESERVED
figure(20)
angles = linspace(0, 2*pi, 720); % 720 is the total number of points
radius = 20;
xCenter = 50;
yCenter = 40;
x = radius * cos(angles) + xCenter; 
y = radius * sin(angles) + yCenter;
% Plot circle.
plot(x, y, 'k-', 'LineWidth', 1.5);
% Plot center.
hold on;
%plot(xCenter, yCenter, 'k+', 'LineWidth', 2, 'MarkerSize', 16);
grid on;
axis equal;
axis off 
%xlabel('X', 'FontSize', 16);
%ylabel('Y', 'FontSize', 16);
title('Preserved', 'Units', 'normalized', 'Position', [0.5, 1.04, 0],'FontName','times','FontSize',14)

% Now get random locations along the circle.
s1 = 19; % Number of random points to get.
randomIndexes = randperm(length(angles), s1)
xRandom = x(randomIndexes);
yRandom = y(randomIndexes);
%plot(xRandom, yRandom, 'ro', 'LineWidth', 2, 'MarkerSize', 5);

% x_all = x([180, 140, 100, 60, 20, 220, 260, 300, 340, 360, 400, 440, 480, 520, 540, 600, 640, 680, 720 ]);
% y_all = y([180, 140, 100, 60, 20, 220, 260, 300, 340, 360, 400, 440, 480, 520, 540, 600, 640, 680, 720  ]);
% plot(x_all , y_all , 'ro', 'LineWidth', 2, 'MarkerSize', 5);
frontocentral_x1 = x([180, 140, 100, 60]);
frontocentral_y1 = y([180, 140, 100, 60]);

frontocentral_x2 = x([  20, 720,680,640,600,540,520]);
frontocentral_y2 = y([ 20, 720,680, 640,600,540,520 ]);
frontocentral_x3 = x( 480);
frontocentral_y3 = y( 480 );
plot(frontocentral_x2  , frontocentral_y2 , 'go', 'LineWidth', 2, 'MarkerSize', 5);
hold on
plot(frontocentral_x1  , frontocentral_y1 , 'go', 'LineWidth', 2, 'MarkerSize', 5);
hold on 
plot(frontocentral_x3  , frontocentral_y3 , 'go', 'LineWidth', 2, 'MarkerSize', 5);
hold on 

plot(x(540)  , y(540) , 'ro', 'LineWidth', 2, 'MarkerSize', 5);

labels1 = {'Fp1', 'Fp2', 'F7','F3'};
labels2 = { 'Fz', 'F4', 'F8', 'T7', 'C3', 'Cz', 'C4'};
labels3= {'T8'};
text(frontocentral_x1, frontocentral_y1, labels1,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left','FontName', 'times')
text(frontocentral_x2, frontocentral_y2, labels2,'VerticalAlignment', 'top', 'HorizontalAlignment', 'left','FontName', 'times')
text(frontocentral_x3, frontocentral_y3, labels3,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right','FontName', 'times')

hold on
parietooccipital_x1 = x([  440, 400, 360]); %13-19
parietooccipital_y1 = y([ 440, 400, 360 ]);
parietooccipital_x2 = x([   340, 300, 260, 220]); %13-19
parietooccipital_y2 = y([ 340, 300, 260, 220 ]);
labels4 = { 'P7','P3','Pz'};
labels5 = {'P4','P8', 'O1','O2'};

plot(parietooccipital_x1  , parietooccipital_y1 , 'bo', 'LineWidth', 2, 'MarkerSize', 5);
hold on 
text(parietooccipital_x1, parietooccipital_y1, labels4,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right','FontName', 'times')
hold on
plot(parietooccipital_x2  , parietooccipital_y2 , 'bo', 'LineWidth', 2, 'MarkerSize', 5);
hold on 
text(parietooccipital_x2, parietooccipital_y2, labels5,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontName', 'times')
hold on
chan1_p_x = [frontocentral_x(:,1) frontocentral_x(:,12); frontocentral_x(:,1) parietooccipital_x(:,1)];
chan1_p_y =[frontocentral_y(:,1) frontocentral_y(:,12); frontocentral_y(:,1) parietooccipital_y(:,1)];
chan2_p_x =[frontocentral_x(:,2) frontocentral_x(:,4);frontocentral_x(:,2) frontocentral_x(:,10);frontocentral_x(:,2) frontocentral_x(:,12);frontocentral_x(:,2) parietooccipital_x(:,1);frontocentral_x(:,2) parietooccipital_x(:,5);frontocentral_x(:,2) parietooccipital_x(:,6) ];
chan2_p_y =[frontocentral_y(:,2) frontocentral_y(:,4);frontocentral_y(:,2) frontocentral_y(:,10);frontocentral_y(:,2) frontocentral_y(:,12);frontocentral_y(:,2) parietooccipital_y(:,1);frontocentral_y(:,2) parietooccipital_y(:,5);frontocentral_y(:,2) parietooccipital_y(:,6) ];
chan3_p_x =[frontocentral_x(:,3) frontocentral_x(:,4);frontocentral_x(:,3) frontocentral_x(:,9);frontocentral_x(:,3) parietooccipital_x(:,1);frontocentral_x(:,3) parietooccipital_x(:,6);frontocentral_x(:,3) parietooccipital_x(:,7) ];
chan3_p_y =[frontocentral_y(:,3) frontocentral_y(:,4);frontocentral_y(:,3) frontocentral_y(:,9);frontocentral_y(:,3) parietooccipital_y(:,1);frontocentral_y(:,3) parietooccipital_y(:,6);frontocentral_y(:,3) parietooccipital_y(:,7) ];
chan4_p_x =[frontocentral_x(:,4) frontocentral_x(:,9);frontocentral_x(:,4) frontocentral_x(:,10);frontocentral_x(:,4) frontocentral_x(:,12);frontocentral_x(:,4) parietooccipital_x(:,1);frontocentral_x(:,4) parietooccipital_x(:,2);frontocentral_x(:,4) parietooccipital_x(:,3);frontocentral_x(:,4) parietooccipital_x(:,4);frontocentral_x(:,4) parietooccipital_x(:,5);frontocentral_x(:,4) parietooccipital_x(:,7) ];
chan4_p_y =[frontocentral_y(:,4) frontocentral_y(:,9);frontocentral_y(:,4) frontocentral_y(:,10);frontocentral_y(:,4) frontocentral_y(:,12);frontocentral_y(:,4) parietooccipital_y(:,1);frontocentral_y(:,4) parietooccipital_y(:,2);frontocentral_y(:,4) parietooccipital_y(:,3);frontocentral_y(:,4) parietooccipital_y(:,4);frontocentral_y(:,4) parietooccipital_y(:,5);frontocentral_y(:,4) parietooccipital_y(:,7) ];
chan5_p_x =[frontocentral_x(:,5) frontocentral_x(:,9);frontocentral_x(:,5) frontocentral_x(:,10);frontocentral_x(:,5) frontocentral_x(:,12);frontocentral_x(:,5) parietooccipital_x(:,1);frontocentral_x(:,5) parietooccipital_x(:,3);frontocentral_x(:,5) parietooccipital_x(:,4);frontocentral_x(:,5) parietooccipital_x(:,5);frontocentral_x(:,5) parietooccipital_x(:,6);frontocentral_x(:,5) parietooccipital_x(:,7) ];
chan5_p_y =[frontocentral_y(:,5) frontocentral_y(:,9);frontocentral_y(:,5) frontocentral_y(:,10);frontocentral_y(:,5) frontocentral_y(:,12);frontocentral_y(:,5) parietooccipital_y(:,1);frontocentral_y(:,5) parietooccipital_y(:,3);frontocentral_y(:,5) parietooccipital_y(:,4);frontocentral_y(:,5) parietooccipital_y(:,5);frontocentral_y(:,5) parietooccipital_y(:,6);frontocentral_y(:,5) parietooccipital_y(:,7) ];
chan6_p_x =[frontocentral_x(:,6) frontocentral_x(:,7);frontocentral_x(:,6) frontocentral_x(:,9);frontocentral_x(:,6) frontocentral_x(:,10);frontocentral_x(:,6) frontocentral_x(:,12);frontocentral_x(:,6) parietooccipital_x(:,1);frontocentral_x(:,6) parietooccipital_x(:,3);frontocentral_x(:,6) parietooccipital_x(:,5);frontocentral_x(:,6) parietooccipital_x(:,6);frontocentral_x(:,6) parietooccipital_x(:,7) ];
chan6_p_y =[frontocentral_y(:,6) frontocentral_y(:,7);frontocentral_y(:,6) frontocentral_y(:,9);frontocentral_y(:,6) frontocentral_y(:,10);frontocentral_y(:,6) frontocentral_y(:,12);frontocentral_y(:,6) parietooccipital_y(:,1);frontocentral_y(:,6) parietooccipital_y(:,3);frontocentral_y(:,6) parietooccipital_y(:,5);frontocentral_y(:,6) parietooccipital_y(:,6);frontocentral_y(:,6) parietooccipital_y(:,7) ];
chan7_p_x = [frontocentral_x(:,7) parietooccipital_x(:,1); frontocentral_x(:,7) parietooccipital_x(:,5)];
chan7_p_y = [frontocentral_y(:,7) parietooccipital_y(:,1); frontocentral_y(:,7) parietooccipital_y(:,5)];
chan8_p_x = [frontocentral_x(:,8) frontocentral_x(:,10); frontocentral_x(:,8) parietooccipital_x(:,3)];
chan8_p_y = [frontocentral_y(:,8) frontocentral_y(:,10); frontocentral_y(:,8) parietooccipital_y(:,3)];
chan9_p_x =[frontocentral_x(:,9) frontocentral_x(:,10);frontocentral_x(:,9) parietooccipital_x(:,1);frontocentral_x(:,9) parietooccipital_x(:,2);frontocentral_x(:,9) parietooccipital_x(:,3);frontocentral_x(:,9) parietooccipital_x(:,4);frontocentral_x(:,9) parietooccipital_x(:,5);frontocentral_x(:,9) parietooccipital_x(:,6);frontocentral_x(:,9) parietooccipital_x(:,7) ];
chan9_p_y =[frontocentral_y(:,9) frontocentral_y(:,10);frontocentral_y(:,9) parietooccipital_y(:,1);frontocentral_y(:,9) parietooccipital_y(:,2);frontocentral_y(:,9) parietooccipital_y(:,3);frontocentral_y(:,9) parietooccipital_y(:,4);frontocentral_y(:,9) parietooccipital_y(:,5);frontocentral_y(:,9) parietooccipital_y(:,6);frontocentral_y(:,9) parietooccipital_y(:,7) ];
chan10_p_x =[frontocentral_x(:,10) frontocentral_x(:,12);frontocentral_x(:,10) parietooccipital_x(:,1);frontocentral_x(:,10) parietooccipital_x(:,2);frontocentral_x(:,10) parietooccipital_x(:,3);frontocentral_x(:,10) parietooccipital_x(:,5);frontocentral_x(:,10) parietooccipital_x(:,6);frontocentral_x(:,10) parietooccipital_x(:,7) ];
chan10_p_y =[frontocentral_y(:,10) frontocentral_y(:,12);frontocentral_y(:,10) parietooccipital_y(:,1);frontocentral_y(:,10) parietooccipital_y(:,2);frontocentral_y(:,10) parietooccipital_y(:,3);frontocentral_y(:,10) parietooccipital_y(:,5);frontocentral_y(:,10) parietooccipital_y(:,6);frontocentral_y(:,10) parietooccipital_y(:,7) ];
chan12_p_x =[frontocentral_x(:,12) parietooccipital_x(:,1);frontocentral_x(:,12) parietooccipital_x(:,3);frontocentral_x(:,12) parietooccipital_x(:,6);frontocentral_x(:,12) parietooccipital_x(:,7) ];
chan12_p_y =[frontocentral_y(:,12) parietooccipital_y(:,1);frontocentral_y(:,12) parietooccipital_y(:,3);frontocentral_y(:,12) parietooccipital_y(:,6);frontocentral_y(:,12) parietooccipital_y(:,7) ];
chan13_p_x = [parietooccipital_x(:,1) parietooccipital_x(:,2); parietooccipital_x(:,1) parietooccipital_x(:,7)];
chan13_p_y = [parietooccipital_y(:,1) parietooccipital_y(:,2); parietooccipital_y(:,1) parietooccipital_y(:,7)];
chan14_p_x = [parietooccipital_x(:,2) parietooccipital_x(:,6); parietooccipital_x(:,2) parietooccipital_x(:,7)];
chan14_p_y = [parietooccipital_y(:,2) parietooccipital_y(:,6); parietooccipital_y(:,2) parietooccipital_y(:,7)];
chan15_p_x = [parietooccipital_x(:,3) parietooccipital_x(:,5);parietooccipital_x(:,3) parietooccipital_x(:,6); parietooccipital_x(:,3) parietooccipital_x(:,7)];
chan15_p_y = [parietooccipital_y(:,3) parietooccipital_y(:,5);parietooccipital_y(:,3) parietooccipital_y(:,6); parietooccipital_y(:,3) parietooccipital_y(:,7)];
chan16_p_y = [parietooccipital_y(:,4) parietooccipital_y(:,7)];
chan16_p_x = [parietooccipital_x(:,4) parietooccipital_x(:,7)];
chan17_p_x = [parietooccipital_x(:,5) parietooccipital_x(:,7)];
chan17_p_y = [parietooccipital_y(:,5) parietooccipital_y(:,7)];


plot(chan1_p_x(1,:),chan1_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan1_p_x(2,:),chan1_p_y(2,:), 'k-', 'LineWidth', 1)
hold on
plot(chan2_p_x(1,:),chan2_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(2,:),chan2_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan2_p_x(3,:),chan2_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(4,:),chan2_p_y(4,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(5,:),chan2_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(6,:),chan2_p_y(6,:), 'k-', 'LineWidth', 1)
hold on
plot(chan3_p_x(1,:),chan3_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(2,:),chan3_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(3,:),chan3_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(4,:),chan3_p_y(4,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(5,:),chan3_p_y(5,:), 'k-', 'LineWidth', 1)
hold on
plot(chan4_p_x(1,:),chan4_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(2,:),chan4_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan4_p_x(3,:),chan4_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(4,:),chan4_p_y(4,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(5,:),chan4_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(6,:),chan4_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(7,:),chan4_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(8,:),chan4_p_y(8,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(9,:),chan4_p_y(9,:), 'k-', 'LineWidth', 1)
hold on
plot(chan5_p_x(1,:),chan5_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(2,:),chan5_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(3,:),chan5_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(4,:),chan5_p_y(4,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(5,:),chan5_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(6,:),chan5_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(7,:),chan5_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(8,:),chan5_p_y(8,:), 'k-', 'LineWidth', 1)
plot(chan5_p_x(9,:),chan5_p_y(9,:), 'k-', 'LineWidth', 1)
hold on
plot(chan6_p_x(1,:),chan6_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(2,:),chan6_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(3,:),chan6_p_y(3,:), 'r-', 'LineWidth', 1)
plot(chan6_p_x(4,:),chan6_p_y(4,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(5,:),chan6_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(6,:),chan6_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(7,:),chan6_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(8,:),chan6_p_y(8,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(9,:),chan6_p_y(9,:), 'k-', 'LineWidth', 1)
hold on
plot(chan7_p_x(1,:),chan7_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan7_p_x(2,:),chan7_p_y(2,:), 'k-', 'LineWidth', 1)
hold on
plot(chan8_p_x(1,:),chan8_p_y(1,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(2,:),chan8_p_y(2,:), 'k-', 'LineWidth', 1)
hold on
plot(chan9_p_x(1,:),chan9_p_y(1,:), 'r-', 'LineWidth', 1)
plot(chan9_p_x(2,:),chan9_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan9_p_x(3,:),chan9_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan9_p_x(4,:),chan9_p_y(4,:), 'k-', 'LineWidth', 1)
plot(chan9_p_x(5,:),chan9_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan9_p_x(6,:),chan9_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan9_p_x(7,:),chan9_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan9_p_x(8,:),chan9_p_y(8,:), 'k-', 'LineWidth', 1)
hold on
plot(chan10_p_x(1,:),chan10_p_y(1,:), 'r-', 'LineWidth', 1)
plot(chan10_p_x(2,:),chan10_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan10_p_x(3,:),chan10_p_y(3,:), 'r-', 'LineWidth', 1)
plot(chan10_p_x(4,:),chan10_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan10_p_x(5,:),chan10_p_y(5,:), 'r-', 'LineWidth', 1)
plot(chan10_p_x(6,:),chan10_p_y(6,:), 'r-', 'LineWidth', 1)
plot(chan10_p_x(7,:),chan10_p_y(7,:), 'r-', 'LineWidth', 1)
hold on
plot(chan12_p_x(1,:),chan12_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan12_p_x(2,:),chan12_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan12_p_x(3,:),chan12_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan12_p_x(4,:),chan12_p_y(4,:), 'k-', 'LineWidth', 1)
hold on
plot(chan13_p_x(1,:),chan13_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan13_p_x(2,:),chan13_p_y(2,:), 'k-', 'LineWidth', 1)
hold on
plot(chan14_p_x(1,:),chan14_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan14_p_x(2,:),chan14_p_y(2,:), 'k-', 'LineWidth', 1)
hold on
plot(chan15_p_x(1,:),chan15_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan15_p_x(2,:),chan15_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan15_p_x(3,:),chan15_p_y(3,:), 'k-', 'LineWidth', 1)
hold on
plot(chan16_p_x(1,:),chan16_p_y(1,:), 'k-', 'LineWidth', 1)
hold on
plot(chan17_p_x(1,:),chan17_p_y(1,:), 'k-', 'LineWidth', 1)




% Plot a circle. UNPRESERVED
figure(21)
angles = linspace(0, 2*pi, 720); % 720 is the total number of points
radius = 20;
xCenter = 50;
yCenter = 40;
x = radius * cos(angles) + xCenter; 
y = radius * sin(angles) + yCenter;
% Plot circle.
plot(x, y, 'k-', 'LineWidth', 1.5);
% Plot center.
hold on;
%plot(xCenter, yCenter, 'k+', 'LineWidth', 2, 'MarkerSize', 16);
grid on;
axis equal;
%xlabel('X', 'FontSize', 16);
%ylabel('Y', 'FontSize', 16);
title('Unpreserved', 'Units', 'normalized', 'Position', [0.5, 1.04, 0],'FontName','times','FontSize',14)
% Now get random locations along the circle.
%plot(xRandom, yRandom, 'ro', 'LineWidth', 2, 'MarkerSize', 5);

% x_all = x([180, 140, 100, 60, 20, 220, 260, 300, 340, 360, 400, 440, 480, 520, 540, 600, 640, 680, 720 ]);
% y_all = y([180, 140, 100, 60, 20, 220, 260, 300, 340, 360, 400, 440, 480, 520, 540, 600, 640, 680, 720  ]);
% plot(x_all , y_all , 'ro', 'LineWidth', 2, 'MarkerSize', 5);

frontocentral_x1 = x([180, 140, 100, 60]);
frontocentral_y1 = y([180, 140, 100, 60]);

frontocentral_x2 = x([  20, 720,680,640,600,540,520]);
frontocentral_y2 = y([ 20, 720,680, 640,600,540,520 ]);
frontocentral_x3 = x( 480);
frontocentral_y3 = y( 480 );
plot(frontocentral_x2  , frontocentral_y2 , 'go', 'LineWidth', 2, 'MarkerSize', 5);
hold on
plot(frontocentral_x1  , frontocentral_y1 , 'go', 'LineWidth', 2, 'MarkerSize', 5);
hold on 
plot(frontocentral_x3  , frontocentral_y3 , 'go', 'LineWidth', 2, 'MarkerSize', 5);
hold on 

plot(x([20, 640,520])  , y([20, 640,520]) , 'ro', 'LineWidth', 2, 'MarkerSize', 5);

labels1 = {'Fp1', 'Fp2', 'F7','F3'};
labels2 = { 'Fz', 'F4', 'F8', 'T7', 'C3', 'Cz', 'C4'};
labels3= {'T8'};
text(frontocentral_x1, frontocentral_y1, labels1,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left','FontName', 'times')
text(frontocentral_x2, frontocentral_y2, labels2,'VerticalAlignment', 'top', 'HorizontalAlignment', 'left','FontName', 'times')
text(frontocentral_x3, frontocentral_y3, labels3,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right','FontName', 'times')

hold on
parietooccipital_x1 = x([  440, 400, 360]); %13-19
parietooccipital_y1 = y([ 440, 400, 360 ]);
parietooccipital_x2 = x([   340, 300, 260, 220]); %13-19
parietooccipital_y2 = y([ 340, 300, 260, 220 ]);
labels4 = { 'P7','P3','Pz'};
labels5 = {'P4','P8', 'O1','O2'};

plot(parietooccipital_x1  , parietooccipital_y1 , 'bo', 'LineWidth', 2, 'MarkerSize', 5);
hold on 
text(parietooccipital_x1, parietooccipital_y1, labels4,'VerticalAlignment', 'top', 'HorizontalAlignment', 'right','FontName', 'times')
hold on
plot(parietooccipital_x2  , parietooccipital_y2 , 'bo', 'LineWidth', 2, 'MarkerSize', 5);
hold on 
text(parietooccipital_x2, parietooccipital_y2, labels5,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right','FontName', 'times')
hold on
chan1_p_x = [frontocentral_x(:,1) frontocentral_x(:,4);frontocentral_x(:,1) frontocentral_x(:,5);frontocentral_x(:,1) frontocentral_x(:,8); frontocentral_x(:,1) frontocentral_x(:,11);frontocentral_x(:,1) parietooccipital_x(:,3);frontocentral_x(:,1) parietooccipital_x(:,5);frontocentral_x(:,1) parietooccipital_x(:,6);frontocentral_x(:,1) parietooccipital_x(:,7)];
chan1_p_y = [frontocentral_y(:,1) frontocentral_y(:,4);frontocentral_y(:,1) frontocentral_y(:,5);frontocentral_y(:,1) frontocentral_y(:,8); frontocentral_y(:,1) frontocentral_y(:,11);frontocentral_y(:,1) parietooccipital_y(:,3);frontocentral_y(:,1) parietooccipital_y(:,5);frontocentral_y(:,1) parietooccipital_y(:,6);frontocentral_y(:,1) parietooccipital_y(:,7)];
chan2_p_x = [frontocentral_x(:,2) frontocentral_x(:,4);frontocentral_x(:,2) frontocentral_x(:,5);frontocentral_x(:,2) frontocentral_x(:,6);frontocentral_x(:,2) frontocentral_x(:,8);frontocentral_x(:,2) frontocentral_x(:,10); frontocentral_x(:,1) frontocentral_x(:,11);frontocentral_x(:,2) parietooccipital_x(:,5);frontocentral_x(:,2) parietooccipital_x(:,6);frontocentral_x(:,2) parietooccipital_x(:,7)];
chan2_p_y = [frontocentral_y(:,2) frontocentral_y(:,4);frontocentral_y(:,2) frontocentral_y(:,5);frontocentral_y(:,2) frontocentral_y(:,6);frontocentral_y(:,2) frontocentral_y(:,8);frontocentral_y(:,2) frontocentral_y(:,10); frontocentral_y(:,1) frontocentral_y(:,11);frontocentral_y(:,2) parietooccipital_y(:,5);frontocentral_y(:,2) parietooccipital_y(:,6);frontocentral_y(:,2) parietooccipital_y(:,7)];
chan3_p_x = [frontocentral_x(:,3) frontocentral_x(:,4);frontocentral_x(:,3) frontocentral_x(:,5);frontocentral_x(:,3) frontocentral_x(:,6);frontocentral_x(:,3) frontocentral_x(:,8);frontocentral_x(:,3) frontocentral_x(:,11);frontocentral_x(:,3) parietooccipital_x(:,1);frontocentral_x(:,3) parietooccipital_x(:,5);frontocentral_x(:,3) parietooccipital_x(:,6);frontocentral_x(:,3) parietooccipital_x(:,7)];
chan3_p_y = [frontocentral_y(:,3) frontocentral_y(:,4);frontocentral_y(:,3) frontocentral_y(:,5);frontocentral_y(:,3) frontocentral_y(:,6);frontocentral_y(:,3) frontocentral_y(:,8);frontocentral_y(:,3) frontocentral_y(:,11);frontocentral_y(:,3) parietooccipital_y(:,1);frontocentral_y(:,3) parietooccipital_y(:,5);frontocentral_y(:,3) parietooccipital_y(:,6);frontocentral_y(:,3) parietooccipital_y(:,7)];
chan4_p_x = [frontocentral_x(:,4) frontocentral_x(:,6);frontocentral_x(:,4) frontocentral_x(:,8);frontocentral_x(:,4) frontocentral_x(:,10);frontocentral_x(:,4) frontocentral_x(:,11);frontocentral_x(:,4) parietooccipital_x(:,1);frontocentral_x(:,4) parietooccipital_x(:,2);frontocentral_x(:,4) parietooccipital_x(:,3);frontocentral_x(:,4) parietooccipital_x(:,4);frontocentral_x(:,4) parietooccipital_x(:,6)];
chan4_p_y = [frontocentral_y(:,4) frontocentral_y(:,6);frontocentral_y(:,4) frontocentral_y(:,8);frontocentral_y(:,4) frontocentral_y(:,10);frontocentral_y(:,4) frontocentral_y(:,11);frontocentral_y(:,4) parietooccipital_y(:,1);frontocentral_y(:,4) parietooccipital_y(:,2);frontocentral_y(:,4) parietooccipital_y(:,3);frontocentral_y(:,4) parietooccipital_y(:,4);frontocentral_y(:,4) parietooccipital_y(:,6)];
chan5_p_x = [frontocentral_x(:,5) frontocentral_x(:,6);frontocentral_x(:,5) frontocentral_x(:,8);frontocentral_x(:,5) frontocentral_x(:,9);frontocentral_x(:,5) frontocentral_x(:,10);frontocentral_x(:,5) frontocentral_x(:,11);frontocentral_x(:,5) parietooccipital_x(:,1);frontocentral_x(:,5) parietooccipital_x(:,2);frontocentral_x(:,5) parietooccipital_x(:,3);frontocentral_x(:,5) parietooccipital_x(:,4);frontocentral_x(:,5) parietooccipital_x(:,6)];
chan5_p_y = [frontocentral_y(:,5) frontocentral_y(:,6);frontocentral_y(:,5) frontocentral_y(:,8);frontocentral_y(:,5) frontocentral_y(:,9);frontocentral_y(:,5) frontocentral_y(:,10);frontocentral_y(:,5) frontocentral_y(:,11);frontocentral_y(:,5) parietooccipital_y(:,1);frontocentral_y(:,5) parietooccipital_y(:,2);frontocentral_y(:,5) parietooccipital_y(:,3);frontocentral_y(:,5) parietooccipital_y(:,4);frontocentral_y(:,5) parietooccipital_y(:,6)];
chan6_p_x = [frontocentral_x(:,6) frontocentral_x(:,8);frontocentral_x(:,6) frontocentral_x(:,9);frontocentral_x(:,6) frontocentral_x(:,10);frontocentral_x(:,6) frontocentral_x(:,11);frontocentral_x(:,6) parietooccipital_x(:,3);frontocentral_x(:,6) parietooccipital_x(:,4);frontocentral_x(:,6) parietooccipital_x(:,6);frontocentral_x(:,6) parietooccipital_x(:,7)];
chan6_p_y = [frontocentral_y(:,6) frontocentral_y(:,8);frontocentral_y(:,6) frontocentral_y(:,9);frontocentral_y(:,6) frontocentral_y(:,10);frontocentral_y(:,6) frontocentral_y(:,11);frontocentral_y(:,6) parietooccipital_y(:,3);frontocentral_y(:,6) parietooccipital_y(:,4);frontocentral_y(:,6) parietooccipital_y(:,6);frontocentral_y(:,6) parietooccipital_y(:,7)];
chan7_p_x = [frontocentral_x(:,7) frontocentral_x(:,11)];
chan7_p_y = [frontocentral_y(:,7) frontocentral_y(:,11)];
chan8_p_x = [frontocentral_x(:,8) frontocentral_x(:,9);frontocentral_x(:,8) frontocentral_x(:,10);frontocentral_x(:,8) frontocentral_x(:,11);frontocentral_x(:,8) parietooccipital_x(:,1);frontocentral_x(:,8) parietooccipital_x(:,3);frontocentral_x(:,8) parietooccipital_x(:,4);frontocentral_x(:,8) parietooccipital_x(:,5);frontocentral_x(:,8) parietooccipital_x(:,6);frontocentral_x(:,8) parietooccipital_x(:,7)];
chan8_p_y = [frontocentral_y(:,8) frontocentral_y(:,9);frontocentral_y(:,8) frontocentral_y(:,10);frontocentral_y(:,8) frontocentral_y(:,11);frontocentral_y(:,8) parietooccipital_y(:,1);frontocentral_y(:,8) parietooccipital_y(:,3);frontocentral_y(:,8) parietooccipital_y(:,4);frontocentral_y(:,8) parietooccipital_y(:,5);frontocentral_y(:,8) parietooccipital_y(:,6);frontocentral_y(:,8) parietooccipital_y(:,7)];
chan10_p_x = [frontocentral_x(:,10) frontocentral_x(:,12);frontocentral_x(:,10) parietooccipital_x(:,2);frontocentral_x(:,10) parietooccipital_x(:,5);frontocentral_x(:,10) parietooccipital_x(:,7)];
chan10_p_y = [frontocentral_y(:,10) frontocentral_y(:,12);frontocentral_y(:,10) parietooccipital_y(:,2);frontocentral_y(:,10) parietooccipital_y(:,5);frontocentral_y(:,10) parietooccipital_y(:,7)];
chan11_p_x = [frontocentral_x(:,11) frontocentral_x(:,12);frontocentral_x(:,11) parietooccipital_x(:,1);frontocentral_x(:,11) parietooccipital_x(:,2);frontocentral_x(:,11) parietooccipital_x(:,3);frontocentral_x(:,11) parietooccipital_x(:,5);frontocentral_x(:,11) parietooccipital_x(:,6);frontocentral_x(:,11) parietooccipital_x(:,7)];
chan11_p_y = [frontocentral_y(:,11) frontocentral_y(:,12);frontocentral_y(:,11) parietooccipital_y(:,1);frontocentral_y(:,11) parietooccipital_y(:,2);frontocentral_y(:,11) parietooccipital_y(:,3);frontocentral_y(:,11) parietooccipital_y(:,5);frontocentral_y(:,11) parietooccipital_y(:,6);frontocentral_y(:,11) parietooccipital_y(:,7)];
chan13_p_x = [parietooccipital_x(:,1) parietooccipital_x(:,3);parietooccipital_x(:,1) parietooccipital_x(:,7)];
chan13_p_y = [parietooccipital_y(:,1) parietooccipital_y(:,3);parietooccipital_y(:,1) parietooccipital_y(:,7)];
chan14_p_x = [parietooccipital_x(:,2) parietooccipital_x(:,3);parietooccipital_x(:,2) parietooccipital_x(:,4);parietooccipital_x(:,2) parietooccipital_x(:,6);parietooccipital_x(:,2) parietooccipital_x(:,7)];
chan14_p_y = [parietooccipital_y(:,2) parietooccipital_y(:,3);parietooccipital_y(:,2) parietooccipital_y(:,4);parietooccipital_y(:,2) parietooccipital_y(:,6);parietooccipital_y(:,2) parietooccipital_y(:,7)];
chan15_p_x = [parietooccipital_x(:,3) parietooccipital_x(:,5);parietooccipital_x(:,3) parietooccipital_x(:,6);parietooccipital_x(:,3) parietooccipital_x(:,7)];
chan15_p_y = [parietooccipital_y(:,3) parietooccipital_y(:,5);parietooccipital_y(:,3) parietooccipital_y(:,6);parietooccipital_y(:,3) parietooccipital_y(:,7)];
chan16_p_x = [parietooccipital_x(:,4) parietooccipital_x(:,7)];
chan16_p_y = [parietooccipital_y(:,4) parietooccipital_y(:,7)];
chan18_p_x = [parietooccipital_x(:,6) parietooccipital_x(:,7)];
chan18_p_y = [parietooccipital_y(:,6) parietooccipital_y(:,7)];


plot(chan1_p_x(1,:),chan1_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan1_p_x(2,:),chan1_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan1_p_x(3,:),chan1_p_y(3,:), 'r-', 'LineWidth', 1)
plot(chan1_p_x(4,:),chan1_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan1_p_x(5,:),chan1_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan1_p_x(6,:),chan1_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan1_p_x(7,:),chan1_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan1_p_x(8,:),chan1_p_y(8,:), 'k-', 'LineWidth', 1)
hold on
plot(chan2_p_x(1,:),chan2_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(2,:),chan2_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan2_p_x(3,:),chan2_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(4,:),chan2_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan2_p_x(5,:),chan2_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(6,:),chan2_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(7,:),chan2_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(8,:),chan2_p_y(8,:), 'k-', 'LineWidth', 1)
plot(chan2_p_x(9,:),chan2_p_y(9,:), 'k-', 'LineWidth', 1)
hold on
plot(chan3_p_x(1,:),chan3_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(2,:),chan3_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan3_p_x(3,:),chan3_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(4,:),chan3_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan3_p_x(5,:),chan3_p_y(5,:), 'r-', 'LineWidth', 1)
plot(chan3_p_x(6,:),chan3_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(7,:),chan3_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(8,:),chan3_p_y(8,:), 'k-', 'LineWidth', 1)
plot(chan3_p_x(9,:),chan3_p_y(9,:), 'k-', 'LineWidth', 1)
hold on
plot(chan4_p_x(1,:),chan4_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(2,:),chan4_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan4_p_x(3,:),chan4_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(4,:),chan4_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan4_p_x(5,:),chan4_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(6,:),chan4_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(7,:),chan4_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(8,:),chan4_p_y(8,:), 'k-', 'LineWidth', 1)
plot(chan4_p_x(9,:),chan4_p_y(9,:), 'k-', 'LineWidth', 1)
hold on
plot(chan5_p_x(1,:),chan5_p_y(1,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(2,:),chan5_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(3,:),chan5_p_y(3,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(4,:),chan5_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(5,:),chan5_p_y(5,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(6,:),chan5_p_y(6,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(7,:),chan5_p_y(7,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(8,:),chan5_p_y(8,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(9,:),chan5_p_y(9,:), 'r-', 'LineWidth', 1)
plot(chan5_p_x(10,:),chan5_p_y(10,:), 'r-', 'LineWidth', 1)
hold on
plot(chan6_p_x(1,:),chan6_p_y(1,:), 'r-', 'LineWidth', 1)
plot(chan6_p_x(2,:),chan6_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(3,:),chan6_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(4,:),chan6_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan6_p_x(5,:),chan6_p_y(5,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(6,:),chan6_p_y(6,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(7,:),chan6_p_y(7,:), 'k-', 'LineWidth', 1)
plot(chan6_p_x(8,:),chan6_p_y(8,:), 'k-', 'LineWidth', 1)
hold on
plot(chan7_p_x(1,:),chan7_p_y(1,:), 'r-', 'LineWidth', 1)
hold on
plot(chan8_p_x(1,:),chan8_p_y(1,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(2,:),chan8_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(3,:),chan8_p_y(3,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(4,:),chan8_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(5,:),chan8_p_y(5,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(6,:),chan8_p_y(6,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(7,:),chan8_p_y(7,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(8,:),chan8_p_y(8,:), 'r-', 'LineWidth', 1)
plot(chan8_p_x(9,:),chan8_p_y(9,:), 'r-', 'LineWidth', 1)
hold on
plot(chan10_p_x(1,:),chan10_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan10_p_x(2,:),chan10_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan10_p_x(3,:),chan10_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan10_p_x(4,:),chan10_p_y(4,:), 'k-', 'LineWidth', 1)
hold on
plot(chan11_p_x(1,:),chan11_p_y(1,:), 'r-', 'LineWidth', 1)
plot(chan11_p_x(2,:),chan11_p_y(2,:), 'r-', 'LineWidth', 1)
plot(chan11_p_x(3,:),chan11_p_y(3,:), 'r-', 'LineWidth', 1)
plot(chan11_p_x(4,:),chan11_p_y(4,:), 'r-', 'LineWidth', 1)
plot(chan11_p_x(5,:),chan11_p_y(5,:), 'r-', 'LineWidth', 1)
plot(chan11_p_x(6,:),chan11_p_y(6,:), 'r-', 'LineWidth', 1)
plot(chan11_p_x(7,:),chan11_p_y(7,:), 'r-', 'LineWidth', 1)
hold on
plot(chan13_p_x(1,:),chan13_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan13_p_x(2,:),chan13_p_y(2,:), 'k-', 'LineWidth', 1)
hold on
plot(chan14_p_x(1,:),chan14_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan14_p_x(2,:),chan14_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan14_p_x(3,:),chan14_p_y(3,:), 'k-', 'LineWidth', 1)
plot(chan14_p_x(4,:),chan14_p_y(4,:), 'k-', 'LineWidth', 1)
hold on
plot(chan15_p_x(1,:),chan15_p_y(1,:), 'k-', 'LineWidth', 1)
plot(chan15_p_x(2,:),chan15_p_y(2,:), 'k-', 'LineWidth', 1)
plot(chan15_p_x(3,:),chan15_p_y(3,:), 'k-', 'LineWidth', 1)
hold on
plot(chan16_p_x(1,:),chan16_p_y(1,:), 'k-', 'LineWidth', 1)
hold on
plot(chan18_p_x(1,:),chan18_p_y(1,:), 'k-', 'LineWidth', 1)
hold on
axis off


%% RT vs PLI 
%from average PLI per patient find one value averaging all channels (but
%zeros)
aver_PLI = []; 
for k = 1:10
    if k ==1
       aver_PLI(k,:) =mean(nonzeros(subj3_PLI));
    elseif k ==2
       aver_PLI(k,:) = mean(nonzeros(subj5_PLI));
    elseif k ==3
       aver_PLI(k,:) = mean(nonzeros(subj6_PLI));
    elseif k ==4
       aver_PLI(k,:) = mean(nonzeros(subj7_PLI));
    elseif k ==5 
       aver_PLI(k,:) = mean(nonzeros(subj8_PLI));
    elseif k ==6
       aver_PLI(k,:) = mean(nonzeros(subj10_PLI));
    elseif k ==7
       aver_PLI(k,:) = mean(nonzeros(subj12_PLI));
    elseif k ==8
       aver_PLI(k,:) = mean(nonzeros(subj13_PLI));
    elseif k ==9
       aver_PLI(k,:) = mean(nonzeros(subj14_PLI));
    elseif k ==10
       aver_PLI(k,:) = mean(nonzeros(subj15_PLI));   
    end
end


% make vector of RT and accuracy per subject 
RTs=[];
Accus = []; 
for j = 1:10
if j ==1 
   subject_nr = 3;
elseif j==2   
   subject_nr = 5; 
elseif j==3  
   subject_nr = 6;
elseif j==4   
   subject_nr = 7;  
elseif j==5   
   subject_nr = 8;    
elseif j==6   
   subject_nr = 10; 
elseif j==7   
   subject_nr = 12;   
elseif j==8   
   subject_nr = 13; 
elseif j==9   
   subject_nr = 14;      
elseif j==10   
   subject_nr = 15;    
end   

if subject_nr <10
   table = readtable(strcat('C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\OSfiles-Tobiifiles-SAGA\tobiinoheader\monkey\trials\output\Output_patient0',string(subject_nr),'.xlsx'));
else
   table = readtable(strcat('C:\Users\vb\Documents\KempHemProject\Eyetracking\EyeTrackingAnalysis\Patients_kemp\OSfiles-Tobiifiles-SAGA\tobiinoheader\monkey\trials\output\Output_patient',string(subject_nr),'.xlsx'));    
end
RT = nonzeros(table.TotalReactionTime_ms_);
RTs(j,:)=mean(RT(RT~=2166.7)); %remove trials with no response 
Accu= table.TotalAccuracy___;
Accus(j,:)= Accu(end);
end


figure(47)
scatter(RTs, aver_PLI,10,'k', 'filled')
title('Rt vs PLI')
xlabel('RT (ms)')
ylabel('PLI')
xlim([600 1500])
ylim([0 1])
hold on
Fit1 = polyfit(RTs,aver_PLI,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit1,RTs);
y1 = plot(RTs,y_est,'color','k','LineWidth',2);


figure(48)
scatter(Accus, aver_PLI,10,'k', 'filled')
title('Correct responses vs PLI')
xlabel('Correct responses (%)')
ylabel('PLI')
xlim([75 105])
ylim([0 1])
hold on
Fit1 = polyfit(Accus,aver_PLI,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
y_est = polyval(Fit1,Accus);
y1 = plot(Accus,y_est,'color','k','LineWidth',2);
[RHO,PVAL] = corr(aver_PLI,Accus,'Type','Pearson');

figure(48) 
eyes_imp=[0,1,0,1,1,0,0,1,0,1];
scatter(eyes_imp, aver_PLI,10,'k', 'filled')

%no clear differences of PLI in two groups
