  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Imports EEG set file of each patient. extract dipoles for each seizure
  % for each patient. Dipoles first allowed both within and outside
  % headmap, and then within headmap only. 
  % 3 output files are saved: 
  % - dipoles within and outside headmap, with smallest RV (residual
  % variance)
  % - index of dipoles outside headmap
  % - dipoles only inside headmap
  %
  %  
  % needs fieldtrip
  % needs standard_bem.mat
  % needs standard_1005.elc
  
  % VBarone Nov, 2022
  % 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; close 

addpath C:\Users\vb\Documents\MATLAB\R2018b\fieldtrip-20220714\
ft_defaults
% load patients 

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

%% define electrodes
elecData = ft_read_sens(EEGfilename);

%% load volume conductor
load('C:\Users\vb\Documents\MATLAB\R2018b\fieldtrip-20220714\template\headmodel\standard_bem.mat') % variable 'vol'

%% load sens from template
elec1005 = ft_read_sens('C:\Users\vb\Documents\MATLAB\R2018b\fieldtrip-20220714\template\electrode\standard_1005.elc');
%%
idx_label = ismember(elec1005.label,eeg_raw.label); %strcmp(elec.label,eeg_raw.label);

elec = [];
elec.chanpos = elec1005.chanpos(idx_label,:);
elec.chantype = elec1005.chantype(idx_label);
elec.chanunit = elec1005.chanunit(idx_label);
elec.elecpos = elec1005.elecpos(idx_label,:);
elec.label = elec1005.label(idx_label);
elec.type = elec1005.type;
elec.unit = elec1005.unit;

%figure, ft_plot_mesh(vol.bnd(1), 'edgecolor', 'none', 'facecolor', 'skin', 'facealpha', 0.3), hold on, 
%ft_plot_sens(elec, 'label', 'yes') %looks ok!
% ft_plot_sens(elec)

%% dipole fitting for the whole trial, 
% look within grid + find dipoles outside head

dip1 = {};
for kk=1:size(trl_startstop,1)
for ii= 1:length(data_trl{1,kk}{1,1}.trial{1})
cfg = [];
cfg.latency = ii/fs; %[0.080 0.110];
cfg.numdipoles = 1;
cfg.gridsearch   = 'yes';                           %only dipole scan
cfg.nonlinear    = 'yes';
cfg.resolution = 10;
cfg.unit = 'mm';
cfg.gridsearch = 'yes';
cfg.headmodel = vol;
cfg.senstype = 'eeg';
cfg.channel = 'all';
source_eeg = ft_dipolefitting(cfg, data_trl{1,kk}{1,1});
dip1{kk,ii} = source_eeg.dip;
end
end
save(strcat('dipfitFT_seizureallasamples_nonlinyesgridsearchyes_subj',string(subject_nr),'.mat'),'dip1')

dippos1 = {};
rvalue1 ={};
for gg=1:size(begin_seizure,2)
  for nn = 1:(end_seizure(:,gg)-begin_seizure(:,gg))+1 
      dippos1{gg}(nn,:) = dip1{gg,nn}.pos;
      rvalue1{gg}(nn,:) = dip1{gg,nn}.rv;
  end
end


% figure(7), ft_plot_mesh(vol.bnd(3), 'edgecolor', 'none', 'facecolor', 'skin', 'facealpha', 0.8), hold on, 
% c=  linspace(0,length(dip)/256,size(dip,2)); %color bar between 0 and length in seconds 
% scatter3(dippos_out(:,1),dippos_out(:,2),dippos_out(:,3), rvvalue_out*100,  c,'filled'), ft_colormap('gray')
% colorbar


%% dipoles outside head
inside_dipoles = {};
for kk=1:size(dip1,1) 
[inside_dip] = ft_inside_headmodel(dippos1{1,kk}, vol); 
inside_dipoles{kk} = inside_dip;
end
save(strcat('insidedip_seizureallasamples_nonlinyesgridsearchyes_subj',string(subject_nr),'.mat'),'inside_dipoles')


%% take dipoles outside headmap only, and rerun dipfit with nonlinear search (so do not take dipoles outside headmap)
data_trl_outside = {};
for nn=1:size(inside_dipoles,2) 
    outside_map_samples = []; 
    outside_map_samples = find(inside_dipoles{1,nn}==0);
   % outside_map_samples = randi(56,1,103)';
    
    data_outside = [];
    for t= 1:size(outside_map_samples,1)
    if t == 1
       data_outside = data_trl{1,nn}{1,1}.trial{1,1}(:, outside_map_samples(t):outside_map_samples(t));
    else 
       data_outside = [data_outside data_trl{1,nn}{1,1}.trial{1,1}(:, outside_map_samples(t):outside_map_samples(t))] ;  
    end
    end
    
    time_out = (0:length(data_outside)-1)/ fs; %in seconds
    data_trl_outside{1,nn} = data_trl{1,nn};
    data_trl_outside{1,nn}{1,1}.trial{1,1} = data_outside;
    data_trl_outside{1,nn}{1,1}.time{1,1} = time_out;
    data_trl_outside{1,nn}{1,1}.hdr.nSamples = size(data_outside,2);

end

dip_out = {};
for kk=1:size(data_trl_outside,2)
for ii= 1:length(data_trl_outside{1,kk}{1,1}.trial{1})
cfg = [];
cfg.latency = ii/fs; %[0.080 0.110];
cfg.numdipoles = 1;
cfg.gridsearch   = 'yes';                           %only dipole scan
cfg.nonlinear    = 'no';
cfg.resolution = 10;
cfg.unit = 'mm';
cfg.gridsearch = 'yes';
cfg.headmodel = vol;
cfg.senstype = 'eeg';
cfg.channel = 'all';
source_eeg = ft_dipolefitting(cfg, data_trl_outside{1,kk}{1,1});
dip_out{kk,ii} = source_eeg.dip;
end
end
save(strcat('dipfitFT_samplesout_nonlinnogridsearchyes_subj',string(subject_nr),'.mat'),'dip_out')

end

