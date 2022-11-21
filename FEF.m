%% load created files from dipole_fit

clear;clc
load('C:\Users\vb\Documents\MATLAB\AbsenceChildren_paper4\scripts\dipfit_map\results_dipfitFT_allpat\pos_unpres.mat');
load('C:\Users\vb\Documents\MATLAB\AbsenceChildren_paper4\scripts\dipfit_map\results_dipfitFT_allpat\pos_pres.mat');

%define matrices with number of samples of each seizure per patient
label_unp_subj=[4353,4865,8052, 5851,0,0,0,0,0,0; 769,1665,1793,1511,1793,2049,1972,2305, 1434,2561; 2817,3649, 3738,5761,3918,8193,5121 0,0,0;  1076,1294, 0,0,0,0,0,0,0,0; 1281, 1281,1537,1281,1537,1281,2049, 0,0,0];
label_p_subj=[1832,1468,0, 0,0,0,0,0,0,0,0,0,0,0; 1537,1844,0, 0,0,0,0,0,0,0,0,0,0,0;3841,999,3841,4353,769,0,0,0,0,0,0,0,0,0;1434,3073,2330,1972,1947,2100,1806,1576,1665,1793,1025,2331,2049,2561;2049,0,0,0,0,0,0,0,0,0,0,0,0,0];

%% determine FEF coordinates from Vernet et al., 2014
FEF1_L = [-32 -2 46]; FEF1_R =[31 -2 47];
FEF2_L= [-35 -18 46]; FEF2_R= [36 -10 47];
FEF3_L= [-30 -7 49]; FEF3_R=[34 -3 47];
FEF4_L=[-37 26 29]; FEF4_R=[37 26 29];
FEF5_L=[-41 12 34]; FEF5_R=[32 10 34];

%from TALAIRACH to MNI
FEF1_L = tal2mni(FEF1_L); FEF1_R = tal2mni(FEF1_R);
FEF2_L = tal2mni(FEF2_L); FEF2_R = tal2mni(FEF2_R);
FEF3_L = tal2mni(FEF3_L); FEF3_R = tal2mni(FEF3_R);
FEF4_L = tal2mni(FEF4_L); FEF4_R = tal2mni(FEF4_R);
FEF5_L = tal2mni(FEF5_L); FEF5_R = tal2mni(FEF5_R);



%% dippos for each subject of the two groups
%unpres
dippos_subj5{1,1} = dippos_unpres(1:label_unp_subj(1,1),:);dippos_subj5{1,2} = dippos_unpres((label_unp_subj(1,1)+1):(label_unp_subj(1,1)+label_unp_subj(1,2)),:);dippos_subj5{1,3} = dippos_unpres((label_unp_subj(1,2)+label_unp_subj(1,1)+1):(label_unp_subj(1,1)+label_unp_subj(1,2)+label_unp_subj(1,3)),:);dippos_subj5{1,4} = dippos_unpres((sum(label_unp_subj(1,1:3))+1):(sum(label_unp_subj(1,1:4))),:);
dippos_subj7{1,1} = dippos_unpres((23121+1):(23121+label_unp_subj(2,1)),:);dippos_subj7{1,2} = dippos_unpres((23890+1):(23890+label_unp_subj(2,2)),:);dippos_subj7{1,3} = dippos_unpres((25555+1):(25555+label_unp_subj(2,3)),:);dippos_subj7{1,4} =dippos_unpres((27348+1):(27348+label_unp_subj(2,4)),:);dippos_subj7{1,5} = dippos_unpres((28859+1):(28859+label_unp_subj(2,5)),:);dippos_subj7{1,6} = dippos_unpres((30652+1):(30652+label_unp_subj(2,6)),:);dippos_subj7{1,7} = dippos_unpres((32701+1):(32701+label_unp_subj(2,7)),:);dippos_subj7{1,8} = dippos_unpres((34673+1):(34673+label_unp_subj(2,8)),:);dippos_subj7{1,9} = dippos_unpres((36978+1):(36978+label_unp_subj(2,9)),:);dippos_subj7{1,10} = dippos_unpres((38412+1):(38412+label_unp_subj(2,10)),:);
dippos_subj8{1,1} = dippos_unpres((40973+1):(40973+label_unp_subj(3,1)),:);dippos_subj8{1,2} = dippos_unpres((43790+1):(43790+label_unp_subj(3,2)),:);dippos_subj8{1,3} = dippos_unpres((47439+1):(47439+label_unp_subj(3,3)),:);dippos_subj8{1,4} = dippos_unpres((51177+1):(51177+label_unp_subj(3,4)),:);dippos_subj8{1,5} = dippos_unpres((56938+1):(56938+label_unp_subj(3,5)),:);dippos_subj8{1,6} = dippos_unpres(( 60856+1):( 60856+label_unp_subj(3,6)),:);dippos_subj8{1,7} = dippos_unpres(( 69049+1):( 69049+label_unp_subj(3,7)),:);
dippos_subj13{1,1} =dippos_unpres(( 74170+1):( 74170+label_unp_subj(4,1)),:);dippos_subj13{1,2} = dippos_unpres(( 75246+1):( 75246+label_unp_subj(4,2)),:);
dippos_subj15{1,1} = dippos_unpres(( 76540+1):( 76540+label_unp_subj(5,1)),:);dippos_subj15{1,2} = dippos_unpres(( 77821+1):( 77821+label_unp_subj(5,2)),:);dippos_subj15{1,3} = dippos_unpres(( 79102+1):( 79102+label_unp_subj(5,3)),:);dippos_subj15{1,4} =dippos_unpres(( 80639+1):( 80639+label_unp_subj(5,4)),:);dippos_subj15{1,5} =dippos_unpres((  81920+1):(  81920+label_unp_subj(5,5)),:);dippos_subj15{1,6} = dippos_unpres((834570+1):( 83457+label_unp_subj(5,6)),:);dippos_subj15{1,7} = dippos_unpres((84738+1):( 84738+label_unp_subj(5,7)),:);

%pres
dippos_subj3{1,1} = dippos_pres(1:label_p_subj(1,1),:);dippos_subj3{1,2} = dippos_pres((label_p_subj(1,1)+1):(label_p_subj(1,1)+label_p_subj(1,2)),:);
dippos_subj6{1,1} = dippos_pres(sum(label_p_subj(1,1:2))+1:sum(label_p_subj(1,1:2))+(label_p_subj(2,1)),:);dippos_subj6{1,2} = dippos_pres(sum(label_p_subj(1,1:2))+(label_p_subj(2,1)+1):sum(label_p_subj(1,1:2))+sum(label_p_subj(2,1:2)),:);
dippos_subj10{1,1} = dippos_pres(6682:6681+label_p_subj(3,1),:);dippos_subj10{1,2} = dippos_pres((6681+label_p_subj(3,1)+1):(6681+label_p_subj(3,1)+label_p_subj(3,2)),:);dippos_subj10{1,3} = dippos_pres((11521+1):(11521+label_p_subj(3,3)),:);dippos_subj10{1,4} = dippos_pres((15362+1):(15362+label_p_subj(3,4)),:);dippos_subj10{1,5} = dippos_pres((19715+1):(19715+label_p_subj(3,5)),:);
dippos_subj12{1,1} = dippos_pres((20484+1:20484+label_p_subj(4,1)),:);dippos_subj12{1,2} = dippos_pres((21918+1):(21918+label_p_subj(4,2)),:);dippos_subj12{1,3} = dippos_pres((24991+1):(24991+label_p_subj(4,3)),:);dippos_subj12{1,4} = dippos_pres((27321+1):(27321+label_p_subj(4,4)),:);dippos_subj12{1,5} = dippos_pres((29293+1):(29293+label_p_subj(4,5)),:);dippos_subj12{1,6} = dippos_pres((31240+1):(31240+label_p_subj(4,6)),:);dippos_subj12{1,7} = dippos_pres((33340+1):(33340+label_p_subj(4,7)),:);dippos_subj12{1,8} = dippos_pres((35146+1):(35146+label_p_subj(4,8)),:);dippos_subj12{1,9} = dippos_pres((36722+1):(36722+label_p_subj(4,9)),:);dippos_subj12{1,10} = dippos_pres((38387+1):(38387+label_p_subj(4,10)),:);dippos_subj12{1,11} = dippos_pres((40180+1):(40180+label_p_subj(4,11)),:);dippos_subj12{1,12} = dippos_pres((41205+1):(41205+label_p_subj(4,12)),:);dippos_subj12{1,13} = dippos_pres((43536+1):(43536+label_p_subj(4,13)),:);dippos_subj12{1,14} = dippos_pres((45585+1):(45585+label_p_subj(4,14)),:);
dippos_subj14{1,1} = dippos_pres((48146+1):(48146+label_p_subj(5,1)),:);

%% make spheres from FEFs
%define sphere with radius 10
[X,Y,Z] = sphere();
r = 10;
X2 = (X * r);
Y2 = (Y * r);
Z2 = (Z * r);


%coordinates to check 
% xp=dippos_subj5{1,1}(:,1);
% yp=dippos_subj5{1,1}(:,2);
% zp = dippos_subj5{1,1}(:,3);

%is within sphere?
% insph = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;

% figure
% mesh(X,Y,Z,'FaceAlpha',0.5)
% hold on
% plot3(xp(insph), yp(insph), zp(insph), '.r')
% plot3(xp(~insph), yp(~insph), zp(~insph), '.b')
% hold off
% axis('equal')
% 
% load('C:\Users\vb\Documents\MATLAB\R2018b\fieldtrip-20220714\template\headmodel\standard_bem.mat') % variable 'vol'
% figure(21), ft_plot_mesh(vol.bnd(3), 'edgecolor', 'none', 'facecolor', 'skin', 'facealpha', 0.3), hold on, 
% % colorbar
% hold on
% hs1 = surf(X2+FEF1_L(1),Y2+FEF1_L(2),Z2+FEF1_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [1 0 0])
% hs1 = surf(X2+FEF1_R(1),Y2+FEF1_R(2),Z2+FEF1_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [1 0 0])
% %hs1 =mesh(X2+FEF1_L(1),Y2+FEF1_L(2),Z2+FEF1_L(3),'FaceAlpha',0.5)
% plot3(xp(insph), yp(insph), zp(insph), '.b','MarkerSize',30)
%surf(X2+5,Y2-5,Z2)  %coord (5,-5,2)

%%
%num of dipoles in FEF unpreserved 
for seiz = 1:length(dippos_subj5) 
    g=1;
    g1=1;
    g2 = 1;
    for t = 1:length(dippos_subj5{1,seiz})
        xp= dippos_subj5{1,seiz}(t,1);
        yp= dippos_subj5{1,seiz}(t,2);
        zp= dippos_subj5{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj5{1,seiz}(g,1) =  t;
           all_FEF_subj5{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj5{1,seiz}(g1,1) =  t;
           all_FEF_subj5{1,seiz}(g2,1)=t;
           g1=g1+1;
           g2=g2+1;
        end    
    end    
end
%percentage of dipoles in FEF
for k = 1:size(right_FEF_subj5,2)
    %LEFT col 1, RIGHT col 2, col 3 all 
           perc_FEF_subj5(k,2) = (length(right_FEF_subj5{1,k})/length(dippos_subj5{1,k}))*100; 
           perc_FEF_subj5(k,1) = (length(left_FEF_subj5{1,k})/length(dippos_subj5{1,k}))*100;
           perc_FEF_subj5(k,3) = (length(all_FEF_subj5{1,k})/length(dippos_subj5{1,k}))*100;
end

for seiz = 1:length(dippos_subj7) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj7{1,seiz})
        xp= dippos_subj7{1,seiz}(t,1);
        yp= dippos_subj7{1,seiz}(t,2);
        zp= dippos_subj7{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj7{1,seiz}(g,1) =  t;
           all_FEF_subj7{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj7{1,seiz}(g1,1) =  t;
           all_FEF_subj7{1,seiz}(g2,1)=t;
           g1=g1+1;
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj7,2)
    %LEFT col 1, RIGHT col 2
      perc_FEF_subj7(k,2) = (length(right_FEF_subj7{1,k})/length(dippos_subj7{1,k}))*100;
      perc_FEF_subj7(k,3) = (length(all_FEF_subj7{1,k})/length(dippos_subj7{1,k}))*100;
      if k<8; perc_FEF_subj7(k,1) = (length(left_FEF_subj7{1,k})/length(dippos_subj7{1,k}))*100;
      else; perc_FEF_subj7(k,1) =0; end
end
for seiz = 1:length(dippos_subj8) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj8{1,seiz})
        xp= dippos_subj8{1,seiz}(t,1);
        yp= dippos_subj8{1,seiz}(t,2);
        zp= dippos_subj8{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj8{1,seiz}(g,1) =  t;
           all_FEF_subj8{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj8{1,seiz}(g1,1) =  t;
           all_FEF_subj8{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj8,2)
    %LEFT col 1, RIGHT col 2
           perc_FEF_subj8(k,2) = (length(right_FEF_subj8{1,k})/length(dippos_subj8{1,k}))*100; 
           perc_FEF_subj8(k,1) = (length(left_FEF_subj8{1,k})/length(dippos_subj8{1,k}))*100;
           perc_FEF_subj8(k,3) = (length(all_FEF_subj8{1,k})/length(dippos_subj8{1,k}))*100;
end
for seiz = 1:length(dippos_subj13) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj13{1,seiz})
        xp= dippos_subj13{1,seiz}(t,1);
        yp= dippos_subj13{1,seiz}(t,2);
        zp= dippos_subj13{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj13{1,seiz}(g,1) =  t;
           all_FEF_subj13{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj13{1,seiz}(g1,1) =  t;
           all_FEF_subj13{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj13,2)
    %LEFT col 1, RIGHT col 2
           perc_FEF_subj13(k,2) = (length(right_FEF_subj13{1,k})/length(dippos_subj13{1,k}))*100; 
           perc_FEF_subj13(k,3) = (length(all_FEF_subj13{1,k})/length(dippos_subj13{1,k}))*100; 
           perc_FEF_subj13(k,1) = 0;          
end
for seiz = 1:length(dippos_subj15) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj15{1,seiz})
        xp= dippos_subj15{1,seiz}(t,1);
        yp= dippos_subj15{1,seiz}(t,2);
        zp= dippos_subj15{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj15{1,seiz}(g,1) =  t;
           all_FEF_subj15{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj15{1,seiz}(g1,1) =  t;
           all_FEF_subj15{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj15,2)
    %LEFT col 1, RIGHT col 2
      perc_FEF_subj15(k,2) = (length(right_FEF_subj15{1,k})/length(dippos_subj15{1,k}))*100;
      perc_FEF_subj15(k,1) = (length(left_FEF_subj15{1,k})/length(dippos_subj15{1,k}))*100;
      perc_FEF_subj15(k,3) = (length(all_FEF_subj15{1,k})/length(dippos_subj15{1,k}))*100;
      
end

FEF_unpres = [perc_FEF_subj5; perc_FEF_subj7; perc_FEF_subj8; perc_FEF_subj13; perc_FEF_subj15];
FEF_unpres(29,:)=[0,0,0];
percentage_left_unpres = mean(nonzeros(FEF_unpres(:,1)));
percentage_right_unpres = mean(nonzeros(FEF_unpres(:,2)));
percentage_all_unpres = mean(nonzeros(FEF_unpres(:,3)));
%%
%num of dipoles in FEF preserved 
for seiz = 1:length(dippos_subj3) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj3{1,seiz})
        xp= dippos_subj3{1,seiz}(t,1);
        yp= dippos_subj3{1,seiz}(t,2);
        zp= dippos_subj3{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj3{1,seiz}(g,1) =  t;
           all_FEF_subj3{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj3{1,seiz}(g1,1) =  t;
           all_FEF_subj3{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj3,2)
    %LEFT col 1, RIGHT col 2
           perc_FEF_subj3(k,2) = (length(right_FEF_subj3{1,k})/length(dippos_subj3{1,k}))*100; 
           perc_FEF_subj3(k,1) = (length(left_FEF_subj3{1,k})/length(dippos_subj3{1,k}))*100;
           perc_FEF_subj3(k,3) = (length(all_FEF_subj3{1,k})/length(dippos_subj3{1,k}))*100;
end
for seiz = 1:length(dippos_subj6) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj6{1,seiz})
        xp= dippos_subj6{1,seiz}(t,1);
        yp= dippos_subj6{1,seiz}(t,2);
        zp= dippos_subj6{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj6{1,seiz}(g,1) =  t;
           all_FEF_subj6{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj6{1,seiz}(g1,1) =  t;
           all_FEF_subj6{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj6,2)
    %LEFT col 1, RIGHT col 2
           perc_FEF_subj6(k,2) = (length(right_FEF_subj6{1,k})/length(dippos_subj6{1,k}))*100;
           perc_FEF_subj6(k,3) = (length(all_FEF_subj6{1,k})/length(dippos_subj6{1,k}))*100;
           if k<2;perc_FEF_subj6(k,1) = (length(left_FEF_subj6{1,k})/length(dippos_subj6{1,k}))*100; 
           else   perc_FEF_subj6(k,1) = 0; end 
end
for seiz = 1:length(dippos_subj10) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj10{1,seiz})
        xp= dippos_subj10{1,seiz}(t,1);
        yp= dippos_subj10{1,seiz}(t,2);
        zp= dippos_subj10{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj10{1,seiz}(g,1) =  t;
           all_FEF_subj10{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj10{1,seiz}(g1,1) =  t;
           all_FEF_subj10{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj10,2)
    %LEFT col 1, RIGHT col 2
           perc_FEF_subj10(k,2) = (length(right_FEF_subj10{1,k})/length(dippos_subj10{1,k}))*100;
           perc_FEF_subj10(k,3) = (length(all_FEF_subj10{1,k})/length(dippos_subj10{1,k}))*100;
           if k<5; perc_FEF_subj10(k,1) = (length(left_FEF_subj10{1,k})/length(dippos_subj10{1,k}))*100;
           else;perc_FEF_subj10(k,1) = 0; end     
end
for seiz = 1:length(dippos_subj12) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj12{1,seiz})
        xp= dippos_subj12{1,seiz}(t,1);
        yp= dippos_subj12{1,seiz}(t,2);
        zp= dippos_subj12{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj12{1,seiz}(g,1) =  t;
           all_FEF_subj12{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj12{1,seiz}(g1,1) =  t;
           all_FEF_subj12{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
for k = 1:size(right_FEF_subj12,2)
    %LEFT col 1, RIGHT col 2
           perc_FEF_subj12(k,2) = (length(right_FEF_subj12{1,k})/length(dippos_subj12{1,k}))*100; 
           perc_FEF_subj12(k,1) = (length(left_FEF_subj12{1,k})/length(dippos_subj12{1,k}))*100; 
           perc_FEF_subj12(k,3) = (length(all_FEF_subj12{1,k})/length(dippos_subj12{1,k}))*100;           
end
for seiz = 1:length(dippos_subj14) 
    g=1;
    g1=1;
    g2=1;
    for t = 1:length(dippos_subj14{1,seiz})
        xp= dippos_subj14{1,seiz}(t,1);
        yp= dippos_subj14{1,seiz}(t,2);
        zp= dippos_subj14{1,seiz}(t,3);
        insphFEF1_R = ((xp-FEF1_R(1)).^2 + (yp-FEF1_R(2)).^2 + (zp-FEF1_R(3)).^2) <= r^2;
        insphFEF1_L = ((xp-FEF1_L(1)).^2 + (yp-FEF1_L(2)).^2 + (zp-FEF1_L(3)).^2) <= r^2;
        insphFEF2_R = ((xp-FEF2_R(1)).^2 + (yp-FEF2_R(2)).^2 + (zp-FEF2_R(3)).^2) <= r^2;
        insphFEF2_L = ((xp-FEF2_L(1)).^2 + (yp-FEF2_L(2)).^2 + (zp-FEF2_L(3)).^2) <= r^2;
        insphFEF3_R = ((xp-FEF3_R(1)).^2 + (yp-FEF3_R(2)).^2 + (zp-FEF3_R(3)).^2) <= r^2;
        insphFEF4_R = ((xp-FEF4_R(1)).^2 + (yp-FEF4_R(2)).^2 + (zp-FEF4_R(3)).^2) <= r^2;
        insphFEF5_R = ((xp-FEF5_R(1)).^2 + (yp-FEF5_R(2)).^2 + (zp-FEF5_R(3)).^2) <= r^2;
        insphFEF3_L = ((xp-FEF3_L(1)).^2 + (yp-FEF3_L(2)).^2 + (zp-FEF3_L(3)).^2) <= r^2;
        insphFEF4_L = ((xp-FEF4_L(1)).^2 + (yp-FEF4_L(2)).^2 + (zp-FEF4_L(3)).^2) <= r^2;
        insphFEF5_L = ((xp-FEF5_L(1)).^2 + (yp-FEF5_L(2)).^2 + (zp-FEF5_L(3)).^2) <= r^2;
        if insphFEF1_R==1 || insphFEF2_R==1 || insphFEF3_R==1 || insphFEF4_R==1 || insphFEF5_R==1
           right_FEF_subj14{1,seiz}(g,1) =  t;
           all_FEF_subj14{1,seiz}(g2,1)=t;
           g=g+1;
           g2=g2+1;
        elseif insphFEF1_L==1 || insphFEF2_L==1 || insphFEF3_L==1 || insphFEF4_L==1 || insphFEF5_L==1   
           left_FEF_subj14{1,seiz}(g1,1) =  t;
           all_FEF_subj14{1,seiz}(g2,1)=t;
           g1=g1+1; 
           g2=g2+1;
        end    
    end    
end
    %LEFT col 1, RIGHT col 2, all col 3
perc_FEF_subj14(1,2) = (length(right_FEF_subj14{1,1})/length(dippos_subj12{1,1}))*100; 
perc_FEF_subj14(1,3) = (length(all_FEF_subj14{1,1})/length(dippos_subj12{1,1}))*100; 
perc_FEF_subj14(1,1) = 0;         


FEF_pres = [perc_FEF_subj3; perc_FEF_subj6; perc_FEF_subj10; perc_FEF_subj12; perc_FEF_subj14];
percentage_left_pres = mean(nonzeros(FEF_pres(:,1)));
percentage_right_pres = mean(nonzeros(FEF_pres(:,2)));
percentage_all_pres = mean(nonzeros(FEF_pres(:,3)));

%% dippos in FEF preserved, to plot
for seiz = 1:length(dippos_subj3)
    %1st raw left, 2nd raw right
    dippos_FEF_subj3{1,seiz}=dippos_subj3{1,seiz}(left_FEF_subj3{1,seiz},:);
    dippos_FEF_subj3{2,seiz}=dippos_subj3{1,seiz}(right_FEF_subj3{1,seiz},:); 
end
dippos_FEF_3 = [ dippos_FEF_subj3{1,1}; dippos_FEF_subj3{1,2}; dippos_FEF_subj3{2,1}; dippos_FEF_subj3{2,2}];

for seiz = 1:length(dippos_subj6)
    %1st raw left, 2nd raw right
    if seiz<2; dippos_FEF_subj6{1,seiz}=dippos_subj6{1,seiz}(left_FEF_subj6{1,seiz},:);
    else;  dippos_FEF_subj6{1,seiz}=[0,0,0]; end
    dippos_FEF_subj6{2,seiz}=dippos_subj6{1,seiz}(right_FEF_subj6{1,seiz},:); 
end
dippos_FEF_6 = [ dippos_FEF_subj6{1,1}; dippos_FEF_subj6{1,2}; dippos_FEF_subj6{2,1}; dippos_FEF_subj6{2,2}];

for seiz = 1:length(dippos_subj10)
    %1st raw left, 2nd raw right
    if seiz<5;dippos_FEF_subj10{1,seiz}=dippos_subj10{1,seiz}(left_FEF_subj10{1,seiz},:);
    else;  dippos_FEF_subj10{1,seiz}=[0,0,0]; end
    dippos_FEF_subj10{2,seiz}=dippos_subj10{1,seiz}(right_FEF_subj10{1,seiz},:); 
end
dippos_FEF_10 = [ dippos_FEF_subj10{1,1}; dippos_FEF_subj10{1,2};dippos_FEF_subj10{1,3}; dippos_FEF_subj10{1,4}; dippos_FEF_subj10{1,5};dippos_FEF_subj10{2,1}; dippos_FEF_subj10{2,2}; dippos_FEF_subj10{2,3}; dippos_FEF_subj10{2,4}; dippos_FEF_subj10{2,5}];

for seiz = 1:length(dippos_subj12)
    %1st raw left, 2nd raw right
    dippos_FEF_subj12{1,seiz}=dippos_subj12{1,seiz}(left_FEF_subj12{1,seiz},:);
   
    dippos_FEF_subj12{2,seiz}=dippos_subj12{1,seiz}(right_FEF_subj12{1,seiz},:); 
end
dippos_FEF_12 = [ dippos_FEF_subj12{1,1}; dippos_FEF_subj12{1,2};dippos_FEF_subj12{1,3}; dippos_FEF_subj12{1,4}; dippos_FEF_subj12{1,5}; dippos_FEF_subj12{1,6};dippos_FEF_subj12{1,7};dippos_FEF_subj12{1,8};dippos_FEF_subj12{1,9};dippos_FEF_subj12{1,10};dippos_FEF_subj12{1,11}; dippos_FEF_subj12{1,12}; dippos_FEF_subj12{1,13};dippos_FEF_subj12{1,14};   dippos_FEF_subj12{2,1}; dippos_FEF_subj12{2,2}; dippos_FEF_subj12{2,3}; dippos_FEF_subj12{2,4}; dippos_FEF_subj12{2,5}; dippos_FEF_subj12{2,6}; dippos_FEF_subj12{2,7}; dippos_FEF_subj12{2,8}; dippos_FEF_subj12{2,9}; dippos_FEF_subj12{2,10}; dippos_FEF_subj12{2,11}; dippos_FEF_subj12{2,12}; dippos_FEF_subj12{2,13}; dippos_FEF_subj12{2,14}];
for seiz = 1:length(dippos_subj14)
    %1st raw left, 2nd raw right
    dippos_FEF_subj14{1,seiz}=[0,0,0];
    dippos_FEF_subj14{2,seiz}=dippos_subj14{1,seiz}(right_FEF_subj14{1,seiz},:); 
end
dippos_FEF_14 = [ dippos_FEF_subj14{1,1};dippos_FEF_subj14{2,1}];

%% dippos in FEF unpr
for seiz = 1:length(dippos_subj5)
    %1st raw left, 2nd raw right
    dippos_FEF_subj5{1,seiz}=dippos_subj5{1,seiz}(left_FEF_subj5{1,seiz},:);
    dippos_FEF_subj5{2,seiz}=dippos_subj5{1,seiz}(right_FEF_subj5{1,seiz},:); 
end
dippos_FEF_5 = [ dippos_FEF_subj5{1,1}; dippos_FEF_subj5{1,2};dippos_FEF_subj5{1,3};dippos_FEF_subj5{1,4}; dippos_FEF_subj5{2,1}; dippos_FEF_subj5{2,2}; dippos_FEF_subj5{2,3}; dippos_FEF_subj5{2,4}];

for seiz = 1:length(dippos_subj7)
    %1st raw left, 2nd raw right
    if seiz<8;dippos_FEF_subj7{1,seiz}=dippos_subj7{1,seiz}(left_FEF_subj7{1,seiz},:);
    else; ;dippos_FEF_subj7{1,seiz}=[0,0,0]; end
    dippos_FEF_subj7{2,seiz}=dippos_subj7{1,seiz}(right_FEF_subj7{1,seiz},:); 
end
dippos_FEF_7 = [ dippos_FEF_subj7{1,1}; dippos_FEF_subj7{1,2};dippos_FEF_subj7{1,3};dippos_FEF_subj7{1,4};dippos_FEF_subj7{1,5};dippos_FEF_subj7{1,6};dippos_FEF_subj7{1,7};dippos_FEF_subj7{1,8};dippos_FEF_subj7{1,9};dippos_FEF_subj7{1,10}; dippos_FEF_subj7{2,1}; dippos_FEF_subj7{2,2}; dippos_FEF_subj7{2,3}; dippos_FEF_subj7{2,4}; dippos_FEF_subj7{2,5}; dippos_FEF_subj7{2,6}; dippos_FEF_subj7{2,7}; dippos_FEF_subj7{2,8}; dippos_FEF_subj7{2,9}; dippos_FEF_subj7{2,10}];

for seiz = 1:length(dippos_subj8)
    %1st raw left, 2nd raw right
    dippos_FEF_subj8{1,seiz}=dippos_subj8{1,seiz}(left_FEF_subj8{1,seiz},:);
    dippos_FEF_subj8{2,seiz}=dippos_subj8{1,seiz}(right_FEF_subj8{1,seiz},:); 
end
dippos_FEF_8 = [ dippos_FEF_subj8{1,1}; dippos_FEF_subj8{1,2};dippos_FEF_subj8{1,3}; dippos_FEF_subj8{1,4}; dippos_FEF_subj8{1,5};dippos_FEF_subj8{1,6};dippos_FEF_subj8{1,7};dippos_FEF_subj8{2,1}; dippos_FEF_subj8{2,2}; dippos_FEF_subj8{2,3}; dippos_FEF_subj8{2,4}; dippos_FEF_subj8{2,5}; dippos_FEF_subj8{2,6}; dippos_FEF_subj8{2,7}];

for seiz = 1:length(dippos_subj13)
    %1st raw left, 2nd raw right
    if seiz<2; dippos_FEF_subj13{1,seiz}=dippos_subj13{1,seiz}(left_FEF_subj13{1,seiz},:);
    else;dippos_FEF_subj13{1,seiz}=[0,0,0]; end
    dippos_FEF_subj13{2,seiz}=dippos_subj13{1,seiz}(right_FEF_subj13{1,seiz},:); 
end
dippos_FEF_13 = [ dippos_FEF_subj13{1,1}; dippos_FEF_subj13{1,2}; dippos_FEF_subj13{2,1}; dippos_FEF_subj13{2,2}];

for seiz = 1:length(dippos_subj15)
    %1st raw left, 2nd raw right
    dippos_FEF_subj15{1,seiz}=dippos_subj15{1,seiz}(left_FEF_subj15{1,seiz},:); 
    dippos_FEF_subj15{2,seiz}=dippos_subj15{1,seiz}(right_FEF_subj15{1,seiz},:); 
end
dippos_FEF_15 = [ dippos_FEF_subj15{1,1};dippos_FEF_subj15{1,2};dippos_FEF_subj15{1,3};dippos_FEF_subj15{1,4};dippos_FEF_subj15{1,5};dippos_FEF_subj15{1,6};dippos_FEF_subj15{1,7};dippos_FEF_subj15{2,1};dippos_FEF_subj15{2,2};dippos_FEF_subj15{2,3};dippos_FEF_subj15{2,4};dippos_FEF_subj15{2,5};dippos_FEF_subj15{2,6};dippos_FEF_subj15{2,7}];


%% plot 

%pres
dippos_FEF_pres=[dippos_FEF_3;dippos_FEF_6; dippos_FEF_10;dippos_FEF_12; dippos_FEF_14];

x=dippos_FEF_pres(:,1);
x=x(x~=0);
y=dippos_FEF_pres(:,2);
y=y(y~=0);
z=dippos_FEF_pres(:,3);
z=z(z~=0);


load('C:\Users\vb\Documents\MATLAB\R2018b\fieldtrip-20220714\template\headmodel\standard_bem.mat') % variable 'vol'
figure(21), ft_plot_mesh(vol.bnd(3), 'edgecolor', 'none', 'facecolor', 'skin', 'facealpha', 0.3), hold on, 
% colorbar
hold on
surf(X2+FEF1_L(1),Y2+FEF1_L(2),Z2+FEF1_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [215,25,28]/256)
surf(X2+FEF1_R(1),Y2+FEF1_R(2),Z2+FEF1_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [215,25,28]/256)
surf(X2+FEF2_L(1),Y2+FEF2_L(2),Z2+FEF2_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [253,174,97]/256)
surf(X2+FEF2_R(1),Y2+FEF2_R(2),Z2+FEF2_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [253,174,97]/256)
surf(X2+FEF3_L(1),Y2+FEF3_L(2),Z2+FEF3_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [255,255,191]/256)
surf(X2+FEF3_R(1),Y2+FEF3_R(2),Z2+FEF3_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [255,255,191]/256)
surf(X2+FEF4_L(1),Y2+FEF4_L(2),Z2+FEF4_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [171,217,233]/256)
surf(X2+FEF4_R(1),Y2+FEF4_R(2),Z2+FEF4_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [171,217,233]/256)
surf(X2+FEF5_L(1),Y2+FEF5_L(2),Z2+FEF5_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [44,123,182]/256)
surf(X2+FEF5_R(1),Y2+FEF5_R(2),Z2+FEF5_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [44,123,182]/256)
%hs1 =mesh(X2+FEF1_L(1),Y2+FEF1_L(2),Z2+FEF1_L(3),'FaceAlpha',0.5)
plot3(x, y, z, 'o','MarkerSize',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','k')
title('Preserved','FontName','times','FontSize',15)
%surf(X2+5,Y2-5,Z2)  %coord (5,-5,2)

%print -depsc FEF_unpr.eps


%unpres
dippos_FEF_unpres=[dippos_FEF_5;dippos_FEF_7; dippos_FEF_8;dippos_FEF_13; dippos_FEF_15];

x=dippos_FEF_unpres(:,1);
x=x(x~=0);
y=dippos_FEF_unpres(:,2);
y([354:356,844],:)=[];
z=dippos_FEF_unpres(:,3);
z=z(z~=0);


load('C:\Users\vb\Documents\MATLAB\R2018b\fieldtrip-20220714\template\headmodel\standard_bem.mat') % variable 'vol'
figure(22), ft_plot_mesh(vol.bnd(3), 'edgecolor', 'none', 'facecolor', 'skin', 'facealpha', 0.3), hold on, 
% colorbar
hold on
surf(X2+FEF1_L(1),Y2+FEF1_L(2),Z2+FEF1_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [215,25,28]/256)
surf(X2+FEF1_R(1),Y2+FEF1_R(2),Z2+FEF1_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [215,25,28]/256)
surf(X2+FEF2_L(1),Y2+FEF2_L(2),Z2+FEF2_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [253,174,97]/256)
surf(X2+FEF2_R(1),Y2+FEF2_R(2),Z2+FEF2_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [253,174,97]/256)
surf(X2+FEF3_L(1),Y2+FEF3_L(2),Z2+FEF3_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [255,255,191]/256)
surf(X2+FEF3_R(1),Y2+FEF3_R(2),Z2+FEF3_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [255,255,191]/256)
surf(X2+FEF4_L(1),Y2+FEF4_L(2),Z2+FEF4_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [171,217,233]/256)
surf(X2+FEF4_R(1),Y2+FEF4_R(2),Z2+FEF4_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [171,217,233]/256)
surf(X2+FEF5_L(1),Y2+FEF5_L(2),Z2+FEF5_L(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [44,123,182]/256)
surf(X2+FEF5_R(1),Y2+FEF5_R(2),Z2+FEF5_R(3),'FaceAlpha',0.5,'EdgeColor','none','FaceColor', [44,123,182]/256)
%hs1 =mesh(X2+FEF1_L(1),Y2+FEF1_L(2),Z2+FEF1_L(3),'FaceAlpha',0.5)
plot3(x, y, z, 'o','MarkerSize',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','k')
title('Unpreserved','FontName','times','FontSize',15)
%surf(X2+5,Y2-5,Z2)  %coord (5,-5,2)

