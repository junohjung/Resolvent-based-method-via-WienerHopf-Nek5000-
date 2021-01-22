close all
clear all
set(0,'defaulttextInterpreter'  ,'latex');
set(0,'defaultlegendInterpreter','latex'); 
set(0,'defaultAxesFontsize'     ,10) ;
set(0,'defaultLegendFontsize'   ,10) ;
set(0,'defaultLineLineWidth'    ,2 ) ;
set(0,'DefaultLineMarkerSize'   ,8 ) ;
%% PertNorm 
forcingDFA=[1:6];
pertNorm_NL=[97.6268,82.5521,30.7217,26.2559,57.8161,100];
pertNorm_L(1:6)=42.3905;
Factor(1:6)=100;
pertNorm_NL=Factor-pertNorm_NL;
pertNorm_L=Factor-pertNorm_L;

%% ProjShape Norm at three targets
ProjNorm_NL=[97.6622,95.5214,52.7707,24.6788,41.6134,100];
ProjNorm_L(1:6)=20.7897;
Factor(1:6)=100;
ProjNorm_NL=Factor-ProjNorm_NL;
ProjNorm_L=Factor-ProjNorm_L;

figure(1)
        plot(forcingDFA,pertNorm_L,'k-o', ...
        forcingDFA,pertNorm_NL,'k-*',...
        forcingDFA,ProjNorm_L,'b--o',...
        forcingDFA,ProjNorm_NL,'b--*');
        xticklabels({'f*1e-3','f*1e-2','f*1e-1','f','f*1e+1','f*1e+2'})
        xlim([1,6])
        ylim([0,100])
        legend('Overall(Li)','Overall(NL)','Target(Li)','Target(NL)')
        xlabel('forcing amplitude');ylabel('perturbation reduction rate(%)');
        title('Linear system vs Non-linear system')
        grid on;
print(gcf,'-dpsc2','pertNormPerformance_LivsNL.eps');

