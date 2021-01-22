close all
clear all
set(0,'defaulttextInterpreter'  ,'latex');
set(0,'defaultlegendInterpreter','latex'); 
set(0,'defaultAxesFontsize'     ,14) ;
set(0,'defaultLegendFontsize'   ,16) ;
set(0,'defaultLineLineWidth'    ,2 ) ;
set(0,'DefaultLineMarkerSize'   ,8 ) ;


folder1 = './'; %folder where the results are stored.
folder2 = '../ControlNL3_uncon/'; %folder where the results are stored.
folder3 = '../ControlNL3_tnc/';
folder4 = '../ControlNL3_nc/';
folder5 = '../ControlNL3_SISO/';
folder6 = '../ControlNL3_SISO/';
folderTarget='../ControlNL3_prooftarget/'
%folder3 = '../fullStep/'; %folder where the results are stored.

%folder = '../fullStep_conT1';
pertNormSource_con = dlmread(sprintf('%s/con1/runNorm.txt',folder1),'',1,0);
ProjShapeSource_con = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder1),'',1,0);
pertNormSource_uncon = dlmread(sprintf('%s/con1/runNorm.txt',folder2),'',1,0);
ProjShapeSource_uncon = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder2),'',1,0);
pertNormSource_tnc = dlmread(sprintf('%s/con1/runNorm.txt',folder3),'',1,0);
ProjShapeSource_tnc = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder2),'',1,0);
pertNormSource_nc = dlmread(sprintf('%s/con1/runNorm.txt',folder4),'',1,0);
ProjShapeSource_nc = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder4),'',1,0);
pertNormSource_SISO_c = dlmread(sprintf('%s/con1_causal/runNorm.txt',folder5),'',1,0);
ProjShapeSource_SISO_c = dlmread(sprintf('%s/con1_causal/ProjShapes_7201.dat',folder5),'',1,0);
pertNormSource_SISO_tnc = dlmread(sprintf('%s/con1/runNorm.txt',folder6),'',1,0);
ProjShapeSource_SISO_tnc = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder6),'',1,0);
pertNormSource_target = dlmread(sprintf('%s/con1/runNorm.txt',folderTarget),'',1,0);
ProjShapeSource_target = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folderTarget),'',1,0);

a=500;
b=1000;
[rows, columns] = size(pertNormSource_con);
for col = 1 : 2
  theSum = 0;
  for row = a : b
    theSum = theSum + pertNormSource_con(row, col);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_con(col) = theSum / (b-a+1);
end
Avg_con=Avg_pertNorm_con(2);

[rows_rnd, columns_rnd] = size(pertNormSource_uncon);
for col_rnd = 1 : 2
  theSum = 0;
  for row_rnd = 833 : 999
    theSum = theSum + pertNormSource_uncon(row_rnd, col_rnd);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_rnd(col_rnd) = theSum / (999-833+1);
end
Avg_rnd=Avg_pertNorm_rnd(2);
% pertNormSource=pertNormSource';
% Avg_pertNorm=mean(pertNormSource,2);
% Avg=Avg_pertNorm(2,1);

% figure(1)
%         plot(pertNormSource_con(:,1),pertNormSource_con(:,2),'r',pertNormSource_uncon(:,1),pertNormSource_uncon(:,2),'b')
%         xlim([500,750])
%         legend('pertNorm(controlled)','pertNorm(uncontrolled)')
%         xlabel('t');ylabel('pertNorm');
%         title('controlled and uncontrolled pertNorm')

figure(2)
        plot(pertNormSource_uncon(1:10000,1),pertNormSource_uncon(1:10000,2),'k', ...
        pertNormSource_SISO_tnc(1:10000,1),pertNormSource_SISO_tnc(1:10000,2),'m--',...
        pertNormSource_SISO_c(1:10000,1),pertNormSource_SISO_c(1:10000,2),'b--',...
        pertNormSource_tnc(1:10000,1),pertNormSource_tnc(1:10000,2),'m',...
        pertNormSource_con(1:10000,1),pertNormSource_con(1:10000,2),'b')
    
        xlim([100,450])
        ylim([0,0.25])
        legend('uncontrolled','TNC(SISO)','Causal(SISO)','TNC(MIMO)','Causal(MIMO)')
        xlabel('t');ylabel('pertNorm');
        title('$Non-linear system$')
        grid on;
print(gcf,'-dpsc2','pertNorm_NL_3MIMO.eps'); 
figure(3)
        plot(ProjShapeSource_uncon(1:10000,1),ProjShapeSource_uncon(1:10000,15),'r',ProjShapeSource_target(1:10000,1),ProjShapeSource_target(1:10000,15),'b')
        xlim([50,150])
        %ylim([0,1])
        legend('uncontrolled','controlled(C)')
        xlabel('t');ylabel('Target');
        title('Non-linear system')
        grid on;        
 print(gcf,'-dpsc2','Target_MIMO_C_xm2.eps');