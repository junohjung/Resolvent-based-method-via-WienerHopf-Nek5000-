close all
clear all
set(0,'defaulttextInterpreter'  ,'latex');
set(0,'defaultlegendInterpreter','latex'); 
set(0,'defaultAxesFontsize'     ,14) ;
set(0,'defaultLegendFontsize'   ,16) ;
set(0,'defaultLineLineWidth'    ,2 ) ;
set(0,'DefaultLineMarkerSize'   ,8 ) ;


% folder1 = '../ControlNL3/'; %folder where the results are stored.
% folder2 = '../ControlNL3_uncon/'; %folder where the results are stored.
% folder3 = '../ControlNL3_tnc/';
% folder4 = '../ControlNL3_nc/';
% folder5 = '../ControlNL3_SISO/';
% folder6 = '../ControlNL3_SISO/';
% folderTarget='../ControlNL3_prooftarget/';
%folder7='./';
folder1='../NL_DFA_1em2';
folder2='../NL_DFA_1em2_uncon';
folder3='../NL_DFA_1em3';
folder4='../NL_DFA_1em3_uncon';
folder5='../NL_DFA_1em4';
folder6='../NL_DFA_1em4_uncon';
folder7='../NL_DFA_1em5';
folder8='../NL_DFA_1em5_uncon';
folder9='../NL_DFA_1em6';
folder10='../NL_DFA_1em6_uncon';



% pertNormSource_con_test0 = dlmread(sprintf('%s/con1_test0/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test0 = dlmread(sprintf('%s/con1_test0/ProjShapes_7201.dat',folder7),'',1,0);
% pertNormSource_con_test0_uncon = dlmread(sprintf('%s/con1_test1_uncon/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test0_uncon = dlmread(sprintf('%s/con1_test1_uncon/ProjShapes_7201.dat',folder7),'',1,0);
% pertNormSource_con_test1 = dlmread(sprintf('%s/con1_test1/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test1 = dlmread(sprintf('%s/con1_test1/ProjShapes_7201.dat',folder7),'',1,0);
% pertNormSource_con_test1_uncon = dlmread(sprintf('%s/con1_test1_uncon/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test1_uncon = dlmread(sprintf('%s/con1_test1_uncon/ProjShapes_7201.dat',folder7),'',1,0);
% 
% pertNormSource_con_test2 = dlmread(sprintf('%s/con1_test2/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test2 = dlmread(sprintf('%s/con1_test2/ProjShapes_7201.dat',folder7),'',1,0);
% pertNormSource_con_test2_uncon = dlmread(sprintf('%s/con1_test2_uncon/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test2_uncon = dlmread(sprintf('%s/con1_test2_uncon/ProjShapes_7201.dat',folder7),'',1,0);
% 
% pertNormSource_con_test3 = dlmread(sprintf('%s/con1_test3/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test3 = dlmread(sprintf('%s/con1_test3/ProjShapes_7201.dat',folder7),'',1,0);
% pertNormSource_con_test3_uncon = dlmread(sprintf('%s/con1_test3_uncon/runNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test3_uncon = dlmread(sprintf('%s/con1_test3_uncon/ProjShapes_7201.dat',folder7),'',1,0);
% 
% pertNormSource_con_test4 = dlmread(sprintf('%s/con1_test4/pertNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test4 = dlmread(sprintf('%s/con1_test4/ProjShapes_01.dat',folder7),'',1,0);
% pertNormSource_con_test4_uncon = dlmread(sprintf('%s/con1_test4_uncon/pertNorm.txt',folder7),'',1,0);
% ProjShapeSource_con_test4_uncon = dlmread(sprintf('%s/con1_test4_uncon/ProjShapes_01.dat',folder7),'',1,0);
% pertNormSource_uncon = dlmread(sprintf('%s/con1/runNorm.txt',folder2),'',1,0);
% ProjShapeSource_uncon = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder2),'',1,0);
% pertNormSource_tnc = dlmread(sprintf('%s/con1/runNorm.txt',folder3),'',1,0);
% ProjShapeSource_tnc = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder2),'',1,0);
% pertNormSource_nc = dlmread(sprintf('%s/con1/runNorm.txt',folder4),'',1,0);
% ProjShapeSource_nc = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder4),'',1,0);
% pertNormSource_SISO_c = dlmread(sprintf('%s/con1_causal/runNorm.txt',folder5),'',1,0);
% ProjShapeSource_SISO_c = dlmread(sprintf('%s/con1_causal/ProjShapes_7201.dat',folder5),'',1,0);
% pertNormSource_SISO_tnc = dlmread(sprintf('%s/con1/runNorm.txt',folder6),'',1,0);
% ProjShapeSource_SISO_tnc = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder6),'',1,0);
% pertNormSource_target = dlmread(sprintf('%s/con1/runNorm.txt',folderTarget),'',1,0);
% ProjShapeSource_target = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folderTarget),'',1,0);
pertNormSource_1e2 = dlmread(sprintf('%s/con1/runNorm.txt',folder1),'',1,0);
ProjShapeSource_1e2 = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder1),'',1,0);
pertNormSource_1e2_uncon = dlmread(sprintf('%s/con1/runNorm.txt',folder2),'',1,0);
ProjShapeSource_1e2_uncon = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder2),'',1,0);
pertNormSource_1e3 = dlmread(sprintf('%s/con1/runNorm.txt',folder3),'',1,0);
ProjShapeSource_1e3 = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder3),'',1,0);
pertNormSource_1e3_uncon = dlmread(sprintf('%s/con1/runNorm.txt',folder4),'',1,0);
ProjShapeSource_1e3_uncon = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder4),'',1,0);
pertNormSource_1e4 = dlmread(sprintf('%s/con1/runNorm.txt',folder5),'',1,0);
ProjShapeSource_1e4 = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder5),'',1,0);
pertNormSource_1e4_uncon = dlmread(sprintf('%s/con1/runNorm.txt',folder6),'',1,0);
ProjShapeSource_1e4_uncon = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder6),'',1,0);
pertNormSource_1e5 = dlmread(sprintf('%s/con1/runNorm.txt',folder7),'',1,0);
ProjShapeSource_1e5 = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder7),'',1,0);
pertNormSource_1e5_uncon = dlmread(sprintf('%s/con1/runNorm.txt',folder8),'',1,0);
ProjShapeSource_1e5_uncon = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder8),'',1,0);
pertNormSource_1e6 = dlmread(sprintf('%s/con1/runNorm.txt',folder9),'',1,0);
ProjShapeSource_1e6 = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder9),'',1,0);
pertNormSource_1e6_uncon = dlmread(sprintf('%s/con1/runNorm.txt',folder10),'',1,0);
ProjShapeSource_1e6_uncon = dlmread(sprintf('%s/con1/ProjShapes_7201.dat',folder10),'',1,0);

% %% controlled pertNorm average test4
% 
% a=1000;
% b=4000;
% [rows, columns] = size(pertNormSource_con_test4);
% for col = 1 : 2
%   theSum = 0;
%   for row = a : b
%     theSum = theSum + pertNormSource_con_test4(row, col);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test4(col) = theSum / (b-a+1);
% end
% Avg_con_test4=Avg_pertNorm_con_test4(2);
% 
% %% uncontrol pertNorm average test4
% [rows_rnd, columns_rnd] = size(pertNormSource_con_test4_uncon);
% for col_rnd = 1 : 2
%   theSum = 0;
%   for row_rnd = a : b
%     theSum = theSum + pertNormSource_con_test4_uncon(row_rnd, col_rnd);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test4_uncon(col_rnd) = theSum / (b-a+1);
% end
% Avg_con_test4_uncon=Avg_pertNorm_con_test4_uncon(2);
% %% performance of pertNorm test4
% 
% Perf_test4=100*(1-(abs(Avg_con_test4_uncon-Avg_con_test4)/abs(Avg_con_test4_uncon)));
% 
% 
% 
% %% Projshapes average test4
% 
% a=1000;
% b=10000;
% [rows, columns] = size(ProjShapeSource_con_test4);
% N_Target=3;
% targetNr1 = 1+11; %add 1 for time column vector
% targetNr2 = 1+12;
% targetNr3 = 1+13;
% %% first target
% theSum = 0;
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test4(row,targetNr1)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test4_uncon(row,targetNr1)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test4_1 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test4_1_uncon = sqrt(theSum_uncon) / (b-a+1);
% %% second target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test4(row,targetNr2)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test4_uncon(row,targetNr2)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test4_2 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test4_2_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% %% third target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test4(row,targetNr3)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test4_uncon(row,targetNr3)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test4_3 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test4_3_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% 
% Avg_Projshapes_con_test4=(Avg_Projshapes_con_test4_1+Avg_Projshapes_con_test4_2+Avg_Projshapes_con_test4_3)/N_Target;
% Avg_Projshapes_con_test4_uncon=(Avg_Projshapes_con_test4_1_uncon+Avg_Projshapes_con_test4_2_uncon+Avg_Projshapes_con_test4_3_uncon)/N_Target;
% 
% Perf_Proj_con_test4=100*(1-(abs(Avg_Projshapes_con_test4_uncon-Avg_Projshapes_con_test4)/abs(Avg_Projshapes_con_test4_uncon)));
% 
% %% controlled pertNorm average test3
% 
% a=10000;
% b=15000;
% [rows, columns] = size(pertNormSource_con_test3);
% for col = 1 : 2
%   theSum = 0;
%   for row = a : b
%     theSum = theSum + pertNormSource_con_test3(row, col);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test3(col) = theSum / (b-a+1);
% end
% Avg_con_test3=Avg_pertNorm_con_test3(2);
% 
% %% uncontrol pertNorm average test3
% [rows_rnd, columns_rnd] = size(pertNormSource_con_test3_uncon);
% for col_rnd = 1 : 2
%   theSum = 0;
%   for row_rnd = a : b
%     theSum = theSum + pertNormSource_con_test3_uncon(row_rnd, col_rnd);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test3_uncon(col_rnd) = theSum / (b-a+1);
% end
% Avg_con_test3_uncon=Avg_pertNorm_con_test3_uncon(2);
% %% performance of pertNorm test3
% 
% Perf_test3=100*(1-(abs(Avg_con_test3_uncon-Avg_con_test3)/abs(Avg_con_test3_uncon)));
% 
% 
% 
% %% Projshapes average test3
% 
% a=10000;
% b=15000;
% [rows, columns] = size(ProjShapeSource_con_test3);
% N_Target=3;
% targetNr1 = 1+11; %add 1 for time column vector
% targetNr2 = 1+12;
% targetNr3 = 1+13;
% %% first target
% theSum = 0;
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test3(row,targetNr1)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test3_uncon(row,targetNr1)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test3_1 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test3_1_uncon = sqrt(theSum_uncon) / (b-a+1);
% %% second target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test3(row,targetNr2)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test3_uncon(row,targetNr2)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test3_2 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test3_2_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% %% third target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test3(row,targetNr3)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test3_uncon(row,targetNr3)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test3_3 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test3_3_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% 
% Avg_Projshapes_con_test3=(Avg_Projshapes_con_test3_1+Avg_Projshapes_con_test3_2+Avg_Projshapes_con_test3_3)/N_Target;
% Avg_Projshapes_con_test3_uncon=(Avg_Projshapes_con_test3_1_uncon+Avg_Projshapes_con_test3_2_uncon+Avg_Projshapes_con_test3_3_uncon)/N_Target;
% 
% Perf_Proj_con_test3=100*(1-(abs(Avg_Projshapes_con_test3_uncon-Avg_Projshapes_con_test3)/abs(Avg_Projshapes_con_test3_uncon)));
% 
% 
% %% controlled pertNorm average test2
% 
% a=10000;
% b=15000;
% [rows, columns] = size(pertNormSource_con_test2);
% for col = 1 : 2
%   theSum = 0;
%   for row = a : b
%     theSum = theSum + pertNormSource_con_test2(row, col);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test2(col) = theSum / (b-a+1);
% end
% Avg_con_test2=Avg_pertNorm_con_test2(2);
% 
% %% uncontrol pertNorm average test2
% [rows_rnd, columns_rnd] = size(pertNormSource_con_test2_uncon);
% for col_rnd = 1 : 2
%   theSum = 0;
%   for row_rnd = a : b
%     theSum = theSum + pertNormSource_con_test2_uncon(row_rnd, col_rnd);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test2_uncon(col_rnd) = theSum / (b-a+1);
% end
% Avg_con_test2_uncon=Avg_pertNorm_con_test2_uncon(2);
% %% performance of pertNorm test2
% 
% Perf_test2=100*(1-(abs(Avg_con_test2_uncon-Avg_con_test2)/abs(Avg_con_test2_uncon)));
% 
% 
% 
% %% Projshapes average test2
% 
% a=10000;
% b=15000;
% [rows, columns] = size(ProjShapeSource_con_test2);
% N_Target=3;
% targetNr1 = 1+11; %add 1 for time column vector
% targetNr2 = 1+12;
% targetNr3 = 1+13;
% %% first target
% theSum = 0;
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test2(row,targetNr1)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test2_uncon(row,targetNr1)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test2_1 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test2_1_uncon = sqrt(theSum_uncon) / (b-a+1);
% %% second target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test2(row,targetNr2)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test2_uncon(row,targetNr2)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test2_2 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test2_2_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% %% third target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test2(row,targetNr3)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test2_uncon(row,targetNr3)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test2_3 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test2_3_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% 
% Avg_Projshapes_con_test2=(Avg_Projshapes_con_test2_1+Avg_Projshapes_con_test2_2+Avg_Projshapes_con_test2_3)/N_Target;
% Avg_Projshapes_con_test2_uncon=(Avg_Projshapes_con_test2_1_uncon+Avg_Projshapes_con_test2_2_uncon+Avg_Projshapes_con_test2_3_uncon)/N_Target;
% 
% Perf_Proj_con_test2=100*(1-(abs(Avg_Projshapes_con_test2_uncon-Avg_Projshapes_con_test2)/abs(Avg_Projshapes_con_test2_uncon)));
% 
% 
% %% controlled pertNorm average test1
% 
% a=10000;
% b=15000;
% [rows, columns] = size(pertNormSource_con_test1);
% for col = 1 : 2
%   theSum = 0;
%   for row = a : b
%     theSum = theSum + pertNormSource_con_test1(row, col);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test1(col) = theSum / (b-a+1);
% end
% Avg_con_test1=Avg_pertNorm_con_test1(2);
% 
% %% uncontrol pertNorm average test1
% [rows_rnd, columns_rnd] = size(pertNormSource_con_test1_uncon);
% for col_rnd = 1 : 2
%   theSum = 0;
%   for row_rnd = a : b
%     theSum = theSum + pertNormSource_con_test1_uncon(row_rnd, col_rnd);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test1_uncon(col_rnd) = theSum / (b-a+1);
% end
% Avg_con_test1_uncon=Avg_pertNorm_con_test1_uncon(2);
% %% performance of pertNorm test1
% 
% Perf_test1=100*(1-(abs(Avg_con_test1_uncon-Avg_con_test1)/abs(Avg_con_test1_uncon)));
% 
% 
% 
% %% Projshapes average test1
% 
% a=10000;
% b=15000;
% [rows, columns] = size(ProjShapeSource_con_test1);
% N_Target=3;
% targetNr1 = 1+11; %add 1 for time column vector
% targetNr2 = 1+12;
% targetNr3 = 1+13;
% %% first target
% theSum = 0;
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test1(row,targetNr1)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test1_uncon(row,targetNr1)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test1_1 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test1_1_uncon = sqrt(theSum_uncon) / (b-a+1);
% %% second target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test1(row,targetNr2)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test1_uncon(row,targetNr2)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test1_2 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test1_2_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% %% third target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test1(row,targetNr3)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test1_uncon(row,targetNr3)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test1_3 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test1_3_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% 
% Avg_Projshapes_con_test1=(Avg_Projshapes_con_test1_1+Avg_Projshapes_con_test1_2+Avg_Projshapes_con_test1_3)/N_Target;
% Avg_Projshapes_con_test1_uncon=(Avg_Projshapes_con_test1_1_uncon+Avg_Projshapes_con_test1_2_uncon+Avg_Projshapes_con_test1_3_uncon)/N_Target;
% 
% Perf_Proj_con_test1=100*(1-(abs(Avg_Projshapes_con_test1_uncon-Avg_Projshapes_con_test1)/abs(Avg_Projshapes_con_test1_uncon)));
% 
% 
% 
% 
% 
% %% controlled pertNorm average test0
% 
% a=10000;
% b=15000;
% [rows, columns] = size(pertNormSource_con_test0);
% for col = 1 : 2
%   theSum = 0;
%   for row = a : b
%     theSum = theSum + pertNormSource_con_test0(row, col);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test0(col) = theSum / (b-a+1);
% end
% Avg_con_test0=Avg_pertNorm_con_test0(2);
% 
% %% uncontrol pertNorm average test0
% [rows_rnd, columns_rnd] = size(pertNormSource_con_test0_uncon);
% for col_rnd = 1 : 2
%   theSum = 0;
%   for row_rnd = a : b
%     theSum = theSum + pertNormSource_con_test0_uncon(row_rnd, col_rnd);
%   end
%   % Now get the mean over all values in this column.
%   Avg_pertNorm_con_test0_uncon(col_rnd) = theSum / (b-a+1);
% end
% Avg_con_test0_uncon=Avg_pertNorm_con_test0_uncon(2);
% %% performance of pertNorm test0
% 
% Perf_test0=100*(1-(abs(Avg_con_test0_uncon-Avg_con_test0)/abs(Avg_con_test0_uncon)));
% 
% 
% 
% %% Projshapes average test0
% 
% a=10000;
% b=15000;
% [rows, columns] = size(ProjShapeSource_con_test0);
% N_Target=3;
% targetNr1 = 1+11; %add 1 for time column vector
% targetNr2 = 1+12;
% targetNr3 = 1+13;
% %% first target
% theSum = 0;
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test0(row,targetNr1)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test0_uncon(row,targetNr1)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test0_1 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test0_1_uncon = sqrt(theSum_uncon) / (b-a+1);
% %% second target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test0(row,targetNr2)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test0_uncon(row,targetNr2)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test0_2 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test0_2_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% %% third target
% theSum = 0;    
% theSum_uncon = 0;
%   for row = a : b
%     theSum = theSum + (ProjShapeSource_con_test0(row,targetNr3)).^2;
%     theSum_uncon = theSum_uncon + (ProjShapeSource_con_test0_uncon(row,targetNr3)).^2;
%   end
%   % Now get the mean over all values in this column.
% Avg_Projshapes_con_test0_3 = sqrt(theSum) / (b-a+1);
% Avg_Projshapes_con_test0_3_uncon = sqrt(theSum_uncon) / (b-a+1);
% 
% 
% Avg_Projshapes_con_test0=(Avg_Projshapes_con_test0_1+Avg_Projshapes_con_test0_2+Avg_Projshapes_con_test0_3)/N_Target;
% Avg_Projshapes_con_test0_uncon=(Avg_Projshapes_con_test0_1_uncon+Avg_Projshapes_con_test0_2_uncon+Avg_Projshapes_con_test0_3_uncon)/N_Target;
% 
% Perf_Proj_con_test0=100*(1-(abs(Avg_Projshapes_con_test0_uncon-Avg_Projshapes_con_test0)/abs(Avg_Projshapes_con_test0_uncon)));



%% controlled pertNorm average

a=10000;
b=15000;
[rows, columns] = size(pertNormSource_1e5);
for col = 1 : 2
  theSum = 0;
  for row = a : b
    theSum = theSum + pertNormSource_1e5(row, col);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e5(col) = theSum / (b-a+1);
end
Avg_con=Avg_pertNorm_1e5(2);

%% uncontrol pertNorm average
[rows_rnd, columns_rnd] = size(pertNormSource_1e5_uncon);
for col_rnd = 1 : 2
  theSum = 0;
  for row_rnd = a : b
    theSum = theSum + pertNormSource_1e5_uncon(row_rnd, col_rnd);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e5_uncon(col_rnd) = theSum / (b-a+1);
end
Avg_rnd=Avg_pertNorm_1e5_uncon(2);
%% performance of pertNorm

Perf_1e5=100*(1-(abs(Avg_rnd-Avg_con)/abs(Avg_rnd)));



%% Projshapes average 1e5

a=10000;
b=15000;
[rows, columns] = size(ProjShapeSource_1e5);
N_Target=3;
targetNr1 = 1+11; %add 1 for time column vector
targetNr2 = 1+12;
targetNr3 = 1+13;
%% first target
theSum = 0;
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e5(row,targetNr1)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e5_uncon(row,targetNr1)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e5_1 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e5_1_uncon = sqrt(theSum_uncon) / (b-a+1);
%% second target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e5(row,targetNr2)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e5_uncon(row,targetNr2)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e5_2 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e5_2_uncon = sqrt(theSum_uncon) / (b-a+1);

%% third target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e5(row,targetNr3)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e5_uncon(row,targetNr3)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e5_3 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e5_3_uncon = sqrt(theSum_uncon) / (b-a+1);


Avg_Projshapes_1e5=(Avg_Projshapes_1e5_1+Avg_Projshapes_1e5_2+Avg_Projshapes_1e5_3)/N_Target;
Avg_Projshapes_1e5_uncon=(Avg_Projshapes_1e5_1_uncon+Avg_Projshapes_1e5_2_uncon+Avg_Projshapes_1e5_3_uncon)/N_Target;

Perf_Proj_1e5=100*(1-(abs(Avg_Projshapes_1e5_uncon-Avg_Projshapes_1e5)/abs(Avg_Projshapes_1e5_uncon)));


%% controlled pertNorm average f=1e-3

a=10000;
b=12000;
[rows, columns] = size(pertNormSource_1e3);
for col = 1 : 2
  theSum = 0;
  for row = a : b
    theSum = theSum + pertNormSource_1e3(row, col);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e3(col) = theSum / (b-a+1);
end
Avg_con=Avg_pertNorm_1e3(2);

%% uncontrol pertNorm average
[rows_rnd, columns_rnd] = size(pertNormSource_1e3_uncon);
for col_rnd = 1 : 2
  theSum = 0;
  for row_rnd = a : b
    theSum = theSum + pertNormSource_1e3_uncon(row_rnd, col_rnd);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e3_uncon(col_rnd) = theSum / (b-a+1);
end
Avg_rnd=Avg_pertNorm_1e3_uncon(2);
%% performance of pertNorm

Perf_1e3=100*(1-(abs(Avg_rnd-Avg_con)/abs(Avg_rnd)));

%% Projshapes average 1e3

a=10000;
b=12000;
[rows, columns] = size(ProjShapeSource_1e3);
N_Target=3;
targetNr1 = 1+11; %add 1 for time column vector
targetNr2 = 1+12;
targetNr3 = 1+13;
%% first target
theSum = 0;
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e3(row,targetNr1)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e3_uncon(row,targetNr1)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e3_1 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e3_1_uncon = sqrt(theSum_uncon) / (b-a+1);
%% second target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e3(row,targetNr2)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e3_uncon(row,targetNr2)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e3_2 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e3_2_uncon = sqrt(theSum_uncon) / (b-a+1);

%% third target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e3(row,targetNr3)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e3_uncon(row,targetNr3)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e3_3 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e3_3_uncon = sqrt(theSum_uncon) / (b-a+1);


Avg_Projshapes_1e3=(Avg_Projshapes_1e3_1+Avg_Projshapes_1e3_2+Avg_Projshapes_1e3_3)/N_Target;
Avg_Projshapes_1e3_uncon=(Avg_Projshapes_1e3_1_uncon+Avg_Projshapes_1e3_2_uncon+Avg_Projshapes_1e3_3_uncon)/N_Target;

Perf_Proj_1e3=100*(1-(abs(Avg_Projshapes_1e3_uncon-Avg_Projshapes_1e3)/abs(Avg_Projshapes_1e3_uncon)));
%% controlled pertNorm average f=1e-4

a=10000;
b=15000;
[rows, columns] = size(pertNormSource_1e4);
for col = 1 : 2
  theSum = 0;
  for row = a : b
    theSum = theSum + pertNormSource_1e4(row, col);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e4(col) = theSum / (b-a+1);
end
Avg_con=Avg_pertNorm_1e4(2);

%% uncontrol pertNorm average
[rows_rnd, columns_rnd] = size(pertNormSource_1e4_uncon);
for col_rnd = 1 : 2
  theSum = 0;
  for row_rnd = a : b
    theSum = theSum + pertNormSource_1e4_uncon(row_rnd, col_rnd);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e4_uncon(col_rnd) = theSum / (b-a+1);
end
Avg_rnd=Avg_pertNorm_1e4_uncon(2);
%% performance of pertNorm

Perf_1e4=100*(1-(abs(Avg_rnd-Avg_con)/abs(Avg_rnd)));
%% Projshapes average 1e4

a=10000;
b=15000;
[rows, columns] = size(ProjShapeSource_1e4);
N_Target=3;
targetNr1 = 1+11; %add 1 for time column vector
targetNr2 = 1+12;
targetNr3 = 1+13;
%% first target
theSum = 0;
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e4(row,targetNr1)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e4_uncon(row,targetNr1)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e4_1 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e4_1_uncon = sqrt(theSum_uncon) / (b-a+1);
%% second target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e4(row,targetNr2)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e4_uncon(row,targetNr2)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e4_2 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e4_2_uncon = sqrt(theSum_uncon) / (b-a+1);

%% third target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e4(row,targetNr3)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e4_uncon(row,targetNr3)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e4_3 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e4_3_uncon = sqrt(theSum_uncon) / (b-a+1);


Avg_Projshapes_1e4=(Avg_Projshapes_1e4_1+Avg_Projshapes_1e4_2+Avg_Projshapes_1e4_3)/N_Target;
Avg_Projshapes_1e4_uncon=(Avg_Projshapes_1e4_1_uncon+Avg_Projshapes_1e4_2_uncon+Avg_Projshapes_1e4_3_uncon)/N_Target;

Perf_Proj_1e4=100*(1-(abs(Avg_Projshapes_1e4_uncon-Avg_Projshapes_1e4)/abs(Avg_Projshapes_1e4_uncon)));

%% controlled pertNorm average f=1e-2

a=10000;
b=15000;
[rows, columns] = size(pertNormSource_1e2);
for col = 1 : 2
  theSum = 0;
  for row = a : b
    theSum = theSum + pertNormSource_1e2(row, col);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e2(col) = theSum / (b-a+1);
end
Avg_con=Avg_pertNorm_1e2(2);

%% uncontrol pertNorm average
[rows_rnd, columns_rnd] = size(pertNormSource_1e2_uncon);
for col_rnd = 1 : 2
  theSum = 0;
  for row_rnd = a : b
    theSum = theSum + pertNormSource_1e2_uncon(row_rnd, col_rnd);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e2_uncon(col_rnd) = theSum / (b-a+1);
end
Avg_rnd=Avg_pertNorm_1e2_uncon(2);
%% performance of pertNorm

Perf_1e2=100*(1-(abs(Avg_rnd-Avg_con)/abs(Avg_rnd)));


%% Projshapes average 1e2

a=10000;
b=15000;
[rows, columns] = size(ProjShapeSource_1e2);
N_Target=3;
targetNr1 = 1+11; %add 1 for time column vector
targetNr2 = 1+12;
targetNr3 = 1+13;
%% first target
theSum = 0;
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e2(row,targetNr1)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e2_uncon(row,targetNr1)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e2_1 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e2_1_uncon = sqrt(theSum_uncon) / (b-a+1);
%% second target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e2(row,targetNr2)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e2_uncon(row,targetNr2)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e2_2 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e2_2_uncon = sqrt(theSum_uncon) / (b-a+1);

%% third target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e2(row,targetNr3)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e2_uncon(row,targetNr3)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e2_3 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e2_3_uncon = sqrt(theSum_uncon) / (b-a+1);


Avg_Projshapes_1e2=(Avg_Projshapes_1e2_1+Avg_Projshapes_1e2_2+Avg_Projshapes_1e2_3)/N_Target;
Avg_Projshapes_1e2_uncon=(Avg_Projshapes_1e2_1_uncon+Avg_Projshapes_1e2_2_uncon+Avg_Projshapes_1e2_3_uncon)/N_Target;

Perf_Proj_1e2=100*(1-(abs(Avg_Projshapes_1e2_uncon-Avg_Projshapes_1e2)/abs(Avg_Projshapes_1e2_uncon)));

%% controlled pertNorm average f=1e-6

a=10000;
b=15000;
[rows, columns] = size(pertNormSource_1e6);
for col = 1 : 2
  theSum = 0;
  for row = a : b
    theSum = theSum + pertNormSource_1e6(row, col);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e6(col) = theSum / (b-a+1);
end
Avg_con=Avg_pertNorm_1e6(2);

%% uncontrol pertNorm average
[rows_rnd, columns_rnd] = size(pertNormSource_1e6_uncon);
for col_rnd = 1 : 2
  theSum = 0;
  for row_rnd = a : b
    theSum = theSum + pertNormSource_1e6_uncon(row_rnd, col_rnd);
  end
  % Now get the mean over all values in this column.
  Avg_pertNorm_1e6_uncon(col_rnd) = theSum / (b-a+1);
end
Avg_rnd=Avg_pertNorm_1e6_uncon(2);
%% performance of pertNorm

Perf_1e6=100*(1-(abs(Avg_rnd-Avg_con)/abs(Avg_rnd)));


%% Projshapes average 1e6

a=10000;
b=15000;
[rows, columns] = size(ProjShapeSource_1e6);
N_Target=3;
targetNr1 = 1+11; %add 1 for time column vector
targetNr2 = 1+12;
targetNr3 = 1+13;
%% first target
theSum = 0;
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e6(row,targetNr1)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e6_uncon(row,targetNr1)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e6_1 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e6_1_uncon = sqrt(theSum_uncon) / (b-a+1);
%% second target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e6(row,targetNr2)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e6_uncon(row,targetNr2)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e6_2 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e6_2_uncon = sqrt(theSum_uncon) / (b-a+1);

%% third target
theSum = 0;    
theSum_uncon = 0;
  for row = a : b
    theSum = theSum + (ProjShapeSource_1e6(row,targetNr3)).^2;
    theSum_uncon = theSum_uncon + (ProjShapeSource_1e6_uncon(row,targetNr3)).^2;
  end
  % Now get the mean over all values in this column.
Avg_Projshapes_1e6_3 = sqrt(theSum) / (b-a+1);
Avg_Projshapes_1e6_3_uncon = sqrt(theSum_uncon) / (b-a+1);


Avg_Projshapes_1e6=(Avg_Projshapes_1e6_1+Avg_Projshapes_1e6_2+Avg_Projshapes_1e6_3)/N_Target;
Avg_Projshapes_1e6_uncon=(Avg_Projshapes_1e6_1_uncon+Avg_Projshapes_1e6_2_uncon+Avg_Projshapes_1e6_3_uncon)/N_Target;

Perf_Proj_1e6=100*(1-(abs(Avg_Projshapes_1e6_uncon-Avg_Projshapes_1e6)/abs(Avg_Projshapes_1e6_uncon)));


% pertNormSource=pertNormSource';
% Avg_pertNorm=mean(pertNormSource,2);
% Avg=Avg_pertNorm(2,1);

% figure(1)
%         plot(pertNormSource_con(:,1),pertNormSource_con(:,2),'r',pertNormSource_uncon(:,1),pertNormSource_uncon(:,2),'b')
%         xlim([500,750])
%         legend('pertNorm(controlled)','pertNorm(uncontrolled)')
%         xlabel('t');ylabel('pertNorm');
%         title('controlled and uncontrolled pertNorm')

% figure(2)
%         plot(pertNormSource_uncon(1:10000,1),pertNormSource_uncon(1:10000,2),'k', ...
%         pertNormSource_SISO_tnc(1:10000,1),pertNormSource_SISO_tnc(1:10000,2),'m--',...
%         pertNormSource_SISO_c(1:10000,1),pertNormSource_SISO_c(1:10000,2),'b--',...
%         pertNormSource_tnc(1:10000,1),pertNormSource_tnc(1:10000,2),'m',...
%         pertNormSource_con(1:10000,1),pertNormSource_con(1:10000,2),'b')
%     
%         xlim([100,450])
%         ylim([0,0.25])
%         legend('uncontrolled','TNC(SISO)','Causal(SISO)','TNC(MIMO)','Causal(MIMO)')
%         xlabel('t');ylabel('pertNorm');
%         title('$Non-linear system$')
%         grid on;
%print(gcf,'-dpsc2','pertNorm_NL_3MIMO.eps'); 
figure(1)
        plot(pertNormSource_1e6_uncon(1:15000,1),pertNormSource_1e6_uncon(1:15000,2),'k', ...
        pertNormSource_1e6(1:15000,1),pertNormSource_1e6(1:15000,2),'r')
    
        %xlim([500,600])
        %ylim([0,0.25])
        legend('uncon(f*1e-2)','Causal(f*1e-2)')
        xlabel('t');ylabel('pertNorm');
        title('$Non-linear system$')
        grid on;
        
figure(2)
        plot(ProjShapeSource_1e6_uncon(1:15000,1),ProjShapeSource_1e6_uncon(1:15000,13),'k',...
            ProjShapeSource_1e6(1:15000,1),ProjShapeSource_1e6(1:15000,13),'b')
        %xlim([500,600])
        %ylim([0,1])
        legend('uncontrolled(f*1e-2)','controlled(C)(f*1e-2)')
        xlabel('t');ylabel('Target');
        title('Non-linear system')
        grid on;   
%print(gcf,'-dpsc2','ProjShape_1e6.eps');

figure(3)
        plot(pertNormSource_1e5_uncon(1:15000,1),pertNormSource_1e5_uncon(1:15000,2),'k', ...
        pertNormSource_1e5(1:15000,1),pertNormSource_1e5(1:15000,2),'r')
    
        %xlim([500,600])
        %ylim([0,0.25])
        legend('uncon(f*1e-1)','Causal(f*1e-1)')
        xlabel('t');ylabel('pertNorm');
        title('$Non-linear system$')
        grid on;
        
figure(4)
        plot(ProjShapeSource_1e5_uncon(1:15000,1),ProjShapeSource_1e5_uncon(1:15000,13),'r',...
            ProjShapeSource_1e5(1:15000,1),ProjShapeSource_1e5(1:15000,13),'b')
        %xlim([500,600])
        %ylim([0,1])
        legend('uncontrolled(f*1e-1)','controlled(C)(f*1e-1)')
        xlabel('t');ylabel('Target');
        title('Non-linear system')
        grid on;      
        
figure(5)
        plot(pertNormSource_1e4_uncon(1:15000,1),pertNormSource_1e4_uncon(1:15000,2),'k', ...
        pertNormSource_1e4(1:15000,1),pertNormSource_1e4(1:15000,2),'r')
    
        %xlim([500,600])
        %ylim([0,0.25])
        legend('uncon(f)','Causal(f)')
        xlabel('t');ylabel('pertNorm');
        title('$Non-linear system$')
        grid on;
        
figure(6)
        plot(ProjShapeSource_1e4_uncon(1:15000,1),ProjShapeSource_1e4_uncon(1:15000,13),'r', ...
            ProjShapeSource_1e4(1:15000,1),ProjShapeSource_1e4(1:15000,13),'b')
        %xlim([500,600])
        %ylim([0,1])
        legend('uncontrolled(f)','controlled(C)(f)')
        xlabel('t');ylabel('Target');
        title('Non-linear system')
        grid on;            
        
figure(7)
        plot(pertNormSource_1e3_uncon(1:15000,1),pertNormSource_1e3_uncon(1:15000,2),'k', ...
        pertNormSource_1e3(1:12000,1),pertNormSource_1e3(1:12000,2),'b')
    
        %xlim([500,600])
        %ylim([0,0.25])
        legend('uncon(f*1e+1)','Causal(f*1e+1)')
        xlabel('t');ylabel('pertNorm');
        title('$Non-linear system$')
        grid on;
        
figure(8)
        plot(ProjShapeSource_1e3_uncon(1:15000,1),ProjShapeSource_1e3_uncon(1:15000,13),'r', ...
            ProjShapeSource_1e3(1:12000,1),ProjShapeSource_1e3(1:12000,13),'b')
        %xlim([500,600])
        %ylim([0,1])
        legend('uncontrolled(f*1e+1)','controlled(C)(f*1e+1)')
        xlabel('t');ylabel('Target');
        title('Non-linear system')
        grid on;        
        

figure(9)
        plot(pertNormSource_1e2_uncon(1:20000,1),pertNormSource_1e2_uncon(1:20000,2),'k', ...
        pertNormSource_1e2(1:40000,1),pertNormSource_1e2(1:40000,2),'b')
    
        %xlim([500,600])
        %ylim([0,0.25])
        legend('uncon(f*1e+2)','Causal(f*1e+2)')
        xlabel('t');ylabel('pertNorm');
        title('$Non-linear system$')
        grid on;
        
figure(10)
        plot(ProjShapeSource_1e2_uncon(1:20000,1),ProjShapeSource_1e2_uncon(1:20000,13),'r', ...
            ProjShapeSource_1e2(1:40000,1),ProjShapeSource_1e2(1:40000,13),'b')
        %xlim([500,600])
        %ylim([0,1])
        legend('uncontrolled(f*1e+2)','controlled(C)(f*1e+2)')
        xlabel('t');ylabel('Target');
        title('Non-linear system')
        grid on;                
 %print(gcf,'-dpsc2','Target_MIMO_C_xm2.eps');
 
 
%  figure(11)
%         plot(pertNormSource_con_test1_uncon(1:15000,1),pertNormSource_con_test1_uncon(1:15000,2),'k', ...
%         pertNormSource_con_test1(1:15000,1),pertNormSource_con_test1(1:15000,2),'b')
%     
%         xlim([500,600])
%         %ylim([0,0.25])
%         legend('uncon(test1)','Causal(test1)')
%         xlabel('t');ylabel('pertNorm');
%         title('$Non-linear system$')
%         grid on;
%         
% figure(12)
%         plot(ProjShapeSource_con_test1_uncon(1:15000,1),ProjShapeSource_con_test1_uncon(1:15000,13),'r', ...
%             ProjShapeSource_con_test1(1:15000,1),ProjShapeSource_con_test1(1:15000,13),'b')
%         xlim([500,600])
%         %ylim([0,1])
%         legend('uncontrolled(test1)','controlled(C)(test1)')
%         xlabel('t');ylabel('Target');
%         title('Non-linear system')
%         grid on;
%         
%         
%   figure(13)
%         plot(pertNormSource_con_test3_uncon(1:10000,1),pertNormSource_con_test3_uncon(1:10000,2),'k', ...
%         pertNormSource_con_test3(1:10000,1),pertNormSource_con_test3(1:10000,2),'b')
%     
%         %xlim([500,600])
%         %ylim([0,0.25])
%         legend('uncon(test3)','Causal(test3)')
%         xlabel('t');ylabel('pertNorm');
%         title('$Non-linear system$')
%         grid on;
%         
% figure(14)
%         plot(ProjShapeSource_con_test3_uncon(1:10000,1),ProjShapeSource_con_test3_uncon(1:10000,8),'r', ...
%             ProjShapeSource_con_test3(1:10000,1),ProjShapeSource_con_test3(1:10000,8),'b')
%         %xlim([500,600])
%         %ylim([0,1])
%         legend('uncontrolled(test3)','controlled(C)(test3)')
%         xlabel('t');ylabel('Target');
%         title('Non-linear system')
%         grid on;       
%   figure(15)
%         plot(pertNormSource_con_test4_uncon(1:4000,1),pertNormSource_con_test4_uncon(1:4000,2),'k', ...
%         pertNormSource_con_test4(1:6000,1),pertNormSource_con_test4(1:6000,2),'b')
%     
%         %xlim([500,600])
%         %ylim([0,0.25])
%         legend('uncon(test4)','Causal(test4)')
%         xlabel('t');ylabel('pertNorm');
%         title('$Non-linear system$')
%         grid on;
%         
% figure(16)
%         plot(ProjShapeSource_con_test4_uncon(1:10000,1),ProjShapeSource_con_test4_uncon(1:10000,13),'r', ...
%             ProjShapeSource_con_test4(1:10000,1),ProjShapeSource_con_test4(1:10000,13),'b')
%         %xlim([500,600])
%         %ylim([0,1])
%         legend('uncontrolled(test4)','controlled(C)(test4)')
%         xlabel('t');ylabel('Target');
%         title('Non-linear system')
%         grid on;       