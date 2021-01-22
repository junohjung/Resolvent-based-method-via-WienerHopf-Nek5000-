%% Script to compute optimal control using Wiener-Hopf formalism, based on a DNS readings
% All frequency dependent matrices and vectors are stored with the
% frequency being the last index.

clear all;
close all;
set(0,'defaulttextInterpreter'  ,'latex');
set(0,'defaultlegendInterpreter','latex'); 
set(0,'defaultAxesFontsize'     ,14) ;
set(0,'defaultLegendFontsize'   ,16) ;
set(0,'defaultLineLineWidth'    ,3 ) ;
set(0,'DefaultLineMarkerSize'   ,8 ) ;


initiClock=tic();
%
computeGains=true;
% Set up some functions to help plotting
sortX = @(x) ifftshift( x-(x>x(end)/2)*x(end)); % Makes a [0,x] array into a [-x/2,x/2].
sortY = @(x) ifftshift( x,2);
pfun = @(x) reshape(x,size(x,1)*size(x,2),size(x,3)); % Function for plot Tensors
getFreqVec = @(t) [((0:length(t)-1) - (0:length(t)-1 > (length(t)-1)/2 )*length(t))/(t(end)-t(1)) ];
%% Defines indexes for sensors, actuators and targets (corresponding to shape index in the Nek code
iy = [ 1:3  ]; ny = length(iy);
ia = [ 4:6 ]; na = length(ia);
iz = [ 11:13  ]; nz = length(iz);

actuatorDataFrom='direct';% adjoint : target and sensor adjoint runs, direct : actuator direct runs
computeTargetsCSD = false;

folder = './'; %folder where the results are stored.

% Discretization used for control Law. Input Tmax and dt. 
Tmax = 800;
dt   = 0.01;
nt   = round(Tmax/dt); 
Tmax = nt*dt;
t    = (0:nt-1)*dt;
ts   = t - nt/2*dt;
df   = 1/t(end);
freq = getFreqVec (t);
freqS= fftshift(freq);

%% Loads Data
[DATA,~]= readData(ts,folder,iy,ia,iz,actuatorDataFrom,computeTargetsCSD);
%%

figure;
    subplot(221)
        plot(DATA.freq,abs(pfun(DATA.estYhat)),'b',DATA.freq,abs(pfun(DATA.estZhat)),'r')
        xlim([-1,1]*1)
        xlabel('$\omega$');ylabel('Sensor/Target');
    subplot(222)
        plot(DATA.t,(pfun(DATA.estY)),'b',DATA.t,(pfun(DATA.estZ)),'r')
        xlim([-10,100])
        xlabel('$t$');ylabel('Sensor/Target');
    subplot(223)
        plot(DATA.freq,abs(pfun(DATA.actYhat)),'b',DATA.freq,abs(pfun(DATA.actZhat)),'r')
        xlim([-1,1]*1)
        xlabel('$\omega$');ylabel('Sensor/Target');
    subplot(224)
        plot(DATA.t,(pfun(DATA.actY)),'b',DATA.t,(pfun(DATA.actZ)),'r')
        xlim([-10,100])
        xlabel('$t$');ylabel('Sensor/Target');

%% Get matrices
HGs = getHGs(DATA,-1e-2,-1e-2,-1i,1e-9); 
%%
figure;
    subplot(211)
        plot(HGs.freq,abs(pfun(HGs.Hhat)),'k',HGs.freq,abs(pfun(HGs.Hminushat)),'b',HGs.freq,abs(pfun (HGs.Hplushat)),'r')
        xlim([-1,1]*2)
        xlabel('$\omega$');ylabel('$|\hat{H}_\pm|$');
    subplot(212)
        yyaxis left
        plot(HGs.t,pfun(HGs.Hminus),'-',HGs.t,pfun(HGs.H),'k-')
        ylabel('$|\hat{H}_-|$')
        yyaxis right
        plot(HGs.t,pfun (HGs.Hplus),'-')
        ylabel('$|\hat{H}_+|$')
        xlim([-1,1]*100)
        xlabel('t');
figure
    subplot(211)
        plot(HGs.freq,abs(pfun(HGs.Ghat)),'k',HGs.freq,abs(pfun(HGs.Gminushat)),'b',HGs.freq,abs(pfun (HGs.Gplushat)),'r:')
        xlim([-1,1]*1)
        xlabel('$\omega$');ylabel('$|\hat{G}_\pm|$');
    subplot(212)
        yyaxis left
        plot(HGs.t,(pfun(HGs.Gminus)),'-',HGs.t,(pfun(HGs.G)),'k-')
        ylabel('$|\hat{G}_-|$');
        yyaxis right
        plot(HGs.t,(pfun (HGs.Gplus)),'-')
        ylabel('$|\hat{G}_+|$');
        xlim([-1,1]*10)
        xlabel('t');
%% Compute Estimation Kernels
Tu= getEstimationKernels(HGs);
%% Plot Estimation Kernels
figure;
    subplot(211)
        plot(Tu.freq,abs(pfun(Tu.nchat)),'-b'); hold on;
        plot(Tu.freq,abs(pfun(Tu.tnchat)),':r');
%         plot(Tu.freq,abs(pfun(Tu.chat)),'--k');
        xlim([-1,1]*2)
        xlabel('$\omega$');ylabel('$|H|$');
    subplot(212)
        plot(Tu.t,(pfun(Tu.nc)),'-b'); hold on;
        plot(Tu.t,(pfun(Tu.tnc)),':r');
%         plot(Tu.t,(pfun(Tu.c)),'--k');
%         xlim([-1,1]*2)
        xlabel('$t$');ylabel('$|H|$');
%% Compute Control Kernels
Gamma = getControlKernels(HGs);
               
figure;
    subplot(211)
        plot(Gamma.freq,abs(pfun(Gamma.nchat)),'b',Gamma.freq,abs(pfun (Gamma.tnchat)),'r:',Gamma.freq,abs(pfun (Gamma.chat)),'k--')
        xlim([-1,1]*2)
        title('Control Kernels');
        xlabel('$\omega$');ylabel('$|\hat{\Gamma}_u|$');
    subplot(212)
        plot(Gamma.t,(pfun(Gamma.nc)),'b',Gamma.t,(pfun(Gamma.tnc)),'r:',Gamma.t,(pfun(Gamma.c)),'k--')
        xlim([-10,50])
        xlabel('t');ylabel('$\Gamma$');
        title('Control Kernels');
        
fig_ker=figure;
fig_ker.Position(3:4)=[1400,400];
        plot(Gamma.t,(pfun(Gamma.c)),'k-',Gamma.t,(pfun(Gamma.nc)),'--r',Gamma.t,(pfun(Gamma.tnc)),'g:'); 
        hold on
%         plot(Gamma.t,(pfun(Gamma.cfb)),'k--',Gamma.t,(pfun(Gamma.tncfb)),'g:')
        xlim([-15,50])
        xlabel('$\tau$');ylabel('$\Gamma$');
        legend('Causal','Non-Causal','TNC')
filename=['Kernel_iy=' strrep(num2str(iy),' ','_') '_ia=' strrep(num2str(ia),' ','_') '_iz=' strrep(num2str(iz),' ','_') ];
saveas(fig_ker,[filename '.fig']);
saveas(fig_ker,[filename '.eps'],'epsc');

        
%% Estimate performance
if computeTargetsCSD
    Czz_nc=zeros(nz,nz,nt);
    Czz_tnc=zeros(nz,nz,nt);
    Czz_c=zeros(nz,nz,nt);
    psd_un=zeros(nt,1);
    psd_nc=zeros(nt,1);
    psd_tnc=zeros(nt,1);
    psd_c=zeros(nt,1);

    for i=1:nt
        Cz1z1 = DATA.tarZhat(:,:,i);
        H=HGs.Hhat(:,:,i);
        h=HGs.hhat(:,:,i);
        G=HGs.Ghat(:,:,i);
        g=HGs.ghat(:,:,i);

        psd_un(i)=trace(Cz1z1);
        gamma = Gamma.nchat(:,:,i);
        Czz_nc(:,:,i) = Cz1z1 + h'*gamma*G*gamma'*h-h'*gamma*g'-g*gamma'*h;
        psd_nc(i)=trace(Czz_nc(:,:,i));

        gamma = Gamma.tnchat(:,:,i);
        Czz_tnc(:,:,i) = Cz1z1 + h'*gamma*G*gamma'*h-h'*gamma*g'-g*gamma'*h;
        psd_tnc(i)=trace(Czz_tnc(:,:,i));

        gamma = Gamma.chat(:,:,i);
        Czz_c(:,:,i) = Cz1z1 + h'*gamma*G*gamma'*h-h'*gamma*g'-g*gamma'*h;
        psd_c(i)=trace(Czz_c(:,:,i));
    end
    
    figrms=figure
        plot(2*pi*freqS,psd_un /(2*pi),'k',...
             2*pi*freqS,psd_c  /(2*pi),'-b', ...
             2*pi*freqS,psd_nc /(2*pi),':r', ...
             2*pi*freqS,psd_tnc/(2*pi),'--g');
         legend('Uncontroled','Causal','Non-Causal','TNC','Location','Southwest');
         xlabel('$\omega$')
         ylabel('PSD')
         xlim([-1,1]*2)
     
    filename=['CSD_iy=' strrep(num2str(iy),' ','_') '_ia=' strrep(num2str(ia),' ','_') '_iz=' strrep(num2str(iz),' ','_') ];
    saveas(figrms,[filename '.fig']);
    saveas(figrms,[filename '.eps'],'epsc');
  
    rms_un =trapz(2*pi*freqS,psd_un /(2*pi));
    rms_c  =trapz(2*pi*freqS,psd_c  /(2*pi));
    rms_nc =trapz(2*pi*freqS,psd_nc /(2*pi));
    rms_tnc=trapz(2*pi*freqS,psd_tnc/(2*pi));

    report= [ ...
    sprintf('---------- Performances ------\n' ) ... 
    sprintf('%9s\t%9s\t%9s\t%9s\n','uncon','causal','non-causal', 'tnc'), ...
    sprintf('%9f\t%9f\t%9f\t%9f\n',[rms_un,rms_c,rms_nc,rms_tnc]), ...
    sprintf('%9f\t%9f\t%9f\t%9f\n',[rms_un,rms_c,rms_nc,rms_tnc]/rms_un*100)];


    filename=['RMS_iy=' strrep(num2str(iy),' ','_') '_ia=' strrep(num2str(ia),' ','_') '_iz=' strrep(num2str(iz),' ','_') ];
     fid = fopen([filename '.txt'],'w');
    fprintf(fid ,'%s',report)
    fprintf('%s\n',report);

end

%% Export Control Law



dt_ext=0.05;
t_ext= 0:dt_ext:50;
CL    = interpTensors(Gamma.t,Gamma.cfb,t_ext);
dlmwrite('ControlLawParams.txt',[ny na length(t_ext) dt_ext],' ');    %%

for i = 1:na
    name = sprintf('ControlLaw%02.0f.txt',i);
    dlmwrite('ControlLawParams.txt',name,'delimiter','','-append');
    CLi= reshape(CL(i,:,:),ny,length(t_ext));
    dlmwrite(name,real(CLi.'),'delimiter',' ','precision','%.6f');
end
%%
disp(['Total time: ' num2str(toc(initiClock)) 's']);

