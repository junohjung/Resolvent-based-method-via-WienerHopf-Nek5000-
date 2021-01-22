function [HGs] = getHGs(DATA,N,P,omega0,tol)
    disp('Computing G,g,H,h matrices, and Wiener-Hopt decmopositions, from data...');clock=tic();
    if ~exist('tol');tol=1e-6;end
    
    ny = size(DATA.estYhat,1);
    na = size(DATA.actYhat,2);
    
    
    nt = length(DATA.t);
    dt = DATA.t(2)-DATA.t(1);
    df = DATA.freq(2)-DATA.freq(1);
     FFT = @(x) ifftshift(ifft( fftshift(    x   ,3),[],3)   ,3)*dt*nt; % fft defined for conventions used;
    IFFT = @(x) ifftshift( fft( fftshift(    x   ,3),[],3)   ,3)*df; % fft defined for conventions used;
    
    % Compute "G"  and "g" 
    
    if length(N)==1 && N(1)<0 
        N=eye(ny)*max(abs(DATA.estYhat(:)))*abs(N);
    end
    if length(P)==1 && P(1)<0 
        P=eye(na)*max(abs(DATA.actZhat(:)))^2*abs(P);
    end
    
    for iw = 1:nt
        HGs.Ghat(:,:,iw) =  (DATA.estYhat(:,:,iw)+DATA.estYhat(:,:,iw)')/2 + N; % removes the non-hermian part of y_esthat
        HGs.ghat(:,:,iw) =  DATA.estZhat(:,:,iw);
    end
    [HGs.Gminushat,HGs.Gplushat]=WienerHopfDec(DATA.freq*2*pi,conj(HGs.Ghat),omega0,[],tol); % reverse factorization
    HGs.Gminushat  = conj(HGs.Gminushat); HGs.Gplushat = conj(HGs.Gplushat);

    % Compute "H"  and "h" 
    for iw = 1:nt
        HGs.hhat(:,:,iw)         = -DATA.actZhat(:,:,iw)';
        HGs.Hhat(:,:,iw)         =  DATA.actZhat(:,:,iw)'*DATA.actZhat(:,:,iw) + P ;
    end    
    [HGs.Hplushat,HGs.Hminushat]=WienerHopfDec(DATA.freq*2*pi,HGs.Hhat,omega0,[],tol);         
    
    % Compute "Ray"  and "h" 
    for iw = 1:nt
        HGs.Ray(:,:,iw)  = DATA.actYhat(:,:,iw);
    end
    
    HGs.G      = IFFT(HGs.Ghat);
    HGs.Gplus  = IFFT(HGs.Gplushat);
    HGs.Gminus = IFFT(HGs.Gminushat);
    HGs.H      = IFFT(HGs.Hhat);
    HGs.Hplus  = IFFT(HGs.Hplushat);
    HGs.Hminus = IFFT(HGs.Hminushat);
    HGs.h      = IFFT(HGs.hhat);
    HGs.g      = IFFT(HGs.ghat);
    
    HGs.t=DATA.t;
    HGs.freq=DATA.freq;

    disp(['Completed in ' num2str(toc(clock)) 's']);
