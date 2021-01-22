function  y =  chilbert(x,overSample)
    if ~exist('overSample'); overSample=10;end
%     disp(['Oversample ' num2str(overSample)]);
    n=length(x);
    xr=zeros(n*overSample,1);
    xi=zeros(n*overSample,1);
    
    y=zeros(size(x));
    for i=1:size(x,1)
        for j=1:size(x,2)
            xr(1:n)= squeeze(real(x(i,j,:)));
            xi(1:n)= squeeze(imag(x(i,j,:)));
            ytmp  = imag(hilbert(squeeze(xr)))+1i*imag(hilbert(squeeze(xi)));
            y(i,j,:) = ytmp(1:n);
        end
    end
end
    