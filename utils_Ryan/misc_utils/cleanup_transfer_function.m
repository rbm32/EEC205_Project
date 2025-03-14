function [H,h] = cleanup_transfer_function(H,t,tmin,tmax)
    h = real(ifft(ifftshift(H)));
    
%     figure(); hold on; grid on;
%     plot(t.*1e9, h,'k');
    
    [~,indmin]=min(abs(t-tmin));
    [~,indmax]=min(abs(t-tmax));
    h(1:indmin)=0;
    h(indmax:end)=0;
    H=fftshift(fft(h));
%     plot(t.*1e9, h,'b');
%     xlabel("Time (ns)");
%     ylabel("Mag");
%     legend("Original","Cleaned");
%     drawnow;
    
    
    