    %% To finish sound 
    T=419*2; T2=419*4; t=[1:.3:T]; t2=[1:.3:T2];
    S2=(sin(t2*pi/T2).^.4+sin(t2*pi/T2).^4).*sin((t2+sin(t2/2)*.2)*1.5)/2;
    soundsc(S2,3*8192);soundsc(S2,3*8192);soundsc(S2,3*8192);

