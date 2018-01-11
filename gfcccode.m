%[data,fs] = audioread('s1_bgaa9a.wav');
data=data(1600:8325);
b=[1 -0.97];a=1;
figure(1);
subplot(211);
fs=25000;
plot((0:6725)/25000,data);
xlabel('time'); ylabel('amplitude');
subplot(212);
y=filter(b,a,data);
plot((0:6725)/25000,y);
xlabel('time'); ylabel('amplitude');
k=0;
for p=1:8   
    dataframe=y(1+k*750:750+k*750);
    datawin=dataframe.*hamming(750);
    figure(2);
    subplot(211);
    plot((0:749)/25000,datawin);
    title('windowed frame - mid part');xlabel('time for windowed frame in seconds');ylabel('magnitude');
    datawinfft=power(abs(fft(datawin,512)),2);
    subplot(212);
    N = length(abs(datawinfft));
    fgrid = fs*(0:(N-1))';
    Xabs = datawinfft(1:floor(N/2));
    fgrid = fgrid(1:floor(N/2));
    plot(fgrid,Xabs);
    xlabel('frequency');ylabel('magnitude');title('fft response');
    pi=3.14;
    numchans=14;
    cfspre = (linspace((21.4*log10(4.37e-3*20+1)),(21.4*log10(4.37e-3*4000+1)),numchans));
    fc=(10.^(cfspre/21.4)-1)/4.37e-3;
    %fc=[100,200,300,400,500,600,700,800,900,1000,1100,1200,1500,2000,2500,3000,4000,5000,5500,6000,6500,7000,7500,8000];
    erb=24.7+(fc/9.26);
    b=1.019*erb;
    o=3;
    a=6;
    t=0:1/fs:717/fs;
    for n=1:14
         g(n,:)=a*power(t,o-1).*exp(-2*pi*b(n)*t).*cos(2*pi*fc(n)*t);
    end
    %gainvalue=[0.0112,0.112,0.237,0.398,0.501,0.631,0.794,0.891,1,1,1,1,1.5,2.51,2.51,2.51,1.78,0.794,0.355,0.282,0.282,0.355,0.224,0.112,0.178];

    figure(3);
    for j=1:14
        gfft(j,:)=abs(fft(g(j,:),512));
        plot(gfft(j,1:256));hold on;
        m(j) = max(abs(gfft(j,:)));
        gfft1(j,:)=abs(fft(g(j,:)./m(j),512));
       plot(gfft1(j,1:256));hold on;
    end
    confftsum=0;
    for i=1:14
        con(i,:)=((Xabs'.*gfft1(i,1:256)));
        %confft(i,:)=abs(fft(con(i,:),512));
        confft=con;
        confftsum=confftsum+(confft(i,1:256));
    end
    gfcc(p,:)=dct(log10(confftsum),14);
    figure(4);
    plot(gfcc(p,:)); hold on;
    k=k+1;
end



