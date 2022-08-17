clear all;
close all;
nPoint=6000;
Fs=1000;
bits=[1 0 1 0 1 1
      1 0 0 1 0 0 
      1 1 0 1 0 1 
      0 1 0 1 1 0];
% bits=[1 1 1 1 1 1
%       0 0 0 0 0 0 
%       1 1 1 1 1 1 
%       0 0 0 0 0 0];
  
divB=2*(reshape(bits,4,6)')-1;

carriers=[];
bw_c=[];
for col = 1:4
    signal=[];
    for row = 1:6
        if divB(row,col)==1
            se=ones(1,10);
        else divB(row,col)=-1;
            se=-ones(1,10);
        end
     signal=[signal se];
    end
    
   t1=1/1000:1/1000:Fs*length(bits)*(1/Fs);
    
   t=linspace(0,length(bits),length(bits)*1000);
    N=length(t);    %Bits size

    Nsb=N/length(signal); % Sample Size for each bit
    disp(Nsb);
    dd=repmat(signal',1,Nsb);
    
    bb=repmat(signal',1,Nsb);
    dw=dd'; 
    dw=dw(:)';
    bw=bb';
    bw=bw(:)';
    w=sin(2*pi*col*t);
    
    bw_c=[bw_c;bw];
    signal=bw.*w; % Modulated Signal
    
carriers=[carriers;signal];
end

subplot(4,1,1);
plot(t,carriers(1,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);hold on;
plot(t1,bw_c(1,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);

subplot(4,1,2);
plot(t,carriers(2,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);hold on;
plot(t1,bw_c(2,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);


subplot(4,1,3);
plot(t,carriers(3,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);hold on;
plot(t1,bw_c(3,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);


subplot(4,1,4);
plot(t,carriers(4,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);hold on;
plot(t1,bw_c(4,:),'lineWidth',1);grid on;ylim([-1.1 1.1]);

figure;
for i = 1:4
    plot(carriers(i,:),'lineWidth',1);grid on;hold on;
end
    plot(sum(carriers),'lineWidth',2);


carXfReal=[];
carXfImag=[];
figure;

for i = 1:4
    carXfReal=[carXfReal; real(ifft((carriers(i,:)),nPoint))/numel(t)];
    carXfImag=[carXfImag; imag(ifft((carriers(i,:)),nPoint))/numel(t)];
    
    F = linspace(-Fs/2 , Fs/2 , numel(t));
    
    
    %subplot(8,1,2*i-1);
    subplot(4,1,i);
    stem(F,ifftshift(carXfReal(i,:)),'o','MarkerSize',5,'lineWidth',1);grid on;xlim([-5 5]);hold on;
    plot(F,ifftshift(carXfReal(i,:)),'lineWidth',1);grid on;xlim([-5 5]);
    %subplot(8,1,i*2);
    %stem(F,carXfImag(i,:),'lineWidth',1);grid on;;xlim([-5 5]);
end
figure;
for i = 1:4
plot(F,ifftshift(carXfReal(i,:))+ifftshift(carXfImag(i,:)));xlim([-5 5]);hold on;
end

  figure;
  plot(F,ifftshift(sum(carXfReal)));xlim([-5 5]);

  
figure;
subplot(2,1,1);
xF_R= (real(ifft(sum(carriers),nPoint)));
plot(F,ifftshift(xF_R),'lineWidth',1);
subplot(2,1,2);
xF_I= (imag(ifft(sum(carriers),nPoint)));
plot(F,ifftshift(xF_I),'lineWidth',1);

figure; 
subplot(2,1,1);
IxF_R= (real((fft(xF_R,nPoint))));
t=linspace(0,6,numel(IxF_R));
plot(t,fftshift(IxF_R),'lineWidth',1);hold on ;
subplot(2,1,2);
IxF_I= (imag(fft(xF_I,nPoint)));
plot(t,fftshift(IxF_I),'lineWidth',1);

xFOut=IxF_R+IxF_I;
%xFOut=-1*xFOut(1690:end);
figure;
for i = 1:4
subplot(4,1,i);
t=linspace(0,6,numel(xFOut));
 

%plot(t,cumtrapz((IxF_R+IxF_I).*cos(2*pi*i*t)),'lineWidth',1);hold on ;
plot(t,-1*fftshift(xFOut),'lineWidth',1);hold on ;

t=linspace(0,6,numel(carriers(i,:)));
%plot(t,carriers(i,:),'lineWidth',1);
plot(t,sum(carriers),'lineWidth',1);


end
