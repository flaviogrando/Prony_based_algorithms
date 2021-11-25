% Teste - Método de Prony
% Flavio L. Grando


clc
clear all
close all

% INICIALIZAÇÕES

f0 = 50;           % Frequência fundamental
N = 256;           % Número de amostras
Fs = 12800;         % Frequência de amostragem    
Ts = 1/Fs;         % Período de amostragem

%h = 10; % numero de componentes de frequência

amort = -3:0.1:3;

for w=1:length(amort)

%amort(w) = w/10;


s = zeros(1,N*50);    % sinais 
% st = zeros(1,N);   % sinal total (soma das multiplas frequências)
% A = zeros(1,N);    % amplitudes
% phi = zeros(1,N);  % angulos
% f = zeros(1,N);    % frequencias
% sigma = zeros(1,N);% fator de amortecimento

% Definindo componentes de frequência harmônica 'h'
% Amplitudes A(.+1)
A(1) = 0;   % nivel cc
A(2) = 1.1;   % fundamental
A(3) = 0;   %
A(4) = 0.0; % 3rd harmonic
A(5) = 0;
A(6) = 0.0;
A(7) = 0;
A(8) = 0.0;

% fases
phi(2) = 0;
phi(4) = 1.8;
phi(6) = 1.16;
phi(8) = 3.91;




% fator de amortecimento
sigma(2) = amort(w);
sigma(3) = 2;
sigma(4) = 2;
sigma(6) = 2;
sigma(8) = 2;

% ------------------------------------------------------------------------
% GERAÇÃO DO SINAL 

for h=1:length(A); 
        
    f(h)=(h-1)*f0;  % componentes de freq harmônicas
    
    for k=0:N*50-1;    
        % equação (1) - sinal
        s(k+1) = s(k+1)+A(h)*exp(sigma(h)*k*Ts)*cos((phi(h)+2*pi*f(h)*k*Ts));
        %s(k) = s(k)+(A(h)/2)*exp(1i*(phi(h)+2*pi*f(h)*k*Ts))+(A(h)/2)*exp(-1i*(phi(h)+2*pi*f(h)*k*Ts)); 
    end
    %st = st + s;
end

% figure(1)
% plot(s)
% grid on

% % ------------------------------------------------------------------------
% % APLICAÇÃO DO MÉTODO DE PRONY
% 
% % %----------prony method------------% 
% % %written by Foad Shiri 
% % %----------parameter set up--------% 
% % close 
% % global order N Z input_signal; 
% % %downsample();  
% % % fs=100; 
% % % T=1/fs; 
% T=Ts; 
% % % t=[5.001761:0.01:9.991761]; 
% % % % input_signal=ias; 
% % order=22; 
%  order=N-1; 
% % % p=load('Results.txt'); 
% % % e=p(:,2); 
% % load 'test2.mat' 
% % input_signal=p'; 
% input_signal = s; 
% % N=length(input_signal);
% t=[0*Ts:Ts:(N-1)*Ts];%.*Ts;
% 
% %--------------step1---------------% 
% X=[]; 
% m=N-order; 
% step=order; 
% for d=1:order 
%     for j=1:m 
%         X(d,j)=input_signal(step-1+j); 
%     end 
%     step=step-1; 
% end 
% X=X'; 
% z=[]; 
% for l=1:(N-order) 
%     z(l,1)=-input_signal(1+order); 
% end 
% theta=pinv(X)*z; 
% %---------------step2---------------% 
% LPM=[1 -theta']; 
% rootz=roots(LPM); 
% Freq=imag(log(rootz))/(2*pi*T); 
% damping_factor=log(abs(rootz))/T; 
% %----------------step3--------------% 
%  Z=[]; 
%  for k=1:N 
%      for m=1:order 
%          Z(k,m)=rootz(m)^(k-1); 
%      end 
%  end 
%  V=[]; 
%  for n=1:N 
%      V(n)=input_signal(n); 
%  end 
%  V=V'; 
%  H=pinv(Z)*V; 
%  phase_rad=angle(H); 
%  amplitude=abs(H); 
% %  result=[Freq  damping_factor  phase_rad  amplitude]; 
% %  disp('  freq   damping_factor   phase_rad   amplitude'); 
% %  disp(result) 
%  %------------display the prony estimation-----------% 
%  y=0; 
%  for n=1:2:(order/1); 
%      y=y+2*amplitude(n)*exp(damping_factor(n)*t).*cos(2*pi*Freq(n)*t+phase_rad(n));
%      
%  end 
 
%  figure, hold on, grid on
%  plot(t,input_signal(1:256),'k') 
%  ylabel('Magnitude(A)') 
%  xlabel('Time(s)') 
%  plot(t,y,':kx') 
%  legend('Current Signal','Prony Estimation'); 
%  %xlim([0 length(y)]); 
 
 %------------------------------------------------------------
 % APLICA FFT
 
Sa = s;
 
 for i=1:50
    
    % índices (janelamento)
    inicio(i) = (i-1)*N+1;
    fim(i) = i*N;
    
    spectro(i,:) = fft( Sa(inicio(i):fim(i)));

    amp_fft(i,:) = 2*abs(spectro(i,:))/N;
    
    amp_ref(i) = Sa(inicio(i));
    
end

% figure, hold on, grid on
% plot(amp_ref, '*-')
% plot(amp_fft(:,2), 'o-')
% legend('Ref','FFT')
% title('amplitudes')


for i=1:50-1
   desvio_amp_fft(i) = (amp_fft(i+1,2) - amp_fft(i,2))/0.02;
   desvio_amp_ref(i) = (amp_ref(i+1) - amp_ref(i))/0.02;
end


% figure, hold on, grid on
% plot(desvio_amp_fft, 'o-')
% plot(desvio_amp_ref, '*-')
% title('Derivada da amplitude'), legend('Ref','FFT')


deriv_ref(w) = desvio_amp_ref(1);
deriv_fft(w) = desvio_amp_fft(1);


end

% figure, hold on, grid on
% plot(-amort, 'o-')
% title('fator de amortecimento')%, legend('Ref','FFT')

figure, hold on, grid on
plot(amort)
plot(deriv_ref)
plot(deriv_fft)
title('Fator de amort e Derivada da amplitude'), legend('Amortecimento','D-Ref','D-FFT')

figure, hold on, grid on
plot(amort - deriv_ref, 'o-')
plot(amort - deriv_fft, '*-')
title('Erros'), legend('D-Ref','D-FFT')


coef = polyfit(amort,deriv_fft,2)