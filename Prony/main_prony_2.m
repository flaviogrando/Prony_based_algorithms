%%% TESTES ESTIMADOR DFT-PRONY

clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');
% amort = -3:0.1:3;
% 
% 
%  for w=1:length(amort)
% 
% amort_vetor = ones(50)*amort(w);   % Vetor de referência para amortecimento

% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 15;      % número de amostras
m = 15;      % passo de janelamento
Fs = 750;   % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 60;    % nível de ruído em dB ( 0 = sem ruído)
T = 0.5;        % tempo total (em segundos)
param = 5;   % amort(w);    % freq. de modulação, por exemplo
type = 1;     % seleciona tipo de teste (ver gerador_sinais.m)

freqs = [50];        % Vetor de frequências

phases = zeros(3,2);
% phases = [0 0;           % vetor de fases
%          -2*pi/3 0;
%          +2*pi/3 0];    
phases(1,1) = deg2rad(0);
phases(2,1) = deg2rad(-120);
phases(3,1) = deg2rad(+120);

amps = [1 0.3;           % vetor de amplitudes
        1 0.3;
        1 0.3]; 

% ------------------------------------------------------------------------------------------------------------------
% GERAÇÃO DO SINAL
[Va, Vb, Vc, t, refs, mod] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param);

% figure, hold on, grid on
% %plot(t,mod*10)
% plot(t,Va, 'o-')
% plot(t,Vb)
% plot(t,Vc)
% legend('Va','Vb','Vc'), title('Sinais')
% xlabel('tempo (s)'), ylabel('amplitude')

%-------------------------------------------------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N*2, m);

% % Plot amostragem (ciclo a ciclo)
% w = 5;   % numero de janelas no plot
% figure
% for i=1:w
%     subplot(w,1,i), hold on, grid on
%     plot(t_seg(i,:),Va_seg(i,:),'o-')
% end

% figure %Plot referências
% subplot(3,1,1), hold on, grid on
% plot(t_seg(:,1), ref_seg(1,:), 'o-')
% xlabel('tempo (s)'), ylabel('freq (Hz)'), title('Referências (segmentada)')
% subplot(3,1,2), hold on, grid on
% plot(t_seg(:,1),ref_seg(2,:), 'o-')
% plot(t_seg(:,1),ref_seg(3,:), 'o-')
% plot(t_seg(:,1),ref_seg(4,:), 'o-')
% xlabel('tempo (s)'), ylabel('amp (pu)')
% subplot(3,1,3), hold on, grid on
% plot(t_seg(:,1),ref_seg(5,:), 'o-')
% plot(t_seg(:,1),ref_seg(6,:), 'o-')
% plot(t_seg(:,1),ref_seg(7,:), 'o-')
% xlabel('tempo (s)'), ylabel('fase (rad)')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL DFT
[Sa, Sb, Sc, Aa, Ab, Ac, phia, phib, phic] = estimador_dft(Va_seg, Vb_seg, Vc_seg, Fs, f0);

granN = (Fs/N);
gradeN = 0:granN:granN*(N-1);
bin = round(granN/f0)+1;       % local da componente fundamental no espectro

% % Plot espectro
% figure
% hold on, grid on
% stem(gradeN,2*(Sa(1,:))/N)
% % stem(abs(Sa))
% % stem(abs(Sa))
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Frequência (Hz)')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL - DFT RECURSIVA

dim = length(Va)-N;

Sar = zeros(1,dim);
Sbr = zeros(1,dim);
Scr = zeros(1,dim);

theta = 2*pi/N;
Sa(1,:) = fft(flip(Va_seg(1,:)));
Sb(1,:) = fft(flip(Vb_seg(1,:)));
Sc(1,:) = fft(flip(Vc_seg(1,:)));

Sar(1) = 2*Sa(1,bin)'/N;
Sbr(1) = 2*Sb(1,bin)'/N;
Scr(1) = 2*Sc(1,bin)'/N;

for  r=1:length(Va)-N;
    
    %if r==1
    %    Sar(r+1) = exp(-1i*theta)*2*Sa(1,bin)'/N + (2/N)*(Va(N+r) - Va(r))*exp(-1i*r*theta);
        %Sar(r) = exp(-1i*theta)*2*Sa(1,bin)'/N + (2/N)*(Va(N+r))*exp(-1i*r*theta);
    %else  
        Sar(r+1) = Sar(r) + (2/N)*(Va(N+r) - Va(r))*exp(-1i*r*theta);
        Sbr(r+1) = Sbr(r) + (2/N)*(Vb(N+r) - Vb(r))*exp(-1i*r*theta);
        Scr(r+1) = Scr(r) + (2/N)*(Vc(N+r) - Vc(r))*exp(-1i*r*theta);
        %Sar(r) = Sar(r-1) + (2/N)*(Va(N+r))*exp(-1i*r*theta);
    %end

end

% figure, hold on, grid on
% plot(abs(Sar))
% plot((Aa(:,bin)))
% 
% figure, hold on, grid on
% plot(angle(Sar), '*-')
% %plot(phia(:,bin), 'o-')

%-------------------------------------------------------------------------------------------------------------------
% DECIMAÇÃO (para algoritmo recursivo)

[num_canais,num_wins] = size(ref_seg);              % tamanho total  
%num_wins = fix(dim/m);      % numero de janelas

% aloca fasores do primeiro ciclo
Aar = ones(num_wins,1)*Aa(1,bin);
Abr = ones(num_wins,1)*Ab(1,bin);
Acr = ones(num_wins,1)*Ac(1,bin);
phiar = ones(num_wins,1)*phia(1,bin);
phibr = ones(num_wins,1)*phib(1,bin);
phicr = ones(num_wins,1)*phic(1,bin);

for n=1:num_wins-1
    Aar(n+1) = abs(Sar(m*n+1));
    Abr(n+1) = abs(Sbr(m*n+1));
    Acr(n+1) = abs(Scr(m*n+1));
    phiar(n+1) = angle(Sar(m*n+1)) + theta;
    phibr(n+1) = angle(Sbr(m*n+1)) + theta;
    phicr(n+1) = angle(Scr(m*n+1)) + theta;
end

% figure, hold on, grid on
% plot(Aar, '*-')
% plot((Aa(:,bin)), 'o-')
% 
% figure, hold on, grid on
% plot(phiar, '*-')
% plot(phia(:,bin), 'o-')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL - MODELO IEEE

Ts = 1/Fs;

Nf = 28; % N = (15 – 1) × 2 = 28


k=-Nf/2:Nf/2;

W = (1 - (2/(Nf+2))*abs(k));

figure
plot(W, 'o-');

Gain = sum(W);

for i=1:num_wins
    
    Sai(i,1) = (1/Gain)*sum(W.*Va_seg(i,1:29).*exp(-1i*k*2*pi*f0*Ts));
    Sbi(i,1) = (1/Gain)*sum(Vb_seg(i,1:29).*exp(-1i*k*2*pi*f0*Ts));
    Sci(i,1) = (1/Gain)*sum(Vc_seg(i,1:29).*exp(-1i*k*2*pi*f0*Ts));
    
end
% for n=1:Nf
%     
%     Sa(n) = 1/G * W * sum( exp(-1i*k*2*pi*f0*Ts ));
% end
phiai = rad2deg(angle(Sai)+0.41855);
Aai = abs(Sai)*2;
% phiai = rad2deg(angle(Sai)+0.445);
% Aai = abs(Sai)*1.054;

erro_amp = (ref_seg(2,2:end) - Aai(1:22,1)')*100;

figure, hold on, grid on
plot(erro_amp)

figure, hold on, grid on
stem(Aai), stem(ref_seg(2,2:end)), title('abs ieee model')



figure, hold on, grid on
stem(phiai), title('angle ieee model')

%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE COMPONENTES SIMÉTRICAS



% COM DADOS FFT
[A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa(:,bin)', Ab(:,bin)', Ac(:,bin)', phia(:,bin)', phib(:,bin)', phic(:,bin)');
[A0r, A1r, A2r, phi0r, phi1r, phi2r] = comp_simetricas(Aar', Abr', Acr', phiar', phibr', phicr');

% COM DADOS DE REFERÊNCIA (para teste do estimador)
%[A0r, A1r, A2r, phi0r, phi1r, phi2r] = comp_simetricas(ref_seg(2,:), ref_seg(3,:), ref_seg(4,:),ref_seg(5,:), ref_seg(6,:), ref_seg(7,:));




% Plot fasores da componente fundamental
figure

subplot(2,3,1), hold on, grid on
plot(A0r, 's-')
plot(A0, 'o-')
title('Amplitudes (pu)')
ylabel('pu'), xlabel('Ciclo')
subplot(2,3,4), hold on, grid on
plot(rad2deg(phi0r), 's-')
plot(rad2deg(phi0), 'o-')

title('Fase (º)')
ylabel('º'), xlabel('Ciclo')

subplot(2,3,2), hold on, grid on
plot(A1r, 's-')
plot(A1, 'o-')
title('Amplitudes (pu)')
ylabel('pu'), xlabel('Ciclo')
subplot(2,3,5), hold on, grid on
plot(rad2deg(phi1r), 's-')
plot(rad2deg(phi1), 'o-')
title('Fase (º)')
ylabel('º'), xlabel('Ciclo')

subplot(2,3,3), hold on, grid on
plot(A2r, 's-')
plot(A2, 'o-')
title('Amplitudes (pu)')
ylabel('pu'), xlabel('Ciclo')
subplot(2,3,6), hold on, grid on
plot(rad2deg(phi2r), 's-')
plot(rad2deg(phi2), 'o-')
title('Fase (º)')
ylabel('º'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva

% COM DADOS DE REFERÊNCIA (teste do estimador)
[freq_final_ref] = estimador_freq_deriv_ang(phi0r, phi1r, phi2r, f0, Fs, ref_seg, m, N);

% COM DADOS FFT
[freq_final_fft] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

figure, hold on, grid on
plot(ref_seg(1,1:end), 'o-')
plot(freq_final_ref, 'o-')
plot(freq_final_fft, 'o-')
title('Frequência (Hz)'), legend('Referência', 'Estimativa (ref)', 'Estimativa (DFT)');
ylabel('Hz'), xlabel('Ciclo')


%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL PRONY

% com freq de ref:
[Apa, Apb, Apc, phipa, phipb, phipc] = estimador_red_prony(Va_seg, Vb_seg, Vc_seg, ref_seg(1,:), Fs, ref_seg);
% com freq estimada (com angulos de ref):
[Apa1, Apb1, Apc1, phipa1, phipb1, phipc1] = estimador_red_prony(Va_seg, Vb_seg, Vc_seg, freq_final_ref, Fs, ref_seg);
% com freq estimada (com angulos da dft):
[Apa2, Apb2, Apc2, phipa2, phipb2, phipc2] = estimador_red_prony(Va_seg, Vb_seg, Vc_seg, freq_final_fft, Fs, ref_seg);

% % Plot componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(Apa, 'o-')
% plot(Apa1, 'o-')
% plot(Apa2, 'o-')
% title('Amplitude (pu) - Prony')
% ylabel('pu'), xlabel('Ciclo')
% 
% subplot(2,1,2), hold on, grid on
% plot(phipa, 'o-')
% plot(phipa1, 'o-')
% plot(phipa2, 'o-')
% title('Fase (rad) - Prony')
% ylabel('rad'), xlabel('Ciclo')

%------------------------------------------------------------------------------
% ESTIMADOR DE AMORTECIMENTO

% COM DADOS DE REFERÊNCIA (teste do estimador)
[amort_final_ref] = estimador_amortec_deriv_amp(ref_seg(2,:), ref_seg(3,:), ref_seg(4,:), f0, Fs, ref_seg, m);

% COM DADOS FFT
[amort_final_fft] = estimador_amortec_deriv_amp(A0, A1, A2, f0, Fs, ref_seg, m, N);


% figure, hold on, grid on
% plot(amort_final_ref, 'o-')
% plot(amort_final_fft, 'o-')
% title('Amortecimento');% legend('Ref', 'DFT');
% ylabel('pu'), xlabel('Ciclo')
% 
% 
% figure, hold on, grid on
% plot((1-exp(amort_final*1/Fs))*100, 'o-')
% title('exp(amort_final*1/Fs)');% legend('Ref', 'DFT');
% ylabel('pu'), xlabel('Ciclo')

% 
% z_ref = exp(amort_final_ref*1/f0);
% z_fft = exp(amort_final_fft*1/f0);
% 
% figure
% subplot(2,1,1), hold on, grid on
% plot(z_ref, 'o-')
% plot(z_fft, 'o-')
% title('exp(amortecimento*1/Fs)'); legend('Ref', 'DFT');
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,1,2), hold on, grid on
% plot(z_ref - z_fft, 'o-')
% title('Erro'); 
% ylabel('pu'), xlabel('Ciclo')

% ----------------------------------------------------------------------------
% ESTIMADOR RM-PRONY (com amortecimento)

% COM DADOS DE REFERÊNCIA (teste do estimador)
[Apa3, Apb3, Apc3, phipa3, phipb3, phipc3] = estimador_mod_red_prony(Va_seg, Vb_seg, Vc_seg, ref_seg(1,:), Fs, ref_seg, amort_final_ref);

% COM DADOS FFT
[Apa4, Apb4, Apc4, phipa4, phipb4, phipc4] = estimador_mod_red_prony(Va_seg, Vb_seg, Vc_seg, freq_final_fft, Fs, ref_seg, amort_final_fft);

% % Plot componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(Apa2, 'o-')
% plot(Apb2, 'o-')
% plot(Apc2, 'o-')
% title('Amplitude (pu) - Prony')
% ylabel('pu'), xlabel('Ciclo')
% 
% subplot(2,1,2), hold on, grid on
% plot(phipa2, 'o-')
% plot(phipb2, 'o-')
% plot(phipc2, 'o-')
% title('Fase (rad) - Prony')
% ylabel('rad'), xlabel('Ciclo')



%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE ERROS

% Seleção automática da componente fundamental
[num_wins,win_size] = size(Aa);
% gran = Fs/win_size/f0;       % granularidade
% bin = round(1/gran)+1;       % local da componente fundamental no espectro

% define indices (recorte dos dados)

i=1;  % inicio
f=floor(T/(1/f0))-1; % fim (numero de ciclos)

% amplitudes
amp_med(1,i:f) = Aa(i:f,bin);    % dft
amp_med(2,i:f) = Apa(i:f,1)';  % com freq de ref
amp_med(3,i:f) = Apa1(i:f,1)'; % com freq estim (com angulos ref)
amp_med(4,i:f) = Apa2(i:f,1)'; % com freq estim (com angulos dft)
amp_med(5,i:f) = Apa3(i:f,1)'; % com freq ref e amort (com freq ref e amp de ref)
amp_med(6,i:f) = Apa4(i:f,1)'; % com freq ref e amort (com ang e amp de dft)
% fases
fase_med(1,i:f) = phia(i:f,bin);   
fase_med(2,i:f) = phipa(i:f,1)';
fase_med(3,i:f) = phipa1(i:f,1)';
fase_med(4,i:f) = phipa2(i:f,1)';
fase_med(5,i:f) = phipa3(i:f,1)';
fase_med(6,i:f) = phipa4(i:f,1)';
% frequência
%freq_med(1,:) = ref_seg(1,:);
freq_med(1,i:f) = freq_final_ref(1,i:f);
freq_med(2,i:f) = freq_final_fft(1,i:f);


[tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, ref_seg(:,i:f));




% % Plot erro frequência
% figure
% subplot(2,1,1),hold on, grid on  
% plot(ref_seg(1,:), 'o-')
% plot(freq_final, 'o-')
% plot(freq_final_prony, 'o-')
% plot(freq_final_ref, '*-')
% title('Frequência (Hz)'), legend('Ref', 'DFT', 'Prony');
% ylabel('Hz'), xlabel('Ciclo')
% subplot(2,1,2),hold on, grid on  
% plot(freq_error(1,:), 'o-')
% plot(freq_error(2,:), 'o-')
% plot(freq_error(3,:), 'o-')
% title('Erro de Frequência (Hz)'), legend('DFT', 'Prony');
% ylabel('Hz'), xlabel('Ciclo')

% Plot fasores da componente fundamental
figure
subplot(2,1,1), hold on, grid on
%plot(ref_seg(2,:),'o-')
plot(t_seg(i:f,1), amp_med(1,:), '*-')
plot(t_seg(i:f,1), amp_med(2,:), 'o-')
%plot(amp_med(3,:), 's-')
plot(t_seg(i:f,1), amp_med(4,:), 'd-')
plot(t_seg(i:f,1), amp_med(5,:), 'x-')
plot(t_seg(i:f,1), amp_med(6,:), '+-')
%title('Amplitudes (pu)')
%legend('DFT','R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
ylabel('Amplitude (pu)'), xlabel('tempo (s)')
subplot(2,1,2), hold on, grid on
% plot(rad2deg(ref_seg(5,:)), 'o-')
% plot(rad2deg(fase_med(1,:)), '*-')
% plot(rad2deg(fase_med(2,:)), 'o-')
% plot(rad2deg(fase_med(3,:)), 's-')
% title('Fase (º)'), legend('Ref','DFT','Prony','Prony2')
%plot((ref_seg(5,:)), 'o-')
plot(t_seg(i:f,1),(fase_med(1,:)), '*-')
plot(t_seg(i:f,1),(fase_med(2,:)), 'o-')
%plot((fase_med(3,:)), 's-')
plot(t_seg(i:f,1),(fase_med(4,:)), 'd-')
plot(t_seg(i:f,1),(fase_med(5,:)), 'x-')
plot(t_seg(i:f,1),(fase_med(6,:)), '+-')
%title('Fase (rad)')
%legend('DFT','R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
ylabel('Fase (rad)'), xlabel('tempo (s)')


% get(gca,'fontname')  % shows you what you are using.
% set(gca,'fontname','times')  % Set it to times

% Plot ERRO de fasores (DFT e Prony)
figure  
subplot(3,1,1), hold on, grid on%
plot(t_seg(i:f,1), tve(1,:), '*-')
plot(t_seg(i:f,1), tve(2,:), 'o-')
%plot(t_seg(i:f,1), tve(3,:), 's-')
plot(t_seg(i:f,1), tve(4,:), 'd-') 
plot(t_seg(i:f,1), tve(5,:), 'x-') % 
plot(t_seg(i:f,1), tve(6,:), '+-') % 
ylabel('TVE (%)'), xlabel('tempo (s)')
%title('TVE (%)') 
%legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
%title('TVE (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
subplot(3,1,2), hold on, grid on
plot(t_seg(i:f,1),amp_error(1,:), '*-')
plot(t_seg(i:f,1),amp_error(2,:), 'o-')
%plot(amp_error(3,:), 's-')
plot(t_seg(i:f,1),amp_error(4,:), 'd-')
plot(t_seg(i:f,1),amp_error(5,:), 'x-')
plot(t_seg(i:f,1),amp_error(6,:), '+-')
%title('Erro de amplitude (%)')
%legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
%title('Erro de amplitude (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
ylabel('Erro de amplitude (%)'), xlabel('tempo (s)')
subplot(3,1,3), hold on, grid on
plot(t_seg(i:f,1),phase_error(1,:), '*-')
plot(t_seg(i:f,1),phase_error(2,:), 'o-')
%plot(phase_error(3,:), 's-')
plot(t_seg(i:f,1),phase_error(4,:), 'd-')
plot(t_seg(i:f,1),phase_error(5,:), 'x-')
plot(t_seg(i:f,1),phase_error(6,:), '+-')
%title('Erro de fase (º)')
legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
%title('Erro de fase (º)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
ylabel('Erro de fase (º)'), xlabel('tempo (s)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:6
   max_tve(i) = max(tve(i,:));
   max_amp_error(i) = max(amp_error(i,:)); 
   max_phi_error(i) = max(phase_error(i,:)); 
end
%