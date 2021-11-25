%%% TESTES ESTIMADOR DFT-PRONY
% teste de numero de amostras
clc
clear all
close all


set(0, 'defaultAxesFontSize',12);
set(0, 'defaultAxesFontName','times');

N_amostras = [8 16 32 64 128 256 512 1024];




 for w=1:length(N_amostras)



% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = N_amostras(w);      % número de amostras
m = N_amostras(w);      % passo de janelamento
Fs = 50*N_amostras(w);   % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 40;    % nível de ruído em dB ( 0 = sem ruído)
T = 10;        % tempo total (em segundos)
param = 1;    % freq. de modulação, por exemplo
type = 0;     % seleciona tipo de teste (ver gerador_sinais.m)

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
% plot(t,Va)
% plot(t,Vb)
% plot(t,Vc)
% legend('Va','Vb','Vc'), title('Sinais')
% xlabel('tempo (s)'), ylabel('amplitude')

%-------------------------------------------------------------------------------------------------------------------
% SEGMENTAÇÃO (JANELAMENTO)
[Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m);

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

% % Plot componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(Aa(:,2), 'o-')
% plot(Ab(:,2), 'o-')
% plot(Ac(:,2), 'o-')
% title('Amplitude (pu) - DFT')
% ylabel('pu'), xlabel('Ciclo')
% 
% subplot(2,1,2), hold on, grid on
% plot(phia(:,2), 'o-')
% plot(phib(:,2), 'o-')
% plot(phic(:,2), 'o-')
% title('Fase (rad) - DFT')
% ylabel('rad'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% CÁLCULO DE COMPONENTES SIMÉTRICAS

% COM DADOS FFT
[A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa, Ab, Ac, phia, phib, phic);

% COM DADOS DE REFERÊNCIA (para teste do estimador)
[A0r, A1r, A2r, phi0r, phi1r, phi2r] = comp_simetricas(ref_seg(2,:), ref_seg(3,:), ref_seg(4,:),ref_seg(5,:), ref_seg(6,:), ref_seg(7,:));




% % Plot fasores da componente fundamental
% figure
% 
% subplot(2,3,1), hold on, grid on
% plot(A0r, 's-')
% plot(A0, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,4), hold on, grid on
% plot(rad2deg(phi0r), 's-')
% plot(rad2deg(phi0), 'o-')
% 
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')
% 
% subplot(2,3,2), hold on, grid on
% plot(A1r, 's-')
% plot(A1, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,5), hold on, grid on
% plot(rad2deg(phi1r), 's-')
% plot(rad2deg(phi1), 'o-')
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')
% 
% subplot(2,3,3), hold on, grid on
% plot(A2r, 's-')
% plot(A2, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,6), hold on, grid on
% plot(rad2deg(phi2r), 's-')
% plot(rad2deg(phi2), 'o-')
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva

% COM DADOS DE REFERÊNCIA (teste do estimador)
[freq_final_ref] = estimador_freq_deriv_ang(phi0r, phi1r, phi2r, f0, Fs, ref_seg, m, N);

% COM DADOS FFT
[freq_final_fft] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N);

% figure, hold on, grid on
% plot(ref_seg(1,1:end), 'o-')
% plot(freq_final_ref, 'o-')
% plot(freq_final_fft, 'o-')
% title('Frequência (Hz)'), legend('Referência', 'Estimativa (ref)', 'Estimativa (DFT)');
% ylabel('Hz'), xlabel('Ciclo')


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

% COM DADOS DE REFERÊNCIA (amplitude de ref)
[Apa3, Apb3, Apc3, phipa3, phipb3, phipc3] = estimador_mod_red_prony(Va_seg, Vb_seg, Vc_seg, ref_seg(1,:), Fs, ref_seg, amort_final_ref);

% COM DADOS FFT
[Apa4, Apb4, Apc4, phipa4, phipb4, phipc4] = estimador_mod_red_prony(Va_seg, Vb_seg, Vc_seg, freq_final_fft, Fs, ref_seg, amort_final_fft);


% COM DADOS de REFERÊNCIA
%[Apa5, Apb5, Apc5, phipa5, phipb5, phipc5] = estimador_mod_red_prony(Va_seg, Vb_seg, Vc_seg, ref_seg(1,:), Fs, ref_seg, amort_ref);

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
gran = Fs/win_size/f0;       % granularidade
bin = round(1/gran)+1;       % local da componente fundamental no espectro

% define indices (recorte dos dados)

i=1;  % inicio
f=floor(T/(1/f0))-1; % fim (numero de ciclos)

% amplitudes
amp_med(1,i:f) = Aa(1,i:f);    % dft
amp_med(2,i:f) = Apa(i:f,1)';  % com freq de ref
amp_med(3,i:f) = Apa1(i:f,1)'; % com freq estim (com angulos ref)
amp_med(4,i:f) = Apa2(i:f,1)'; % com freq estim (com angulos dft)
amp_med(5,i:f) = Apa3(i:f,1)'; % com freq ref e amort (com amp de ref)
amp_med(6,i:f) = Apa4(i:f,1)'; % com freq estim e amort estim (com ang e amp de dft)
%amp_med(7,i:f) = Apa5(i:f,1)'; % com freq ref e amort ref
% fases
fase_med(1,i:f) = phia(1,i:f);   
fase_med(2,i:f) = phipa(i:f,1)';
fase_med(3,i:f) = phipa1(i:f,1)';
fase_med(4,i:f) = phipa2(i:f,1)';
fase_med(5,i:f) = phipa3(i:f,1)';
fase_med(6,i:f) = phipa4(i:f,1)';
%fase_med(7,i:f) = phipa5(i:f,1)';
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

% % Plot fasores da componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(ref_seg(2,:),'o-')
% plot(amp_med(1,:), '*-')
% plot(amp_med(2,:), 'o-')
% %plot(amp_med(3,:), 's-')
% plot(amp_med(4,:), 's-')
% plot(amp_med(5,:), 'd-')
% title('Amplitudes (pu)'),legend('Ref','DFT','R-Prony f_{ref}', 'R-Prony (ref)', 'R-Prony (dft)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,1,2), hold on, grid on
% % plot(rad2deg(ref_seg(5,:)), 'o-')
% % plot(rad2deg(fase_med(1,:)), '*-')
% % plot(rad2deg(fase_med(2,:)), 'o-')
% % plot(rad2deg(fase_med(3,:)), 's-')
% % title('Fase (º)'), legend('Ref','DFT','Prony','Prony2')
% plot((ref_seg(5,:)), 'o-')
% plot((fase_med(1,:)), '*-')
% plot((fase_med(2,:)), 'o-')
% %plot((fase_med(3,:)), 's-')
% plot((fase_med(4,:)), 's-')
% plot((fase_med(5,:)), 'd-')
% title('Fase (rad)'), legend('Ref','DFT','R-Prony f_{ref}', 'R-Prony (ref)', 'R-Prony (dft)')
% ylabel('rad'), xlabel('Ciclo')


% 
% % Plot ERRO de fasores (DFT e Prony)
% figure  
% subplot(3,1,1), hold on, grid on%
% plot(t_seg(i:f,1), tve(1,:), '*-')
% plot(t_seg(i:f,1), tve(2,:), 'o-')
% %plot(t_seg(i:f,1), tve(3,:), 's-')
% plot(t_seg(i:f,1), tve(4,:), 'd-') 
% plot(t_seg(i:f,1), tve(5,:), 'x-') % 
% plot(t_seg(i:f,1), tve(6,:), '+-') % 
% ylabel('TVE (%)'), xlabel('tempo (s)')
% %title('TVE (%)') 
% %legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
% %title('TVE (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
% subplot(3,1,2), hold on, grid on
% plot(t_seg(i:f,1),amp_error(1,:), '*-')
% plot(t_seg(i:f,1),amp_error(2,:), 'o-')
% %plot(amp_error(3,:), 's-')
% plot(t_seg(i:f,1),amp_error(4,:), 'd-')
% plot(t_seg(i:f,1),amp_error(5,:), 'x-')
% plot(t_seg(i:f,1),amp_error(6,:), '+-')
% %title('Erro de amplitude (%)')
% %legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
% %title('Erro de amplitude (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
% ylabel('Erro de amplitude (%)'), xlabel('tempo (s)')
% subplot(3,1,3), hold on, grid on
% plot(t_seg(i:f,1),phase_error(1,:), '*-')
% plot(t_seg(i:f,1),phase_error(2,:), 'o-')
% %plot(phase_error(3,:), 's-')
% plot(t_seg(i:f,1),phase_error(4,:), 'd-')
% plot(t_seg(i:f,1),phase_error(5,:), 'x-')
% plot(t_seg(i:f,1),phase_error(6,:), '+-')
% %title('Erro de fase (º)')
% legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
% %title('Erro de fase (º)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
% ylabel('Erro de fase (º)'), xlabel('tempo (s)')

% ERROS PARA DFT
for i=1:length(tve(1,:)) 
    tve_dft(:,w) = tve(1,:);         % 1 = DFT
    tve_r_prony_ref(:,w) = tve(2,:); % 2 = Prony (com freq ref)
    tve_r_prony(:,w) = tve(4,:);     % 4 = R_prony (com freq estim)
    tve_mr_prony_ref(:,w) = tve(5,:);% 5 = MR_prony (com dados ref)
    tve_mr_prony(:,w) = tve(6,:);    % 6 = MR_prony (com freq e alfa estim)
end



 end

% % 
% i=1;  % inicio
% f=length(ruido(1,:)); % fim (numero de ciclos)
% 
% % Plot ERRO de fasores (DFT e Prony)
% figure  
% subplot(3,1,1), hold on, grid on%
% plot(ruido(i:f), ruido_tve(1,:), '*-')
% plot(ruido(i:f), ruido_tve(2,:), 'o-')
% %plot(amort(i:f), tve(3,:), 's-')
% plot(ruido(i:f), ruido_tve(4,:), 'd-') 
% plot(ruido(i:f), ruido_tve(5,:), 'x-') % 
% plot(ruido(i:f), ruido_tve(6,:), '+-') % 
% %plot(ruido(i:f), ruido_tve(7,:), 's-') % 
% ylabel('TVE (%)'), xlabel('fator de amortecimento (\alpha)')
% %title('TVE (%)') 
% %legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
% %title('TVE (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
% subplot(3,1,2), hold on, grid on
% plot(ruido(i:f), ruido_amp_error(1,:), '*-')
% plot(ruido(i:f), ruido_amp_error(2,:), 'o-')
% %plot(amp_error(3,:), 's-')
% plot(ruido(i:f), ruido_amp_error(4,:), 'd-')
% plot(ruido(i:f), ruido_amp_error(5,:), 'x-')
% plot(ruido(i:f), ruido_amp_error(6,:), '+-')
% %plot(ruido(i:f), ruido_amp_error(7,:), 's-')
% %title('Erro de amplitude (%)')
% %legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)')
% %title('Erro de amplitude (%)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
% ylabel('Erro de amplitude (%)'), xlabel('fator de amortecimento (\alpha)')
% subplot(3,1,3), hold on, grid on
% plot(ruido(i:f), ruido_phase_error(1,:), '*-')
% plot(ruido(i:f), ruido_phase_error(2,:), 'o-')
% %plot(amort(i:f), amort_phase_error(3,:), 's-')
% plot(ruido(i:f), ruido_phase_error(4,:), 'd-')
% plot(ruido(i:f), ruido_phase_error(5,:), 'x-')
% plot(ruido(i:f), ruido_phase_error(6,:), '+-')
% %plot(ruido(i:f), ruido_phase_error(7,:), 's-')
% %title('Erro de fase (º)')
% legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (A-ref)', 'MR-Prony (A-dft)', 'MR-Prony (\alpha-ref))')
% %title('Erro de fase (º)'), legend('R-Prony (ref)', 'R-Prony (dft)', 'MR-Prony (ref)', 'MR-Prony (dft)') 
% ylabel('Erro de fase (º)'), xlabel('fator de amortecimento (\alpha)')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:length(ruido)
%    max_tve(i) = max(ruido_tve(i,:));
%    max_amp_error(i) = max(ruido_amp_error(i,:)); 
%    max_phi_error(i) = max(ruido_phase_error(i,:)); 
% end

figure
subplot(3,1,1), hold on, grid on%
boxplot(tve_dft, N_amostras)
title('DFT'), xlabel('N'), ylabel('TVE (%)'), legend('DFT')
subplot(3,1,2), hold on, grid on%
%boxplot(tve_r_prony_ref, N_amostras)
boxplot(tve_r_prony, N_amostras)
title('R-Prony (dft)'), xlabel('N'), ylabel('TVE (%)'), legend('R-Prony')
subplot(3,1,3), hold on, grid on%
%boxplot(tve_mr_prony_ref, N_amostras)
boxplot(tve_mr_prony, N_amostras)
title('MR-Prony (dft)'), xlabel('N'), ylabel('TVE (%)'), legend('MR-Prony')


figure, hold on, grid on%
plot(max(tve_dft), '*-')
plot(max(tve_r_prony_ref), 'o-')
plot(max(tve_r_prony), 'd-')
plot(max(tve_mr_prony_ref),'x-')
plot(max(tve_mr_prony),'+-')
xlabel('N'), ylabel('TVE (%)'), legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)','RM-Prony (ref)', 'RM-Prony (dft)')

figure, hold on, grid on%
plot(mean(tve_dft), '*-')
plot(mean(tve_r_prony_ref), 'o-')
plot(mean(tve_r_prony), 'd-')
plot(mean(tve_mr_prony_ref),'x-')
plot(mean(tve_mr_prony),'+-')
xlabel('N'), ylabel('TVE (%)'), legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)','RM-Prony (ref)', 'RM-Prony (dft)')

figure, hold on, grid on%
plot(mode(tve_dft), '*-')
plot(mode(tve_r_prony_ref), 'o-')
plot(mode(tve_r_prony), 'd-')
plot(mode(tve_mr_prony_ref),'x-')
plot(mode(tve_mr_prony),'+-')
xlabel('N'), ylabel('TVE (%)'), legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)','RM-Prony (ref)', 'RM-Prony (dft)')

figure, hold on, grid on%
plot(median(tve_dft), '*-')
plot(median(tve_r_prony_ref), 'o-')
plot(median(tve_r_prony), 'd-')
plot(median(tve_mr_prony_ref),'x-')
plot(median(tve_mr_prony),'+-')
xlabel('N'), ylabel('TVE (%)'), legend('DFT', 'R-Prony (ref)', 'R-Prony (dft)','RM-Prony (ref)', 'RM-Prony (dft)')
%subplot(3,1,2), hold on, grid on%
%boxplot(tve_r_prony_ref, N_amostras)
%plot(mean(tve_r_prony), N_amostras)
%xlabel('N'), ylabel('TVE (%)'), legend('R-Prony')
%subplot(3,1,3), hold on, grid on%
%boxplot(tve_mr_prony_ref, N_amostras)
%plot(mean(tve_mr_prony), N_amostras)
%xlabel('N'), ylabel('TVE (%)'), legend('RM-Prony')

% Monta vetores para o boxplot logo abaixo
% comparativo entre cada método para um SNR específico
a = [tve_dft(:,1), tve_r_prony_ref(:,1), tve_r_prony(:,1), tve_mr_prony_ref(:,1), tve_mr_prony(:,1)];
b = [tve_dft(:,3), tve_r_prony_ref(:,3), tve_r_prony(:,3), tve_mr_prony_ref(:,3), tve_mr_prony(:,3)];
c = [tve_dft(:,5), tve_r_prony_ref(:,5), tve_r_prony(:,5), tve_mr_prony_ref(:,5), tve_mr_prony(:,5)];
d = [tve_dft(:,7), tve_r_prony_ref(:,7), tve_r_prony(:,7), tve_mr_prony_ref(:,7), tve_mr_prony(:,7)];


figure
subplot(2,2,1), hold on, grid on
boxplot(a, 'Labels',{'a','b', 'c','d', 'e'})
ylabel('TVE (%)')
title(['SNR = ', num2str(N_amostras(1))])
subplot(2,2,2), hold on, grid on
boxplot(b, 'Labels',{'a','b', 'c','d', 'e'})
ylabel('TVE (%)')
title(['SNR = ', num2str(N_amostras(3))])
subplot(2,2,3), hold on, grid on
boxplot(c, 'Labels',{'a','b', 'c','d', 'e'})
ylabel('TVE (%)')
title(['SNR = ', num2str(N_amostras(5))])
subplot(2,2,4), hold on, grid on
boxplot(d, 'Labels',{'a','b', 'c','d', 'e'})
ylabel('TVE (%)')
title(['SNR = ', num2str(N_amostras(7))])