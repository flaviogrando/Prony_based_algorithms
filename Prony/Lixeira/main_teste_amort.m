%%% TESTES ESTIMADOR DFT-PRONY

clc
clear all
close all



amort = -3:0.1:3;


 for w=1:length(amort)

amort_vetor = ones(50)*amort(w);   % Vetor de referência para amortecimento

% ---------------------------------------------------------------------
% INICIALIZAÇÕES

N = 256;      % número de amostras
m = 256;      % passo de janelamento
Fs = 12800;   % Taxa de amostragem
f0 = 50;      % freq. fundamental (teórica)
noise = 0;   % nível de ruído em dB ( 0 = sem ruído)
T = 1;       % tempo total (em segundos)
param = amort(w);    % freq. de modulação, por exemplo
type = 4;     % seleciona tipo de teste (ver gerador_sinais.m)

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
[Va, Vb, Vc, t, refs] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param);

% figure, hold on, grid on
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
% xlabel('tempo (s)'), ylabel('amp (pu)')
% subplot(3,1,3), hold on, grid on
% plot(t_seg(:,1),ref_seg(5,:), 'o-')
% xlabel('tempo (s)'), ylabel('fase (rad)')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL DFT
[Sa, Sb, Sc, Aa, Ab, Ac, phia, phib, phic] = estimador_dft(Va_seg, Vb_seg, Vc_seg);

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
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva
% COM DADOS DFT

[freq_final] = estimador_freq_deriv_ang_fft(Sa, Sb, Sc, f0, Fs, ref_seg, m);

% figure, hold on, grid on
% plot(ref_seg(1,2:end), 'o-')
% plot(freq_final, 'o-')
% title('Frequência (Hz)'), legend('Ref', 'DFT');
% ylabel('Hz'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR FASORIAL PRONY

[Ha, Hb, Hc, Apa, Apb, Apc, phipa, phipb, phipc] = estimador_red_prony(Va_seg, Vb_seg, Vc_seg, freq_final, Fs, ref_seg);

% % Plot componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(Apa, 'o-')
% plot(Apb, 'o-')
% plot(Apc, 'o-')
% title('Amplitude (pu) - Prony')
% ylabel('pu'), xlabel('Ciclo')
% 
% subplot(2,1,2), hold on, grid on
% plot(phipa, 'o-')
% plot(phipb, 'o-')
% plot(phipc, 'o-')
% title('Fase (rad) - Prony')
% ylabel('rad'), xlabel('Ciclo')


%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva
% COM DADOS PRONY

[freq_final_prony] = estimador_freq_deriv_ang_prony(Ha, Hb, Hc, f0, Fs, ref_seg, m);

% figure, hold on, grid on
% plot(ref_seg(1,2:end), 'o-')
% plot(freq_final, 'o-')
% title('Frequência (Hz)'), legend('Ref', 'DFT');
% ylabel('Hz'), xlabel('Ciclo')

%--------------------------------------------------------------------------------------------------------------------
% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo do fasor de seq. positiva
% COM DADOS DE REFERÊNCIA (teste do estimador)

[freq_final_ref] = estimador_freq_deriv_ang_ref(f0, Fs, ref_seg, m, N);

% figure, hold on, grid on
% plot(ref_seg(1,2:end), 'o-')
% plot(freq_final, 'o-')
% title('Frequência (Hz)'), legend('Ref', 'DFT');
% ylabel('Hz'), xlabel('Ciclo')

%----------------------------------------------------------------------------
% CÁLCULO DE COMPONENTES SIMÉTRICAS
% Seleção automática da componente fundamental
[num_wins,win_size] = size(Aa);
gran = Fs/win_size/f0;       % granularidade
bin = round(1/gran)+1;       % local da componente fundamental no espectro


% Transformação ABC -> 012
alfa = exp(2*pi*1j/3); % operador alfa
abc012 = (1/3)*[1 1 1; 1 alfa alfa^2; 1 alfa^2 alfa]; % Matriz de transf.

% amplitudes
amp_med(1,:) = Aa(:,bin)';   
amp_med(2,:) = Ab(:,bin)';
amp_med(3,:) = Ac(:,bin)';
% fases
fase_med(1,:) = phia(:,bin)';  
fase_med(2,:) = phib(:,bin)';  
fase_med(3,:) = phic(:,bin)';  

for j=1:3 % Conversão polar para retangular   
    %[x,y] = pol2cart(ref_seg(j+4,:),ref_seg(j+1,:)); % dados de referência
    [x,y] = pol2cart(fase_med(j,:),amp_med(j,:));     % dados da FFT
    fasor_abc(j,:) = x+1i*y;
end

fasor_012 = abc012*fasor_abc;    % Calcula componentes simétricas


[phi0,A0] = cart2pol(real(fasor_012(1,:)),imag(fasor_012(1,:)));
[phi1,A1] = cart2pol(real(fasor_012(2,:)),imag(fasor_012(2,:)));
[phi2,A2] = cart2pol(real(fasor_012(3,:)),imag(fasor_012(3,:)));


% % Plot fasores da componente fundamental
% figure
% 
% subplot(2,3,1), hold on, grid on
% plot(A0, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,4), hold on, grid on
% plot(rad2deg(phi0), 'o-')
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')
% 
% subplot(2,3,2), hold on, grid on
% plot(A1, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,5), hold on, grid on
% plot(rad2deg(phi1), 'o-')
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')
% 
% subplot(2,3,3), hold on, grid on
% plot(A2, 'o-')
% title('Amplitudes (pu)')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,3,6), hold on, grid on
% plot(rad2deg(phi2), 'o-')
% title('Fase (º)')
% ylabel('º'), xlabel('Ciclo')

%------------------------------------------------------------------------------
% ESTIMADOR DE AMORTECIMENTO
% com dados de referência

[amort_final_ref] = estimador_amortec_deriv_amp_ref(Sa, Sb, Sc, f0, Fs, ref_seg, m);

%amort_final_ref = amort_final_ref;
% figure, hold on, grid on
% plot(amort_final_ref, 'o-')
% title('Amortecimento');% legend('Ref', 'DFT');
% ylabel('pu'), xlabel('Ciclo')
% 
% 
% figure, hold on, grid on
% plot((1-exp(amort_final*1/Fs))*100, 'o-')
% title('exp(amort_final*1/Fs)');% legend('Ref', 'DFT');
% ylabel('pu'), xlabel('Ciclo')

%------------------------------------------------------------------------------
% ESTIMADOR DE AMORTECIMENTO
% com dados da fft

[amort_final_fft] = estimador_amortec_deriv_amp_fft(Sa, Sb, Sc, f0, Fs, ref_seg, m);

% figure, hold on, grid on
% plot(amort_final_ref, 'o-')
% plot(amort_final_fft, 'o-')
% title('Amortecimento'); legend('Ref', 'DFT');
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
% ESTIMADOR PRONY COM AMORTECIMENTO

[Ha2, Hb2, Hc2, Apa2, Apb2, Apc2, phipa2, phipb2, phipc2] = estimador_red_prony_2(Va_seg, Vb_seg, Vc_seg, freq_final, Fs, ref_seg, amort_final_ref);

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

% % Seleção automática da componente fundamental
% [num_wins,win_size] = size(Aa);
% gran = Fs/win_size/f0;       % granularidade
% bin = round(1/gran)+1;       % local da componente fundamental no espectro

% seleciona dado da comp. fund. a partir do ciclo 2 (da fase X)
% amplitudes
amp_med(1,:) = Aa(:,bin)';   
amp_med(2,:) = Apa(:,1)';
amp_med(3,:) = Apa2(:,1)';
% fases
fase_med(1,:) = phia(:,bin)';  
fase_med(2,:) = phipa(:,1)';
fase_med(3,:) = phipa2(:,1)';
% frequência
freq_med(1,:) = freq_final;
freq_med(2,:) = freq_final_prony;
freq_med(3,:) = freq_final_ref;

[tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, ref_seg);




% % % Plot erro frequência
% % figure
% % subplot(2,1,1),hold on, grid on  
% % plot(ref_seg(1,:), 'o-')
% % plot(freq_final, 'o-')
% % plot(freq_final_prony, 'o-')
% % plot(freq_final_ref, '*-')
% % title('Frequência (Hz)'), legend('Ref', 'DFT', 'Prony');
% % ylabel('Hz'), xlabel('Ciclo')
% % subplot(2,1,2),hold on, grid on  
% % plot(freq_error(1,:), 'o-')
% % plot(freq_error(2,:), 'o-')
% % plot(freq_error(3,:), 'o-')
% % title('Erro de Frequência (Hz)'), legend('DFT', 'Prony');
% % ylabel('Hz'), xlabel('Ciclo')
% 
% % Plot fasores da componente fundamental
% figure
% subplot(2,1,1), hold on, grid on
% plot(ref_seg(3,:), 'o-')
% plot(amp_med(1,:), '*-')
% plot(amp_med(2,:), 'o-')
% plot(amp_med(3,:), 's-')
% title('Amplitudes (pu)'),legend('Ref','DFT','Prony','Prony2')
% ylabel('pu'), xlabel('Ciclo')
% subplot(2,1,2), hold on, grid on
% % plot(rad2deg(ref_seg(5,:)), 'o-')
% % plot(rad2deg(fase_med(1,:)), '*-')
% % plot(rad2deg(fase_med(2,:)), 'o-')
% % plot(rad2deg(fase_med(3,:)), 's-')
% % title('Fase (º)'), legend('Ref','DFT','Prony','Prony2')
% plot((ref_seg(6,:)), 'o-')
% plot((fase_med(1,:)), '*-')
% plot((fase_med(2,:)), 'o-')
% plot((fase_med(3,:)), 's-')
% title('Fase (rad)'), legend('Ref','DFT','Prony','Prony2')
% ylabel('rad'), xlabel('Ciclo')
% 
% % Plot ERRO de fasores (DFT e Prony)
% figure  
% subplot(3,1,1), hold on, grid on%
% plot(tve(1,:), '*-')
% plot(tve(2,:), 'o-')
% plot(tve(3,:), 's-')
% ylabel('%'), xlabel('Ciclo')
% title('TVE (%)'), legend('DFT','Prony','Prony2')
% subplot(3,1,2), hold on, grid on
% plot(amp_error(1,:), '*-')
% plot(amp_error(2,:), 'o-')
% plot(amp_error(3,:), 's-')
% title('Erro de amplitude (%)'), legend('DFT','Prony', 'Prony2')
% ylabel('%'), xlabel('Ciclo')
% subplot(3,1,3), hold on, grid on
% plot(phase_error(1,:), '*-')
% plot(phase_error(2,:), 'o-')
% plot(phase_error(3,:), 's-')
% title('Erro de fase (º)'), legend('DFT','Prony', 'Prony2') %'DFT', 
% ylabel('deg'), xlabel('Ciclo')



% 
for i=1:3
    amort_tve(i,w) = tve(i,1);
    amort_amp_error(i,w) = amp_error(i,1);
    amort_phase_error(i,w) = phase_error(i,1);
end

for i=1:3
    amort_amp_med(i,w) = amp_med(i,1);
    amort_fase_med(i,w) = fase_med(i,1);
end

amort_ref_seg(1,w) = ref_seg(2,1);
amort_ref_seg(2,w) = ref_seg(5,1);

end

% Plot fasores da componente fundamental
figure
subplot(2,1,1), hold on, grid on
plot(amort_ref_seg(1,:), 'o-')
plot(amort_amp_med(1,:), '*-')
plot(amort_amp_med(2,:), 'o-')
plot(amort_amp_med(3,:), 's-')
title('Amplitudes (pu)'),legend('Ref','DFT','Prony','Prony2')
ylabel('pu'), xlabel('Ciclo')
subplot(2,1,2), hold on, grid on
% plot(rad2deg(ref_seg(5,:)), 'o-') plot(rad2deg(fase_med(1,:)), '*-')
% plot(rad2deg(fase_med(2,:)), 'o-') plot(rad2deg(fase_med(3,:)), 's-')
% title('Fase (º)'), legend('Ref','DFT','Prony','Prony2')
plot(amort_ref_seg(2,:), 'o-')
plot(amort_fase_med(1,:), '*-')
plot(amort_fase_med(2,:), 'o-')
plot(amort_fase_med(3,:), 's-')
title('Fase (rad)'), legend('Ref','DFT','Prony','Prony2')
ylabel('rad'), xlabel('Ciclo')

% Plot ERRO de fasores (DFT e Prony)
figure  
subplot(3,1,1), hold on, grid on%
plot(amort_tve(1,:), '*-')
plot(amort_tve(2,:), 'o-')
plot(amort_tve(3,:), 's-')
ylabel('%'), xlabel('Ciclo')
title('TVE (%)'), legend('DFT','Prony','Prony2')
subplot(3,1,2), hold on, grid on
plot(amort_amp_error(1,:), '*-')
plot(amort_amp_error(2,:), 'o-')
plot(amort_amp_error(3,:), 's-')
title('Erro de amplitude (%)'), legend('DFT','Prony', 'Prony2')
ylabel('%'), xlabel('Ciclo')
subplot(3,1,3), hold on, grid on
plot(amort_phase_error(1,:), '*-')
plot(amort_phase_error(2,:), 'o-')
plot(amort_phase_error(3,:), 's-')
title('Erro de fase (º)'), legend('DFT','Prony', 'Prony2') %'DFT', 
ylabel('deg'), xlabel('Ciclo')