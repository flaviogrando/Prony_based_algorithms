%%% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo da sequência positiva

% ENTRADA COM DADOS DE ÂNGULO

% Sf_x - Espectro de frequência (usado apenas a comp. fundamental)
% f0 - frequência fundamental (teórica)
% Fs - Frequência de amostragem


function [amort_final] = estimador_amortec_deriv_amp(A0, A1, A2, f0, Fs, ref_seg, m, N)

num_wins = length(A0);

desvio_amp = zeros(3,num_wins-1);
amort_measured = zeros(3,num_wins-1);
amort_final = zeros(1,num_wins);% ***** cria vetor c/ dados de freq nominal *******


dt = m/Fs;                % passo da janela (s)


% Cria ângulo constante de referência - para janelamento 'm' < N 
% (desvio angular dif de zero mesmo com freq nominal)
Ang_ref = (m/256)*2*pi;



for i=1:num_wins-1
    
%     dec_log = log( A1(i)/A1(i+1));
%     damp_fact(i) = 1/(sqrt(1 + ((2*pi)/dec_log)^2 ));
    
    desvio_amp(1,i) = A0(i+1) - A0(i);
    desvio_amp(2,i) = A1(i+1) - A1(i); 
    desvio_amp(3,i) = A2(i+1) - A2(i);


    for j=1:3
        amort_measured(j,i) = desvio_amp(j,i)/dt; 
    end  
    
    amort_final(i) = amort_measured(2,i); % '2' seleciona freq. da seq. pos.
end

end