%%% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo da sequência positiva

% ENTRADA COM DADOS DE ÂNGULO

% phix - ângulos 
% f0 - frequência nominal
% Fs - Frequência de amostragem


function [freq_final] = estimador_freq_deriv_ang(phi0, phi1, phi2, f0, Fs, ref_seg, m, N)

num_wins = length(phi0);

desvio_fase = zeros(3,num_wins-1);
freq_measured = ones(3,num_wins-1);
freq_final = ones(1,num_wins)*ref_seg(1,1);% ***** cria vetor c/ dados de freq de ref *******


dt = m/Fs;                % passo da janela (s)

% Cria ângulo constante de referência - para janelamento 'm' < N 
% (desvio angular dif de zero mesmo com freq nominal)
num_cic = f0*N/Fs;
Ang_ref = (m/N)*2*pi*num_cic;


for i=1:num_wins-1;
    
    desvio_fase(1,i) = phi0(i+1) - phi0(i) - Ang_ref;
    desvio_fase(2,i) = phi1(i+1) - phi1(i) - Ang_ref;
    desvio_fase(3,i) = phi2(i+1) - phi2(i) - Ang_ref;

    % Corrige desvio ângular maior que 180º
    for j=1:3
        % Enquanto houver desvios >ou< que +/-180º executa correção
        while desvio_fase(j,i)>pi || desvio_fase(j,i)<-pi
            if desvio_fase(j,i)>pi
                desvio_fase(j,i) = desvio_fase(j,i)-2*pi;
            elseif desvio_fase(j,i)<-pi
                desvio_fase(j,i) = desvio_fase(j,i)+2*pi;
            end
        end
        % Finalmente, calcula freq. com desvios corrigidos
        freq_measured(j,i) = f0 + desvio_fase(j,i)/(2*pi*dt); 
    end         
    freq_final(i+1) = freq_measured(2,i); % '2' seleciona freq. da seq. pos.
end

end