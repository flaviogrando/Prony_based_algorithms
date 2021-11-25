%%% ESTIMADOR DE FREQUÊNCIA - Derivada do ângulo da sequência positiva

% ENTRADA COM DADOS DE ÂNGULO

% Sf_x - Espectro de frequência (usado apenas a comp. fundamental)
% f0 - frequência fundamental (teórica)
% Fs - Frequência de amostragem


function [freq_final] = estimador_freq_deriv_ang_fft(Sa, Sb, Sc, f0, Fs, ref_seg, m)

[num_wins,win_size] = size(Sa);

desvio_fase = zeros(3,num_wins-1);
freq_measured = ones(3,num_wins-1);
freq_final = ones(1,num_wins)*ref_seg(1,1);% ***** cria vetor c/ dados de freq nominal *******

% Transformação ABC -> 012
alfa = exp(2*pi*1j/3); % operador alfa
abc012 = (1/3)*[1 1 1; 1 alfa alfa^2; 1 alfa^2 alfa]; % Matriz de transf.

% Seleção automática da componente fundamental
gran = Fs/win_size/f0;    % granularidade
bin = round(1/gran)+1;           % local da componente fundamental no espectro
dt = m/Fs;                % passo da janela (s)
% seleciona fasor da componente fundamental
fasor_abc(1,:) = Sa(:,bin);
fasor_abc(2,:) = Sb(:,bin);
fasor_abc(3,:) = Sc(:,bin);

fasor_012 = abc012*fasor_abc; % Calcula componentes simétricas

% Cria ângulo constante de referência - para janelamento 'm' < N 
% (desvio angular dif de zero mesmo com freq nominal)
Ang_ref = (m/win_size)*2*pi;


for i=1:num_wins-1;
    %desvio_fase(i) = phia_seg(i+1,bin)-phia_seg(i,bin);
    desvio_fase(1,i) = (angle(fasor_012(1,i+1)) - angle(fasor_012(1,i))) - Ang_ref;
    desvio_fase(2,i) = (angle(fasor_012(2,i+1)) - angle(fasor_012(2,i))) - Ang_ref;
    desvio_fase(3,i) = (angle(fasor_012(3,i+1)) - angle(fasor_012(3,i))) - Ang_ref;

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