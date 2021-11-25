%%%%%%% CALCULA ERROS %%%%%%%

% Vetores de medições, 1 linha por canal
% amp_med - amplitudes
% fase_med - fases em radianos
% freq_med - frequências

% ref_seg - Referências, obtidas no segmentador


function [tve, amp_error, phase_error, freq_error] = calcula_erros(amp_med, fase_med, freq_med, ref_seg)

% ERROS DE FREQUÊNCIA
[num_canais, num_med] = size(freq_med);     % identifica tamanho dos dados de freq 
freq_error = zeros(num_canais,num_med);  % cria vetor de erro de freq

for i=1:num_canais % calcula erro de freq (1ª linha do ref_seg)
    freq_error(i,:) = ref_seg(1,:) - freq_med(i,:);
end

% ERROS FASORIAIS
[num_canais, num_med] = size(amp_med);   % identifica tamanho dos dados fasoriais
amp_error = zeros(num_canais,num_med);   % cria vetor de erro de amplitude
phase_error = zeros(num_canais,num_med); % cria vetor de erro de fase
tve = zeros(num_canais,num_med);

% Converte para coord. cart. para calculo do TVE
[Xre,Xie] = pol2cart(fase_med,amp_med);        % Medições
[Xr,Xi] = pol2cart(ref_seg(5,:),ref_seg(2,:)); % Referências

for i=1:num_canais  % calcula erro a 
    
    amp_error(i,:) = ( (ref_seg(2,:) - amp_med(i,:) )./ref_seg(2,:))*100;
    
    phase_error(i,:) = rad2deg(ref_seg(5,:) - fase_med(i,:));
    
    % Corrige desvios de fase maior que 180º
    for j=1:length(phase_error(i,:))
        if phase_error(i,j)>180
            phase_error(i,j) = phase_error(i,j)-360;
        elseif phase_error(i,j)<-180
            phase_error(i,j) = phase_error(i,j)+360;
        end
    end

    
    for j=1:length(Xr) 
        %tve(i,j) = 100*sqrt(((Xr(1,j)-Xre(i,j))^2 + (Xi(1,j)-Xie(i,j))^2) /(Xr(1,j)^2 + Xi(1,j)^2));
        tve(i,j) = 100*((Xr(1,j)-Xre(i,j)) + (Xi(1,j)-Xie(i,j))) /(Xr(1,j) + Xi(1,j));
    end
end


end
