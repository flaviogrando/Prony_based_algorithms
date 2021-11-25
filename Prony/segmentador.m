%%%% SEGMENTAÇÃO DOS SINAIS (JANELAMENTO) %%%%%

% N - tamanho da janela (numero de amostras)
% m - deslocamento da janela (numero de amostras)
% Va, Vb, Vc - sinais trifásicos
% t - vetor de tempo
% refs - vetor de referências

function [Va_seg, Vb_seg, Vc_seg, t_seg, ref_seg] = segmentador(Va, Vb, Vc, t, refs, N, m)

if m == 0
    m=1;
end

[num_canais,dim] = size(refs);              % tamanho total  
num_wins = fix(dim/m);      % numero de janelas

% correção para o final do janelamento
limite = fix(N/m);

% Va_seg = zeros(num_wins - limite ,N);
% Vb_seg = zeros(num_wins - limite ,N);
% Vc_seg = zeros(num_wins - limite ,N);
% t_seg = zeros(num_wins - limite ,N);

ref_seg = zeros(num_canais,num_wins - limite);

%%%%%%%
%w = flattopwin(N);
% figure
% stem(Va_seg(1,:).*w')
%%%%

%%% Segmentação do sinal (janelamento)
for i=1:num_wins - limite;
    inicio = (i-1)*m+1;               % indice inicial da janela
    fim = inicio + N-1;               % indice final da janela
    Va_seg(i,:) = Va(inicio:fim);     % segmentação(uma janela/linha)
    Vb_seg(i,:) = Vb(inicio:fim);     
    Vc_seg(i,:) = Vc(inicio:fim);
    t_seg(i,:) = t(inicio:fim);      % vetor de tempos
    
    % segmentação das referências (relativo ao instante inicial da janela)
    for j=1:num_canais;
        ref_seg(j,i) = refs(j,inicio);
    end
end
    

end