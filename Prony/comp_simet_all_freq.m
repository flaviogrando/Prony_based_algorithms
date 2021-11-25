% CÁLCULO DE COMPONENTES SIMÉTRICAS
% Seleção automática da componente fundamental


function [S0, S1, S2] = comp_simet_all_freq(Sa, Sb, Sc)


[num_wins,win_size] = size(Sa);

% gran = Fs/win_size/f0;       % granularidade
% bin = round(1/gran)+1;       % local da componente fundamental no espectro


% Transformação ABC -> 012
alfa = exp(2*pi*1j/3); % operador alfa
abc012 = (1/3)*[1 1 1; 1 alfa alfa^2; 1 alfa^2 alfa]; % Matriz de transf.


for l=1:num_wins;
    for c=1:win_size;
        fasor_012(:,c) = abc012*[Sa(l,c); Sb(l,c); Sc(l,c)];
        S0(l,c) = fasor_012(1,c);
        S1(l,c) = fasor_012(2,c);
        S2(l,c) = fasor_012(3,c);
    end  
end

%fasor_012 = abc012*fasor_abc;    % Calcula componentes simétricas
%fasor_012 = abc012*[Sa; Sb; Sc];

%
end