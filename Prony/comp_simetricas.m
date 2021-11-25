% CÁLCULO DE COMPONENTES SIMÉTRICAS
% Seleção automática da componente fundamental


function [A0, A1, A2, phi0, phi1, phi2] = comp_simetricas(Aa, Ab, Ac, phia, phib, phic)


%[num_wins,win_size] = size(Aa);

% gran = Fs/win_size/f0;       % granularidade
% bin = round(1/gran)+1;       % local da componente fundamental no espectro


% Transformação ABC -> 012
alfa = exp(2*pi*1j/3); % operador alfa
abc012 = (1/3)*[1 1 1; 1 alfa alfa^2; 1 alfa^2 alfa]; % Matriz de transf.

% amplitudes (fundamental)
amp_med(1,:) = Aa(1,:);   
amp_med(2,:) = Ab(1,:);  
amp_med(3,:) = Ac(1,:);  
% fases (fundamental)
fase_med(1,:) = phia(1,:);    
fase_med(2,:) = phib(1,:);    
fase_med(3,:) = phic(1,:);    

for j=1:3 % Conversão polar para retangular   
    %[x,y] = pol2cart(ref_seg(j+4,:),ref_seg(j+1,:)); % dados de referência
    [x,y] = pol2cart(fase_med(j,:),amp_med(j,:));     % dados da FFT
    fasor_abc(j,:) = x+1i*y;
end

fasor_012 = abc012*fasor_abc;    % Calcula componentes simétricas


[phi0,A0] = cart2pol(real(fasor_012(1,:)),imag(fasor_012(1,:)));
[phi1,A1] = cart2pol(real(fasor_012(2,:)),imag(fasor_012(2,:)));
[phi2,A2] = cart2pol(real(fasor_012(3,:)),imag(fasor_012(3,:)));

end