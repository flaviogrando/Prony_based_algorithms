%%%% GERADOR DE SINAIS

% -----------------   SAÍDAS  ---------------- %
% Va, Vb, Vc - Sinais senoidais a taxa de amostragem Fs e período total T.
% t - vetor de tempos (relativo a amostragem)
% refs - vetor de referências (freq, amplitudes e fases) - por amostra
% mod - modulação (utilizado para plotar em comparação com o sinal)

% -----------------   ENTRADAS  ---------------- %
% Exemplo:
% freqs = [50 150];  % Vetor de frequências
% phases = [0 0];    % vetor de fases
% amps = [1 0.3];    % vetor de amplitudes
%
% Fs = 12800;   % Taxa de amostragem
% T = 2;        % tempo total (em segundos)
% noise = 40;   % nível de ruído em dB
% type = 1;     % seleciona tipo de teste
% param = 2;    % parâmetro do teste (freq. de modulação, por exemplo)
%
% Tipos de teste
% 0 - Estado estacionário
% 1 - Mod em amplitude
% 2 - Mod em fase
% 3 - Rampa em freq
% 4 - amortecimento
% 5 - degrau em amplitude
% 6 - degrau em fase
% 7 - Extra: degrau em frequência (implementar)
% 8 - Extra: mod em fase e amplitude

function [Va, Vb, Vc, t, refs, mod] = gerador_sinais(f0, freqs, phases, amps,Fs, T,noise, type, param)

mod = 0;
Va = 0;
Vb = 0;
Vc = 0;

% Parâmetros default -----------------------------------------------
ka = 0;                 % amplitude da modulação em fase
kx = 0;                 % amplitude da modulação em amplitude
R = 0;                  % Rampa em frequência em Hz/s
alfa = 0;               % Fator de amortecimento
sa = 0;                 % Degrau em fase
sx = 0;                 % Degrau em amplitude
% ------------------------------------------------------------------

Ts = 1/Fs;              % Período de amostragem
w = 2*pi*freqs;         % freq. do sinal (rad/s)
t = 0:Ts:T-Ts;          % vetor de tempos

% Modulação
fm = param;               % frequência de modulação
wm = 2*pi*fm;             % freq. modulação (rad/s)

switch type
    case 0
        %disp('Teste: Estado estacionário')
    case 1
        kx = 0.1;
        mod = kx*cos(wm*t);
        disp('Teste: Modulação em amplitude')
    case 2
        ka = 0.1;
        mod = ka*cos(wm*t);
        disp('Teste: Modulação em fase')
    case 3
        R = param;
        disp('Teste: Rampa em frequência')
    case 4
        alfa = param;
        disp('Teste: Amortecimento')
    case 5
        sx1 = zeros(1,floor(length(t)/2));
        sx2 = amps(1,1)*param*ones(1,floor(length(t)/2));
        sx = [sx1 sx2];
        disp('Teste: degrau em amplitude')
    case 6
        sa1 = zeros(1,floor(length(t)/2));
        sa2 = param*ones(1,floor(length(t)/2));
        sa = [sa1 sa2];
        disp('Teste: degrau em fase')
    case 7
        disp('Teste: degrau em frequência (NÃO IMPLEMENTADO)')
    case 8
        kx = 0.1;
        ka = 0.1;
        mod = kx*cos(wm*t); % para plotar
        disp('Teste: modulação em fase e amplitude')
    otherwise
        disp('Tipo de teste inexistente')
end


for i=1:length(freqs)      % itera cada componente de freq
    
    % modulações
    mod_a = ka*cos(wm*t - pi) + pi*R*(t.^2);
    mod_b = ka*cos(wm*t - pi) + pi*R*(t.^2);
    mod_c = ka*cos(wm*t - pi) + pi*R*(t.^2);
    
    % argumento do sinal
    arg_a(i,:) = w(i)*t + phases(1,i) + mod_a + sa;
    arg_b(i,:) = w(i)*t + phases(2,i) + mod_b + sa;
    arg_c(i,:) = w(i)*t + phases(3,i) + mod_c + sa;

    % amplitudes do sinal
    Amp_a(i,:) = amps(1,i)*(exp(alfa*t)) + kx*cos(wm*t) + sx;
    Amp_b(i,:) = amps(2,i)*(exp(alfa*t)) + kx*cos(wm*t) + sx;
    Amp_c(i,:) = amps(3,i)*(exp(alfa*t)) + kx*cos(wm*t) + sx;

    % gera o sinal
    Sa = Amp_a(i,:).*cos(arg_a(i,:));                  % Tensão fase a
    Sb = Amp_b(i,:).*cos(arg_b(i,:));                  % Tensão fase b
    Sc = Amp_c(i,:).*cos(arg_c(i,:));                  % Tensão fase c
    
    % soma componentes de frequência
    Va = Va + Sa;
    Vb = Vb + Sb;
    Vc = Vc + Sc;
end

% insere ruído (ou não)
if noise ~= 0
    Va = awgn(Va,noise,'measured');
    Vb = awgn(Vb,noise,'measured');
    Vc = awgn(Vc,noise,'measured');
end

% Gera referências (da comp. fundamental)
f_ref = freqs(1) - ka*fm*sin(wm*t - pi) + R*t; 

amp_ref(1,:) = Amp_a(1,:);
amp_ref(2,:) = Amp_b(1,:);
amp_ref(3,:) = Amp_c(1,:);

phi_ref(1,:) = arg_a(1,:);
phi_ref(2,:) = arg_b(1,:);
phi_ref(3,:) = arg_c(1,:);


% Corrige fase maior que 180º
for j=1:3
    i=1;
    while i <= length(t)
        if phi_ref(j,i)>pi
            phi_ref(j,i) = phi_ref(j,i)-2*pi;
        else
            i=i+1;
        end  
    end
end




% % Calcula referências de fase para freq. fora do nominal
% for i=1:length(t);
%     %phi_ref(i+1) = phi_ref(i) + 2*pi*Ts*(f_ref(i+1)-f0);
%     % Corrige desvios de fase maior que 180º
%     if phi_ref(i)>pi
%         phi_ref(i) = phi_ref(i)-2*pi;
%     elseif phi_ref(i)<-pi
%         phi_ref(i) = phi_ref(i)+2*pi;
%     end
% end


refs = [f_ref; amp_ref; phi_ref;];

end


