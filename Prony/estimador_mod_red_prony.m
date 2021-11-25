%%%% ESTIMADOR FASORIAL PRONY


function [Apa, Apb, Apc, phipa, phipb, phipc] = estimador_mod_red_prony(Va_seg, Vb_seg, Vc_seg, freq_final, Fs, ref_seg, amort_final)


alfa = amort_final;
Ts = 1/Fs;
T = Ts;%Ts;

[num_wins,win_size] = size(Va_seg);
% Declarando vetores
Apa = zeros(num_wins,2);
Apb = zeros(num_wins,2);
Apc = zeros(num_wins,2);
phipa = zeros(num_wins,2);
phipb = zeros(num_wins,2);
phipc = zeros(num_wins,2);

 z_mod= exp(2*alfa*T);

for i=1:num_wins  % inicia a partir da 2nd janela devida a med de 
                         % frequência (por desvio angular)
                         
    % freq_final(i);  % freq. estimada
    % ref_seg(1,i);   % freq de referência
    
    f = freq_final(i); %ref_seg(1,i); 
    
    order = length(f)*2;  % numero de comp. de freq. no sinal
    
     
    %--- step - reduced prony ---
    re_z = sqrt( z_mod(i)/(1+(tan(2*pi*Ts*f)).^2) );
    im_z = sqrt( (z_mod(i)*(tan(2*pi*Ts*f).^2)) / (1+(tan(2*pi*Ts*f)).^2) );
   %  re_z = sqrt( 1/(1+(tan(2*pi*Ts*f)).^2) );
   %  im_z = sqrt( 1*(tan(2*pi*Ts*f).^2)) / (1+(tan(2*pi*Ts*f)).^2) );

    z = re_z + im_z*1i;
    z = [z conj(z)];

    % %----------------step3--------------% 
    Z=zeros(win_size,order); 
    for k=1:win_size 
        for m=1:order 
            Z(k,m)=z(m)^(k-1); %Z(k,m)=rootz(m)^(k-1);
        end 
    end 

    Va = Va_seg(i,:);  % sinal de entrada
    Vb = Vb_seg(i,:);  % sinal de entrada
    Vc = Vc_seg(i,:);  % sinal de entrada
    
    % spectro - nas colunas 
    Ha(:,i)=pinv(Z)*Va';
    Hb(:,i)=pinv(Z)*Vb'; 
    Hc(:,i)=pinv(Z)*Vc'; 
    
    c(i) = 64*(1-z_mod(i));
    
    % seleciona amplitude da fundamental
    Apa(i,:) = ( abs(Ha(:,i))*2 )';% + c(i);
    Apb(i,:) = ( abs(Hb(:,i))*2 )';
    Apc(i,:) = ( abs(Hc(:,i))*2 )';

    % seleciona fase da fundamental.
    phipa(i,:) = (angle(Ha(:,i)))';%+ 10*(1-z_mod(i));
    phipb(i,:) = (angle(Hb(:,i)))'; 
    phipc(i,:) = (angle(Hc(:,i)))'; 
    
end