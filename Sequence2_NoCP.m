close all
clear all
N = 64;                                 % number of sub-carriers
%pos = [0:63];                          % case 1
pos = [5:15];                          %case2
Nu = length(pos);                       % number of useful sub-carriers
N_symb_OFDM = 1000;                     % number of symbols per sub-carrier
N_symb = Nu * N_symb_OFDM;              % total number of symbols = Nu*Nombre symboles OFDM



M = 4; % QPSK
EbN0dB = [100]; % Eb/No range


    
     % Channel impulse response
    hc = [1 0 0.5 0.25 zeros(1, N-4)];   %case1 
    Hc = fft(hc);
    trans = conj(Hc.');
    vc = ones(1, N_symb_OFDM);
    B = trans * vc;
    
    % Equalizer coefficients
    equalizer_ML = B;
    
    % Modulation
    bits = 2 * randi(2, 1, 2 * N_symb) - 3;      % +1/-1
    symb = bits(1:2:end) + 1i * bits(2:2:end);   % mapping bits => symboles
    entree_ifft = zeros(N, N_symb_OFDM); 
    
    for ii = 1:Nu 
        entree_ifft(pos(ii) + 1, :) = symb(ii:Nu:end);    % Correct indexing for subcarriers 5-15
    end
    
     % OFDM signal, matrix form, no CP
    sortie_ifft = ifft(entree_ifft); 
    signal_ofdmnoCP = reshape(sortie_ifft, 1, N* N_symb_OFDM);
    

    % Plot frequency response
    figure;
    [hhfreq, w] = freqz(hc);
    plot(w/(2*pi)*N, 20*log10(abs(hhfreq)/max(abs(hhfreq))));
    grid on;
    xlabel('f/Rs');
    title('Frequency Response of the Channel');

    
    % Signal PSD
    [pxx, f] = pwelch(signal_ofdmnoCP(:), 1024, 512, 1024, N);
    figure
    plot(f, 10 * log10(pxx / max(abs(pxx))))
    grid on
    xlabel ('f/Rs')
    ylabel('dB')
    title (strcat('OFDM PSD , Nu= : ',num2str(Nu),' carriers'))
    legend(strcat('positions = ', num2str(pos)))
    
    % Apply channel
    signal_recu = filter(hc, 1, signal_ofdmnoCP(:));

    % Plot PSD at receiver input
    figure;
    [pxx, f] = pwelch(signal_recu, 1024, 512, 1024, N);
    plot(f, 10*log10(pxx/max(abs(pxx))));
    grid on;
    xlabel('f/Rs');
    ylabel('dB');
    title('OFDM PSD at Receiver Input');
        
    % Demodulation
    mat_recu = reshape(signal_recu, N, N_symb_OFDM);    % vector => matrix
    mat_demod = fft(mat_recu);                                    % FFT (demodulation) 
    mat_demod_eq = mat_demod.*B;
   
    % Useful subcarriers recovery
    mat_demod_utile = mat_demod_eq(pos+1,:);
    mat_demod_utilenoeq = mat_demod(pos+1,:);
    % Apply equalization
    matrice_apres_egalisation = zeros(length(pos), N_symb_OFDM);
    for ii = 1:length(pos)
        matrice_apres_egalisation(ii, :) = mat_demod_utile(ii, :);% * equalizer_ML(pos(ii));
    end

    figure
    hold on
    plot(real(mat_demod_utilenoeq(ii, :)), imag(mat_demod_utilenoeq(ii, :)), '*')
    grid on
    title('Constellation without CP no equalization')    % no Equalization plot without CP
    xlabel('Real part')
    ylabel('Imaginary part')
    
    figure
    hold on
    plot(real(matrice_apres_egalisation(ii, :)), imag(matrice_apres_egalisation(ii, :)), '*')
    grid on
    title('Constellation without CP with equalization')    %Equalization plot without CP
    xlabel('Real part')
    ylabel('Imaginary part')
 

 % Plot constellations
   % ii = 15;              %subcarrier 15 for scatter diagram
    ii = 11;               %change to for case 2
    
    figure
    hold on
    plot(real(mat_demod_utilenoeq(ii, :)), imag(mat_demod_utilenoeq(ii, :)), '*')
    grid on
    title('Constellation for subcarrier 15 without equalizer, no noise')
    xlabel('Real part')
    ylabel('Imaginary part')

    figure
    hold on
    plot(real(matrice_apres_egalisation(ii, :)), imag(matrice_apres_egalisation(ii, :)), '*')
    grid on
    title('constellation of subcarrier 15 without CP with equalization, no noise')
    xlabel('real part')
    ylabel('imaginary part')

