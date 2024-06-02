% Date: 25/04/2024
% Experiment 1:
% Q) Generate OFDM Waveform, with no: of subcarriers = 4096. 
% Perform  equalisation with perfect CSI and imperfect CSI and compare the BER performance.
close all;clearvars;clc;


nfft = 4096; 
nsym = 14;
L = 6;    % number of channel taps 
ncp = L-1;  % length of cyclic prefix
SNRdB = -5:10;
SNR  = 10.^(SNRdB/10);

err_perf = zeros(length(SNRdB),1);
err_imperf = zeros(length(SNRdB),1);
for k=1:length(SNRdB)

    % Generate I and Q bits for QPSK
    BitsI = randi([0,1],nfft,nsym);
    BitsQ = randi([0,1],nfft,nsym);

    % Generate the QPSK transmit symbols
    InSym = (1/sqrt(2))*((1-(2*BitsI)) + 1i*(1-(2*BitsQ)));

    % Selecting a column in Time frequency grid
    for i=1:nsym
        % Input symbols to the IFFT block
        InIFFT = InSym(:,i);
        % Output of the IFFT block
        OutFFT = sqrt(nfft)*ifft(InIFFT,nfft);
        % P to S conversion
        PtoSout = OutFFT.';
        % Add cyclic prefix
        AddCP = [PtoSout(nfft-ncp+1:nfft),PtoSout];

        %% Channel 
        h_Ltaps = 1/sqrt(2)*(randn(1,L) + 1i*randn(1,L));
        RxSamCP = conv(h_Ltaps,AddCP);

        %% Receiver
        ChNoise = sqrt(1/2)*(randn(1,L+nfft+ncp-1) + 1i*randn(1,L+nfft+ncp-1));
        RxSamCP = sqrt(SNR(k)/2)*RxSamCP + ChNoise;

        % Remove cyclic prefix
        RemoveCP = RxSamCP(ncp+1:ncp+nfft);

        % OFDM reception
        RxSym = sqrt(1/nfft)*fft(RemoveCP,nfft);

        % Generate channel frequency response( We perform perfect CSI
        % equalisation using this)

        hFFT_perfect = fft(h_Ltaps,nfft);
        pilot_freq   = 4;
        h_imperfect  = calc_imperfectCSI(RxSym,InIFFT,pilot_freq);

 
        % Equalization(ZF)
        DecSymZF_perf = RxSym./hFFT_perfect;
        DecSymZF_imperf = RxSym./h_imperfect;

        % Decoded I and Q bits
        DecBitsI_ZF = (real(DecSymZF_perf)<0).';
        DecBitsQ_ZF = (imag(DecSymZF_perf)<0).';
        DecBitsI_imperf = (real(DecSymZF_imperf)<0).';
        DecBitsQ_imperf = (imag(DecSymZF_imperf)<0).';

        % Calculate BER
        err_perf(k)   = err_perf(k)   + sum(DecBitsI_ZF ~= BitsI(:,i)) + sum(DecBitsQ_ZF ~= BitsQ(:,i));
        err_imperf(k) = err_imperf(k) + sum(DecBitsI_imperf ~= BitsI(:,i)) + sum(DecBitsQ_imperf ~= BitsQ(:,i));

    end
end
BER_perf = err_perf/(numel(BitsQ)+numel(BitsI));
BER_imperf = err_imperf/(numel(BitsQ)+numel(BitsI));
semilogy(SNRdB, BER_perf,'bo-','linewidth',2.0,'MarkerFaceColor','b','MarkerSize',7.5);
hold on;
semilogy(SNRdB, BER_imperf,'ko-','linewidth',2.0,'MarkerFaceColor','k','MarkerSize',7.5);
grid on;
title('BER for OFDM');
legend('OFDM BER(perfectCSI)','OFDM BER(imperfectCSI)');
xlabel('SNR (dB)');
ylabel('BER');



function  h_imperfect  = calc_imperfectCSI(RxSym,InIFFT,pilot_freq)
    
    pilot_idx = 1:pilot_freq:length(RxSym);
    RxSym = RxSym(:);
    InIFFT = InIFFT(:);
    pilots = InIFFT(pilot_idx);
    h_imperfect = RxSym(pilot_idx)./pilots;
    h_imperfect = (repelem(h_imperfect,pilot_freq)).';

end

