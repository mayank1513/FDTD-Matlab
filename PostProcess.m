% computing frequency domain outputs
NFFT=2^nextpow2(number_of_time_steps);
frequencies_ = 1/(2*dt)*linspace(0,1,NFFT/2+1);

[X] = fft(ProbValues,NFFT)/number_of_time_steps;    %time_to_frequency_domain(x, dt, frequency_array, time_shift);
fft_results = 2*abs(X(1:NFFT/2+1,:)).^2; % as only peaks are important we have sqared the signal

P = sum(fft_results,2);
if exist('MaxDetactablePeak','var')
else
    MaxDetactablePeak=1e-7;
end

MinimumPower=MaxDetactablePeak*max(P(2:end));
[pks, locs]=findpeaks(P(2:end));
Peaks=frequencies_(locs(P(locs)>MinimumPower)+1)*a/c; % normalized frequency

fprintf('Calculation is done for %g set of wave vector(k) out of %g \n rInd= %g out of %g  plasmaFreqInd= %g out of %g', kInd,length(kx),rInd,length(r_by_a),plasmaFreqInd,length(omega_p_Array))
figure(2)
plot(kInd-1,Peaks(1:min(length(Peaks),30)),'o','LineWidth',1)
Pks = Peaks(1:min(length(Peaks),30)); Pks(end+1:30)=10;
BandSq(kInd,:)= Pks;
hold on
axis([0 length(kx)-1 0 2])
drawnow;