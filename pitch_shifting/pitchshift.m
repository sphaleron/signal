## Usage: output = pitchshift(signal, step, N)
##
## Phase vocoder implementation.
## Shifts the frequency content by 'step' semitones, by taking windowed
## Fourier transform of length 'N' (currently hamming window is used),
## and stepping the window N/4 points between iterations.
##
## Note: this handles now the whole signal at once. Real-time implementation
## should work differently.

## At sampling frequency Fs the time between successive samples is 1/Fs,
## and the time difference between frames is ta = hopa / Fs. As a first step
## we want to scale the length of the frame (and interval between frames) so
## that the interval becomes ts = s*ta
## When using window (and DFT) length N, the true angular frequency of bin 'k'
## is 2*pi*k*Fs (for k < N/2).


function output = pitchshift(signal, steps, N)
  warning ("off", "Octave:broadcast");

  s    = 2**(steps/12);
  % Some constants:
  H = fix(N/2);

  hopa = fix(N/4);
  hops = s*hopa;
  % Pad with half window length on both sides.
  pad = H;
  
  work_signal = [zeros(pad, 1); signal(:); zeros(pad, 1)];

  % Hann window. Normalization needs some thinking, seems to work quite well.
  window = hanning(N);
  % window = window / scale

  % N + (nf-1)*hopa <= len -> nf <= (len - N)/hopa + 1 
  nframes = fix((length(work_signal) - N)/hopa) + 1;
  data = zeros(N, nframes);
  for k = 0:(nframes - 1)
    data(:, k + 1) = work_signal(1 + k*hopa : N + k*hopa) .* window;
  endfor
  
  sgm = fft(data);

  power = abs(sgm);
  phase = arg(sgm);

  % Phase difference between successive frames
  % Add also the original phase in the beginning
  dphase = diff(phase, 1, 2);
  dphase = [phase(:, 1) dphase];

  % Frequency bin central frequency
  % We are only interested in the first H+1 bins here, rest will be fixed
  % by symmetry below. TODO: only compute for the interesting bins, not all
  center_freq = 0:N-1;
  center_freq = center_freq' / N;
  % Subtract the "natural" phase difference due to bin center frequency
  dphase -= 2*pi*hopa*center_freq;
  % Normalize the phase difference and restore the natural component
  dphase = mod(dphase + pi, 2*pi) - pi;
  dphase += 2*pi*hopa*center_freq;

  % Scale the phase differences by the time scaling s and sum:
  phase = cumsum(s*dphase, 2);
  sgm = power .* exp(i*phase);

  % Explicitly restore conjugate symmetry before inverse DFT:
  sgm(H+2:end, :) = flipud(conj(sgm(2:N-H, :)));

  % Inverse windowed DFT
  data = window .* real(ifft(sgm));

  % Synthesize the time dilated signal from separate frames.

  % index of the last frames last element: N + (nframes-1)*hops
  work_signal = zeros(round((nframes-1)*hops) + N, 1);
  start_time = 1.0;
  for k = 1:nframes
    idx = round(start_time);
    work_signal(idx:idx + N - 1) += data(:, k);
    start_time += hops;
  endfor
  % Relevant data starts from work_signal(s*pad + 1)
  work_signal = work_signal(round(s*pad) + 1:end);

  output = interp1(0:length(work_signal) - 1, work_signal, s*(0:length(signal) - 1), 'linear');
  output = reshape(output, size(signal));

endfunction
