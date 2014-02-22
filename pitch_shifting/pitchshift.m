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
## and the time difference between frames is ta = inc / Fs. As a first step
## we want to scale the length of the frame (and interval between frames) so
## that the interval becomes ts = s*ta
## When using window (and DFT) length N, the true angular frequency of bin 'k'
## is 2*pi*k*Fs (for k < N/2, the rest should really be interpreted as negative frequencies).


function output = pitchshift(signal, steps, N)
  warning ("off", "Octave:broadcast");

  s   = 2**(steps/12);
  hopa = floor(N/4);
  hops = s*hopa;
  % Pad with zeros on both sides
  pad = fix(N/2);
  work_signal = [zeros(pad, 1); signal(:); zeros(N, 1)];

  % Hamming window. Normalization is nontrivial, as we hop frames.
  window = hamming(N);
  window = window / sqrt(sumsq(window(1:hopa:end)));

  nframes = fix((length(work_signal) - N)/hopa);
  data = zeros(N, nframes);
  for k = 0:(nframes - 1)
    data(:, k + 1) = work_signal(1 + k*hopa:N + k*hopa) .* window;
  endfor
  
  sgm = fft(data);

  power = abs(sgm);
  phase = arg(sgm);

  % Add a column of zeros before taking the differences.
  phase = [zeros(size(phase, 1), 1), phase];
  % Phase difference between successive frames
  dphase = diff(phase, 1, 2);
  center_freq = mod((0:N-1)' + fix(N/2) , N) - fix(N/2);
  center_freq = center_freq / N;

  % Subtract the "natural" phase difference due to bin center frequency
  dphase -= 2*pi*hopa*center_freq;
  % Normalize the phase and restore the natural component
  dphase = fmod(dphase + pi, 2*pi) - pi;
  dphase += 2*pi*hopa*center_freq;
  % During time ts the phase should then increase by w*ts:
  phase = cumsum(s * dphase, 2);
  sgm = power .* exp(i*phase);

  % Inverse windowed DFT
  data = real(ifft(sgm)) .* window;

  % Synthesize the time dilated signal from separate frames.
  % Required length is  nframes*hops + N
  work_signal = zeros(fix(s*length(signal)) + 2*N, 1);
  start_time = 1;
  for k = 1:nframes
    idx = fix(start_time);
    work_signal(idx:idx + N - 1) += data(:, k);
    start_time += hops;
  endfor
  
  % Relevant data starts from work_signal(w+1)
  output = work_signal;

%  output = reshape(output(1:length(signal)), size(signal));

endfunction