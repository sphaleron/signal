## Usage: output = pitchshift(signal, step, w)
##
## Phase vocoder implementation.
## Shifts the frequency content by 'step' semitones, by taking windowed
## Fourier transform of length '2w+1' (currently hamming window is used),
## and stepping the window w/2 points between iterations.
##
## Note: this handles now the whole signal at once. Real-time implementation
## should work differently.

## At sampling frequency Fs the time between successive samples is 1/Fs,
## and the time difference between frames is ta = inc / Fs. As a first step
## we want to scale the length of the frame (and interval between frames) so
## that the interval becomes ts = s*ta
## When using window (and DFT) length N, the true angular frequency of bin 'k'
## is 2*pi*k*Fs (for k < N/2, the rest should be interpreted as negative frequncies,
## but that does not matter here).


function output = pitchshift(signal, steps, w)
  warning ("off", "Octave:broadcast");

  N   = 2*w + 1;
  s   = 2**(steps/12);
  hopa = floor(w/4);
  hops = s*hopa;
  % Pad on the left with zeros.
  % For simplicity, just add enough zeros in the end
  work_signal = [zeros(w, 1); signal(:); zeros(N, 1)];

  % Hamming window. Normalization is nontrivial, as we hop frames.
  window = hamming(N);

  nframes = fix(length(signal)/hopa);
  data = zeros(N, nframes);
  for k = 0:(nframes - 1)
    data(:, k + 1) = work_signal(1 + k*hopa:1 + k*hopa + 2*w) .* window;
  endfor
  
  sgm = fft(data);

  power = abs(sgm);
  phase = arg(sgm);

  % Add a column of zeros before taking the differences.
  phase = [zeros(size(phase, 1), 1), phase];
  % Phase difference between successive frames
  dphase = diff(phase, 1, 2);

  % Subtract the "natural" phase difference due to bin center frequency (not necessary?)
  dphase -= 2*pi*hopa*(0:N-1)';
  % Normalize the phase and restore the natural component
  dphase = fmod(dphase + pi, 2*pi) - pi;
  dphase += 2*pi*hopa*(0:N-1)';
  % During time ts the phase should then increase by w*ts:
  phase = cumsum(s * dphase, 2);
  sgm = power .* exp(i*phase);

  % Inverse windowed DFT
  data = real(ifft(sgm)) .* window;

  % Synthesize the time dilated signal from separate frames.
  work_signal = zeros(fix(s*length(signal)) + N, 1);
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