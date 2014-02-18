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
  N   = 2*w + 1;
  inc = floor(w/4);
  s   = 2**(steps/12);
  % Pad on the left with zeros.
  % For simplicity, just add enough zeros in the end
  work_signal = [zeros(w, 1); signal(:); zeros(2*w, 1)];

  % Hamming window, normalized to one
  window = hamming(N);
  window = window / sumsq(window);

  nframes = fix(length(signal)/inc);
  data = zeros(N, nframes);
  for idx = 0:(nframes - 1)
    data(:, idx + 1) = work_signal(1 + idx*inc:1 + idx*inc + 2*w) .* window;
  endfor
  
  sgm = fft(data);

  power = abs(sgm);
  phase = arg(sgm);

  % Add a column of zeros before taking the differences.
  phase = [zeros(size(phase, 1), 1), phase];
  % Phase difference between successive frames
  dphase = diff(phase, 1, 2);

  % Subtract the "natural" phase difference due to bin center frequency (not necessary?)
  dphase -= 2*pi*inc*(0:N-1)';
  % Normalize the phase and restore the natural component
  dphase = fmod(dphase + pi, 2*pi) - pi;
  dphase += 2*pi*inc*(0:N-1)';
  % During time ts the phase should then increase by w*ts:
  dphase = cumsum(s * dphase, 2);
  sgm = power .* exp(i*dphase);

  % Inverse windowed DFT
  data = ifft(sgm) .* window;
  % Sum all contributions to get back to the original frame.

%  output = reshape(output(1:length(signal)), size(signal));

endfunction