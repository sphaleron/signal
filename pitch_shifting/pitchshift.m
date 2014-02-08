## Usage: output = pitchshift(signal, step, N)
##
## Phase vocoder implementation, using octave signal processing package.
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
## is 2*pi*k*Fs (for k < N/2, the rest should be interpreted as negative frequncies,
## but that does not matter here).

function output = pitchshift(signal, steps, N)
  inc = floor(N/4);
  s   = 2**(steps/12);
#   signal = 
  [sgm, pars] = stft(signal, N, inc, N/2, 2);
  power = abs(sgm);
  % Wrapped phase
  phase = arg(sgm);
  % Add a column of zeros before taking the differences.
  phase = prepad(phase, size(phase, 2) + 1, 0, 2);
  % Phase difference between successive frames
  dphase = phase - shift(phase, 1, 2);
  dphase = dphase(:, 2:end);

  % Subtract the "natural" phase difference due to bin center frequency (not necessary?)
  dphase -= 2*pi*inc*(0:N-1)';
  % Normalize the phase and restore the natural component
  dphase = fmod(dphase + pi, 2*pi) - pi;
  dphase += 2*pi*inc*(0:N-1)';
  % During time ts the phase should then increase by w*ts:
  dphase = cumsum(s * dphase, 2);
  sgm = power .* exp(i*dphase);

  % Inverse windowed DFT
  output = synthesis(sgm, pars);


endfunction