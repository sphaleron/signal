## Usage: output = pitchshift(signal, step, N)
##
## Phase vocoder implementation, using octave signal processing package.
## Shifts the frequency content by 'step' semitones, by taking windowed
## Fourier transform of length 'N' (currently hamming window is used),
## and stepping the window N/4 points between iterations.
##
## Note: this handles now the whole signal at once. Real-time implementation
## should be quite different.


function output = pitchshift(signal, steps, N)
  inc = floor(N/4);
  s   = 2**(steps/12);
  [sgm, pars] = stft(signal, N, inc, N, 2);
  power = abs(sgm);
  % Unwrap the phase to get correct phase shifts
  phase = unwrap(arg(sgm), pi, 2);
  % imagesc(power)
  % imagesc(phase)
  dphase = phase - shift(phase, 1, 2);
  % imagesc(abs(dphase))
  
endfunction