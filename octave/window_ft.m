## Usage: wft = window_ft(signal, w, type, pad)
##
## Computes the windowed Fourier transform of the signal,
## with window length 2*w+1 and given window type.
## By default the signal is assumed non-periodic and padded with zeros,
## setting the last argument to false will create a cyclic DFT instead.
## 
## The valid window types are:
##  'gaussian', 'rect', 'hamming', 'hann', 'blackman'

function sgm = window_ft(sig, w, name = "gaussian", pad = true)
  type = validatestring(name, {"gaussian", "rectangle", "hamming", "hanning", "blackman"});
  win  = make_window(w, type);
  n    = length(sig);
  sgm  = zeros(n, n);

  % This is now the transform defined in the book, of the same length as the signal.
  % We always take DFT of the same interval, with the window highlighting different
  % sections of it.
  if pad == true
    sig_work = [zeros(1,n), sig, zeros(1,n)];
  else
    sig_work = [sig, sig, sig];
  endif
  
  % Another (faster?) way to do this would be a sliding DFT, try that later.
  for t = 1:n
    windowed = zeros(1, 3*n);
    windowed(n+t-w:n+t+w) = win .* sig_work(n+t-w:n+t+w);
    sgm(:, t) = fft(windowed(n+1:2*n));
  endfor
endfunction




% These are from the WTSP book
function win = make_window(w, type)
  t   = linspace(-1/2, 1/2, 2*w + 1);
  switch (type)
    case "gaussian"
      win = exp(-18*t.*t);
    case "hanning"
      win = (cos(pi*t)).^2;
    case "hamming"
      win = 0.54 + 0.46*cos(2*pi*t);
    case "rectangle"
      win = ones(1, 2*w + 1);
    case "blackman"
      win = 0.42 + 0.5*cos(2*pi*t) + 0.08*cos(4*pi*t);
  endswitch

  % Normalize to one
  norm = sqrt(sumsq(win));
  win = win / norm;
endfunction