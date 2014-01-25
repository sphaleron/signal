## Usage: wft = window_ft(signal, w, type, pad)
##
## Computes the windowed Fourier transform of the signal,
## with window length 2*w+1 and given window type.
## By default the signal is assumed non-periodic and padded with zeros,
## setting the last argument to false will generate a cyclic DFT instead.
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

  % Shift the window center to the beginning of the signal. If padding is enabled,
  % pad the window with zeros so that it does not leak to the other end.
  % grow to signal length (as we work with the whole signal length DFTs all the time)
  win = postpad(win, n);
  if pad == true
    win = [win, zeros(1, w)];
  endif
  win = shift(win, -(w+1));
  
  % We could also generate transforms of the window length via sliding DFT, try that later.
  for t = 1:n
    win = shift(win, 1);
    sgm(:, t) = fft(sig .* win(1:n));
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