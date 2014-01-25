## Usage: signal = chirps(N, eps)
##
## A collection of reference signals from the WTSP book.
## Adapted from Wavelab, cleaning up a bit and adding a possibility
## of additive noise (of size eps).


function signal = chirps(N, eps = 0)
  t  = linspace(0, 1, N);
  % linear chirp
  sig1 = cos( (10*pi*t).^2 );

  % quadratic chirp
  sig2 = cos( 30*(pi*(1-t)).^3 );

  % gaussian blobs
  ix = linspace(-20, 20, 2*N + 1);
  g  = exp(-4*ix.^2);
  sig3 = g(N/2 + 1: N/2 + N) .* cos(50*pi*t);
  sig4 = g(N/8 + 1: N/8 + N) .* cos(350*pi*t);

  signal = sig1 + sig2 + sig3 + sig4;
  if (eps > 1e-6)
    signal += eps*randn(size(t));
  endif

  % Smoothly vanish on the boundaries.
  % Compared w/ Wavelab: (1 + sin(-pi/2 + x))/2 = (1 - sin(pi/2 - x))/2 = (1 - cos(x))/2
  %                       = sin^2(x/2)
  t    = linspace(0, pi/2, N/8);
  ramp = sin(t).^2;
  env  = [ramp, ones(1, N - 2*length(t)), reverse(ramp)];
  signal = env .* signal;;
endfunction
