## Usage: signal = chirps(N, eps)
##
## Creates a signal of length N consisting of two linear chirps,
## whose instantaneous frequencies vary between 1 and 100 
## and cross halfway through the signal.
##
## The signal components in continuous time t = 0..1 are:
## c1 = sin ( t + (99/2)*t^2 )
## c2 = sin ( 100*t - (99/2)*t^2 )
##
## The optional parameter eps gives the amplitude of additive white noise.

function signal = chirps(N, eps = 0)
  t  = linspace(0, 1, N);
  c1 = sin( t + 49.5*t.*t);
  c2 = sin( 100*t + 49.5*t.*t);
  signal = 1.1*c1 + 0.9*c2 + eps*randn(size(t));
endfunction

