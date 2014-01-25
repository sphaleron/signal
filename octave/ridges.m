## Usage: signal = ridges(sgm, method, cutoff)
##
## Computes the instantaneous frequencies from the ridges
## of the given spectrogram / scalogram. Method can be either
## 'abs' (locally maximum absolute value) or 'phase' (stationary
## point of the phase). Cutoff limits the number of possible maxima
## to given fraction of the highest value / smallest derivative,
## and nneig defines the number of points on both sides belonging
## to the same local neigborhood.


function rgm = ridges(sgm, method = "abs", nneig = 2, cutoff = 0.05)
  type = validatestring(method, {"absolute", "phase"});

  data  = abs(sgm);
  limit = cutoff * max(data(:));

  switch (type)
    case "absolute"
      rgm   = zeros(size(data));
      % Local maxima
      for k = -nneig:nneig
        rgm = max(rgm, shift(data, k));
      endfor
      rgm = (data > limit) .* (rgm == data);
    case "phase"
      % We still require a clearly nonzero argument
      idx   = data > limit;
      phase = arg(sgm);
      % TBC: this seems pretty hard, typically the phase varies so much.
      % Some kind of unwrapping and smoothing required first?
  endswitch

endfunction