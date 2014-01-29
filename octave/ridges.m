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

  % Perhaps we should use octave's findpeaks instead of own hacks?

  switch (type)
    case "absolute"
      rgm   = zeros(size(data));
      % Replace every point by its max neighbor
      for k = -nneig:nneig
        rgm = max(rgm, shift(data, k));
      endfor
      % Large enough and local maximum
      rgm = (data > limit) .* (rgm == data);
    case "phase"
      phase = arg(sgm);
      rgm   = ones(size(data));
      % Plotting the two-sided difference reveals a slow linear trend near the ridges,
      % which could then be pinpointed at the zero crossings of the linear lines.
      % Everywhere else the difference oscillates wildly.
      % (So much that there is not much point doing any unwrapping over the whole column.)
      dphase = abs(shift(phase, 1, 2) - shift(phase, -1, 2));
      % It is probably best if we first get rid of most of the oscillations
      % by only considering high amplitude regions:
      dphase(find(data < limit)) = nan;

      % Replace every point by its min neighbor
      for k = -nneig:nneig
        rgm = min(rgm, shift(dphase, k));
      endfor
      % Large enough and local maximum
      rgm = (data > limit) .* (rgm == dphase);
  endswitch

endfunction



function phase = uwphase(sgm)
  
endfunction