function thresh(x, sorh, t)
  y = zeros(x)
  if sorh == "s"
    tmp = abs(x) - t
    tmp = (tmp + abs(tmp)) / 2
    y   = sign(x).*tmp
  end
  if sorh == "h"
    n = length(x)
    for i=1:n
      if abs(x[i]) <= t
        x[i] = 0.
      end
    end
    y = copy(x)
  end
  return y
end
