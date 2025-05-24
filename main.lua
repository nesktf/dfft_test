local bit = require("bit")
local cmplx = require("cmplx")

local function copy_array(arr)
  local out = {}
  for _, v in ipairs(arr) do
    table.insert(out, cmplx(v.real, v.imag))
  end
  return out
end

local function sample_signal(signal, max, min, n_samples, use_cmplx)
  use_cmplx = use_cmplx or false
  max = max or 1
  n_samples = n_samples or 64

  -- Samples in the range [min, max)
  local dt = (max-min)/n_samples
  local samples = {}
  for i = min, n_samples-1, 1 do
    table.insert(samples, cmplx(signal(i*dt)))
  end
  return samples
end

local function dft_naive(samples, inv)
  inv = inv or false
  local n = #samples
  local omega = (inv and 1 or -1)*2*math.pi/n
  local scale = (inv and 1/n or 1)

  local copy = copy_array(samples)
  for k = 0, n-1 do
    local sum = cmplx(0, 0)
    for i = 0, n-1 do
      sum = sum + copy[i+1]*cmplx.expi(omega*k*i)
    end
    samples[k+1] = sum*scale
  end
end

local function dft_fast_ct(samples, inv)
  inv = inv or false
  local n = #samples -- Assumes that n is a power of 2
  if (n == 1) then
    return
  end

  local half = n/2
  local evens = {}
  local odds = {}
  for i = 0, half-1 do
    table.insert(odds, samples[2*i + 1])
    table.insert(evens, samples[2*i + 2])
  end
  dft_fast_ct(evens, inv)
  dft_fast_ct(odds, inv)

  local ang = (inv and -1 or 1)*2*math.pi/n
  local scale = (inv and .5 or 1)
  for k = 0, half-1 do
    local w = cmplx.expi(ang*k)
    samples[k+1] = (evens[k+1] + w*odds[k+1])*scale
    samples[k+half+1] = (evens[k+1] - w*odds[k+1])*scale
  end
end

local function dft_fast_ct_inplace(samples, method, inv)
  -- https://cp-algorithms.com/algebra/fft.html
  inv = inv or false
  method = method or true

  local n = #samples -- Assumes that n is a power of 2
  if (method) then
    local function bit_reverse(num, exp)
      local res = 0
      for i = 0, exp-1 do
        if (bit.band(num, bit.lshift(1, i)) > 0) then
          res = bit.bor(res, bit.lshift(1, (exp-1-i)))
        end
      end
      return res
    end

    local lg_n = 0
    while (bit.lshift(1, lg_n) < n) do
      lg_n = lg_n + 1
    end

    for i = 0, n-1 do
      if (i < bit_reverse(i, lg_n)) then
        local idx_a = i+1
        local idx_b = bit_reverse(i, lg_n)+1

        local a = samples[idx_a]
        samples[idx_a] = samples[idx_b]
        samples[idx_b] = a
      end
    end
  else
    local i = 1
    local j = 0
    while(i < n) do
      local curr_bit = bit.rshift(n, 1)
      while (bit.band(j, curr_bit) > 0) do
        j = bit.bxor(j, curr_bit)
        curr_bit = bit.rshift(curr_bit, 1)
      end
      j = bit.bxor(j, curr_bit)

      if (i < j) then
        local a = samples[i+1]
        samples[i+1] = samples[j+1]
        samples[j+1] = a
      end

      i = i + 1
    end
  end

  local len = 2
  while (len <= n) do
    local ang = cmplx.expi((inv and 1 or -1)*2*math.pi/len)
    for i = 1, n, len do
      local w = cmplx(1)
      for j = 0, len/2-1 do
        local u = samples[i+j]
        local v = samples[i+j+len/2]*w
        samples[i+j] = u+v
        samples[i+j+len/2] = u-v
        w = w * ang
      end
    end
    len = len * 2
  end

  if (inv) then
    local scale = 1/n
    for i,_ in ipairs(samples) do
      samples[i] = samples[i]*scale
    end
  end
end

local function print_samples(samples)
  for i, v in ipairs(samples) do
    print(string.format("- x[%d] = (%.2f, %.2f) [%.2f]", i-1, v.real, v.imag, v:length()))
  end
end

do
  local max = 2*math.pi
  local min = 0
  local n_samples = 32
  print(string.format("Samples for sin(t) in [%.2f, %.2f) (%d samples)", min, max, n_samples))
  local samples = sample_signal(math.sin, 2*math.pi, 0, 32)
  print_samples(samples)

  print("\nNaive DFT")
  local naive_dft = copy_array(samples)
  dft_naive(naive_dft)
  print_samples(naive_dft)

  print("\nFast DFT recursive")
  local fast_dft_recursive = copy_array(samples)
  dft_fast_ct(fast_dft_recursive)
  print_samples(fast_dft_recursive)

  print("\nFast DFT inplace")
  local fast_dft_inplace = copy_array(samples)
  dft_fast_ct_inplace(fast_dft_inplace)
  print_samples(fast_dft_inplace)
end
