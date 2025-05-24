local cmplx = {
  sqlength = function(self)
    return self.real*self.real + self.imag*self.imag
  end,
  length = function(self)
    return math.sqrt(self:sqlength())
  end,
  normalize = function(self)
    local len = self:length()
    self.real = self.real/len
    self.imag = self.imag/len
  end,
  print = function(self)
    print(self.real, self.imag)
  end,
}

local function new(real, imag)
  local out = setmetatable({}, cmplx._mt)
  out.real = real or 0
  out.imag = imag or 0
  return out
end

local function polar(r, theta)
  local out = { real = r*math.cos(theta), imag = r*math.sin(theta) }
  return setmetatable(out, cmplx._mt)
end


local function exp(c)
  local out = setmetatable({}, cmplx._mt)
  local s = math.exp(c.real)
  out.real = s*math.cos(c.imag)
  out.imag = s*math.sin(c.imag)
  return out
end

local function exps(real, imag)
  real = real or 0
  imag = imag or 0
  local out = setmetatable({}, cmplx._mt)
  local s = math.exp(real)
  out.real = s*math.cos(imag)
  out.imag = s*math.sin(imag)
  return out
end

local function expi(real)
  real = real or 0
  local out = setmetatable({}, cmplx._mt)
  out.real = math.cos(real)
  out.imag = math.sin(real)
  return out
end

local function cast_input(a, b)
  if (type(a) ~= "table" or getmetatable(a) ~= cmplx._mt) then
    a = new(tonumber(a))
  end
  if (type(b) ~= "table" or getmetatable(b) ~= cmplx._mt) then
    b = new(tonumber(b))
  end
  return a, b
end

cmplx._mt = {
  __index = cmplx,
  __add = function(a, b)
    a, b = cast_input(a, b)
    return new(a.real+b.real, a.imag+b.imag)
  end,
  __sub = function(a, b)
    a, b = cast_input(a, b)
    return new(a.real-b.real, a.imag-b.imag)
  end,
  __mul = function(a, b)
    a, b = cast_input(a, b)
    return new(a.real*b.real - a.imag*b.imag, a.real*b.imag + a.imag*b.real)
  end,
  __div = function(a, b)
    a, b = cast_input(a, b)
    local d = b:sqlength()
    return new((a.real*b.real + a.imag*b.imag)/d, (a.imag*b.real - a.real*b.imag)/d)
  end,
}

return setmetatable({
  _base = cmplx,
  new = new,
  polar = polar,
  exp = exp,
  exps = exps,
  expi = expi,
},{
  __call = function(_, real, imag) return new(real, imag) end,
})
