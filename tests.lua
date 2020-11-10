--
-- Basic tests
--
nx = tonumber(args[2])
NX_IN = tonumber(args[3])
NY_IN = tonumber(args[4])
nstep = tonumber(args[5])
-- vskip = math.floor(nx/200)
vskip = 1

pond = {
  init = function(x,y) return 1, 0, 0 end,
  out = "pond.out",
  nx = nx,
  NX_IN = NX_IN,
  NY_IN = NY_IN,
  nstep = nstep,
  vskip = vskip
}

river = {
  init = function(x,y) return 1, 1, 0 end,
  out = "river.out",
  nx = nx,
  NX_IN = NX_IN,
  NY_IN = NY_IN,
  nstep = nstep,
  vskip = vskip
}

dam = {
  init = function(x,y)
    if (x-1)*(x-1) + (y-1)*(y-1) < 0.25 then
      return 1.5, 0, 0
    else
      return 1, 0, 0
    end
  end,
  out = "dam_break.out",
  nx = nx,
  NX_IN = NX_IN,
  NY_IN = NY_IN,
  nstep = nstep,
  vskip = vskip
}

wave = {
  init = function(x,y)
    return 1.0 + 0.2 * math.sin(math.pi * x), 1, 0
  end,
  out = "wave.out",
  frames = 100,
  nx = nx,
  NX_IN = NX_IN,
  NY_IN = NY_IN,
  nstep = nstep,
  vskip = vskip
}

simulate(_G[args[1]])
