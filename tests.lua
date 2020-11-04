--
-- Basic tests
--
nx = tonumber(args[2])
block_nx = tonumber(args[3])
block_ny = tonumber(args[4])
vskip = math.floor(nx/200)

pond = {
  init = function(x,y) return 1, 0, 0 end,
  out = "pond.out",
  nx = nx,
  block_nx = block_nx,
  block_ny = block_ny,
  vskip = vskip
}

river = {
  init = function(x,y) return 1, 1, 0 end,
  out = "river.out",
  nx = nx,
  block_nx = block_nx,
  block_ny = block_ny,
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
  block_nx = block_nx,
  block_ny = block_ny,
  vskip = vskip
}

wave = {
  init = function(x,y)
    return 1.0 + 0.2 * math.sin(math.pi * x), 1, 0
  end,
  out = "wave.out",
  frames = 100,
  nx = nx,
  block_nx = block_nx,
  block_ny = block_ny,
  vskip = vskip
}

simulate(_G[args[1]])
