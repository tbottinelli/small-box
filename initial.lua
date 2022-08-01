#!/usr/bin/env halmd
--
-- Copyright © 2011-2014 Felix Höfling
-- Copyright © 2010-2012 Peter Colberg
--
-- This file is part of HALMD.
--
-- HALMD is free software: you can redistribute it and/or modify
-- it under the terms of the GNU Lesser General Public License as
-- published by the Free Software Foundation, either version 3 of
-- the License, or (at your option) any later version.
--
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU Lesser General Public License for more details.
--
-- You should have received a copy of the GNU Lesser General
-- Public License along with this program.  If not, see
-- <http://www.gnu.org/licenses/>.
--

-- grab modules
local log = halmd.io.log
local mdsim = halmd.mdsim
local numeric = halmd.numeric
local observables = halmd.observables
local random = halmd.random
local writers = halmd.io.writers
local readers = halmd.io.readers
local utility = halmd.utility
local random1 = math.random
math.randomseed(os.time())
-- Setup and run simulation
--
function main(args)
    -- total number of particles from sum of particles per species
    local nspecies = #args.particles
    local nfluid = numeric.sum(args.particles)
   -- derive edge lengths from number of particles, density and edge ratios
    local file = readers.h5md({path = args.input})
    local reader , sample = observables.phase_space.reader({file = file, location = {"particles", "pore"}, fields = {"position"}})

    reader:read_at_step(0)
    local npore = assert(sample.nparticle)
    print(npore)

    local edges = mdsim.box.reader({file = file, location = {"particles", "pore"}, fields = {"box"}})
    local box = mdsim.box({edges = edges})

    local eps = 0.01
    -- create simulation domain with periodic boundary conditions

    -- create system state
    local particle = {
        fluid = mdsim.particle({dimension = 3, particles = nfluid, species = 2, label = "fluid"})
      , pore = mdsim.particle({dimension = 3, particles = npore, species = 2, label = "pore"})
    }
    local all_group = {}
    for label, p in pairs(particle) do
    	all_group[label] = mdsim.particle_groups.all({particle = p, label = label})
    end
    
    local phase_space = {}
    for label, g in pairs(all_group) do
	phase_space[label] = observables.phase_space({box = box, group = g})
    end
    phase_space["pore"]:set(sample)

    -- set particle species
    local species = {}
    --for i = 1, nfluid do table.insert(species,0) end
    for i = 1, npore do table.insert(species, 1) end
    particle["pore"].data["species"] = species

    -- calculate how to fill the box to stay inside and not overcrowd
    margin = 1
    offset = (args.slab*edges[1][1])/2 + margin
    fill_box = (edges[1][1]/2 - offset)*2/edges[1][1]
    print(fill_box)
   
    -- set initial particle positions sequentially on an fcc lattice
    local lattice = mdsim.positions.lattice({box = box, particle = particle["fluid"], slab = {fill_box,1,1}})
    lattice:set()
    local positions = particle["fluid"].data["position"]

    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle["fluid"]
      , group = all_group["fluid"]
      , temperature = args.temperature
    })
    boltzmann:set()
    local velocities = particle["fluid"].data["velocity"]
    
    --oss lua tables start with index 1 so this is for position x
    --this moves the fluid away from the pore particles
    for i = 1,nfluid do
        if (positions[i][1] < 0) then 
            positions[i][1] = positions[i][1] - ((args.slab*edges[1][1]/2) + margin)
        end
        if (positions[i][1] >= 0) then
            positions[i][1] = positions[i][1] + ((args.slab*edges[1][1]/2) + margin)
        end
    end
    particle["fluid"].data["position"] = positions
    particle["fluid"].data["velocity"] = velocities

    local potential = mdsim.potentials.pair.lennard_jones({species = 2
      , epsilon = 1
      , sigma = 1
    })

    cut = math.pow(2, 1 /6) 
    -- smoothing at potential cutoff
    potential = potential:truncate({"smooth_r4",
       cutoff = {
	       {2.5, cut}
	     , {cut, cut}
     }
      , h = 0.005
    })

    local binning = {
	    fluid = mdsim.binning({
		    box = box
		  , particle = particle["fluid"]
		 -- , skin = 0.2
		  , r_cut = potential.r_cut
	  })
	  , pore = mdsim.binning({
		    box = box
		  , particle = particle["pore"]
		 -- , skin = 0.2
		  , r_cut = potential.r_cut
		  , occupancy = 0.5 * particle["pore"].nparticle / particle["fluid"].nparticle
	  })
  }

    -- compute forces
    for label, p2 in pairs(particle) do
	   local neighbour = mdsim.neighbour({
		   box = box
		 , particle = {particle["fluid"], p2}
		 , r_cut = potential.r_cut
		 , binning = {binning["fluid"], binning[label]}
	 })
	   mdsim.forces.pair({
		   box = box
		 , particle = { particle["fluid"], p2 }
		 , potential = potential
		 , neighbour = neighbour
	   })
    end

    -- add velocity-Verlet integrator with Andersen thermostat (NVT)
    local integrator = mdsim.integrators.verlet_nvt_andersen({
        box = box, group = all_group["fluid"]
      , particle = particle, timestep = args.timestep, temperature = args.temperature, rate = 10
    })

    --****************************************
    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("pore_result.h5"):format(args.output), overwrite = args.overwrite})

    -- sample phase space
    -- write trajectory of particle groups to H5MD file
    halmd.observables.phase_space({box = box, group = all_group["fluid"]})
       :writer({file = file, fields = {"position", "species", "image", "velocity", "mass", "force"}, every = steps})
    halmd.observables.phase_space({box = box, group = all_group["pore"]})
       :writer({file = file, fields = {"position", "species"}, every = steps+1})
   
    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({
        steps = steps, first = 10, interval = 900, sample = 60
    })

    --write density modes to H5MD file
    local kmax = (201 + 1) / 2 * (2 * math.pi / edges[1][1])
    local density_wavevector = observables.utility.wavevector({box = box, wavenumber = {kmax}, filter = {1, 0, 0}, dense = true})
    local density_mode = observables.density_mode({group = all_group["fluid"], wavevector = density_wavevector})

    density_mode:writer({file = file, every = 300})


    -- run remaining part of the simulation in NVE ensemble
    -- to prepare for the NVE production run
    observables.sampler:run(steps)

    -- log profiler results
    utility.profiler:profile()
end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "initial", help = "prefix of output files"})
    parser:add_argument("input", {type = 'string', action = function(args, key, value)
	    readers.h5md.check(value)
	    args[key] = value
	    end, help = "input file h5"})
    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})
    parser:add_argument("particles", {type = "vector", dtype = "integer", default = {3000}, help = "number of particles"})
    parser:add_argument("density", {type = "number", default = 0.5, help = "particle number density"})
    parser:add_argument("ratios", {type = "vector", dtype = "number", action = function(args, key, value)
        if #value ~= 2 and #value ~= 3 then
            error(("box ratios has invalid dimension '%d' "):format(#value), 0)
        end
        args[key] = value
    end, default = {1, 1, 1}, help = "relative aspect ratios of simulation box"})
    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})
    parser:add_argument("temperature", {type = "number", default = 1.5, help = "target temperature"})
    parser:add_argument("rate", {type = "number", default = 4, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default =10 , help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.005, help = "integration time step"})
    parser:add_argument("slab",{type = 'number', default= 0.465})
    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 10, help = "for state variables"})
end
