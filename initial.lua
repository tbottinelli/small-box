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
    local nobstacle = assert(sample.nparticle)
    print(nobstacle)

    local length = {}
        length[1] = 40
        length[2] = 10
        length[3] = 25

    local dimension = 3
    local eps = 0.01
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nfluid+nobstacle, species = 2})
    local obstacle = mdsim.particle({particles = nobstacle, dimension = dimension, species = 1})
    local obst_group = mdsim.particle_groups.all({particle = obstacle, label = "allobst"})
    local obst_phasespace = observables.phase_space({box = box, group = obst_group})
    obst_phasespace:set(sample)
    local obst_pos = obstacle.data["position"]

    -- set particle species
    local species = {}
    local positions = {}
    for i = 1, nfluid do table.insert(species,0) end
    for i = nfluid + 1, nfluid+nobstacle do table.insert(species, 1) end
    particle.data["species"] = species

    --local groups = { mdsim.particle_groups.id_range({ particle = particle, range = {1,nfluid}, label = 'fluid'}),
    --                 mdsim.particle_groups.id_range({ particle = particle, range = {nfluid + 1,nfluid+ nobstacle}, label = 'obst'})
--	     }

   -- local particle_ = {
   --   fluid = groups[1]:to_particle({label = 'fluid'})
   -- , obst = groups[2]:to_particle({label = 'obst'})
   -- }
    
   -- particle = particle_; particle_ = nil

    -- set initial particle positions sequentially on an fcc lattice
    local lattice = mdsim.positions.lattice({box = box, particle = particle, slab = {0.6,1,1}})
    lattice:set()
    positions = particle.data["position"]

    -- select all particles
    local allbox = mdsim.geometries.cuboid({lowest_corner = {-length[1]/2 - eps ,-length[2]/2,-length[3]/2}, length = {length[1]+2*eps,length[2],length[3]}})
    local fluid_group = mdsim.particle_groups.region_species({ particle = particle, box = box, species=0,selection="included",label = "fluid",geometry=allbox})
    local all_fluid_group = mdsim.particle_groups.id_range({particle = particle, range={1,nfluid},label="allfluid"})
    local all_obst_group = mdsim.particle_groups.id_range({particle = particle, range={nfluid+1,nfluid+nobstacle},label="allobst"})
    local all_group = mdsim.particle_groups.all({particle = particle, label="all"})
    -- set initial particle velocities
    local boltzmann = mdsim.velocities.boltzmann({
        particle = particle
      , group = all_group
      , temperature = args.temperature
    })
    boltzmann:set()
    local velocities = particle.data["velocity"]
    
    --oss lua tables start with index 1 so this is for position x
    for i = 1,nfluid do
       -- if (positions[i][1] < 0) then 
         --   positions[i][1] = positions[i][1] - ((args.slab*length[1]/2) +4)
       -- end
       -- if (positions[i][1] >= 0) then
         --   positions[i][1] = positions[i][1] + ((args.slab*length[1]/2) +4)
       -- end
       velocities[i]= {1,1,1}
    end
    for i = nfluid+1,nfluid+nobstacle do
        velocities[i] = {0,0,0}
        positions[i] = obst_pos[i - nfluid]
        --positions[i][1] = positions[i][1] + 5
    end
    particle.data["position"] = positions
    particle.data["velocity"] = velocities

    local potential = mdsim.potentials.pair.lennard_jones({species = 2,
        epsilon = 1
      , sigma = 1
    })
    ct = 2.5 --math.pow(2, 1 /6) 
    -- smoothing at potential cutoff
    potential = potential:truncate({"smooth_r4",
       cutoff = ct
      , h = 0.005
    })

    -- compute forces
    local force = mdsim.forces.pair({
        box = box
      , particle = particle
      , potential = potential
    })
   -- local force1 = mdsim.forces.pair({box = box, particle = {all_fluid_group, all_fluid_group}, potential = potential})
   -- local force2 = mdsim.forces.pair({box = box, particle = {obstacle_group, obstacle_group} , potential = potential})


    -- add velocity-Verlet integrator with Andersen thermostat (NVT)
    local integrator = mdsim.integrators.verlet_nvt_andersen({
        box = box, group = fluid_group
      , particle = particle, timestep = args.timestep, temperature = args.temperature, rate = 10
    })

    --****************************************
    -- convert integration time to number of steps
    local steps = math.ceil(args.time / args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("pore_result.h5"):format(args.output), overwrite = args.overwrite})

    -- sample phase space
    local phase_space = halmd.observables.phase_space({box = box, group = all_group})

    -- write trajectory of particle groups to H5MD file
    phase_space:writer({file = file, fields = {"position", "species", "image", "velocity", "mass", "force"}, every = steps})

    -- sample initial state
    observables.sampler:sample()

    -- estimate remaining runtime
    local runtime = observables.runtime_estimate({
        steps = steps, first = 10, interval = 900, sample = 60
    })

    --write density modes to H5MD file
    local kmax = (201 + 1) / 2 * (2 * math.pi / box.length[1])
    local density_wavevector = observables.utility.wavevector({box = box, wavenumber = {kmax}, filter = {1, 0, 0}, dense = true})
    local density_mode = observables.density_mode({group = all_fluid_group, wavevector = density_wavevector})

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
    parser:add_argument('input', {type = 'string', action = function(args, key, value)
	    readers.h5md.check(value)
	    args[key] = value
	    end, help = 'input file h5'})
    parser:add_argument("overwrite", {type = "boolean", default = true, help = "overwrite output file"})
    parser:add_argument("particles", {type = "vector", dtype = "integer", default = {2000}, help = "number of particles"})
    parser:add_argument("density", {type = "number", default = 0.35, help = "particle number density"})
    parser:add_argument("ratios", {type = "vector", dtype = "number", action = function(args, key, value)
        if #value ~= 2 and #value ~= 3 then
            error(("box ratios has invalid dimension '%d'"):format(#value), 0)
        end
        args[key] = value
    end, default = {1, 1, 1}, help = "relative aspect ratios of simulation box"})
    parser:add_argument("masses", {type = "vector", dtype = "number", default = {1}, help = "particle masses"})
    parser:add_argument("initial-temperature", {type = "number", default = 1.5, help = "initial temperature"})
    parser:add_argument("temperature", {type = "number", default = 1, help = "target temperature"})
    parser:add_argument("rate", {type = "number", default = 4, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default =10 , help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.005, help = "integration time step"})
    parser:add_argument('slab',{type = 'number', default= 0.465})

    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 10, help = "for state variables"})
end
