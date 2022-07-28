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
    -- total number of particles from sum of particles
    local length = {}
    	length[1]=40
	length[2]=10
	length[3]=25

    local dimension = 3
    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({length = length})

    local nparticle = 2600 -- math.ceil( (args.slab*length[1]*length[2]*length[3]*4)/math.pow(1.748, 3) )
   

    -- create system state
    local particle = mdsim.particle({dimension = dimension, particles = nparticle})

    -- set initial particle positions sequentially on an fcc lattice
    local lattice = mdsim.positions.lattice({box = box, particle = particle, slab = {args.slab,1,1} })
    lattice:set()
    position = particle.data["position"]

    --create the slot
    local pore = mdsim.geometries.cuboid({ lowest_corner = {-args.slab*length[1]/2, -args.pore_diameter/2, -args.pore_diameter/2 }, length = { args.slab*length[1], args.pore_diameter, args.pore_diameter } } )
    local pore_group = mdsim.particle_groups.region_species({particle = particle, box = box, species = 0, geometry = pore, selection = 'excluded', label = 'pore'})
  

    --steps
    local steps = math.ceil(args.time/args.timestep)

    -- H5MD file writer
    local file = writers.h5md({path = ("obstacle_slab.h5"):format(args.output), overwrite = true})


    -- select all particles
    local all_group = mdsim.particle_groups.all({particle = particle, label = "all"})

    --sample phase space
    local phase_space = observables.phase_space({box=box, group = pore_group , every = steps})

    --write positions in h5 file
    --if steps > 0 then
    phase_space:writer({file=file, fields={'position'}})
    --end

    --write observables to the file, find msv
   -- local interval = args.sampling.trajectory
   -- local msv = observables.thermodynamics({box=box, group=all_group})
   -- msv:writer({file=file, fileds={'msv'}, every = interval})


    observables.sampler:sample()
    
    --run simulation
    local integrator = mdsim.integrators.verlet_nvt_andersen({
    box=box
    , group = pore_group
    , particle = particle
    , timestep = args.timestep
    , temperature = args.temperature
    , rate = args.rate})
   
    observables.sampler:run(steps)




end

--
-- Parse command-line arguments.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "obstacle_slab", help = "prefix of output files"})
    parser:add_argument("temperature", {type = "number", default = 0, help = "target temperature"})
    parser:add_argument("rate", {type = "number", default = 4, help = "heat bath collision rate"})
    parser:add_argument("time", {type = "number", default =10 , help = "integration time"})
    parser:add_argument("timestep", {type = "number", default = 0.005, help = "integration time step"})
    parser:add_argument('slab', {type = 'number', default = 0.465, help = 'box fraction occupied'})
    parser:add_argument('pore_diameter', {type = 'number', default = 5, help = 'pore diameter'})

end
