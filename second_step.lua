--
-- Copyright © 2019 Abbas Gholami, Roya Ebrahimi
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
local writers = halmd.io.writers
local readers = halmd.io.readers
local utility = halmd.observables.utility
local random = halmd.random
local random1 = math.random
math.randomseed(os.time())

function main(args)
    --parameters
    local timestep = 0.002
    local temperature = 1.0
    local steps = math.ceil(args.time / timestep)
    local equibliration_steps = 0.5 * steps --math.min(50000, math.floor(steps/2))
    local nknots = {args.nknots, 2, 2}
    local Source_size = 3
    local rest_size = 40 - Source_size
    local eps = 0.001
    local source_velocity = {0.25,0,0}

    --open H5MD file for reading and read from last step the fluid and pore parameters
    local file_read = readers.h5md({path = args.input})
    local samples = {}
    local reader, fluid_sample = observables.phase_space.reader({file = file_read, location = {"particles", "fluid"}, fields = {"position", "velocity", "species", "mass"}})
    samples["fluid"] = fluid_sample
    reader:read_at_step(-1)
    local nfluid = fluid_sample.nparticle
    
    local reader, pore_sample = observables.phase_space.reader({file = file_read, location = {"particles", "pore"}, fields = {"position", "species"}})
    samples["pore"] = pore_sample
    reader:read_at_step(-1)
    local npore = pore_sample.nparticle
   
    -- determine system parameters 
    local nspecies = 2
    local dimension = 3

    -- read edge vectors of simulation domain from file
    local edges = mdsim.box.reader({file = file_read, location = {"particles", "fluid"}})

    -- create simulation domain with periodic boundary conditions
    local box = mdsim.box({edges = edges})

    local particle = {
	    fluid = mdsim.particle({dimension = dimension, particles = nfluid, species = nspecies, label = "fluid"})
	  , pore = mdsim.particle({dimension = dimension, particles = npore, species = nspecies, label = "pore"})
  }

    --all_group
    local all_group = {}
    for label, p in pairs(particle) do
	    all_group[label] = mdsim.particle_groups.all({particle = p , label = label})
    end

    local phase_space = {}
    for label, g in pairs(all_group) do
	    phase_space[label] = observables.phase_space({box = box, group = g})
    end
    phase_space["fluid"]:set(fluid_sample)
    phase_space["pore"]:set(pore_sample)

    local source_group = mdsim.particle_groups.region({
    	particle = particle["fluid"]
   	, selection = "included"
    	, geometry = mdsim.geometries.cuboid({lowest_corner = {-edges[1][1]/2 -eps , -edges[2][2]/2, -edges[3][3]/2}, length = {Source_size + eps, edges[2][2] , edges[3][3]}})
	, box = box
	, label = "source"
	})

    local rest_group = mdsim.particle_groups.region({
	particle = particle["fluid"]
    	, selection = "included"
	, geometry = mdsim.geometries.cuboid({lowest_corner = {-edges[1][1]/2 + Source_size , -edges[2][2]/2, -edges[3][3]/2}, length = {eps + rest_size, edges[2][2] , edges[3][3]}})
	, box = box
	, label = "rest"})

    local potential = mdsim.potentials.pair.lennard_jones({
       epsilon = {
            {1, 1} -- onAT_fromAT, onAT_fromObst
          , {1, 0} -- onObst_fromAT, onObst_fromObst
        }
      , sigma = {
            {1, 1.5}
          , {1.5, 2}
       }
    })

    local cut = math.pow(2,1/6)
    potential = potential:truncate({"smooth_r4", 
           cutoff = {
            {2.5,cut}
          , {cut,cut}
        }
          , h = 0.005
        })

    -- compute forces
    local binning = {
	    fluid = mdsim.binning({
		    box = box
		  , particle = particle["fluid"]
--		  , skin = 1
		  , r_cut = 2.5
	  })
	   , pore = mdsim.binning({
	            box = box
		  , particle = particle["pore"]
--		  , skin = 1
		  , r_cut = 2.5
	  })
   }

	

    for label, p2 in pairs(particle) do
	    local nieghbour = mdsim.neighbour({
		    box = box
		  , particle = { particle["fluid"], p2 }
		  , r_cut = potential.r_cut
		  , binning = { binning["fluid"], binning[label] }
	  })


	    mdsim.forces.pair({
		    box = box
	          , particle = {particle["fluid"], p2} 
		  , potential = potential
		  , neighbour = neighbour
	  })
	end
    
    
    observables.sampler:sample()
    log.info(("Closing H5MD file %s"):format(file_read.path))
    file_read:close()


    --integrator
    local source_integrator = mdsim.integrators.verlet_nvt_boltzmann_value({group = source_group, mean_velocity = source_velocity, box = box, timestep = timestep, temperature = temperature, rate = 8})
    local rest_integrator = mdsim.integrators.verlet({group = rest_group, box = box, timestep = timestep})
    -- H5MD file writer
    local file_write = writers.h5md({path = args.output, mode = "truncate", overwrite = true})

    -- Sample macroscopic state variables.
    local msv
    local interval = 300
    if interval > 0 then
        msv = observables.thermodynamics({box = box, group = all_group["fluid"]})
        msv:writer({file = file_write, every = interval})
    end
    local runtime = observables.runtime_estimate({
            steps=steps, first=10,interval=900,sample=60
    })

    -- run half of the simulation
    observables.sampler:run(equibliration_steps) 

    --writer

    phase_space["fluid"]:writer({file = file_write, fields = {"position", "velocity", "species", "force", "image","mass"}, every = 300})
    phase_space["pore"]:writer({file = file_write, fields = {"position", "velocity", "species", "force", "image","mass"}, every = 300})

    -- Create 16 slabs
    local geometry = {}
    local group = {}
    local msv = {}

    for i= 0,15 do
        geometry[i] = mdsim.geometries.cuboid({lowest_corner = {(-edges[1][1]/2+i*2.5),-edges[2][2]/2, -edges[3][3]/2}, length = {2.5,edges[2][2], edges[3][3]}})

        group[i] = mdsim.particle_groups.region_species({particle = particle["fluid"] ,species = 0, geometry = geometry[i], selection = "included", box = box, label = string.format("%s%d",'region' , i)})
    end

    for i= 0,15 do
        msv[i] = observables.thermodynamics({box = box, group = group[i]})
    end

    for j = 0,15 do
        msv[j]:writer({
        file = file_write
            , fields = {
                  "potential_energy", "pressure", "temperature" , "total_force"
                , "internal_energy", "kinetic_energy","virial", "center_of_mass_velocity" , "heat_flux"
              }
            , every = 300
          })
      end
    --write density modes to H5MD file
    local kmax = (nknots[1] + 1) / 2 * (2 * math.pi / box.length[1])
    local wavevector = observables.utility.wavevector({box = box, wavenumber = {kmax}, filter = {1, 0, 0}, dense = true})
    local density_mode = observables.density_mode({group = all_group["fluid"], wavevector = wavevector})
    density_mode:writer({file = file_write, every = 300})
    local current_density_mode = observables.current_density_mode({group = all_group["fluid"], wavevector = wavevector})
    current_density_mode:writer({file = file_write, every = 300})
    local kinetic_energy_density_mode = observables.kinetic_energy_density_mode({group = all_group["fluid"], wavevector = wavevector})
    kinetic_energy_density_mode:writer({file = file_write, every = 300})

--average run
    -- run rest of the simulation
    observables.sampler:run(steps - equibliration_steps)
    
    -- log profiler results
    halmd.utility.profiler:profile()


end

-- Define command line parser.
--
function define_args(parser)
    parser:add_argument("output,o", {type = "string", action = parser.action.substitute_date_time,
        default = "result_out_%Y%m%d_%H%M%S", help = "prefix of output files"})
    parser:add_argument("overwrite", {type = "boolean", default = false, help = "overwrite output file"})

    parser:add_argument("input", {type = "string", required = true, action = function(args, key, value)
        readers.h5md.check(value)
        args[key] = value
    end, help = "H5MD input file"})

    parser:add_argument("time", {type = "number", default = 1000, help = "integration time"})
    parser:add_argument("temperature", {type = "number", default = 1.5, help = "temperature of heat bath"})
    parser:add_argument("nknots", {type = "number", default = 201, help = "number of knots"})

    parser:add_argument("smoothing", {type = "number", default = 0.005, help = "cutoff smoothing parameter"})
    parser:add_argument("cutoff", {type = "float32", default = math.pow(2, 1 / 6), help = "potential cutoff radius"})
    parser:add_argument("force_threshold", {type = "float32", default = 50, help = "force_threshold"})
    local sampling = parser:add_argument_group("sampling", {help = "sampling intervals (0: disabled)"})
    sampling:add_argument("trajectory", {type = "integer", help = "for trajectory"})
    sampling:add_argument("state-vars", {type = "integer", default = 1000, help = "for state variables"})
    sampling:add_argument("structure", {type = "integer", default = 1000, help = "for density modes, static structure factor"})
    sampling:add_argument("correlation", {type = "integer", default = 100, help = "for correlation functions"})
    sampling:add_argument("average", {type = "integer", help = "output averages of given number of samples"})

    local wavevector = parser:add_argument_group("wavevector", {help = "wavevector shells in reciprocal space"})
    observables.utility.wavevector.add_options(wavevector, {tolerance = 0.01, max_count = 7})
    observables.utility.semilog_grid.add_options(wavevector, {maximum = 15, decimation = 0})
end
