import csv
from numpy.random import binomial
import datetime
from multiprocessing import Pool
import pickle
from utils import *



class IBMmap:
    """
    Container class for all hexagons on simulation map
    Has a few methods that should probably end up in
    simulation manager or hexes
    """
    def __init__(self):
        self.hexes = {}
        self.neighbor_manager = HexNeighbors()

    def from_pickle(self, path):
        """
        Load 'beagle_suitability.map' as simulation map
        Hex class has been modified since pickling
        so need to translate to new format hexes.
        #TODO fresh pickle
        :param path: path to pickled map
        :return:
        """
        with open(path, 'rb') as hex_pickle:
            tmp_hexes = pickle.load(hex_pickle)
        for key,hex in tmp_hexes.items():

            new_hex = Hex(self,key)
            new_hex.suitability = hex.properties['suitability']

            new_hex.set_quality()
            self.hexes[key] = new_hex

        for hex in self.hexes.values():
            hex.is_occupied = 0
            hex.fono = 0
            hex.set_fon()

    def get_neighbors(self,hex,end_distance,start_distance=0):
        return (self.hexes.get(hex) for hex in self.neighbor_manager.get_radius_coords(start_distance=start_distance,
                                                                                       origin=hex.axial_coords,
                                                                                       end_distance=end_distance))
    def get_neighbors_ring(self,hex,distance):
        return [self.hexes.get(hex) for hex in self.neighbor_manager.get_ring_coords(distance, hex.axial_coords)]

    def save_state(self,path):
        with open(path,'a') as f:
            writer = csv.writer(f)
            for key in (hex for hex in self.hexes if self.hexes[hex].is_occupied==1):
                writer.writerow([key])

class Hex:
    """
    Single hexagon class
    Each penguin should have a reference to a single hexagon.
    Hexagon stores nesting quality
    """
    def __init__(self,grid,axial_coords):
        self.grid = grid
        self.axial_coords = axial_coords
        self.fon = []
        self.fono = 0
        self.quality = 0

    def __str__(self):
        return "[Hex {},{},{}]".format(*self.axial_coords)

    def set_quality(self):
        """
        update hexagons quality - call whenever neighborhood changes
        :return: None
        """
        p = self.suitability + 1.15 * self.fono
        self.quality = np.exp(p) / (1 + np.exp(p))

    def get_neighbors(self, end_distance, start_distance=0):
        """
        Get neighbors of a hexagon in range between start and end
        :param end_distance: int
        :param start_distance: int
        :return: list of hexagons
        """
        return self.grid.get_neighbors(self, end_distance=end_distance, start_distance=start_distance)

    def get_neighbors_ring(self,distance):
        """
        Get neighbors of hexagon at single distance
        :param distance: int
        :return: list of hexagons
        """
        return self.grid.get_neighbors_ring(self,distance)

    def set_fon(self):
        """
        Get first order neighborhood (radius 1). This is used frequently so we
        cache it. Sets self.fon to list of hexagons
        :return:
        """
        self.fon = [hex for hex in self.grid.get_neighbors(self, 1) if hex is not None]

    def remove_neighbor(self):
        """
        Reduce count of occupied nearest neighbors by 1
        :return:  None
        """
        self.fono -= 1

    def add_neighbor(self):
        """
        Inrease count of occupied nearest neighbords by 1
        :return: None
        """
        self.fono += 1

    def unoccupied(self):
        """
        Set hexagon as unoccupied and update neighbors
        :return: None
        """
        self.is_occupied = 0
        for hex in self.fon:
            hex.remove_neighbor()
            hex.set_quality()

    def occupied(self):
        """
        Set hexagon as occupied and update neighbors
        :return: None
        """
        self.is_occupied = 1
        for hex in self.fon:
            hex.add_neighbor()
            hex.set_quality()

class Agent:
    """
    Class for an 'agent' representing a single penguin
    """
    __slots__ = ['hex', 'age', 'history', 'fidelity_threshold'] # Specifying slots reduces memory usage

    def __init__(self,hex=None):
        self.hex = hex
        if self.hex is not None:
            self.hex.occupied()
        self.age = 0
        self.history = []
        self.fidelity_threshold = 0.5

    @property
    def fidelity(self):
        # stay if successfull in t-1.. leave otherwise
        # Can be changed to alter nest site fidelity decision
        # should return truth (1/True/Objects...) for stay,
        # false for leave (False or None)
        # Self.history contains reproductive success in nest
        # return self.history[-1]
        return sum(self.history[-5:]) > 0
        # return sum(self.history[-3:]) > 1


    def get_fledge_probability(self):
        # parameters of reproductive success
        # Modify this function to modify probability of successful
        # reproduction
        alpha = -5.239
        beta = 0.928
        prob = inv_logit(alpha + beta * self.age+1)
        prob *= self.hex.quality
        return prob

    def fledged(self,prob = None):
        """
        Check for reproductive success
        Returns boolean of success, and updates
        reproductive history
        :param prob:
        :return:
        """
        if prob is None:
            fledging_prob = self.get_fledge_probability()
        else:
            fledging_prob = prob
        fledging_prob = min(fledging_prob,1)
        if binomial(1,fledging_prob) == 1:
            self.history.append(1)
            return True
        else:
            self.history.append(0)
            return False

    def survived(self,rates):
        """
        Check for survival. Returns boolean
        :param rates:
        :return:
        """
        # Juvenile or adult survival
        if self.age > 3:
            prob = rates[1]
        else:
            prob = rates[0]
        if np.random.binomial(1,prob):
            return True
        else:
            return False

    def n_neighbors(self,n):
        """
        Get number of occupied neighbors in ring at distance n
        :param n:
        :return:
        """
        return sum(1 for x in self.hex.get_neighbors_ring(n) if x is not None and x.is_occupied == 1)

    def get_move_options(self,target,distance,sort_func=None,extend=True):
        """
        Get list of available moves for a given starting point.
        Change this to change potential moves
        :param target: starting point hex
        :param distance: initial search distance
        :param sort_func: function used to sort options
        :param extend: should search radius be extended if initial search fails
        :return: list of avialable hexagons
        """
        options = [x for x in target.get_neighbors(distance) if x is not None and x.is_occupied == 0]
        while extend and len(options) < 1:
            new_distance = distance * 2
            options = [x for x in target.get_neighbors(end_distance=new_distance, start_distance=distance) if x is not None and x.is_occupied == 0]
            distance = new_distance
        if sort_func is not None:
            return sorted(options, key=sort_func)
        else:
            return options

    def move(self, hex):
        """
        Move agent to new hexagon
        :param hex: destination hexagon
        :return: None
        """
        # If current has nest, set nest location to unoccupied
        if self.hex is not None:
            self.hex.unoccupied()
        # Set nest site to new hexagon
        self.hex = hex
        # Update occupancy of new hexagon
        self.hex.occupied()
        self.history = []

    def die(self):
        # Set nest of dead penguin to unoccupied
        self.hex.unoccupied()
        self.hex = None

    def inc_age(self):
        # Increment age
        self.age += 1

class SimulationManager:
    def __init__(self,params):


        # General - Parameter setup
        self.params = params

        self.output_path = params['output_path']
        self.log_file = os.path.join(self.output_path,"log.csv")
        self.state_file = os.path.join(self.output_path, "{}_state.csv")
        self.timesteps = params['timesteps']
        self.save_state_frequency = params.get('save_state_frequency',1)
        self.save_stats_frequency = params.get('save_stats_frequency',1)

        self.cluster_hist_breaks = [10**x for x in range(6)]
        self.iterfp = IncrementOutputFile("frame",6,ext='.png')

        # Rates - Parameter setup
        # Setting these up for future use,
        # values are overwritten later
        self.adult_survival_prob = 0
        self.survival_rates = [0,0]

        self.kernel_size = params['kernel_size']

        # Map creation
        self.setup_map()

        self.vital_var = params['vital_var']

        # Agent list setup
        self.agents = []
        self.juveniles = []

        # Attributes
        self.time = 0
        self.pop_size = 0

    def setup_map(self):
        self.hexmap = IBMmap()
        self.hexmap.from_pickle('beagle_suitability.map')
        print("Loaded {} hexes".format(len(self.hexmap.hexes)))

    def seed(self,n,key=None):

        # Pick seed nest on basis of quality
        t1 = datetime.datetime.now()
        origin = self.hexmap.hexes[select_weighted(self.hexmap.hexes,lambda x: x.quality**2)]
        t2 = datetime.datetime.now() - t1
        print("found random dict item in {}".format(t2.total_seconds()))

        # Create founder penguin
        # Starting age 4 means adult survival is used
        # and chances of pop going extinct in first few steps
        # are lower
        parent = Agent(origin)
        parent.age = 4
        self.agents.append(parent)
        self.pop_size += 1

        # Create 'children' from founder, but set age to
        # 4 so that adult survival is used.
        for child in range(n):
            fledgling = Agent()
            fledgling.age = 4
            self._fledgling_move(fledgling, parent.hex)
            self.agents.append(fledgling)
            self.pop_size += 1

        print("seeded with {} agents".format(n))

    def set_survival_probability(self,var):
        alpha_adult = 10.818
        beta_adult = 2.641

        # Calculate params of beta distribution
        # which give required mean and variance
        # for juvenile survival
        alpha_juv,beta_juv = estBetaParams(0.7061055,var)
        self.survival_rates[1] = np.random.beta(alpha_adult,beta_adult)
        # Density dependent juvenile survival
        k = 800000 # Proportional to carrying capactiy
        self.survival_rates[0] = ((k-self.pop_size)/k) * np.random.beta(alpha_juv,beta_juv)


    def timestep(self):
        """
        Single time step of the model
        - Survive
        - Move
        - Update survival probabilities for next timestep
        - Reproduce
        - write stats if requested (modify frequency in pararms, defaults to every timestep).
        :return: Bool: population greater than 0
        """
        self.survive()
        self.move()
        self.set_survival_probability(self.vital_var)
        self.reproduce()
        # self.save_frame()
        if self.time % self.save_state_frequency == 0:
            self.save_state(self.state_file.format(self.time))
        if self.time % self.save_stats_frequency == 0:
            self.write_stats()
        # self.connected_component()
        self.time += 1

        if self.pop_size > 0:
            return True
        else:
            return False

    def run(self,timesteps = None):
        if timesteps is not None:
             self.timesteps = timesteps
        for step in range(self.timesteps):
            if not self.timestep():
                break

    def survive(self):
        """
        Go through agents and test for survival.
        If survived, increment age.
        If died, remove agent, decrement population size
        :return:
        """
        for ind, agent in enumerate(self.agents):
            if not agent.survived(self.survival_rates):
                        # agent.die sets nest to unoccupied
                        agent.die()
                        # remove dead individual, this pop may be a bottleneck
                        self.agents.pop(ind)
                        self.pop_size -= 1
                        # delete agent
                        # this is inefficient, also needs gc call
                        del agent
            else:
                agent.inc_age()

    def reproduce(self,prob=None):
        """
        Go through agents and test for reproduction.
        If successful, create offspring and find its
        new nest location. Increment population size.
        :param prob:
        :return:
        """
        for parent in self.agents:
            if parent.fledged(prob):
                fledgling = Agent()
                self._fledgling_move(fledgling,parent.hex)
                self.agents.append(fledgling)
                self.pop_size += 1

    def move(self):
        """
        Test each agent for fidelity, and move
        all that return false
        :return:
        """
        for agent in self.agents:
            if not agent.fidelity:
                options = agent.get_move_options(agent.hex, self.kernel_size, None, extend=True)
                target = random36.choices(population=options,weights=[x.quality**2 for x in options])
                agent.move(target[0])

    def _fledgling_move(self,fledgling,parent_hex):
        """
        Select an initial nest site for a fledgling, and move it
        :param fledgling: agent
        :param parent_hex: location of parent nest
        :return: None
        """
        options = fledgling.get_move_options(parent_hex, self.kernel_size, None, extend=True)
        target = random36.choices(population=options, weights=[x.quality**2 for x in options])
        fledgling.move(target[0])

    def get_average_age(self):
        """
        Get average ages of all living penguins
        :return: float: average age
        """
        return np.mean([agent.age for agent in self.agents])

    def get_average_survival(self):
        """
        Get mean of survival rates being used
        :return: float: average survival
        """
        return np.mean(self.survival_rates)

    def get_average_repro(self):
        """
        Get mean reproductive probability
        :return: float: mean reproductive probability
        """
        return np.mean([agent.get_fledge_probability() for agent in self.agents])

    def get_average_neighbors(self,radius):
        """
        Get average number of occupied nests at radius
        Discrete version of PCF
        :param radius: int
        :return: float: mean number of occupied nests
        """
        return np.mean([agent.n_neighbors(radius) for agent in self.agents])

    def connected_component(self):
        """
        find all connected components of occupied nests (ie. subcolonies)
        :return: number of clusters, histogram of cluster sizes
        """
        t1 = datetime.datetime.now()
        nodes = set(x.hex for x in self.agents)
        result = []
        while nodes:
            node = nodes.pop()
            # This set will contain the next group of nodes connected to each other.
            group = {node}
            # Build a queue with this node in it.
            queue = [node]
            # Iterate the queue.
            # When it's empty, we finished visiting a group of connected nodes.
            while queue:
                # Consume the next item from the queue.
                node = queue.pop(0)
                # Fetch the neighbors.
                neighbors = set(x for x in node.fon if x.is_occupied == 1)
                # Remove the neighbors we already visited.
                neighbors.difference_update(group)
                # Remove the remaining nodes from the global set.
                nodes.difference_update(neighbors)
                # Add them to the group of connected nodes.
                group.update(neighbors)
                # Add them to the queue, so we visit them in the next iterations.
                queue.extend(neighbors)

            # Add the group to the list of groups.
            result.append(len(group))
        td = datetime.datetime.now() - t1
        print("calculated {} connected components in {} seconds".format(len(result),td.total_seconds()))
        return len(result), np.histogram(result, self.cluster_hist_breaks)[0]

    def create_logs(self):
        """
        Create file to store stats at each iteration,
        and writes headers for csv
        :return:
        """
        print("creating logs...")
        with open(self.log_file,'w') as log:
            writer = csv.writer(log)
            writer.writerow(['population',
                             'avg_age',
                             'avg_surv',
                             'avg_repro',
                             # 'avg_neighbors_1',
                             # 'avg_neighbors_2',
                             # 'avg_neighbors_3',
                             # 'avg_neighbors_4',
                             # 'avg_neighbors_5',
                             # 'avg_neighbors_6',
                             # 'avg_neighbors_7',
                             # 'avg_neighbors_8',
                             'number_of_clusters',
                             'clusters_10e1',
                             'clusters_10e2',
                             'clusters_10e3',
                             'clusters_10e4',
                             'clusters_10e5'])
        print("Logs created @ {}".format(self.log_file))

    # Write current simulation statistics to file
    def write_stats(self):
        """
        Collect statistics on simulation and write to log file
        :return:
        """
        with open(self.log_file,'a') as output:
            writer = csv.writer(output)
            n_comps,comp_size = self.connected_component() # Calculate number of connected components (sub-colonies)
            writer.writerow([self.pop_size,
                   self.get_average_age(),
                   self.get_average_survival(),
                    # Nearest neighbor logging disabled for speed
                    # Use c++ tool to calculate nearest neighbors after runs
                    # or uncomment line below to calculate in python (slower)
                   # self.get_average_repro()] + [self.get_average_neighbors(r) for r in range(0,16)] +
                   self.get_average_repro()] +
                    [n_comps,",".join(map(str,comp_size))])

    # Save axial coordinates of occupied nests to state file
    def save_state(self,path):
        with open(path,'a') as f:
            writer = csv.writer(f)
            for agent in self.agents:
                writer.writerow([agent.hex.axial_coords]+[agent.age]) # added age



def simulate(params):
    """
    Wrapped function to run simulation.
    Creates a simulation manager, sets up some
    parameters, logs setup actions, seed sim,
    run
    :param params: dictionary of model run parameters
    :return: None
    """
    sys.setswitchinterval(1000) # slight performance boost, but not essential if it causes issues

    # extract parameters from argument, with defaults
    # if arguments aren't in params dict
    kernel_size = params.get('kernel_size',10)
    vital_var = params.get('vital_var',0.006)
    timesteps = params.get('timesteps',1500)
    path = params['path'] # no default path

    # Create log of parameters being used
    with open('{}/params.txt'.format(path), 'w') as log:
        print(params, file=log)
        print(kernel_size, file=log)

    # Run the simulation, logging any errors that
    # occur to 'error.txt' inside the output subdirectory
    with open('{}/error.txt'.format(path), 'w') as log:
        try:
            params['kernel_size'] = kernel_size
            params['vital_var'] = vital_var
            params['output_path'] = path
            params['timesteps'] = timesteps
            print("Setting output to {}".format(path),file=log)

            # Create simulation manager with specified parameters
            log.write("establishing simulation\n")
            sim = SimulationManager(params)
            log.write('simulation established\n')

            # Seed simulation with initial nests
            log.write("seeding simulation\n")
            sim.seed(500)

            # Setup output files
            sim.create_logs()

            # Run model for set number of timesteps
            sim.run(timesteps)

        # Catch any errors, write them to log file, then reraise error
        except Exception as e:
                log.write('an error occured\n')
                log.write(str(e))
                raise

if __name__ == '__main__':
    import argparse

    # Arguments passed to script from stdin
    # these are unnamed arguments
    # eg...
    #  model.py 0.5 output 1
    # runs the model with variance 0.5 in repro rates,
    # outputs to the subdirectory 'output'
    # and uses 1 core
    parser = argparse.ArgumentParser(description='Hexagonal IBM model')
    parser.add_argument('variance',help='variance on repro rates')
    parser.add_argument('output', help='output folder')
    parser.add_argument('ncores',help='number of cores')
    args = parser.parse_args()


    # setup unique output directory generator
    try:
        # op is a function which returns a new unique subdirectory
        # of args.output each time its called
        op = UniqueOutputDirectory(args.output)
    except Exception as e:
        with open('logfile.txt', 'w') as log:
            print(e, file=log)
            raise

    # Logging info
    with open('logfile.txt','w') as log:
        print("created outputdir at {}".format(op),file=log)

    # Parameters are passed to the model as a dictionary
    # 'Jobs' is a list of parameters to run the model one
    # 'path' is set as the result of op() which generates
    # a unique directory each time its called
    jobs = []
    # Currently set up to run scenarios of varying kernel size
    # for 10 reps. To change the model thats run, fill 'jobs'
    # with dictionaries of parameters.
    # Currently params can be kernel_size, vital_var, path,
    
    #for rep in range(10):
        #for k in range(10):
            #params = {'path':op(),'kernel_size':(1+k),'vital_var':float(args.variance)}
            #jobs.append(params)
    jobs = [{'path':op(),'timesteps':50}] 
    # Logging info
    with open('logfile.txt', 'a') as log:
        print(jobs, file=log)
        print("created {} tasks".format(len(jobs)), file=log)

    pool = Pool(processes=int(args.ncores))  # start worker processes
    results = pool.map(simulate, jobs) # Run one instance of simulate on each worker process




