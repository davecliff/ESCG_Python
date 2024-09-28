# -*- coding: utf-8 -*-
#
# RPSLS_ESCG.py: Rock-Paper-Scissors-Lizard-Spock (RPSLS) Evolutionary Spatial Cyclic Game (ESCG).
#
#
# Copyright (c) 2024, Dave Cliff
#
#
# ------------------------
#
# MIT Open-Source License:
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
# associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial
# portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ------------------------

# This started as an implementation and extension of the kind of ablated RPSLS models explored in the paper
# "Species coexistence in spatial cyclic game of five species"
# by Linwu Zhong, Liming Zhang, Haihong Li, Qionglin Dai & Junzhong Yang,
# Published in the journal "Chaos, Solitons and Fractals" 156 (2022) 111806.
# Referred to in comments below as ZZLDY22.

# The code incorporates the "Orginal Elementary Step" (OES) as used by Zhong et al. (2022) and also
# the "Revised Elementary Step" (RES) introduced in this paper:
# Cliff, D. (2024a) "Never Mind The No-Ops: Faster and Less Volatile Simulation Modelling of Co-Evolutionary Species
# Interactions via Spatial Cyclic Games". Accepted for presentation/publication at the 36th European Modeling and
# Simulation Symposium (EMSS2024); Tenerife, Spain; 18th-20th September 2024.
# Available at SSRN: https://ssrn.com/abstract=4883174

# The code also implements parameterised circulant networks as reported in this paper:
# Cliff, D. (2024b) "Tournament versus Circulant: On Simulating 7-Species Evolutionary Spatial Cyclic Games with
# Ablated Predator-Prey Networks as Models of Biodiversity". Accepted for presentation/publication in Proceedings of
# the 36th European Modeling and Simulation Symposium (EMSS2024). Tenerife, Spain, 18-20 September 2024.
# Available at SSRN: https://ssrn.com/abstract=4961889

# Some suggestions on getting results with less compute...
# ZZLDY22 use N=200^2 but their Figure 4 shows comparable results should be possible at N=50^2 or N=75^2.
# ZZLDY22 run experiments for 100000 Monte Carlo Steps (MCSs) but their Fig 2 shows "asymptotic" state once MCS=50000
# ZZLDY22 show all M-response data for M swept over four orders of magnitude [10E-7, 10E-3] but...
# ZZLDY22 Figures 3, & 5-7 all show approx steady-state over M in [10E-4, 10E-3], so could possibly end runs at ~10E-4.
# ZZLDY22 Figs 3, 5, & 6 show approx steady-state for M in [10E-7, 10E-6] so could possibly start runs at ~10E-6.
# ZZLDY22 Figs 3, 5, 6 & 7 show 5 data points per power-of-ten range, but could use fewer (e.g. 4 or 3 or 2).
#
# The code here has IID repetition in the outer loop, and M-sweep on the inner loop, to give whole M-sweep data first

import sys
import os
import math
import random
import numpy as np
# import graphviz

if __name__ == "__main__":

    # key hyperparameters from ZZLDY2022
    L = 200   # side-length of square grid of cells
    N = L * L           # number of cells
    mu = 1.0            # probability of state transition in ZZLDY22's Equations 1 and 2 (competition)
    sigma = 1.0         # probability of state transition in ZZLDY22's Equation 3 (reproduction)
    M_min = 1.00E-7     # lowerbound on M range  -- ZZDLY22 used 1E-7 for all figures except 10E-8 for Fig 2 (heatmap)
    M_max = 1.00E-3     # upperbound on M range  -- ZZDLY22 used 1E-3
    n_mcs = 100000      # number of "Monte Carlo Steps" in any one experiment -- ZZDLY22 used 100000
    bc_period = False   # if True use periodic boundary condition (toroidal wraparound); otherwise no-flux boundary
    ns_oom = 4          # Number of Steps per Order Of Magnitude: how many data-points per power-of-10 range for M sweep
    n_generations = N   # ZZLDY unclear on this
    nbrhood_count = 4   # Should be either 4 (von Neumann) or 8 (Moore)

    n_expts = 64        # number of IID experiments to run for any one value of M
    ablation_mcs = 0    # ablation occurs at the start of this MCS
    
    n_species = 5       # how many species
    n_ablations = 0     # how many links in the dominance network to randomly select and delete
    cx_offsets = []     # connection offsets for the dominance-network constructor: leave blank for tournament graph
    if cx_offsets == []:
        for s in range(1, int(math.floor(n_species/2))+1):
            cx_offsets.append(s)

    verbose = False          # print messages giving a running commentary on progress, or not
    use_OES = False          # which version of elementary step to use: OES or RES
    draw_networks = False    # draw network diagrams?
    write_domnets = True     # if True then write the dominance networks (original and then ablated) for each experiment
    write_grids = True       # if True then record the grid at every g_write_freq MCsteps
    g_write_freq = 5000      # frequency of writing grids
    empty_cell = -1         # constant value to mark empty cell in grid

    if use_OES:
        es_str = 'OES'
    else:
        es_str = 'RES'
    es_str += '_NS%02d_NA%d_L%d_MCS%.2e' % (n_species, n_ablations, L, n_mcs)
    exp_id = es_str + '_M{:.2e}_{:.2e}'.format(M_min, M_max)  # header string for folder name


    def populate_grid(sl, ns, vrbs):
        """
        Populate a square LxL grid/lattice of cells
        :param sl: side-length of the grid
        :param ns: number of species
        :param vrbs: verbosity
        :return: the grid, as a numpy array
        """
        g = np.zeros((sl, sl))
        for x in range(sl):
            for y in range(sl):
                g[x, y] = np.random.randint(0, ns)
        if vrbs:
            print(g)
        return g


    def write_grid(mcs, g, f):
        """
        Write the grid of cells as csv
        :param mcs: the current Monte-Carlo-step
        :param g: the grid to be written
        :param f: the file to be written to
        :return: <nothing>>
        """
        sl = np.shape(g)[0]
        f.write('MCS=,{}\n'.format(mcs))
        for gridrow in range(sl):
            for gridcol in range(sl):
                f.write('%d,' % g[gridrow, gridcol])
            f.write('\n')
        f.write('\n')

    def draw_net(net, fns):
        """
        Use graphviz to draw network diagram
        :param net: the network to draw
        :param fns: the file-name string, to be written to
        :return: <nothing>>
        """
        n_sp = len(net)

        fns += '.gv'
        # g = graphviz.Digraph(filename=fns, format='png', engine='circo') doesn't preserve node orderings
        # which looks bad on smaller networks but is not a problem for bigger ones
        # so for smaller nets instead calculate a set of x,y positions and use engine='neato'

        if n_sp > 21:
            g = graphviz.Digraph(filename=fns, format='png', engine='circo')
        else:
            g = graphviz.Digraph(filename=fns, format='png', engine='neato')

        radius = 50
        radial = (2 * math.pi) / n_sp
        for sp in range(n_sp):
            label = str(sp)
            x = math.cos(sp * radial) * radius
            y = math.sin(sp * radial) * radius
            pos_str = '%f,%f' % (x, y)
            g.node(label, shape='circle', pos=pos_str)

        for sp in net:
            sp_label = str(sp[0])
            for dominated_sp in sp[1]:
                dom_label = str(dominated_sp)
                g.edge(sp_label, dom_label, len='5')
        g.render(view=True)

    
    def species_density(ns, g, vrbs):
        """
        calculate the density of each species in the grid
        :param ns: number of species
        :param g: the grid
        :param vrbs: verbosity
        :return: list of densities for each species
        """
        densities = np.zeros(ns)
        sl = np.shape(g)[0]
        i_count = 0
        # count absolute numbers
        for x in range(sl):
            for y in range(sl):
                i = int(g[x, y])
                if i >= 0:
                    # nonempty cell
                    dns = densities[i]
                    densities[i] = dns + 1
                    i_count += 1
        # convert absolute counts to percentages, and count number of nonzero species
        percentages = np.zeros(ns)
        species_count = 0
        for sp in range(ns):
            percentages[sp] = densities[sp] / i_count
            if densities[sp] > 0:
                species_count += 1
        if vrbs:
            print('percentages={}, s_count={}'.format(percentages, species_count))

        return percentages, species_count

    
    def dominance_network(ns, cxspec, vrbs):
        """
        Construct the full dominance network, then optionally delete some links from it.
        :param ns: number of species
        :param cxspec: connectivity specification -- list of ringplus offsets
        :param vrbs: verbosity
        :return: the dominance network: for each species, list of other species that it dominates
        """
        # create full dominance network using list of ringplus offsets
        if vrbs:
            print("Creating full dominance net for {} species with connection offsets ={}".format(n_species, cxspec))
        network = np.zeros((ns, ns))
        for sp in range(ns):
            for cx_offset in cxspec:
                dominated = (sp + cx_offset) % ns
                network[sp, dominated] = 1
                if vrbs:
                    print('Species {} dominates species {}'.format(sp, dominated))
        if vrbs:
            print("dominance_network:")
            print(network)

        return network

  
    def dnet_ablate(na, net, vrbs):
        """
        Randomly ablate the dominance network.
        :param na: number of dominance links to ablate from the network
        :param net: the network to ablate
        :param vrbs: verbosity
        :return: the ablated dominance network
        """
        if na > 0:
            # to generate comparable results, be consistent in choosing which species gets ablation
            # always start the ablation at species 0
            for sp_del in range(na):
                # Todo: Come back and tidy this up, make it more general.
                #  Currently written just for a single random ablation from species 0
                dominator_species = 0
                n_dominated = np.count_nonzero(net[dominator_species, :])
                dom_count = 0
                to_be_ablated = np.random.randint(1, n_dominated + 1)
                for dominated_index in range(net.shape[0]):
                    if net[dominator_species, dominated_index] == 1:
                        dom_count += 1
                    if dom_count == to_be_ablated:
                        net[dominator_species, dominated_index] = 0
                        to_be_ablated = -1
                        if vrbs:
                            print('species {} dominating species {} ablated'.format(dominator_species, dominated_index))
        if vrbs:
            print("ablated dominance_network:")
            print(net)

        return net

    
    def write_dnet_adjmat(dn, ofile):
        """
        Write adjacency matrix for a dominance network to a csv output file
        :param dn: the network to write the adjacency matrix for
        :param ofile: the file to write to
        :return: <nothing>"
        """
        n_species = dn.shape[0]
        for species in range(n_species):
            for dominated_opponent in range(n_species):
                ofile.write('%d,' % dn[species, dominated_opponent])
            ofile.write('\n')
        ofile.write('\n')

    def generate_m_vals(min_m, max_m, nsteps_oom, logstep, vrbs):
        """
        Generate a list of the M values that experiments will be run at,
        This is only called once
        :param min_m: lowerbound on m values
        :param max_m: upperbound on m values
        :param nsteps_oom: number of discrete steps in m-value per order of magnitude
        :param logstep: Boolean; if True then steps are spaced equally on log scale; otherwise equally on linear scale.
        :param vrbs: verbosity
        :return: sorted list of M values
        """
        m_vals = []
        if vrbs:
            print('generate_m_vals({:e}, {:e}, {}, logstep={})'.format(min_m, max_m, nsteps_oom, logstep))
        
        if logstep:
            
            # find the order of magnitude that Min_M is in
            # -- NB expts don't necessarily start at M=M_min; M_Min is lower bound
            this_M_oom = math.pow(10, math.floor(math.log(M_min, 10)))
            next_M_oom = this_M_oom * 10
            
            # find the first m_val step that is greater than min_m
            oom_step = 0
            m_val = min_m - 1
            while m_val < min_m:
                m_val = math.pow(10, math.log(this_M_oom, 10) + (oom_step / nsteps_oom))
                oom_step += 1
            
            while m_val <= max_m:
                
                if vrbs:
                    print('OOM:[{:.2e},{:.2e}]'.format(this_M_oom, next_M_oom))

                m_vals.append(m_val)
                oom_step += 1
                
                if m_val >= next_M_oom or oom_step >= ns_oom:
                    m_val = next_M_oom
                    this_M_oom = next_M_oom
                    next_M_oom = this_M_oom * 10
                    oom_step = 0
                
                m_val = math.pow(10, math.log(this_M_oom, 10) + (oom_step / nsteps_oom))
                    
        if vrbs:
            print('M_vals=[', end='')
            for mv in m_vals:
                print('%.2e, ' % mv, end='')
            print(']')
            
        return m_vals


    def rand_nbr(sidelen, ix, iy, nbr_offsets, max_nbr_offset_index, pbc=False):
        """
        Random choice of a neighbor of the cell at (ix,iy) in squaregrid of sidelength sl.
        This gets called millions of times, so it is written for fast execution, rather than elegance
        :param sidelen: sidelength of square grid
        :param ix: x-coord of individual in grid to choose neighbor for
        :param iy: y-coord of individual in grid to choose neighbor for
        :param nbr_offsets: numpy array 4x2 or 8x2 of x,y offsets for either von Neumann or Moore nbrhoods
        :param max_nbr_offset_index: highest index in the nbr_offsets array -- either 3 (von Neumann) or 7 (Moore)
        :param pbc: if True, periodic boundary condition (toroidal wraparound); if False, no-flux (cutoff)
        :return: (nx, ny) coords of randmly chosen neighbor
        """
            
        nx, ny = ix, iy
         
        while nx == ix and ny == iy:
            # randomly choose integer index into offset array
            offset_index = np.random.randint(0, max_nbr_offset_index + 1)
            nbr_oset = nbr_offsets[offset_index, :]
            
            if pbc:
                nx = (ix + nbr_oset[0]) % sidelen
                ny = (iy + nbr_oset[1]) % sidelen
            else:
                nx = min(max(0, (ix + nbr_oset[0])), sidelen - 1)
                ny = min(max(0, (iy + nbr_oset[1])), sidelen - 1)
                    
        return nx, ny


    def original_elementary_step(g, domnet, n_empty, p_mu, p_sigma, p_epsilon, periodic_boundary, nbr_offsets,
                                 max_offset_i, vrbs):
        """
        The original elementary step, as used by many authors
        :param g: the grid (lattice)
        :param domnet: the dominance network
        :param n_empty: count of how many cells are empty in g
        :param p_mu: mu
        :param p_sigma: sigma
        :param p_epsilon: epsilon
        :param periodic_boundary: are boundary conditions periodic or no-flux?
        :param nbr_offsets: the neeighborhood function
        :param max_offset_i: maximum osset index
        :param vrbs: verbosity
        :return:
            g: updated grid
            n_empty: updated n_empty count
        """
        
        sl = g.shape[0]
        # randomly select an individual
        N = sl * sl
        i_rand = np.random.randint(0, N-1)
        i_row = int(math.floor(i_rand/sl))
        i_col = i_rand - (i_row * sl)
        
        # randomly select a neighbour to be the opponent of the individual
        o_row, o_col = rand_nbr(sl, i_row, i_col, nbr_offsets, max_offset_i, periodic_boundary)

        # initially assume that the opponent will lose
        loser_row = o_row
        loser_col = o_col
        is_loser = False
        v_str = ''

        sum_prs = p_mu + p_sigma + p_epsilon
        mu_d_sum_prs = p_mu / sum_prs
        sigma_d_sum_prs = p_sigma / sum_prs
        
        rv_action = np.random.random()
        
        # do they compete?
        if rv_action <= mu_d_sum_prs:
            i_species = int(g[i_row, i_col])
            o_species = int(g[o_row, o_col])
            if i_species != empty_cell and o_species != empty_cell:
                # v_str += 'i=g[{}][{}]={} o=g[{}][{}]={} '.format(i_row, i_col, i_species, o_row, o_col, o_species)
                i_dominates = int(domnet[i_species, o_species])
                o_dominates = int(domnet[o_species, i_species])
                if i_dominates == 1:
                    # the individual wins
                    v_str += 'i win; '
                    is_loser = True
                    loser_row, loser_col = o_row, o_col
                    winner = i_species
                elif o_dominates == 1:
                    # the opponent wins
                    v_str += 'o win; '
                    is_loser = True
                    loser_row, loser_col = i_row, i_col
                    winner = o_species
                else:
                    # no change
                    v_str += 'no change; '
                    is_loser = False
                    winner = None
                if is_loser:
                    # delete the loser
                    g[loser_row, loser_col] = empty_cell
                    n_empty += 1
            else:
                v_str += 'noop: one or both cells are empty'
                pass

        # reproduce?
        elif mu_d_sum_prs < rv_action and rv_action <= mu_d_sum_prs + sigma_d_sum_prs:
            i_species = int(g[i_row, i_col])
            o_species = int(g[o_row, o_col])
            if i_species == empty_cell and o_species != empty_cell:
                g[i_row, i_col] = o_species
                v_str += 'o reproduced into [{},{}]'.format(i_row, i_col)
                n_empty += -1
            if o_species == empty_cell and i_species != empty_cell:
                g[o_row, o_col] = i_species
                v_str += 'i reproduced into [{},{}]'.format(o_row, o_col)
                n_empty += -1
            if i_species != empty_cell and o_species != empty_cell:
                v_str += 'noop: neither cell is empty'
                pass
            if i_species == empty_cell and o_species == empty_cell:
                v_str += 'noop: both cells are empty'
                pass

        # move?
        else:
            # swap positions
            g[i_row, i_col], g[o_row, o_col] = g[o_row, o_col], g[i_row, i_col]
            v_str += 'moved.'
            
        if vrbs:
            print(v_str)
            
        return g, n_empty


    def revised_elementary_step(g, domnet, p_mu, p_sigma, p_epsilon, periodic_boundary, nbr_offsets, max_offset_i, vrbs):
        """
        The revised elementary step, introduced in Cliff (2024)
        :param g: the grid (lattice)
        :param domnet: the dominance network
        :param p_mu: mu
        :param p_sigma: sigma
        :param p_epsilon: epsilon
        :param periodic_boundary: are boundary conditions periodic or no-flux?
        :param nbr_offsets: the neeighborhood function
        :param max_offset_i: maximum osset index
        :param vrbs: verbosity
        :return:
            g: updated grid
        """

        sl = g.shape[0]
        
        # randomly select an individual
        N = sl * sl
        i_rand = np.random.randint(0, N-1)
        i_row = int(math.floor(i_rand/sl))
        i_col = i_rand - (i_row * sl)
        
        # randomly select a neighbour to be the opponent of the individual
        o_row, o_col = rand_nbr(sl, i_row, i_col, nbr_offsets, max_offset_i, periodic_boundary)

        # initially assume that the opponent will lose
        loser_row = o_row
        loser_col = o_col
        is_loser = False
        v_str = ''
        
        # do they compete?
        rv = 0
        if p_mu < 1:
            # for efficiency: only call random() if we really need to
            rv = np.random.random()
        if p_mu == 1 or rv < p_mu:
            i_species = int(g[i_row, i_col])
            o_species = int(g[o_row, o_col])
            # v_str = 'i=g[{}][{}]={} o=g[{}][{}]={} '.format(i_row, i_col, i_species, o_row, o_col, o_species)
            i_dominates = int(domnet[i_species, o_species])
            o_dominates = int(domnet[o_species, i_species])
            if i_dominates == 1:
                # the individual wins
                # v_str += 'i win; '
                is_loser = True
                loser_row, loser_col = o_row, o_col
                winner = i_species
            elif o_dominates == 1:
                # the opponent wins
                # v_str += 'o win; '
                is_loser = True
                loser_row, loser_col = i_row, i_col
                winner = o_species
            else:
                # no change
                # v_str += 'no change; '
                is_loser = False
                winner = None

            if is_loser:
                # someone lost, delete them and maybe also replace them with the winner
                g[loser_row, loser_col] = empty_cell
                rv = 0
                if p_sigma < 1.0:
                    rv = np.random.random()
                if p_sigma == 1.0 or rv < p_sigma:
                    g[loser_row, loser_col] = winner
                    # v_str += 'g[{}][{}]={}. '.format(loser_row, loser_col, winner)
                else:
                    # v_str += '. '
                    pass

        # end of compete/replace
        # does this individual move?
        rv = 0
        if p_epsilon < 1.0:
            rv = np.random.random()
        if p_epsilon == 1.0 or rv < p_epsilon:
            # swap positions
            g[i_row, i_col], g[o_row, o_col] = g[o_row, o_col], g[i_row, i_col]
            # v_str += 'Moved.'

        if vrbs:
            print(v_str)

        return g


    # pre-flight sanity checks
    if n_species % 2 == 0 or n_species < 3:
        # n_species is an even number or <3
        sys.exit('FAIL: n_species=%d but should be an odd number greater than or equal to 3' % n_species)
    for offset in cx_offsets:
        if offset >= n_species or offset < 1 or type(offset) is not int:
            sys.exit('FAIL: all cx_offset values must be integers between 1 and (n_species-1)')
    if len(cx_offsets) > len(set(cx_offsets)):
        # duplicated values in cx_offsets
        sys.exit('FAIL: all cx_offset values must be unique')
    if M_min >= M_max:
        sys.exit('FAIL: M_min must be less than M_max')
        
    # folder creation
    if not os.path.isdir(exp_id):
        os.mkdir(exp_id)
        
    # create the Numpy arrays for the (x,y) offsets for the nbrhood function
    if nbrhood_count == 4:
        # von Neumann nbrhood
        nbr_offsets = np.array([[-1, 0], [+1, 0], [0, -1], [0, +1]])
        max_offset_i = 3
    elif nbrhood_count == 8:
        # Moore nbrhood
        nbr_offsets = np.array([[-1, -1], [-1, 0], [-1, +1], [0, -1], [0, +1], [+1, -1], [+1, 0], [+1, +1]])
        max_offset_i = 7
    else:
        sys.exit('FAIL: bad nbrhood_count value={}'.format(nbrhood_count))

    # main IID experiment loop
    for exp in range(n_expts):

        # generate the set of M values that we'll run experiments for
        M_vals_list = generate_m_vals(M_min, M_max, ns_oom, True, True)
        M_vals = np.array(M_vals_list)
        
        # create the full (un-ablated) dominance network, used repeatedly for this M-sweep experiment
        full_dom_net = dominance_network(n_species, cx_offsets, verbose)
        # ablate it
        ablated_dom_net = np.copy(full_dom_net)
        ablated_dom_net = dnet_ablate(n_ablations, ablated_dom_net, verbose)

        # draw network diagram, but only if number of nodes is 30 or less (otherwise diagram is too crowded)
        net_fname_str = 'Net_exp%04d' % exp
        if exp == 0 and verbose and n_species < 30 and draw_networks:
            draw_net(ablated_dom_net, net_fname_str)

        # populate the grid with species, used repeatedly for this M-sweep experiment
        start_grid = populate_grid(L, n_species, False)
        n_empty = 0
        
        if use_OES:
            m_max = n_mcs
        else:
            m_max = 1.0/(2*L*L)
        
        for m_index in range(M_vals.shape[0]):

            M = M_vals[m_index]
            
            if M > m_max:
                continue

            # copy the starting grid arrangement into the grid for this experiment
            grid = np.copy(start_grid)

            fname_str = '%s/NA%d_MCS%.2e_N%.2e_NS%d_M%.2e_%04d' % (exp_id, n_ablations, n_mcs, N, n_species, M, exp)
            outfile = open(fname_str + '_densities.csv', 'w')
            # write header line for CSV file
            header_str = 'MCS,n_species,n_empty,'
            for species in range(n_species):
                header_str += "F(n_s(t)=%d)" % species
                if species < (n_species - 1):
                    header_str += ","
            outfile.write(header_str + '\n')

            if write_grids:
                gridfile = open(fname_str + '_grids.csv', 'w')

            if verbose:
                print('Writing lattice grids to %s_grids.csv' % fname_str)

            epsilon = 2 * M * N

            if verbose or True:
                print('M=%e, epsilon=%f' % (M, epsilon))

            if write_domnets:
                # record the full and ablated dominance networks
                # each as a binary adjacency matrix
                domnetfile = open(fname_str + '_domnets.csv', 'w')
                if verbose:
                    print('Writing dominance networks to %s_domnets.csv' % fname_str)
                write_dnet_adjmat(full_dom_net, domnetfile)
                write_dnet_adjmat(ablated_dom_net, domnetfile)
                domnetfile.close()

            # run the specified number of Monte Carlo Steps
            dom_net = np.copy(full_dom_net)

            for mcstep in range(int(n_mcs)):
                # one Monte Carlo Step (MCS)
                
                if mcstep == ablation_mcs:
                    dom_net = np.copy(ablated_dom_net)
                
                if verbose:
                    for r in range(L):
                        print('MCS=,%d, grid row %5d=, %s' % (mcstep, r, grid[r]))

                # count how many individuals of each species (density) and how many nonzero species (s_count)
                density, s_count = species_density(n_species, grid, False)
                    
                # write a line to the output file
                outstr = ('{},{},{},'.format(mcstep + 1, s_count, n_empty))
                for d in range(density.shape[0]):
                    outstr += '{:.4e}'.format(float(density[d]))
                    if d < density.shape[0] - 1:
                        outstr += ','
                outstr += '\n'
                outfile.write(outstr)

                if write_grids and mcstep % g_write_freq == 0:
                    write_grid(mcstep, grid, gridfile)

                if verbose:
                    print(outstr)

                # one Monte Carlo Step is one "generation" of "elementary steps"
                # NB random choice of individual to update on each timestep is not the same as shuffling the population
                for e_step in range(n_generations):
                    if use_OES:
                        grid, n_empty = original_elementary_step(grid, dom_net, n_empty, mu, sigma, epsilon, bc_period,
                                                                 nbr_offsets, max_offset_i, verbose)
                    else:
                        revised_elementary_step(grid, dom_net, mu, sigma, epsilon, bc_period, nbr_offsets,
                                                max_offset_i, verbose)

            # end of Monte Carlo Step loop
            outfile.close()
            del grid

        # end of M-sweep loop
    # end of IID experiments loop
