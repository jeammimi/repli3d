
import argparse
import hoomd
from hoomd import data, init, md, group, dump, deprecated, analyze, comm
from scipy.spatial.distance import cdist
import os
import numpy as np
import pandas as pd
#from yeast_nucleosomal_resolution.tools import fl
#from hyperopt import fmin, tpe, hp, STATUS_OK


from repli3d.replication import chromosome


def initialize_snap(syst, params={}):

    syst["np"] = sum(syst["len_polymers"]) * (1+syst["excess_DNA"]) + syst["ndiff"]

    syst["nbond"] = sum(syst["len_polymers"])-len(syst["len_polymers"])

    if syst["cohesin_real"]:
        syst["np"] += syst["n_cohesin"] * syst["particule_per_cohesin"]
        syst["nbond"] += syst["n_cohesin"] * syst["particule_per_cohesin"]

    syst["Rf"] = (sum(syst["len_polymers"])*0.5**3/syst['density'])**0.33
    print("Radius of the cell", syst["Rf"])
    R = syst["Rf"]+1
    snapshot = data.make_snapshot(N=syst["np"], box=data.boxdim(L=2 * R), bond_types=['polymer'])

    syst["bond_list"] = ['DNA', 'cohesin', "fFactor", "weak"]

    syst["plist"] = ["DNA", "uDNA", "pDNA", "fFactor", 'cohesine', "rDNA", "fDNA", "tDNA"]
    # uDNA unreplicated DNA
    # rDNA replicated DNA
    # fDNA freshly replicated DNA
    # tDNA template replicated DNA

    return snapshot, syst


def create_snapshot(snapshot, syst, seed=False):

    snapshot.bonds.types = syst["bond_list"]
    snapshot.bonds.resize(syst["nbond"])
    snapshot.particles.types = syst["plist"]

    ##################################################
    # Define particle type and positions
    ref_particle = 0
    ref_bond = 0
    for iX_len, X_len in enumerate(syst["len_polymers"]):

        # Particle type
        for _ in range(X_len):
            snapshot.particles.typeid[ref_particle] = syst["plist"].index("DNA")
            ref_particle += 1

        # Particle position
        maxi = 2*syst["Rf"]
        if seed:
            np.random.seed(0)
        while maxi > syst["Rf"]:
            p0 = np.array([0.0, 0, 0])
            for i in range(X_len):
                snapshot.particles.position[i] = p0
                p0 += 0.2*(1-2*np.random.rand(3))

            Cm = np.mean(snapshot.particles.position, axis=0)
            # print(Cm)
            for i in range(X_len):
                snapshot.particles.position[i] -= Cm[:]
            maxi = np.max(np.abs(snapshot.particles.position))
            print("Maximum size", maxi, sum(syst["len_polymers"])**0.5)

        # Define bonds
        start = sum(syst["len_polymers"][:iX_len])
        for i in range(X_len-1):
            snapshot.bonds.group[ref_bond] = [start+i, start+i+1]
            snapshot.bonds.typeid[ref_bond] = syst["bond_list"].index('DNA')  # polymer_A
            ref_bond += 1

    ##################################################

    for i in range(syst["excess_DNA"]*sum(syst["len_polymers"])):
        snapshot.particles.typeid[ref_particle] = syst["plist"].index("uDNA")
        ref_particle += 1

    for i in range(syst["ndiff"]):
        snapshot.particles.typeid[ref_particle] = syst["plist"].index("fFactor")
        snapshot.particles.position[ref_particle] = random_sampling_sphere(syst["Rf"]*0.95)
        ref_particle += 1

    if syst["cohesin_real"]:
        for i in range(2*syst["n_cohesin"] * syst["particule_per_cohesin"]):
            snapshot.particles.typeid[ref_particle] = syst["plist"].index("cohesine")
            ref_particle += 1

    if syst["cohesin_real"]:
        print("Todo")

    return snapshot


def random_sampling_sphere(R):
    phi = np.random.rand()*2*np.pi
    costheta = 2*np.random.rand() - 1
    u = np.random.rand()

    theta = np.arccos(costheta)
    r = R * u**0.33
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x, y, z])


def simulate(syst, n_steps, data_folder="./repli", params={}, seed=False):
    import time as Time
    global t0
    t0 = Time.time()

    def time(where):
        global t0
        print(where, "elapsed %.1f" % (Time.time()-t0))
        t0 = Time.time()

    stretch = syst["stretch"]
    verbose = syst["verbose"]

    data_folder = os.path.join(data_folder)

    os.makedirs(data_folder, exist_ok=True)

    print(data_folder)

    time("Start")
    snapshot, syst = initialize_snap(syst)
    time("Initialize")
    length_steps = syst["length_steps"]  # 50000

    if comm.get_rank() == 0:
        snapshot = create_snapshot(snapshot, syst, seed=seed)

    snapshot.broadcast()

    system = init.read_snapshot(snapshot)

    bond = md.bond.harmonic(name="mybond")
    bond.bond_coeff.set(syst["bond_list"], k=100.0, r0=0.84)
    bond.bond_coeff.set("weak", k=10.0, r0=0.84)

    nl = md.nlist.cell()
    #nl = md.nlist.tree()

    sc=0.5
    gauss = md.pair.gauss(r_cut=3.0*sc, nlist=nl)
    gauss.pair_coeff.set(syst["plist"], syst["plist"], epsilon=4.0, sigma=1.0 * sc)
    gauss.pair_coeff.set("fDNA", syst["plist"], epsilon=0.5, sigma=.5 * sc)
    gauss.pair_coeff.set("uDNA", syst["plist"], epsilon=0, sigma=1.0 * sc)

    ##################################################
    # wall
    sphere = md.wall.group()
    r_extrap = 0.5
    r0 = 0.5

    sphere.add_sphere(r=syst["Rf"], origin=(0.0, 0.0, 0.0), inside=True)
    wall_force_slj = md.wall.lj(sphere, r_cut=1.12)
    wall_force_slj.force_coeff.set(syst["plist"], epsilon=1, sigma=0.5 + r0 * 1.0,
                                   r_cut=1.12 * (0.5 + r0 * 1.0), mode="shift", r_extrap=0.5 + r0 * r_extrap)

    all = group.all()
    period = length_steps

    if stretch:
        period = 1000
    gsd = dump.gsd(group=all, filename=os.path.join(data_folder, 'poly.gsd'),
                   period=period, overwrite=True, dynamic=["attribute", "topology"], phase=0)

    ##################################################
    # Run the simulation

    sim_dt = 0.01

    snp = system
    md.integrate.mode_standard(dt=sim_dt)
    if seed:
        seed = 0
    else:
        seed = np.random.randint(10000)

    method = md.integrate.langevin(group=all, kT=1, seed=seed, dscale=False)


    group_hic = all  # group.tags(name="hic", tag_min=0, tag_max=Nparticule)

    r_hic = []
    Oml = []
    Nel = []

    cdists = []

    print(len(snapshot.particles.typeid))
    Free_firing_factor = [i for i in range(
        syst["np"]) if snapshot.particles.typeid[i] == syst["plist"].index("fFactor")]
    Unrep = [i for i in range(syst["np"]) if snapshot.particles.typeid[i]
             == syst["plist"].index("uDNA")]

    g_firing_factor = group.tag_list("free", Free_firing_factor)
    g_unrep = group.type("uDNA", update=True)

    if stretch:
        md.force.constant(fx=-2.0, fy=0, fz=0, group=group.tag_list("0", [0]))
        md.force.constant(fx=2.0, fy=0, fz=0, group=group.tag_list(
            "0", [syst["len_polymers"][0]-1]))

    print("Firing factors", Free_firing_factor)
    #hoomd.run(length_steps*10, profile=False,quiet=True)

    #from repli3d.replication import replicator
    global iname
    iname = 0

    def gname(name):
        global iname
        iname += 1
        return name+str(iname)
    time("End define all")
    l_ch = []

    for i, X_len in enumerate(syst["len_polymers"]):
        start = sum(syst["len_polymers"][:i])
        rand = list(set(np.random.choice(
            range(start, X_len+start), int(syst["ori_density"]*X_len))))
        if "origin_position" in syst:
            Potential_ori = syst["origin_position"]
        else:
            Potential_ori = rand
        print("Potential ori", Potential_ori)
        l_ch.append(chromosome(start, start+X_len-1,
                               Potential_ori, attached=syst["attached"], verbose=syst["verbose"],
                               verbose_replicon=syst["verbose"]))

    hoomd.run(syst["equi"], profile=False, quiet=True)
    time("Start loop")
    for i in range(n_steps):
        time("Start run")
        hoomd.run(length_steps, profile=False, quiet=True)
        #snapshot = system.take_snapshot(all=True)
        # snapshot.broadcast()
        time("End run")

        p_ori_tag = []
        for X in l_ch:
            p_ori_tag.extend(X.l_ori)

        p_ori = np.array([system.particles.get(ptag).position for ptag in p_ori_tag])

        free = np.array([p.position for p in g_firing_factor])
        free_tag = [p.tag for p in g_firing_factor]

        if len(p_ori) > 0 and len(free) > 0:
            D1 = cdist(p_ori, free)
            if stretch:
                d = 4
            else:
                d = 2
            D1[D1 > d] = 0
            used_ori = []
            used_free = []
            for ifree, free_number in enumerate(free_tag):
                for ipori, pori_number in enumerate(p_ori_tag):

                    # insert the replicon correctly in the list
                    if D1[ipori, ifree] > 1e-7 and pori_number not in used_ori and free_number not in used_free:

                        if verbose:
                            print("Activate", pori_number)
                        for X in l_ch:
                            if pori_number in X.l_ori:
                                X.add_forks(free_number, pori_number, system, time=i)
                        # Remove ori
                        # Remove firing factor
                        g_firing_factor = group.difference(
                            gname("delta2"), g_firing_factor, group.tag_list(gname("rand"), [free_number]))

                        used_ori.append(pori_number)
                        used_free.append(free_number)

                        continue
        if verbose:
            for iX, X in enumerate(l_ch):
                print("Chromosome %i" % iX)
                print("List replicator")
                for repl in X.l_replicator:
                    print(repl.left, repl.right)
                print("End list")
        time("Ori association")
        hoomd.run(length_steps, profile=True, quiet=False)
        time("Run length %i" % i)
        # get list of position where to add material

        if i % 15 == 0:
            for X in l_ch:
                Frees = X.propagate(system, g_unrep, time=i)
                g_firing_factor = group.union(
                    gname("delta3"), g_firing_factor, group.tag_list(gname("rand"), Frees))
            time("Propagation ")
            if verbose:
                print("Firing", [p.tag for p in g_firing_factor])

        for iX,X in enumerate(l_ch):
            np.save(data_folder+"/rep_%i.npy" % iX,X.rfd)

"""
def test_16(attached=True, nch=16):
    syst = {}
    syst["len_polymers"] = [200]*nch
    syst["excess_DNA"] = 2
    syst["ndiff"] = 50
    syst["n_cohesin"] = 40
    syst["particule_per_cohesin"] = 10
    syst["density"] = 0.05
    syst["cohesin_real"] = False
    syst["ori_density"] = 20/syst["len_polymers"][0]
    syst["stretch"] = True
    syst["verbose"] = True
    syst["attached"] = attached
    simulate(syst, 400, data_folder="./repli")
"""


def test_one(attached=True, nch=16, args={}):
    syst = {}
    syst["len_polymers"] = [args["lengthp"]]
    syst["excess_DNA"] = 1
    syst["ndiff"] = 50
    syst["n_cohesin"] = 40
    syst["particule_per_cohesin"] = 10
    syst["density"] = args["density"]
    syst["cohesin_real"] = False
    syst["ori_density"] = 0.2
    syst["stretch"] = False
    syst["verbose"] = True
    syst["attached"] = attached
    syst["length_steps"] = 10000
    syst["equi"] = 100000
    Xp = pd.read_csv("data//nn_K562_from_all.csv", sep="\t")
    # skip 5 Mg
    start = int(5000/5)
    d3p = Xp["signalValue"][start:start+syst["len_polymers"][0]]
    d3p /= np.sum(d3p)
    noise = np.ones_like(d3p)*np.sum(d3p)*0.1/len(d3p)
    d3p += noise
    syst["origin_position"] = list(set(np.random.choice(
        range(len(d3p)), p=d3p/np.sum(d3p), size=int(len(d3p)*5/20))))
    syst["origin_position"].sort()
    print(syst["origin_position"])

    simulate(syst, args["nsteps"], data_folder=args["root"])




def test_one_ori(attached=True, nch=1, args={}):
    syst = {}
    syst["len_polymers"] = [5000]
    syst["excess_DNA"] = 1
    syst["ndiff"] = 200
    syst["n_cohesin"] = 40
    syst["particule_per_cohesin"] = 10
    syst["density"] = args["density"]
    syst["cohesin_real"] = False
    syst["ori_density"] = 0.2
    syst["stretch"] = False
    syst["verbose"] = True
    syst["attached"] = attached
    syst["length_steps"] = 10000
    syst["equi"] = 20000
    #Xp = pd.read_csv("data//nn_K562_from_all.csv", sep="\t")
    # skip 5 Mg
    start = int(5000/5)

    syst["origin_position"] = [500]
    syst["origin_position"].sort()
    print(syst["origin_position"])

    simulate(syst, args["nsteps"], data_folder=args["root"])


    #ch4
    # YDL089w 296820..298274
    # Tetos1 298000
    # YDL088C 298417..300003

    #ARS environ 325

    # YDL055c 355674..356759
    # Lac 359000
    # YDL054c 359825..361285


    # YDL020c 415113..416708
    # Tet os 2 417000
    # YDL019c 417663..421514



if __name__ == "__main__":
    import hoomd

    parser = argparse.ArgumentParser()
    parser.add_argument('--nch', type=int, default=16)
    parser.add_argument('--cpu', action="store_true")
    parser.add_argument('--debug', action="store_true")
    parser.add_argument('--attached', action="store_true")
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--nsteps', type=int, default=2000)
    parser.add_argument('--lengthp', type=int, default=10000)
    parser.add_argument('--density', type=float, default=0.05)

    parser.add_argument('--root', type=str, default="./repli")

    args = parser.parse_args()
    nch = args.nch
    if args.cpu:
        init_cmd = "--mode cpu"
    else:
        init_cmd = "--mode=gpu --gpu=%i " % args.gpu

    if args.debug:
        init_cmd += " --notice-level=10"
    hoomd.context.initialize(init_cmd)
    print(nch)
    my_dict = args.__dict__
    print(my_dict)

    #test_one(attached=False, nch=nch, args=args.__dict__)
    test_one_ori(attached=args.attached, nch=nch, args=args.__dict__)
