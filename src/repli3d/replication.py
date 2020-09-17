from scipy.spatial.distance import cdist
import os
import numpy as np


class replicator():
    """
    class that deals with bonds and merging to replicate the system
    """

    def __init__(self, free_number, pori_number, pmin, pmax, system,
                 attached=True, check=True, verbose=False):
        self.num = free_number  # bead number of diffusing element
        self.left = pori_number  # bead number where it is attached
        self.right = pori_number  # bead number where it is attached
        self.attached = attached  # Are the two forks attached
        self.pmin = pmin        # Bead min of the polymer
        self.pmax = pmax        # Bead max of the polymer
        self.to_merge = None  # in case the same bead is replicated store the fork here

        self.check = check
        self.verbose = verbose

        self.bond_dna_dna_tag = None
        self.one_free = False  # is one extrmity free

        self.cl = None  # connection replicon template
        self.cr = None
        self.internal_bonds = []

        # create bond firing DNA
        self.bond_firing_dna_tag = system.bonds.add("DNA", self.num, self.left)
        self.which = "l"
        self.replicated = None  # to store the last bead replicated

        self.replicon = []

    def merge(self, replicator, system):

        verbose1 = self.verbose
        if len(replicator.replicon) == 0:
            system.bonds.remove(replicator.bond_firing_dna_tag)
            return replicator.num

        if verbose1:
            print("merging")
            print(self)
            print(replicator)

        # first clean replicator:
        if replicator.cl is not None:
            system.bonds.remove(replicator.cl)

        if replicator.bond_dna_dna_tag is not None:
            system.bonds.remove(replicator.bond_dna_dna_tag)

        if self.bond_dna_dna_tag is not None:
            system.bonds.remove(self.bond_dna_dna_tag)
            self.bond_dna_dna_tag = None

        self.unweak("r", system)
        replicator.unweak("l", system)

        # add bond between replicon:
        if self.replicon != [] and replicator.replicon != []:
            self.internal_bonds.append(system.bonds.add(
                "DNA", self.replicon[-1], replicator.replicon[0]))

        # transform weak particle left and right:

        self.replicon += replicator.replicon
        self.internal_bonds += replicator.internal_bonds
        if self.cr is not None:
            system.bonds.remove(self.cr)
        self.cr = replicator.cr
        # Should check if None
        self.right = replicator.right

        if verbose1:
            print("one_f", self.one_free, replicator.one_free)

        self.one_free = self.one_free | replicator.one_free

        if self.attached and (not self.one_free):
            if verbose1:
                print("adding bond")
            self.bond_dna_dna_tag = system.bonds.add("DNA", self.left, self.right)
        else:
            if verbose1:
                print("not adding",)
        # then release diffusing:,

        system.bonds.remove(replicator.bond_firing_dna_tag)

        if self.check:
            self.checks(system)

        return replicator.num

    def unweak(self, side, system):
        # transform weak into strong before merging
        if self.verbose:
            print("Unweak")
        if len(self.replicon) > 1:
            if side == "r":
                p = system.particles.get(self.replicon[-1])
            else:
                p = system.particles.get(self.replicon[0])
            assert p.type == "fDNA"
            p.type = "rDNA"

        if len(self.replicon) >= 2:
            if side == "r":
                system.bonds.remove(self.internal_bonds[-1])
                self.internal_bonds[-1] = system.bonds.add("DNA",
                                                           self.replicon[-1], self.replicon[-2])
            else:
                system.bonds.remove(self.internal_bonds[0])
                self.internal_bonds[0] = system.bonds.add("DNA", self.replicon[0], self.replicon[1])

    def __repr__(self):
        r = "l %i, size %i, right %i\n" % (self.left, len(self.replicon), self.right)
        r += "Onef %s %s" % (str(self.one_free), str(self.bond_dna_dna_tag))
        return r

    def checks(self, system):
        if len(self.replicon) > 0:
            # Check right and left are freshDNA
            p = system.particles.get(self.replicon[0])
            try:
                assert p.type == "fDNA"
            except:
                print(self)
                print("Error fresh", p.type, self.replicon[0])
                raise

            p = system.particles.get(self.replicon[-1])
            try:
                assert p.type == "fDNA"
            except:
                print(self)
                print("Error fresh", p.type)
                raise
        if len(self.replicon) > 2:
            for p in self.replicon[1:-1]:
                pi = system.particles.get(p)
                try:
                    assert(pi.type == "rDNA")
                except:
                    print("Problem freh DNA in the middle")
                    print(self)
                    print(p, pi.type)
                    raise

        if self.bond_dna_dna_tag is not None:
            # Check that right and left are attached
            b = system.bonds.get(self.bond_dna_dna_tag)
            try:
                assert((b.a, b.b) == (self.left, self.right)
                       or (b.a, b.b) == (self.right, self.left))
            except:
                print(self.bond_dna_dna_tag)
                print(self)
                print(b)
                raise

    def propagate(self, free_bead, system):
        # possible verification

        self.replicated = None
        if self.check:
            self.checks(system)

        verbose1 = self.verbose
        if verbose1:
            print("Propagate")
            print(self)
        # one bead by one bead
        if self.which == "l":
            if self.left == self.pmin:
                self.which = "r"
                # Detach template from new:
                if self.cl is not None:
                    system.bonds.remove(self.cl)
                    self.cl = None
                if self.attached and self.bond_dna_dna_tag is not None:
                    system.bonds.remove(self.bond_dna_dna_tag)
                    self.bond_dna_dna_tag = None

                self.one_free = True

                return
            self.left -= 1
            self.replicated = self.left

            self.replicon.insert(0, free_bead)
        else:
            if self.right == self.pmax:
                self.which = "l"
                if self.cr is not None:
                    if verbose1:
                        print("removing right")
                    system.bonds.remove(self.cr)
                    self.cr = None
                if self.attached and self.bond_dna_dna_tag is not None:
                    system.bonds.remove(self.bond_dna_dna_tag)
                    if verbose1:
                        print("removing bond")
                    self.bond_dna_dna_tag = None
                if verbose1:
                    print("Free")
                self.one_free = True
                return

            self.right += 1
            self.replicated = self.right
            self.replicon.append(free_bead)

        #print("Propagate replicon %i size,l %i r %i" %(len(self.replicon),self.left,self.right))

        # change tag of particle repliacted
        #print("Need to change tag")
        p = system.particles.get(free_bead)
        # print(p)
        p.type = "fDNA"
        #system.particles[free_bead].typeid  = syst["plist"].index("DNA")

        system.bonds.remove(self.bond_firing_dna_tag)
        self.bond_firing_dna_tag = system.bonds.add("DNA", self.num, self.left)

        if self.attached and not self.one_free:
            if self.bond_dna_dna_tag is not None:
                system.bonds.remove(self.bond_dna_dna_tag)
                # print("DNA",self.left,self.right)

            self.bond_dna_dna_tag = system.bonds.add("DNA", self.left, self.right)
            # for i in range(10,len(system.bonds)):
            #    print(system.bonds[i])

        # deal with replicon

        if len(self.replicon) == 1:
            self.cl = system.bonds.add("DNA", self.replicon[0], self.left)
            p = system.particles.get(self.left)
            #assert p.type == "DNA"
            p.type = "tDNA"
            self.cr = system.bonds.add("DNA", self.replicon[0], self.right)
            p = system.particles.get(self.right)
            #assert p.type == "DNA"
            p.type = "tDNA"

        else:
            if self.which == "l":
                system.bonds.remove(self.cl)
                self.cl = system.bonds.add("weak", self.replicon[0], self.left)

                self.internal_bonds.insert(0, system.bonds.add(
                    "weak", self.replicon[0], self.replicon[1]))
                if len(self.internal_bonds) > 2:
                    system.bonds.remove(self.internal_bonds[1])
                    self.internal_bonds[1] = system.bonds.add(
                        "DNA", self.replicon[1], self.replicon[2])
                    p = system.particles.get(self.replicon[1])
                    p.type = "rDNA"
                if len(self.internal_bonds) == 2:
                    p = system.particles.get(self.replicon[1])
                    p.type = "rDNA"

                p = system.particles.get(self.left)
                #assert p.type == "DNA"
                p.type = "tDNA"

                # reinforce right
            if self.which == "r":
                system.bonds.remove(self.cr)
                self.cr = system.bonds.add("weak", self.replicon[-1], self.right)
                self.internal_bonds.append(system.bonds.add(
                    "weak", self.replicon[-2], self.replicon[-1]))
                if len(self.internal_bonds) > 2:
                    system.bonds.remove(self.internal_bonds[-2])
                    self.internal_bonds[-2] = system.bonds.add("DNA",
                                                               self.replicon[-2], self.replicon[-3])
                    p = system.particles.get(self.replicon[-2])
                    p.type = "rDNA"

                if len(self.internal_bonds) == 2:
                    p = system.particles.get(self.replicon[1])
                    p.type = "rDNA"

                p = system.particles.get(self.right)
                #assert p.type == "DNA"
                p.type = "tDNA"

        if self.which == "l":
            self.which = "r"
        else:
            self.which = "l"

        print(self)
        # for r in self.replicon:
        #    print(r,system.particles.get(r).type)

    def to_extend(self):
        if self.replicon == []:
            return self.left
        if self.which == "l":
            return self.left
        else:
            return self.right


class chromosome():
    def __init__(self, start, end, l_ori, verbose=False, verbose_replicon=False, attached=True):
        self.start = start
        self.end = end
        self.l_replicator = []
        self.l_ori = l_ori
        self.verbose = verbose
        self.verbose_replicon = verbose_replicon
        self.attached = attached
        self.rfd = np.zeros(self.end+1-self.start)
        for ori in self.l_ori:
            assert(self.start <= ori <= self.end)

    def add_forks(self, free_number, pori_number, system, time):

        assert(pori_number in self.l_ori)

        replicon = replicator(free_number, pori_number, self.start,
                              self.end, system, verbose=self.verbose_replicon,
                              attached=self.attached)
        # Position the replicator in a ordered way:
        if self.l_replicator == []:
            self.l_replicator.append(replicon)

        elif len(self.l_replicator) >= 2:
            found = False
            for ip, (r1, r2) in enumerate(zip(self.l_replicator[:-1], self.l_replicator[1:])):
                if r1.right < pori_number and r2.left > pori_number:
                    found = True
                    self.l_replicator.insert(ip+1, replicon)
            if not found:
                if self.l_replicator[0].left > pori_number:
                    self.l_replicator.insert(0, replicon)
                else:
                    self.l_replicator.append(replicon)
        else:
            if self.l_replicator[0].left > pori_number:
                self.l_replicator.insert(0, replicon)
            else:
                self.l_replicator.append(replicon)

        self.rfd[pori_number] = np.random.choice([time, -time])
        self.l_ori.remove(pori_number)

        self.checks()

    def checks(self):
        if len(self.l_replicator) >= 2:
            for ip, (r1, r2) in enumerate(zip(self.l_replicator[:-1], self.l_replicator[1:])):
                try:
                    assert(r1.right <= r2.left)
                except:
                    for repl in self.l_replicator:
                        print(repl.left, repl.right)
                    raise

    def propagate(self, system, g_unrep, time):

        print(self.rfd)
        verbose = self.verbose
        if verbose:
            print("Propagating")

        # Building list to extend,
        # Not replicating same or adjacent DNA
        to_extend = []
        pop = []
        for ir, r in enumerate(self.l_replicator):
            e = r.to_extend()
            if e is not None:
                if verbose:
                    print("extension", e)
                if e not in to_extend and e-1 not in to_extend:
                    to_extend.append(e)
                else:
                    self.l_replicator[ir-1].to_merge = r
                    pop.append(r)
        if verbose:
            print("list extension", to_extend)

        # remove cross replicator
        for r in pop:
            self.l_replicator.remove(r)

        remove = []
        for e in to_extend:
            if e in self.l_ori:
                remove.append(e)
        if remove != []:
            if verbose:
                print("Passivation")
            for r in remove:
                self.l_ori.remove(r)

        if verbose:
            print("N replicon", len(self.l_replicator), [ri.left for ri in self.l_replicator])

        to_extend_p = np.array([system.particles[i].position for i in to_extend])
        # get list of unreplicated DNA
        g_unrep.force_update()
        unrep = np.array([p.position for p in g_unrep])
        tags_unrep = [p.tag for p in g_unrep]

        if len(to_extend_p) > 0 and len(unrep) > 0:
            D2 = cdist(to_extend_p, unrep)
            selected = np.argmin(D2, axis=1)
            if verbose:
                print(selected)
            if len(set(selected)) == len(selected):
                for r, s in zip(self.l_replicator, selected):
                    if verbose:
                        print("Propagate", r.left)

                    r.propagate(tags_unrep[s], system)
                    # Extend rfd
                    if r.replicated is not None:
                        if self.rfd[r.replicated] == 0:
                            if r.which == "l":
                                self.rfd[r.replicated] = time  # Inversed as already swhiched
                            else:
                                self.rfd[r.replicated] = -time
                        else:
                            print("Reassiging RFD", r.replicated, r.which)

                    to_extend.pop(0)
            else:
                # in case of overlap between available DNA
                propagated = np.zeros_like(to_extend, dtype=np.bool)
                while np.sum(propagated) != len(propagated):
                    to_extend_p = np.array([system.particles[i].position for i in to_extend])
                    # get list of unreplicated DNA
                    g_unrep.force_update()
                    unrep = np.array([p.position for p in g_unrep])
                    tags_unrep = [p.tag for p in g_unrep]
                    D2 = cdist(to_extend_p, unrep)
                    selected = np.argmin(D2, axis=1)
                    used = []
                    for ir, (r, s, p) in enumerate(zip(self.l_replicator, selected, propagated)):
                        if (not p) and (s not in used):
                            if verbose:
                                print("propagate", r.left)
                            used.append(s)

                            r.propagate(tags_unrep[s], system)
                            # Extend rfd
                            if r.replicated is not None:
                                if self.rfd[r.replicated] == 0:
                                    if r.which == "l":
                                        # Inversed as already swhiched
                                        self.rfd[r.replicated] = time
                                    else:
                                        self.rfd[r.replicated] = -time
                                else:
                                    print("Reassiging RFD", r.replicated, r.which)

                            propagated[ir] = True

        Frees = []
        for r in self.l_replicator:
            if r.to_merge != None:

                free = r.merge(r.to_merge, system)
                Frees.append(free)
                r.to_merge = None

        # Check for passivation again

        remove = []
        for r in self.l_replicator:
            if r.left in self.l_ori:
                remove.append(r.left)
            if r.right in self.l_ori:
                remove.append(r.right)
        remove = list(set(remove))
        if remove != []:
            if verbose:
                print("Passivation")
            for r in remove:
                self.l_ori.remove(r)

        # check if two replicons are attached to the same bead
        merged = True
        while merged:
            merged = False
            if len(self.l_replicator) >= 2:
                for index, (r1, r2) in enumerate(zip(self.l_replicator[:-1], self.l_replicator[1:])):
                    if r1.right == r2.left:
                        # l_replicator.remove(r1)
                        if verbose:
                            print("merging", r1.right, r2.left)
                        self.l_replicator.remove(r2)
                        free = self.l_replicator[index].merge(r2, system)
                        Frees.append(free)
                        merged = True
                        break
        return Frees
