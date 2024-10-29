import numpy as np
from constants import *
import matplotlib.pyplot as plt
from collections import deque


def S(x, y, L):
    """Salinity field

    Args:
        x (float): x position
        y (float): y position
        L (int): dimension of study (m)
        sl (float): slope parameter from 0 to 1

    Returns:
        float: salinity field
    """
    xy_inf_L = x + y <= L
    sl_01 = 0 <= sl and sl < 1
    result = 0
    if sl_01 and xy_inf_L:
        a = L*(1-sl)/(x+y)
        b = 1/(1-sl)
        result = max(36/a + 36*(1-b), 0)
    elif sl == 1 and xy_inf_L:
        # 0 or 72
        result = 0
    elif sl_01 and not xy_inf_L:
        a = L*(1-sl)/(x+y)
        b = 1/(1-sl)
        return min(36/a + 36*(1-b), 1)
    return result

def I(x, y, L):
    """inundation field

    Args:
        x (float): x position
        y (float): y position
        L (int): dimension of study (m)
        sl (float): slope parameter from 0 to 1

    Returns:
        float: inundation field
    """
    xy_inf_L = x + y <= L
    sl_01 = 0 <= sl and sl < 1
    result = 0
    if sl_01 and xy_inf_L:
        a = L*(1-sl)/(x+y)
        b = 1/(1-sl)
        result = max(0.5/a + 0.5*(1-b), 0)
    elif sl == 1 and xy_inf_L:
        # 0 or 1
        result = 0
    elif sl_01 and not xy_inf_L:
        a = L*(1-sl)/(x+y)
        b = 1/(1-sl)
        result = min(0.5/a + 0.5*(1-b), 1)
    return result

def sigma(tree, L):
    """salinity response of a given tree

    Args:
        tree (Tree): tree
        L (int): dimension of study (m)

    Returns:
        float: salinity response
    """
    x = tree.x
    y = tree.y
    return 1/(1 + np.exp(ds*(Ui-S(x, y, L))))

def eta(tree, L):
    """inundation response of a given tree

    Args:
        tree (Tree): tree
        L (int): dimension of study (m)

    Returns:
        float: inundation response
    """
    x = tree.x
    y = tree.y
    return 1 - I(x, y, L)

def C(tree):
    """competition factor of a given tree

    Args:
        tree (Tree): tree

    Returns:
        float: competition factor
    """
    FA = tree.FA
    return 1 - 2*FA

def H(dbh):
    """Height in cm

    Args:
        dbh (float): diameter at breast height (cm)

    Returns:
        float: height in (cm)
    """
    return 137 + b2*dbh - b3*dbh**2

# maximum height (cm)
H_max = H(D_max)

def cold_immunity(dbh, gen):
    """cold immunity of a given tree dbh (cm) and it's cold resistance

    Args:
        dbh (float): diameter at breast height (cm)
        gen (float): cold resistance

    Returns:
        bool: True if the tree has cold resistance
    """
    if dbh < sapp_dim:
        return np.random.rand() < gen*0.01
    elif dbh < mature_dim:
        return np.random.rand() < max(gen*0.1, 0.7)
    else:
        return np.random.rand() < max(gen*0.15, 0.95)

def f(tree):
    """Compute tree growth in cm (1 year interval)

    Args:
        tree (Tree): tree

    Returns:
        Tree: tree after growth
    """
    D = tree.dbh

    num = G * D * (1-D*H(D)/D_max/H_max)
    denom = 274 + 3*b2*D - 4*b3*D**2
    solve = num / denom * sigma(tree, L) * C(tree) * eta(tree, L)
    
    tree.add_memory(solve)
    if np.random.rand() < epsilon:
        tree.kill(source="natural")


    if np.mean(tree.memory) < half_D_over_age:
        # the tree dies
        tree.kill(source="negative")
        # nothing happens if memory is positive
    else:
        # tree grows
        tree.dbh = D + solve if solve > 0 else D
        tree.update_info()
    return tree

def random_seeds(x, y, gen, dist=5, dist_min=0, num_points=5):
    """Plant seeds around a tree at position (x, y)

    Args:
        x (float): tree x position
        y (float): tree y position
        gen (float): tree cold resistance
        dist (float, optional): distance around tree where seeds can be planted. Defaults to 5.
        dist_min (float, optional): distance around tree where seeds can't be planted. Defaults to 0.
        num_points (int, optional): number of seeds. Defaults to 5.

    Returns:
        ndarray[Tree]: array of seeds
    """
    distances = np.random.uniform(dist_min, dist, np.random.randint(0, 6))
    angles = np.random.uniform(0, 2 * np.pi, num_points)
    
    coordinates = [(x + d * np.cos(angle), y + d * np.sin(angle)) for d, angle in zip(distances, angles)]
    # return np.array([])
    return np.array([Tree(x_new, y_new, dbh=np.random.uniform(sap_dim, 1.5*sap_dim), gen=gen) for x_new, y_new in coordinates if x_new > 0 and y_new > 0 and x_new < L and y_new < L])

def can_establish(tree, forest):
    """check if a given seed can establish in a forest

    Args:
        tree (Tree): seed
        forest (Forest): forest

    Returns:
        bool: True if seed can establish
    """
    force = sum([FA_kn(tree, tree_n, 2) for tree_n in forest.trees])
    if force/(np.pi * tree.R**2) > 0.5:
        return False
    return True

class Forest:
    """Mangrove forest
    """
    def __init__(self, seed=0, n_cold=2) -> None:
        """

        Args:
            seed (int, optional): seed number of simulation. Defaults to 0.
            n_cold (int, optional): number of cold events per year. Defaults to 2.
        """
        self.trees = np.array([])
        self.AGB = []
        self.BGB = []
        self.seed = seed
        self.n_cold = n_cold
        # np.random.seed(seed)
    
    def total_AGB(self):
        """Above ground biomass of the forest (tons)

        Returns:
            float: above ground biomass (tons)
        """
        result = 0
        for tree in self.trees:
            result += tree.AGB()
        return result

    def total_BGB(self):
        """Below ground biomass of the forest (tons)

        Returns:
            float: below ground biomass (tons)
        """
        result = 0
        for tree in self.trees:
            result += tree.BGB()
        return result

    def cold_event(self, prob_death=0.11):
        """simulating a cold event

        Args:
            prob_death (float, optional): proportion of deaths caused by cold. Defaults to 0.11.
        """
        new_trees = []
        for tree in self.trees:
            if np.random.rand() < prob_death and not cold_immunity(tree.dbh, tree.gen):
                tree.kill(source="cold")
            else:
                if tree.dbh < sap_dim:
                    tree.gen += 0.5
                elif tree.dbh < sapp_dim:
                    tree.gen += 0.1
                else:
                    tree.gen += 1/30
                new_trees.append(tree)
        self.trees = np.array(new_trees)
    

    def update(self, years, nStrides=2):
        """Simulate forest growth over a number of years

        Args:
            years (int): number of years to simulate
            nStrides (int, optional): numerical parameter for competition computation. Can be increased for more precise computation. Defaults to 2.
        """
        self.AGB.append(self.total_AGB())
        self.BGB.append(self.total_BGB())
        self.plot()
        for i in range(1, years+1):
            print("year", i)
            self.update_step(nStrides, self.n_cold)
            if i%5==0:
                self.plot(annee=i)
        self.plot_BIOM()
        for t in self.trees:
            if t.dbh > 2:
                print(t)
            

    
    def initialize(self, nb_trees, gen=1):
        """Plant seeds in the forest

        Args:
            nb_trees (int): number of seeds to plant
            gen (float, optional): cold resistance. Defaults to 1.
        """
        X = []
        Y = []   
        while len(X) < nb_trees:
            # Generate random x and y values in the range [0, L]          
            x = np.random.uniform(0, L)
            y = np.random.uniform(0, L)
            # Check if the point satisfies x + y <= L
            if x + y <= L:
                X.append(x)
                Y.append(y)
        # Convert lists to numpy arrays
        X = np.array(X)
        Y = np.array(Y)

        for k in range(nb_trees):
            P = np.array([Tree(X[k], Y[k], dbh=np.random.uniform(sap_dim, 1.5*sap_dim), gen=gen)])
            # P = np.array([Tree(X[k], Y[k], dbh=1.37)])
            self.trees = np.concatenate((self.trees, P))
        
        print(f"nb alive after ini {self.nb_alive()}")

    def nb_alive(self):
        """number of trees alive in the forest

        Returns:
            int: trees alive
        """
        return self.trees.shape[0]

    def plot(self, seed=True, sap=True, mat=True, annee=0):
        """Plot the simulation

        Args:
            seed (bool, optional): True if seeds are represented. Defaults to True.
            sap (bool, optional): True if saplings are represented. Defaults to True.
            mat (bool, optional): True if mature trees are represented. Defaults to True.
            annee (int, optional): moment (year) of the simulation. Defaults to 0.
        """
        # Create an array of x values
        x = np.linspace(0, L, 100)
        # Calculate the corresponding y values
        y = L - x
        # Plot the area where x + y <= L
        plt.fill_between(x, y, color='lightgreen', alpha=0.7)
        plt.fill_between(x, y, L, color='lightblue', alpha=0.7)
        if seed:
            X_seed = [tree.x for tree in self.trees if tree.dbh < sapp_dim ]
            Y_seed = [tree.y for tree in self.trees if tree.dbh < sapp_dim ]
            Z_seed = [tree.get_crown_area()*2 for tree in self.trees if tree.dbh < sapp_dim ]
            plt.scatter(X_seed, Y_seed, Z_seed, c='purple', label=f"Seed ({len(X_seed)})")
        if sap:
            X_sapp = [tree.x for tree in self.trees if tree.dbh >= sapp_dim and tree.dbh < mature_dim]
            Y_sapp = [tree.y for tree in self.trees if tree.dbh >= sapp_dim and tree.dbh < mature_dim]
            Z_sapp = [tree.get_crown_area()*3 for tree in self.trees if tree.dbh >= sapp_dim and tree.dbh < mature_dim]
            plt.scatter(X_sapp, Y_sapp, Z_sapp,c ='r', label=f"Sapling ({len(X_sapp)})")
        if mat:
            X_mat = [tree.x for tree in self.trees if tree.dbh >= mature_dim]
            Y_mat = [tree.y for tree in self.trees if tree.dbh >= mature_dim]
            Z_mat = [tree.get_crown_area()*4 for tree in self.trees if tree.dbh >= mature_dim]
            plt.scatter(X_mat, Y_mat, Z_mat, c='blue', label=f"Mature ({len(X_mat)})")

        plt.xlim(0, L)
        plt.ylim(0, L)
        plt.xlabel("Longueur (m)")
        plt.ylabel("Largeur (m)")
        plt.title(f"Année {annee} ({self.n_cold} épisode(s) de froid)" if annee != 0 else f"Initialisation")
        plt.legend()
        plt.savefig(f"images\\seed_{self.seed}\\cold_{self.n_cold}\\seed_{self.seed}_cold_{self.n_cold}_annee_{annee}.png")
        plt.close()

    def plot_BIOM(self):
        """Plot total biomass of the forest over the simulation
        """
        AGB = [val/(L*L/10000)*2 for val in self.AGB]
        BGB = [val/(L*L/10000)*2 for val in self.BGB]
        BIOM = np.array([sum(x) for x in zip(AGB, BGB)])
        np.save(f"BIOM_{self.n_cold}_sans.npy", BIOM)

        BIOM_zero = np.load(f"BIOM_0.npy") # reference with 0 cold events

        ll = len(AGB)
        plt.plot(list(range(ll)), AGB, label="AGB", c="b")
        plt.plot(list(range(ll)), BGB, label="BGB", c="r")
        plt.plot(list(range(ll)), BIOM, label="BIOM", c="g")
        plt.plot(list(range(ll)), BIOM_zero[:ll], label="BIOM Zero", c="g", linestyle="dashed")

        plt.xlabel("Années")
        plt.ylabel("Carbone t/ha")

        plt.title(f"BIOM en fonction du temps ({self.n_cold} épisode(s) de froid)")
        
        plt.annotate("%0.2f"%BIOM[-1], xy=(1, BIOM[-1]),xytext=(8, 0), 
             xycoords=('axes fraction', 'data'), textcoords='offset points')
        plt.annotate("%0.2f"%BIOM_zero[ll-1], xy=(1, BIOM_zero[ll-1]),xytext=(8, 0), 
             xycoords=('axes fraction', 'data'), textcoords='offset points')
        
        plt.axhline(y=BIOM_zero[ll-1], color="g", linestyle="dashed")
        plt.axhline(y=BIOM[-1], color='y', linestyle="dashed")

        plt.axhline(y=80, color="m", linestyle='-.', label="France")

        plt.legend()
        plt.savefig(f"images\\seed_{self.seed}\\cold_{self.n_cold}\\seed_{self.seed}_cold_{self.n_cold}_BIOM.png")
        plt.close()

    def update_step(self, nStrides, n_cold):
        """simulation of a single year

        Args:
            nStrides (int): number of strides for competition computation
            n_cold (int): number of cold events
        """
        # kill trees due to cold
        for _ in range(n_cold):
            self.cold_event()
        print(f"alive after cold {self.nb_alive()}")

        self.trees = np.array([tree for tree in self.trees if tree.dbh > 0])

        CalculateFA(self.trees, nStrides)

        f_vect = np.vectorize(f)
        self.trees = f_vect(self.trees)

        self.trees = np.array([tree for tree in self.trees if tree.dbh > 0])
        self.AGB.append(self.total_AGB()) 
        self.BGB.append(self.total_BGB())

        # new random trees
        new_trees = np.array([])
        for tree in self.trees:
            if tree.dbh >= mature_dim:
                new_seeds = random_seeds(tree.x, tree.y, tree.gen+ 1, 
                                         dist=5+tree.get_crown_diameter(), 
                                         dist_min=tree.get_crown_diameter())
                new_trees = np.concatenate((new_trees, new_seeds))

        real_new_trees = np.array([trees for trees in new_trees if can_establish(trees, self)])
        self.trees = np.concatenate((self.trees, real_new_trees))
        if self.nb_alive() < pop:
            L = [t.gen for t in self.trees if t.dbh>mature_dim]
            L.append(1)
            # add seeds
            self.initialize(pop-self.nb_alive(), gen=max(L))
        print(f"nb after year {self.nb_alive()}")


class Tree:
    """1 tree
    """
    def __init__(self, x, y, dbh=sap_dim, FA=0.0, gen=1):
        """

        Args:
            x (float): x position of tree (m)
            y (float): y position of tree (m)
            dbh (float, optional): diamenter at breast height (cm). Defaults to sap_dim.
            FA (float, optional): field effect. Defaults to 0.0.
            gen (float, optional): cold resistance. Defaults to 1.
        """
        self.x = x
        self.y = y
        self.dbh = dbh # half diameter at breast height
        self.R = a * np.sqrt(self.dbh/2)  # radius of ZOI (m)
        self.FA = FA  # field effect
        self.gen = gen
        self.memory = deque(maxlen=5) # memory for tree survival

    def __str__(self):
        return f"x {self.x:.2f}, y {self.y:.2f}, dbh {self.dbh:.2f}cm, R {self.R:.2f}m, FA {self.FA:.2f}, gen {self.gen:.2f} years, AGB {self.AGB():.2f}"

    def add_memory(self, delta):
        """add memory to tree

        Args:
            delta (float): growth in 1 year
        """
        self.memory.append(delta)

    def update_info(self):
        self.R = a * np.sqrt(self.dbh)

    def get_crown_diameter(self):
        return 0.222 * self.dbh**0.654
    
    def get_crown_area(self):
        return np.pi * self.get_crown_diameter()**2

    def AGB(self):
        """Computes AGB

        Returns:
            float: AGB (tons)
        """
        result = 0
        if self.dbh >= 1 and self.dbh < 4:
            result += 0.0002004 * self.dbh**2.1
        elif self.dbh >= 4:
            result += 0.00014 * self.dbh**2.4
        return result

    def BGB(self):
        """Computes BGB

        Returns:
            float: BGB (tons)
        """
        result = 0
        h_meters = H(self.dbh)/100
        AGB = self.AGB()
        if h_meters <= 1:
            result += 0.75*AGB
        elif h_meters <= 2.5:
            result += 0.89*AGB
        else:
            result += 0.98*AGB
        return result


    def is_alive(self):
        return self.dbh > 0

    def kill(self, source="natural"):
        """kill the tree

        Args:
            source (str, optional): cause of death. Defaults to "natural".
        """
        # print(f'{source} dbh {self.dbh} x+y {self.x + self.y} L {L}, salt : {sigma(self, L)}, inundation : {eta(self, L)}, compet : {self.FA}, gen : {self.gen}')
        # print(self)
        self.x = 0.1
        self.y = 0.1
        self.dbh = 0
        self.R = 0

    def FON_PCR(self, phi, r0, r1):
        """compute FON integral from r0 to r1 and from -phi to phi

        Args:
            phi (float):    
            r0 (float): 
            r1 (float): 

        Returns:
            float: FON integral
        """
        return 2 * phi * self.integrate(r0, r1)
    
    def partition(self, r0, r1):
        """partition [r0, r1] with respect to FON conditions

        Args:
            r0 (float): 
            r1 (float): 

        Returns:
            tuple[float]: partition
        """
        return (min(r0, self.R), min(r1, self.R))

    def integrate(self, r0, r1):
        """integrate FON from r0 to r1

        Args:
            r0 (float): 
            r1 (float): 

        Returns:
            float: 
        """

        c = 1 
        x1, x2 = self.partition(r0, r1)
        # print(r0, r1, self.dbh)
        rbh = self.dbh/2/100
        # print(result, self)
        return max((np.exp(-c*(x1-rbh)) - np.exp(-c*(x2-rbh))) / c, F_min)
        

def FA_kn(treeK, treeN, nStrides):
    """Calculate the integral over the overlapping area of two trees.

    Args:
        treeK (Tree): influenced tree
        treeN (Tree): influencer tree
        nStrides (int): computing parameter

    Returns:
        float: integral
    """
    result = 0.0  # suppose they do not overlap
    # distance between the streaming positions (m)
    A = np.sqrt((treeK.x-treeN.x)**2 + (treeK.y-treeN.y)**2)  
    

    # no overlap
    if A >= treeN.R + treeK.R: 
        return result  

    # first we compute the max. diameter of the overlapping area
    r0 = A - treeK.R  # ‘left’ side of overlapping area

    # a whole concentric circle (grayed in Fig. A1) is overlapped
    if r0 < 0.0:
        r0 = -r0
        result = treeN.FON_PCR(np.pi, 0.0, r0)

        # treeN is completely overlapped by treeK
        if r0 >= treeN.R:
            return result  # return the exact result to the caller

    r1 = min(treeN.R, A + treeK.R)  # ‘right’ side of overlapping area
    dr = (r1 - r0) / nStrides  # division into nStrides pitch circle rings
    r = r0 + dr / 2  # centre radius of the first PCR

    for i in range(nStrides):
        # Pythagoras theorem to calculate the intersection coordinates
        x = (A**2 + r**2 - treeK.R**2) / (2 * (A+1e-7))
        y = np.sqrt(r**2 - x**2)
        
        # calculate the angle phi
        phi = np.arctan2(y, x)
        
        # add the integral over this pitch circle ring
        result += treeN.FON_PCR(phi, r - dr / 2, r + dr / 2)
        
        # go to the next pitch circle ring
        r += dr

    return result  # return the result to the caller


def CalculateFA(trees, nStrides):
    """summation of the overlapping FON of all trees

    Args:
        trees (): array of trees
        nStrides (int): computing parameter
    """
    N = len(trees)
    
    # Initialize FA for all trees
    for n in range(N):
        trees[n].FA = 0.0
    
    # Calculate FA for each pair of trees
    for k in range(N):
        for n in range(k + 1, N):
            trees[k].FA += FA_kn(trees[k], trees[n], nStrides)
            trees[n].FA += FA_kn(trees[n], trees[k], nStrides)
        
        # Normalize by the ZOI
        trees[k].FA = trees[k].FA / (np.pi * trees[k].R**2)

