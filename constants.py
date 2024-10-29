# dimension of study (m)
L = 100
# initial population (int)
pop = 100

# years of study
years = 75
# seedling dim (cm)
sap_dim = 0.5
# sappling dim (cm)
sapp_dim = 2.5
# mature dim (cm)
mature_dim = 5

# Physical Properties
# slope parameter [0, 1]
sl = 0.03
# effects on the salinity gradient [-0.75, 0.25]
ds = -0.2

# INFO FOR A.germinans
# scaling factor
a = 10
# minimum value of the FON
F_min = 0.1
# grow constant
G = 162
# maximum dbh (cm)
D_max = 50
# maximum age (year)
age_max = 300
# ratio
half_D_over_age = 0.5*D_max/age_max
# constant in height to dbh relationship
b2 = 48.04
# constant in height to dbh relationship
b3 = 0.172
# constant for salt effect on growth
d = -0.18
# salt effect on growth (g/kg)
Ui = 72.0
# additional yearly mortality
epsilon = 1/600
