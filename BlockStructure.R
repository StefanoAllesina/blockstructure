##############################################################################################
# Code accompanying the manuscript
# "Modularity and stability in ecological communities"
# by Jacopo Grilli & Stefano Allesina
##############################################################################################
# v. 0.1 --- June 2015
# For any question/comment/bug, please contact 
# Stefano Allesina sallesina@uchicago.edu
##############################################################################################
# LICENSE:
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################################
# USE:
# This file shows how to use of the programs provided in the files "BuildMatrices.R" 
# and "Predictions.R". 
#
# There are three parts:
# Part 1: How to build the matrices
# Part 2: How to predict the eigenvalues
# Part 3: Examples taken from Figure 8 of the manuscript.
##############################################################################################

# code for building matrices
source("BuildMatrices.R")
# code for predicting eigenvalues and plotting
source("Predictions.R")

##############################################################################################
# 1. HOW TO BUILD THE MATRICES
##############################################################################################
# The file "BuildMatrices.R" contains the function BuildMatrices, which builds all the matrices
mats <- BuildMatrices(S = 1000,          # size of the system
                      connectance = 0.1, # overall connectance (C in the manuscript)
                      a = 0.4,           # proportion of species in the smaller subsystem (\alpha in the manuscript)
                      mu = 0.5,          # mean of the coefficients in W
                      sigma = 1,         # standard deviation coefficients in W
                      rho = -0.5,        # correlation between pairs in W
                      Q = 0.48,          # modularity
                      Distr = "Normal",  # either "Normal" (as in the manuscript), or "FourCorner" (introduced in Allesina et al., Nature Communications 2015)
                      Eyeball = FALSE)   # if TRUE, build a food web using the modified cascade model, as explained in the manuscript
# see the structure of the return object
str(mats)

# mats$parameters contains all the parameters
print(mats$parameters)

# mats$membership the subsystem membership (\gamma in the manuscript)
print(mats$membership[1:10])

# mats$matrices is a list of the three matrices M, A, and B
print(mats$matrices$M[1:10, 1:10])

##############################################################################################
# 2. PREDICTING THE EIGENVALUES
##############################################################################################
# In the case of random ecosystems (i.e., Eyeball == FALSE), the manuscript introduces 
# analytical predictions for the rightmost eigenvalue (Re(\lambda_{M, 1})) when:
# 1. the subsystems have equal size (i.e., a = 1/2)
# 2. the network is perfectly modular (i.e., Cb = 0)
# 3. the network is perfectly bipartite (i.e., Cw = 0)

# The function Predictions in the file "Predictions.R" calculates the eigenvalues of M
# and (if it's one of the cases treated in the manuscript) the corresponding analytical prediction
# By default, the function also plots the predictions
pr <- Predictions(mats, plot = TRUE)

# The return object contains the observed as well as the predicted real part of the leading eigenvalue:
print(paste("Re(lambda_{M, 1}) Observed:", pr$ReL1.observed, "Predicted:", pr$ReL1.predicted))

# The object also stores the rightmost eigenvalue on the bulk
print(pr$Re.bulk)
# and rightmost outlier
print(pr$Re.outlier)

# The object also stores the eigenvalues of M
plot(pr$eigenvalues.M)

# For the corresponding unstructured case, repeat the same procedure, but setting Q = 0

##############################################################################################
# 3. EXAMPLES
##############################################################################################
# Interesting cases (those in Figure 8 in the supplementary information)
set.seed(10)

# Figure 8, Equal Size a)
EqualSize.a <- BuildMatrices(S = 1000, 
                            connectance = 0.1, 
                            a = 0.5, 
                            mu = 0.5, 
                            sigma = 1, 
                            rho = -0.5, 
                            Q = -0.2, 
                            Distr = "Normal", 
                            Eyeball = FALSE)
pr <- Predictions(EqualSize.a)
# We can change the distribution to a discrete one (See Allesina et al. Nature Communications, 2015)
# and the results are qualitatively the same
EqualSize.a1 <- BuildMatrices(S = 1000, 
                              connectance = 0.1, 
                              a = 0.5, 
                              mu = 0.5, 
                              sigma = 1, 
                              rho = -0.5, 
                              Q = -0.2, 
                              Distr = "FourCorner", 
                              Eyeball = FALSE)
pr <- Predictions(EqualSize.a1)

# Figure 8, Equal Size b)
EqualSize.b <- BuildMatrices(S = 1000, 
                            connectance = 0.2, 
                            a = 0.5, 
                            mu = -0.5, 
                            sigma = 1, 
                            rho = 0.5, 
                            Q = 0, 
                            Distr = "Normal", 
                            Eyeball = FALSE)
pr <- Predictions(EqualSize.b)

# Figure 8, Equal Size c)
EqualSize.c <- BuildMatrices(S = 1000, 
                            connectance = 0.4, 
                            a = 0.5, 
                            mu = 0.5, 
                            sigma = 1, 
                            rho = 0, 
                            Q = 0.2, 
                            Distr = "Normal", 
                            Eyeball = FALSE)
pr <- Predictions(EqualSize.c)

# Now perfectly modular structures
# Figure 8, Perfectly Modular a)
Modular.a <- BuildMatrices(S = 1000, 
                          connectance = 0.1, 
                          a = 0.4, 
                          mu = 1.0, 
                          sigma = 1, 
                          rho = -0.5, 
                          Q = 0.48, 
                          Distr = "Normal", 
                          Eyeball = FALSE)
pr <- Predictions(Modular.a)

# Figure 8, Perfectly Modular b)
Modular.b <- BuildMatrices(S = 1000, 
                          connectance = 0.3, 
                          a = 0.3, 
                          mu = -1.0, 
                          sigma = 1, 
                          rho = -0.25, 
                          Q = 0.42, 
                          Distr = "Normal", 
                          Eyeball = FALSE)
pr <- Predictions(Modular.b)

# Figure 8, Perfectly Modular c)
Modular.c <- BuildMatrices(S = 1000, 
                          connectance = 0.4, 
                          a = 0.1, 
                          mu = 0.5, 
                          sigma = 1, 
                          rho = 0.35, 
                          Q = 0.18, 
                          Distr = "Normal", 
                          Eyeball = FALSE)
pr <- Predictions(Modular.c)


# Finally, perfectly bipartite structures
# Figure 8, Perfectly Bipartite a)
Bipartite.a <- BuildMatrices(S = 1000, 
                             connectance = 0.1, 
                             a = 0.1, 
                             mu = -0.5, 
                             sigma = 1, 
                             rho = -0.5, 
                             Q = -0.82, 
                             Distr = "Normal", 
                             Eyeball = FALSE)
pr <- Predictions(Bipartite.a)

# Figure 8, Perfectly Bipartite b)
Bipartite.b <- BuildMatrices(S = 1000, 
                             connectance = 0.25, 
                             a = 0.25, 
                             mu = -0.5, 
                             sigma = 1, 
                             rho = 0.5, 
                             Q = -0.625, 
                             Distr = "Normal", 
                             Eyeball = FALSE)
pr <- Predictions(Bipartite.b)

# Figure 8, Perfectly Bipartite c)
Bipartite.c <- BuildMatrices(S = 1000, 
                             connectance = 0.2, 
                             a = 0.3,
                             mu = 0.25, 
                             sigma = 1, 
                             rho = -0.25, 
                             Q = -0.58, 
                             Distr = "Normal", 
                             Eyeball = FALSE)
pr <- Predictions(Bipartite.c)
