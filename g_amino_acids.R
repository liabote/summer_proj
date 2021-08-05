install.packages("GillespieSSA")
library(GillespieSSA)

#stochastic simulations for first 5 reactions in CO2 fixation pathway

#reaction constants
#only C1/3/6 involved in first 5 rxns
#values set based on relative frequency of rxn category in full pathway 
#ie. C1 is more common than C6, wc is more common than C3
parms <- c(k1=2.53, k2= 2.53, k3=1, k4=1.86, k6=1.13) 

#initial concentrations of reactants and products involved 
#concentrations arbitrarily assigned 
x0 <- c(CO2 = 5000, NADPH = 5000, formate = 0, NADP = 0,
        tetrahyd = 5000, ATP = 5000, Nformyl = 0, ADP = 0, P = 0, 
        H = 5000, methenyltet = 0, H2O = 5000, 
        methylenetet = 0, red_FeS = 5000, methyltet = 0, ox_FeS = 5000,
        NADH = 5000, NAD = 0, CoI_FeS = 5000, metCoIII_FeS = 5000, 
        CO = 0, coA = 5000, acetyl_coA = 0, 
        lipo_lys = 5000, dihyd_lys = 5000, NH4 = 5000, 
        amino_dihyd = 5000, gly = 0, oxo_glut = 5000, 
        glu = 0, pyruvate = 5000, ala = 0,
        aceto_lac = 5000, dihyd_metbut = 5000, 
        met_oxobut = 5000, val = 0, 
        HCO3 = 5000, oxaloac = 0, malate = 0, 
        fumarate = 0, menaquinol = 5000, 
        succinate = 0, menaquinone = 5000, 
        succ_coA = 5000, asp = 0)

#rate equations for forward rxns 
#rate constant assigned based on EC number 

      
a <- #carbon fixation 
      c("k1*CO2*NADPH", 
       "k6*tetrahyd*ATP*formate", 
       "k3*Nformyl*H", 
       "k1*methenyltet*NADPH", 
       "k1*methylenetet*red_FeS^2*H^2",
       "k1*methylenetet*NADH*H", 
       "k1*methylenetet*ox_FeS^2*NADH^2",
       "k2*CoI_FeS*methyltet*H",
       "k1*CO2*red_FeS^2*H^2",
       "k2*CO*metCoIII_FeS*coA", 
       
      #reductive glycine synthesis
       "k1*lipo_lys*NADH*H",
       "k2*dihyd_lys*methylenetet*NH4",
       "k1*amino_dihyd*CO2",
       
      #glutamate synthesis 
       "k1*oxo_glut*NH4*NADPH*H",
       
      #alanine synthesis 
       "k2*glu*pyruvate",
      
      #aliphatic amino acid synthesis
       "k2*pyruvate^2*H",
       "k1*aceto_lac*NADPH*H",
       "k4*dihyd_metbut",
       "k2*glu*met_oxobut",
       
      #reductive incomplete TCA
       "k1*acetyl_coA*CO2*red_FeS^2*H",
       "k6*pyruvate*HCO3*ATP",
       "k1*oxaloac*NADH*H",
       "k4*malate",
       "k1*fumarate*menaquinol",
       "k6*succinate*ATP*coA",
       "k1*succ_coA*CO2*red_FeS^2*H",
      
      #aspartate synthesis
       "k2*oxaloac*glu")

#state change matrix
#stoichiometric ratios of reactants and products 
#rows: reactants/products, columns: reactions
mat <- read.csv("matrix.csv", header=F)
mat <- as.matrix(mat[-1])


out <- ssa(x0, a, mat, parms, tf=.0005)
ssa.plot(out, show.legend =T)

data <- out[["data"]]

par(mfrow = (c(2,3)))
plot(data[,1], data[,29], main ="glycine", pch=20)
plot(data[,1], data[,31], main ="glutamine", pch=20)
plot(data[,1], data[,33], main ="alanine", pch=20)
plot(data[,1], data[,37], main ="valine", pch=20)
plot(data[,1], data[,46], main ="aspartate", pch=20)