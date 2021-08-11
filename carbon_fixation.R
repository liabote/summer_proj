#install.packages("GillespieSSA")
library(GillespieSSA)

#stochastic simulations for first 5 reactions in CO2 fixation pathway

#reaction constants
#only C1/3/6 involved in first 5 rxns
#values set based on relative frequency of rxn category in full pathway 
#ie. C1 is more common than C6, wc is more common than C3
parms <- c(k1=2.53, k2= 2.53, k3=1.0, k4=1.86, k6=1.13) 


#intro to code how it works 
#what ive done ec numbers proportional to rates for 5 aa 
#graphs and what they do 

nmol <- 5000
#initial concentrations of reactants and products involved 
#concentrations arbitrarily assigned 

x <- c(CO2 = nmol, NADPH = nmol, formate = 0, NADP = 0,
        tetrahyd = nmol, ATP = nmol, Nformyl = 0, ADP = 0, P = 0, 
        H = nmol, methenyltet = 0, H2O = nmol, 
        methylenetet = 0, red_FeS = nmol, methyltet = 0, ox_FeS = nmol,
        NADH = nmol, NAD = 0, CoI_FeS = nmol, metCoIII_FeS = 0,
        CO = 0, coA = nmol, acetyl_coA = 0)


        #lipo_lys = nmol, dihyd_lys = nmol, NH4 = nmol,
        #amino_dihyd = 0, gly = 0)
        #pyruvate = nmol, ala = 0,
        #aceto_lac = nmol, dihyd_metbut = nmol, 
        #met_oxobut = nmol, val = 0, 
        #HCO3 = nmol, oxaloac = 0, malate = 0, 
        #fumarate = 0, menaquinol = nmol, 
        #succinate = 0, menaquinone = nmol, 
        #succ_coA = nmol, oxo_glut = nmol, 
       #glu = 0, asp = 0)

a <- #carbon fixation 
  c("k1*CO2*NADPH", 
    "k6*tetrahyd*ATP*formate", 
    "k3*Nformyl*H", 
    "k1*methenyltet*NADPH", 
    "k1*methylenetet*NADH*H", 
    "k2*CoI_FeS*methyltet*H",
    "k1*CO2*red_FeS^2*H^2",
    "k2*CO*metCoIII_FeS*coA",
    "formate*NADP",
    "Nformyl*ADP*P",
    "methenyltet*H2O",
    "methylenetet*NADP",
    "methyltet*NAD",
    "metCoIII_FeS*tetrahyd",
    "CO*ox_FeS^2*H2O",
    "acetyl_coA*CoI_FeS*H")
    

    #reductive glycine synthesis
    #"k1*lipo_lys*NADH*H",
    #"k2*dihyd_lys*methylenetet*NH4",
    #"k1*amino_dihyd*CO2")
    
    #alanine synthesis 
    #"k2*NH4*NADPH*H*pyruvate",
    
    #aliphatic amino acid synthesis -> valine
    #"k2*pyruvate^2*H",
    #"k1*aceto_lac*NADPH*H",
    #"k4*dihyd_metbut",
    #"k2*NH4*NADPH*H*met_oxobut",
    
    #reductive incomplete TCA
    #"k1*acetyl_coA*CO2*red_FeS^2*H",
    #"k6*pyruvate*HCO3*ATP",
    #"k1*oxaloac*NADH*H",
    #"k4*malate",
    #"k1*fumarate*menaquinol",
    #"k6*succinate*ATP*coA",
    #"k1*succ_coA*CO2*red_FeS^2*H",
    
    
    #glutamate synthesis 
    #"k1*oxo_glut*NH4*NADPH*H",
    
    #aspartate synthesis
    #"k2*oxaloac*NH4*NADPH*H")

#state change matrix
#stoichiometric ratios of reactants and products 
#rows: reactants/products, columns: reactions
#mat_t <- read.csv("matrix_trial.csv", header=T)
#mat_t <- as.matrix(mat_t[-1])


mat <- read.csv("Book2.csv", header=F)
mat <- as.matrix(mat[-1])


out <- ssa(x, a, mat, parms, tf=0.00001)

#ssa.plot(out, show.legend =T)


data <- out[["data"]]

par(mfrow = (c(5,5)))

plot(data[,1], data[,2], main ="CO2", pch=20)
plot(data[,1], data[,3], main ="NADPH", pch=20)
plot(data[,1], data[,4], main ="formate", pch=20)
plot(data[,1], data[,5], main ="NADP+", pch=20)
plot(data[,1], data[,6], main ="tetrahydrofolate", pch=20)
plot(data[,1], data[,7], main ="ATP", pch=20)
plot(data[,1], data[,8], main ="N10-formyltetrahydrofolate", pch=20)
plot(data[,1], data[,9], main ="ADP", pch=20)
plot(data[,1], data[,10], main ="P", pch=20)
plot(data[,1], data[,11], main ="H+", pch=20)
plot(data[,1], data[,12], main ="methenyltet", pch=20)
plot(data[,1], data[,13], main ="H2O", pch=20)
plot(data[,1], data[,14], main ="methylenetet", pch=20)

plot(data[,1], data[,15], main ="red FeS", pch=20)
plot(data[,1], data[,16], main ="methyltet", pch=20)
plot(data[,1], data[,17], main ="ox FeS", pch=20)

plot(data[,1], data[,18], main ="NADH", pch=20)
plot(data[,1], data[,19], main ="NAD+", pch=20)
plot(data[,1], data[,20], main ="corrinoid I FeS", pch=20)
plot(data[,1], data[,21], main ="methyl co III", pch=20)
plot(data[,1], data[,22], main ="CO", pch=20)
plot(data[,1], data[,23], main ="coA", pch=20)
plot(data[,1], data[,24], main ="acetyl coA", pch=20)

#plot(data[,1], data[,25], main ="lipolyl-lysine", pch=20)
#plot(data[,1], data[,26], main ="dihydrolipolyl-", pch=20)
#plot(data[,1], data[,27], main ="NH4", pch=20)
#plot(data[,1], data[,28], main ="aminomethyldihydro-", pch=20)
#plot(data[,1], data[,29], main ="glycine", pch=20)


#plot(data[,1], data[,11], main ="H+", pch=20)

#plot(data[,1], data[,29], main ="glycine", pch=20)
#plot(data[,1], data[,45], main ="glutamate", pch=20)
#plot(data[,1], data[,31], main ="alanine", pch=20)
#plot(data[,1], data[,35], main ="valine", pch=20)
#plot(data[,1], data[,46], main ="aspartate", pch=20)
#plot(data[,1], data[,"oxo_glut"], main ="oxoglutarate", pch=20)
#plot(data[,1], data[,"malate"], main ="malate", pch=20)

