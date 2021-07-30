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
        pyruvate = 0, HCO3 = 5000, oxaloacetate = 0, 
        S_malate = 0, fumarate = 5000, 
        menaquinol = 5000, succinate = 0, menaquinone = 5000, 
        succinyl_coA = 5000, oxoglut = 0, 
        acetolac = 0, dihyd_metbut = 0, met_oxobut = 5000,
        glu = 0, val = 0, isop_malate = 0, isop_maleate = 0, 
        isop_3_malate = 0, isop_oxosucc = 0, met_42_oxopent = 0, 
        leu = 0, thr = 0, aminobut = 5000, iminobut = 0, 
        oxobut = 0, NH3 = 0, acetohydbut = 5000,
        dihyd_metpent = 0, met_32_oxopent = 0, ile = 0)

#rate equations for forward rxns 
#rate constant assigned based on EC number 
a <- c("k1*CO2*NADPH", 
       "k6*tetrahyd*ATP*formate", 
       "k3*Nformyl*H", 
       "k1*methenyltet*NADPH", 
       "k1*methylenetet*red_FeS^2*H^2",
       "k1*methylenetet*NADH*H", 
       "k1*methylenetet*ox_FeS^2*NADH^2",
       "k2*CoI_FeS*methyltet*H",
       "k1*CO2*red_FeS^2*H^2",
       "k2*CO*metCoIII_FeS*coA", 
       "k1*acetyl_coA*CO2*red_FeS^2*H",
       "k6*pyruvate*HCO3*ATP",
       "k1*oxaloacetate*NADH*H",
       "k4*S_malate",
       "k1*fumarate*menaquinol",
       "k6*succinate*ATP*coA",
       "k1*succinyl_coA*CO2*red_FeS^2*H", 
       "k2*pyruvate^2*H", 
       "k1*acetohydbut*NADPH*H", 
       "k4*dihyd_metpent",
       "k2*glu*met_oxobut",
       "k2*met_oxobut*acetyl_coA*H2O",
       "isop_malate*NAD", 
       "isop_maleate*NAD",
       "isop_malate*NAD",
       "isop_oxosucc*H", 
       "k2*glu*met_42_oxopent",
       "thr",
       "aminobut",
       "k3*iminobut*H2O*H",
       "k2*pyruvate*oxobut",
       "k1*acetolac*NADPH*H",
       "k4*dihyd_metbut",
       "k2*glu*met_32_oxopent")


#state change matrix
#stoichiometric ratios of reactants and products 
#rows: reactants/products, columns: reactions
mat <- read.csv("matrix.csv", header=F)
mat <- as.matrix(mat[-1])

out <- ssa(x0, a, mat, parms, tf=.0005)
ssa.plot(out, show.legend =T)

#################################################################################

x <- data.frame(out$data[,1])
y <- data.frame(out$data[,4]) 
y
out_2 <- ssa(x0, a, nu, parms, tf=10)
ssa.plot(out_2, plot.from=4, plot.to=4)

par(mfrow=c(3,3))
ssa.plot(out,  plot.from=2, plot.to=2) #CO2 
ssa.plot(out,  plot.from=4, plot.to=4) #formate 
ssa.plot(out,  plot.from=7, plot.to=7) #Nformyltetrahydrofolate 
ssa.plot(out,  plot.from=12, plot.to=12) #methenyltetrahydrofolate 
ssa.plot(out,  plot.from=14, plot.to=14) #methylenetetrahydrofolate 
ssa.plot(out,  plot.from=16, plot.to=16) #methyltetrahydrofolate 
ssa.plot(out,  plot.from=21, plot.to=21) #methyl Corrinoid(III)-FeS 
ssa.plot(out,  plot.from=22, plot.to=22) #carbon monoxide 
ssa.plot(out,  plot.from=24, plot.to=24) #acetyl coA 

library(ggplot2)
g <- ggplot(out, ggplot2::aes_string(x,y)) + geom_point()  + geom_smooth(method="lm") 
plot(g)