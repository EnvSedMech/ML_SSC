## Demo for ESM group
options(digits = 10)

SSC_list2 <- readRDS("SSC_list.rds")
Q_h <- readRDS("Q_h_list.rds")
dgo_list <- readRDS("dgo_list.rds")
Q_v <- Q_h[[2]]
h_v <- Q_h[[3]]

###
Qin <- Q_v[1]
hinitial <- h_v[1]

###
SSC <- SSC_list2[[1]]
dgo <-dgo_list[[1]]


## Main script
### Input Data
S <- c(0.0039, 0.0039, 0.0039, 0.0040, 0.0031, 0.0044, 0.0044, 0.0022)
S <- rev(S)
## Grid parameters
L <- c(720, 720, 720, 600, 740, 520, 1040,2080)
L <- rev(L)
dx <- 20

width <- c(25, 25, 25, 35, 34, 39, 39, 39)
width <- rev(width) 
Ls <- sum(L)
nc = (Ls/dx)+2
nL = L/dx
numS <- length(S)

## River Channel characteristics
b <- rep(0,nc)
b[-c(1,nc)] <- rep(width,((L/dx)))
b[1] <- b[2]
b[nc] <- b[nc-1]

nf <- c(0.024, 0.024, 0.024,	0.059,	0.046,	0.051, 0.051, 0.051)
nf <- rev(nf)
n <- rep(0,nc)
n[-c(1,nc)] <- rep(nf,((L/dx)))
n[1] <- n[2]
n[nc] <- n[nc-1]
nf = n

## Flow characterisitics



pw <- 1000
kvis <- 1e-6

## Input sediment characteristics
ps = 2650
k = 9

## time marching parameters
iter_max = 1000000
CFL = 0.01
dts = 60

### Boundary conditions
hin <- 0
hout <- 0


## Simulation options
Xi_cr = 0.1
g = 9.81

# Riverbed discretization
cd <- function(L,S,dx){
  Ls <- sum(L)
  nc = (Ls/dx)+2
  nL = L/dx
  numS <- length(S)
  
  zd = 0
  z = zd
  
  zcd = zd - (S[1]*(dx/2))
  zc = zcd
  
  for(i in 1: numS){
    for( j in 1:nL[i]){
      zu = zd + S[i]*dx
      
      zcu = zcd + S[i]*dx
      zc = c(zc,zcu)
      zcd = zcu
      
      z = c(z,zu)
      zd = zu
    }
  }
  return(as.numeric(zc))
}

xx <- cd(L,S,dx)
xxx <- rep(0,nc)
xxx[-nc] <- xx
xxx[nc] <- xxx[nc-1] + (S[as.numeric(length(S))]*dx)
zc  = xxx

## initial riverbed profile
zco = zc
zco
# Initialize variables
## Hydraulic sub model
h = rep(0,nc)
z = rep(0,nc)
A = rep(0,nc)
P = rep(0,nc)
R = rep(0,nc)
Q = rep(0,nc)
V = rep(0,nc)
Q[nc] = Qin
Sf <- rep(0,nc)

## Sediment transport sub model
Qb = rep(0,nc)
Qt = rep(0,nc)
Qs = rep(0,nc)
Ct = rep(0,nc)
C = rep(0,nc)

Qb_k = matrix(0,k,nc)
Qt_k = matrix(0,k,nc)
Qs_k = matrix(0,k,nc)
Ct_k = matrix(0,k,nc)
C_k = matrix(0,k,nc)
ws_k = matrix(0,k,nc)

## adaptive numerical scheme
Psi = rep(0,nc)
Psi[c(1,nc)] <- 0
## initialize variables for functions

## 2) Deposition non - uniform
## Preallocating variables
alpha <- rep(0,nc)
Cb_k  = matrix(0,k,nc)
D_k = matrix(0,k,nc)

### 3) Entrainment non - uniform
## Preallocating variables

ws_star_k = matrix(0,k,nc)
Cb_star_k = matrix(0,k,nc)
E_k = matrix(0,k,nc)


### 4) bedload non uniform
# Prellocating variables
Tau_star_k = matrix(0,k,nc)
qb_star_k = matrix(0,k,nc)
Tau_star_c_k = rep(0.047,nc)
Tau_star_c_k <- matrix(Tau_star_c_k,k,nc)
# porosity of the riverbed
lambdao = rep(0,nc)

# density of the riverbed surface layer
po = rep(1000, nc)

# Fractional settling velocity
wso_k <- matrix(0,k,nc)

## grain size
library(pracma)
# Limits and representative grain sizes considered for each sediment class
#dk_limits <- c(16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.062)
dk_limits <- rev(c(0.125,0.5,1,4,8,16,32,64))
dk_limits <- sort(dk_limits, decreasing = F)

## sediment classes
k = 9
###Input sediment supply
pk1 <- c(0.0,0.16,0.34,0.40,0.10,0.00,0.00,0.00,0.00)
pk1 <- matrix(pk1,k,nc)

pk <- c(0.0,0.0,0.0,0.14,0.30,0.28,0.22,0.08,0.08)
pk = matrix(pk,k,nc)

sum <- matrix(0,k,nc)
pk_limits <- matrix(0,8,nc)

# Grain size diameters in m
#dk <- c(0.062, 0.0935, 0.1875, 0.375, 0.75, 1.5, 3, 6, 12, 16)
dk <- c(0.125, 0.375, 0.75,2,6,12,24,48,64)

#dk <- c(0.062, 0.1874, 1,2,3,4,8,16,32,64)
dk = matrix(dk/1000,k,nc)
d50 = 10/1000
d90 = 48/1000

d50 = rep(d50,nc)
d90 = rep(d90,nc)

## Input flow with sediment characteristics
Qtin <- rep(0,nc)
Css <- 0
Qtin[nc] <- Css*Qin
Qtin_k <- t(Qtin*t(pk1))
#bay <- c(0.00,0.01,0.02,0.04,0.06,0.08,0.1,0.3)
#Qtin[nc] <- Q_v[1]*bay[1]
#Qtin[nc]
#Qtin_k <- t(Qtin*t(pk1))

# Set initial conditions for boundary condition
# before iter 1
h = rep(hinitial,nc)
z[-c(1,nc)] = h[-c(1,nc)] + zc[-c(1,nc)]
A[-c(1,nc)] = h[-c(1,nc)]*b[-c(1,nc)]
P[-c(1,nc)] = b[-c(1,nc)] + (2*h[-c(1,nc)])
R[-c(1,nc)] = A[-c(1,nc)]/P[-c(1,nc)]
V[-c(1,nc)] = Q[-c(1,nc)]/A[-c(1,nc)]
Sf[-c(1,nc)] = ((nf[-c(1,nc)]*Q[-c(1,nc)])^2)/((R[-c(1,nc)]^(4/3))*(A[-c(1,nc)]^2))

## Submerged specific gravity of sediment particles
Rs = (ps-pw)/pw
d50_xc <- d50

# porosity of the riverbed
lambdao[-c(1,nc)] = 0.13+0.21/((d50_xc[-c(1,nc)]*1000+0.002)^0.21)
lambdao[1] = lambdao[2]
lambdao[nc] = lambdao[nc-1]

# density of the riverbed surface layer
pw  = 1000
ps = 2650
po[-c(1,nc)] = pw*lambdao[-c(1,nc)]+ps*(1-lambdao[-c(1,nc)])
po[1] = po[2]
po[nc] = po[nc-1]

# Fractional settling velocity
rs2 <- matrix(Rs,k,nc)
wso_k[1:k,1:nc] = ((13.95*kvis/dk)^2+(1.09*rs2*g*dk))^(0.5)-(13.95*kvis/dk)
ws_k <- wso_k
alpha <- rep(0 ,nc)
E <- rep(0,nc)
D <- rep(0,nc)

## Simulation variables
## Dummy variables
a1 <- seq(nc-1,2,-1)
b1 <- seq(nc,3,-1)
c1 <- seq(nc-2,1,-1)
Hsl_k = matrix(0,k,nc)
AC_k = matrix(0,k,nc)
dAb_k = matrix(0,k,nc)
dQb_kdx = matrix(0,k,nc)
L2_norm_h <- 1000
L2_norm_Q <- 1000
Hc = rep(0,nc-1)
w1 = rep(0,nc-1)
w2 = rep(0,nc-1)
dzdx_down = rep(0,nc-1)
dzdx_up = rep(0,nc-1)
Source_cwf = rep(0,nc-1)
Source_slf = rep(0,nc-1)
Hm = rep(0,nc-1)

# friction slope
CFL <- 0.05

### while Loop 
iter = 0
t = 0
count_ts = 1
dts = 60

## intitialize model output variables
hit <- matrix(0,count_ts,nc)
zit <- matrix(0,count_ts,nc)
Qit <- matrix(0,count_ts,nc)
Vit <- matrix(0,count_ts,nc)
Qtit <- matrix(0,count_ts,nc)
Qbit <- matrix(0,count_ts,nc)
Qsit <- matrix(0,count_ts,nc)
Ctit <- matrix(0,count_ts,nc)
Cit <- matrix(0,count_ts,nc)

## Morphological evolution model
Eit <- matrix(0,count_ts,nc)
Dit <- matrix(0,count_ts,nc)
zcit <- matrix(0,count_ts,nc)
pkit = pk

## Adaptive sub model
Xiit <- matrix(0,count_ts,nc)


## while loop
p <- rep(1000,nc)
p[-c(1,nc)] = pw*(1-Ct[-c(1,nc)])+ps*(Ct[-c(1,nc)])

dQbdx <- rep(0,nc)
dCtdx <- rep(0,nc)
dAb <- rep(0,nc)
dzc <- rep(0,nc)
delta_m = rep(1,nc)
A_m = rep(0,nc)
pkn = matrix(0,k,nc)
pk_star = matrix(0,k,nc)
pk_sub = pk1
d50_xc <- d50
d90_xc <- d90
Fr <- rep(0,nc)
Xi <- rep(0,nc)
Xi_cr <- rep(0.1, nc)



Num_model2 <- function(Css){
  
  Qtin[nc] <- Css*Qin
  Qtin_k <- t(Qtin*t(pk1))
  ## While loop for simulation across space and time:
  
  
  while (L2_norm_h > 1e-2 | L2_norm_Q > 1e-2){
    if(iter>1000000){
      break
    }
    
    if(is.nan(L2_norm_h) | is.nan(L2_norm_Q)){
      break
    }
    
    
    ## iteration counter
    iter  = iter + 1
    
    ##time step
    CFL <- 0.05
    dt = CFL* min(dx/(V+((g*A/b)^0.5)))
    
    ## simulation time
    t = t+dt
    ## Set boundary conditions
    ## Upstream boundary
    Q[nc] = Qin
    z[nc] = z[nc-1]
    h[nc] = z[nc] - zc[nc]
    A[nc] = h[nc]*b[nc]
    P[nc] = b[nc] + (2*h[nc])
    R[nc] = A[nc]/P[nc]
    V[nc] = Q[nc]/A[nc]
    
    ## sediment transport sub-model
    Sf[nc] <- ((nf[nc]*Q[nc])^2)/((R[nc]^(4/3))*(A[nc]^2))
    
    ## Dimensionless boundary stress
    Tau_star_k[, nc] = t((R[nc]*Sf[nc])*t((1/Rs*dk[,nc])))
    
    ## Dimensionless fractional bedload transport rate per unit channel width
    qb_star_k = matrix(0,k,nc)
    
    # Preallocating variables
    for(j in 1:k)
    {
      if(Tau_star_k[j,nc] >= Tau_star_c_k[j,nc])
      {
        qb_star_k[j,nc] = (17*(Tau_star_k[j,nc]-Tau_star_c_k[j,nc]))*(((Tau_star_k[j,nc])^0.5)-(((Tau_star_c_k[j,nc])^0.5)))
      }
      else{qb_star_k[j,nc] == 0}
    }
    ## Fractional bedload transport rate
    Qb_k[, nc] = t(b[nc] *t(qb_star_k[, nc] * pk[, nc] * (dk[,nc]) * ((Rs*g*dk[,nc])^0.5)))
    
    ## Total bedload transport rate
    Qb = colSums(Qb_k[,])
    Qtin = colSums(Qtin_k[,])
    
    if(Qtin[nc] == 0){
      ## Total and fractional total load transport rate and concentration of total sediment load
      Qt[nc] = Qb[nc]
      Qt_k[1:k, nc] = Qb_k[1:k, nc]
      
      Ct[nc] = Qt[nc]/Q[nc]
      Ct_k[, nc] = t((1/Q[nc]) *(t(Qt_k[, nc])))
    }else if(Qtin[nc] <= Qb[nc])
    {
      ## Total and fractional suspended load transport rate and concentration of suspended sediment
      Qs[nc] = Qtin[nc]
      Qs_k[, nc] = Qtin_k[,nc]
      
      C[nc] = Qs[nc]/Q[nc]
      C_k[, nc] = t((1/Q[nc]) *(t(Qs_k[, nc])))
      
      ## Total and fractional total load transport rate and concentration of total load sediment
      Qt[nc] = Qs[nc]+Qb[nc]
      Qt_k[, nc] = Qs_k[, nc] + Qb_k[, nc]
      
      Ct[nc] = Qt[nc]/Q[nc]
      Ct_k[, nc] = t((1/Q[nc]) *(t(Qt_k[, nc])))
    }else
    {
      ## Total and fractional total load transport rate and concentration of suspended sediment
      Qt[nc] = Qtin[nc]
      Qt_k[, nc] = Qtin_k[,nc]
      
      Ct[nc] = Qt[nc]/Q[nc]
      Ct_k[, nc] = t((1/Q[nc]) *(t(Qt_k[, nc])))
      
      ## Total and fractional suspended load transport rate and concentration of suspended sediment
      Qs[nc] = Qt[nc] - Qb[nc]
      for(i in 1:k)
      {
        if(Qt_k[i,nc] > 0)
        {
          Qs_k[i, nc] = Qt_k[i, nc] - Qb_k[i, nc]
        }
        else
        {
          Qs_k[i,nc] = 0
        }
      }
      C[nc] = Qs[nc]/Q[nc]
      C_k[1:k, nc] = t((1/Q[nc])*(t(Qs_k[1:k, nc])))
    }
    ## Downstream Boundary
    ## Hydraulic sub-model
    z[1] = 2*z[2] - z[3]
    h[1] = z[1] - zc[1]
    A[1] = h[1] * b[1]
    P[1] = b[1] + (2*h[1])
    R[1] = A[1]/P[1]
    Q[1] = 2*Q[2] - Q[3]
    V[1] = Q[1]/A[1]
    
    ## sediment transport sub-model
    Qt[1] = Qt[2]
    Qt_k[1:k,1] = Qt_k[1:k,2]
    Qb[1] = Qb[2]
    Qb_k[1:k, 1] = Qb_k[1:k, 2]
    Qs[1] = Qs[2]
    Qs_k[1:k, 1] = Qs_k[1:k, 2]
    Ct[1] = Ct[2]
    Ct_k[1:k, 1] = Ct_k[1:k, 2]
    C[1] = C[2]
    C_k[1:k, 1] = C_k[1:k, 2]
    
    ## Update variables at internal points
    ## hydraulic sub-model
    h[2:nc-1] = h[2:nc-1]
    z[2:nc-1] = z[2:nc-1]
    A[2:nc-1] = A[2:nc-1]
    P[2:nc-1] = P[2:nc-1]
    R[2:nc-1] = R[2:nc-1]
    
    ## Sediment transport sub model
    Qt[2:nc-1] = Qt[2:nc-1]
    Qt_k[1:k, 2:nc-1] = Qt_k[1:k, 2:nc-1]
    Qb[2:nc-1] = Qb[2:nc-1]
    Qb_k[1:k, 2:nc-1] = Qb_k[1:k, 2:nc-1]
    Qs[2:nc-1] = Qs[2:nc-1]
    Qs_k[1:k, 2:nc-1] = Qs_k[1:k, 2:nc-1]
    Ct[2:nc-1] = Ct[2:nc-1]
    Ct_k[1:k, 2:nc-1] = Ct_k[1:k, 2:nc-1]
    C[2:nc-1] = C[2:nc-1]
    C_k[1:k, 2:nc-1] = C_k[1:k, 2:nc-1]
    
    ## Save friction slope
    Sf[-c(1,nc)] = ((nf[-c(1,nc)]*Q[-c(1,nc)])^2)/((R[-c(1,nc)]^(4/3))*(A[-c(1,nc)]^2))
    
    ## iteration 1
    if( iter ==1){
      Tau_star_k[1:k,1:nc] <- 0
      qb_star_k[1:k, 1:nc] <- 0
      
      # Fractional bedload transport rate
      for(j in 1:k)
      {
        # Dimensionless boundary shear stress (1)
        Tau_star_k[j,2:nc-1] = t((R[2:nc-1]*Sf[2:nc-1])/(t((Rs*dk[j,2:nc-1]))))
        # Dimensionless fractional bedload transport rate (1)
        for(i in 2:nc-1){
          if(Tau_star_k[j,i] >= Tau_star_c_k[j,i]){
            qb_star_k[j,i] = (17*(Tau_star_k[j,i]-Tau_star_c_k[j,i]))*(((Tau_star_k[j,i])^0.5)-(((Tau_star_c_k[j,i])^0.5)))
          }
          else{qb_star_k[j,i] = 0}
        }
        # Fractional bedload transport rate (m^3/s)
        Qb_k[j,2:nc-1] = qb_star_k[j,2:nc-1]*pk[j,2:nc-1]*dk[j,2:nc-1]*((Rs*g*dk[j,2:nc-1])^0.5)
      }
      
      Qb_k[,-c(1,nc)] <- t(b[-c(1,nc)]*t(Qb_k[,-c(1,nc)]))
      # Total bedload transport rate (m^3/s)
      Qb[2:nc-1] = colSums(Qb_k[,2:nc-1])
      
      # Boundary conditions
      Qb[nc] = Qb[nc]
      Qb[1] = Qb[1]
      Qb_k[1:k,nc] = Qb_k[1:k,nc]
      Qb_k[1:k,1] = Qb_k[1:k,1]
      
      ## Total and fractional total load transport rate (m3/s) and concentration of total load sediment
      Qt[-c(1,nc)] <- Qb[-c(1,nc)] + Qs[-c(1,nc)]
      Qt_k[, -c(1,nc)] <- Qb_k[,-c(1,nc)] + Qs_k[,-c(1,nc)]
      for(i in 2:nc-1){
        if (Q[i] == 0){
          Ct[i] = 0
        }else
        {
          Ct[i] = Qt[i]/Q[i]
        }
      }
      
      for(i in 2:nc-1){
        if (Q[i] == 0)
        {
          Ct_k[1:k,i] <- 0
        }
        else
        {
          Ct_k[1:k,i] = t((1/Q[i])*t(Qt_k[1:k,i]))}
        
      }
      
      ## Density of water-sediment mixture
      p[-c(1,nc)] = (pw*(1-Ct[-c(1,nc)]))+(ps*(Ct[-c(1,nc)]))
      
      ## Fractional settling velocity od sediment particles in sediment laden water
      for(j in 1:k){
        ws_k[j,-c(1,nc)] <- wso_k[j,-c(1,nc)]
      }
      
      ## Total and fractional entrainment rate per unit channel width
      ws_star_k <- matrix(0,k,nc)
      Cb_star_k <- matrix(0,k,nc)
      E_k <- matrix(0,k,nc)
      for(j in 1:k){
        
        ## Dimensionless settling velocity of sediment particles
        ws_star_k[j,2:nc-1] = ws_k[j,2:nc-1]/ ((Rs*g*dk[j,2:nc-1])^0.5)
        
        ## Equilibrium near bed suspended sediment concentration
        for(i in 2:nc-1){
          if(Tau_star_k[j,i] >= Tau_star_c_k[j,i]){
            Cb_star_k[j,i] = pk[j,i]*(0.001*(Tau_star_k[j,i]/ws_star_k[j,i])^2)
          }
          else {Cb_star_k[j,i] = 0}
        }
        ## Fractional entrainment rate per unit channel width
        E_k[j,2:nc-1] = Cb_star_k[j,2:nc-1] * ws_k[j,2:nc-1]
      }
      E[2:nc-1] <- colSums(E_k[,2:nc-1])
      
      ##Total and fractional deposition rate per unit channel width
      Cb_k <- matrix(0,k,nc)
      D_k <- matrix(0,k,nc)
      kk <- (1-lambdao)/(C)
      kl <- rep(2,nc)
      alpha = pmin(kl,kk)
      for(i in 2:nc-1)
      {
        if(C[i] > 0.2)
        {
          Cb_k[1:k,i] = C_k[1:k,i]
        }
        else{
          Cb_k[1:k,i] = t(alpha[i] * t(C_k[1:k,i]))
        } 
        
      }
      D_k[,2:nc-1] <- Cb_k[,2:nc-1] * ws_k[,2:nc-1] 
      D[-c(1,nc)] <- colSums(D_k[,-c(1,nc)])
      
      ## Total and fractional bedload transport rate gradient
      dQbdx[-c(1,nc)] = (Qb[-c(1,nc)]-Qb[-c(1,2)])/dx
      dQb_kdx[,-c(1,nc)] = (Qb_k[,-c(1,nc)] - Qb_k[,-c(1,2)])/dx
      
      ## Concentration of total transport gradient
      dCtdx[-c(1,nc)] = (Ct[-c(1,nc)] - Ct[-c(1,2)])/dx
    }
    
    ## varaiables at time step n
    # Hydraulic sub model
    hn = h 
    zn = z 
    An = A 
    Pn = P 
    Rn = R 
    Vn = V 
    Qn = Q
    
    # Sediment transport sub-model
    Qtn = Qt
    Qt_kn = Qt_k
    Qbn = Qb
    Qb_kn = Qb_k
    Qsn = Qs
    Qs_kn = Qs_k
    Ctn = Ct
    Ct_kn = Ct_k
    Cn = C
    C_kn = C_k
    
    # Morphological evolution sub model
    En = E
    E_kn = E_k
    Dn = D
    D_kn = D_k
    zcn = zc
    pkn = pk
    
    ## Gradients
    dQbndx = dQbdx
    dQb_kndx = dQb_kdx
    dCtndx = dCtdx
    
    ## Variables at initial time
    if (iter == 1){
      ts = dts
      count_ts = 1
      
      ## Hydraulic sub-model
      hit[count_ts,] = h
      zit[count_ts,] = z
      Qit[count_ts,] = Q
      Vit[count_ts,] = V
      
      ## Sediment transport sub model
      Qtit[count_ts,] = Qt
      Qbit[count_ts,] = Qb
      Qsit[count_ts,] = Qs
      Ctit[count_ts,] = Ct
      Cit[count_ts,] = C
      
      ## Morphological evolution model
      Eit[count_ts,] = E
      Dit[count_ts,] = D
      zcit[count_ts,] = zc
      pkit = pk
      
      ## Adaptive sub model
      Xiit[count_ts,] = Xi
    }  
    ## Check En values
    ## Hydraulic sub model
    Hc[a1] <- (-((Qn[a1]-Qn[b1])/dx))+(Psi[a1]*((1/(1-lambdao[a1]))*(dQbndx[a1]+b[a1]*(En[a1]-Dn[a1]))))
    A[a1] = An[a1]+(dt*Hc[a1])
    
    # Update variables at time step n+1
    # Water depth (m)
    h[-c(1,nc)] = A[-c(1,nc)]/b[-c(1,nc)]
    
    # Water surface elevation (m)
    z[-c(1,nc)] = zc[-c(1,nc)]+h[-c(1,nc)]
    
    # Wetted perimeter (m)
    P[-c(1,nc)] = b[-c(1,nc)]+2*h[-c(1,nc)]
    
    # Hydraulic Radius (m)
    R[-c(1,nc)] = A[-c(1,nc)]/P[-c(1,nc)]
    
    ########################### Momentum Equation ##########################
    ## Weighting factors based on the Courant (Cr) number (1)
    
    w1[a1] = 1-((dt/dx)*(Vn[c1]+Vn[a1])/2)
    w2[a1] = (dt/dx)*((Vn[a1]+Vn[b1])/2)
    
    ## Water surface elevation gradients (1)
    dzdx_down[a1] = (z[c1]-z[a1])/dx
    dzdx_up[a1] = (z[a1]-z[b1])/dx
    
    ## Source terms for clear-water flow
    Source_cwf[a1] = -(g*A[a1]*((w1[a1]*dzdx_down[a1])+(w2[a1]*dzdx_up[a1])))-(g*An[a1]*Sf[a1])
    
    ## Source terms for sediment laden-flow
    Source_slf[a1] <- -(g*(0.5*b[a1]*((hn[a1])^2))*((ps-pw)/p[a1])*dCtndx[a1])-(Vn[a1]*((po[a1]-p[a1])/p[a1])*((1/(1-lambdao[a1]))*(dQbndx[a1]+b[a1]*(En[a1]-Dn[a1]))))
    
    ## Discharge at time step n+1 (m^3/s)
    Hm[a1] = -((((Qn[a1]^2)/An[a1])-((Qn[b1]^2)/An[b1]))/dx)+Source_cwf[a1]+(Psi[a1]*Source_slf[a1])
    
    ### simuating Q
    Q[a1] = Qn[a1]+(dt*Hm[a1])
    
    ## Update variables at time step n+1
    ## Flow velocity (m/s)
    V[-c(1,nc)] = Q[-c(1,nc)]/A[-c(1,nc)]
    
    ## Sediment transport sub model
    Hsl_k[,a1] <- (-((t(Qn[a1]*t(C_kn[1:k, a1]))) - (t(Qn[b1]*t(C_kn[1:k,b1]))))/dx) + (t(b[a1]*t((E_kn[1:k,a1] - D_kn[1:k,a1]))))
    AC_k[1:k,a1] <- (t(An[a1]*t(C_kn[1:k,a1]))) + (dt*Hsl_k[1:k, a1])
    
    # Update variables at time step n+1
    # Fractional concentration of suspended sediment (1)
    C_k[1:k,2:nc-1] = t((1/A[2:nc-1]) * t(AC_k[1:k,2:nc-1]))
    
    # Total concentration of suspended sediment (1)
    C[-c(1,nc)] = colSums(C_k[,-c(1,nc)])
    
    # Total and fractional suspended load transport rate (m^3/s)
    Qs[-c(1,nc)] = C[-c(1,nc)]*Q[-c(1,nc)]
    Qs_k[1:k, 2:nc-1] = t(Qs[2:nc-1]*t(pk[1:k, 2:nc-1]))
    
    # Total and fractional bedload transport rate (m^3/s)
    Tau_star_k[1:k,1:nc] <- 0
    qb_star_k[1:k, 1:nc] <- 0
    
    # Fractional bedload transport rate
    for(j in 1:k)
    {
      # Dimensionless boundary shear stress (1)
      Tau_star_k[j,2:nc-1] = t((R[2:nc-1]*Sf[2:nc-1])/(t((Rs*dk[j,2:nc-1]))))
      # Dimensionless fractional bedload transport rate (1)
      for(i in 2:nc-1){
        if(Tau_star_k[j,i] >= Tau_star_c_k[j,i]){
          qb_star_k[j,i] = (17*(Tau_star_k[j,i]-Tau_star_c_k[j,i]))*(((Tau_star_k[j,i])^0.5)-(((Tau_star_c_k[j,i])^0.5)))
        }
        else{qb_star_k[j,i] = 0}
      }
      # Fractional bedload transport rate (m^3/s)
      Qb_k[j,2:nc-1] = qb_star_k[j,2:nc-1]*pk[j,2:nc-1]*dk[j,2:nc-1]*((Rs*g*dk[j,2:nc-1])^0.5)
    }
    
    Qb_k[,-c(1,nc)] <- t(b[-c(1,nc)]*t(Qb_k[,-c(1,nc)]))
    # Total bedload transport rate (m^3/s)
    Qb[2:nc-1] = colSums(Qb_k[,2:nc-1])
    
    # Boundary conditions
    Qb[nc] = Qb[nc]
    Qb[1] = Qb[1]
    Qb_k[1:k,nc] = Qb_k[1:k,nc]
    Qb_k[1:k,1] = Qb_k[1:k,1]
    
    # Total and fractional total load transport rate (m^3/s)& concentration of total load sediment (1)
    Qt[-c(1,nc)] = Qb[-c(1,nc)]+Qs[-c(1,nc)]
    Qt_k[,-c(1,nc)] = Qb_k[,-c(1,nc)]+Qs_k[,-c(1,nc)]
    
    for(i in 2:nc-1){
      if(Q[i] == 0){
        Ct[i] = 0}
      else{Ct[i] = Qt[i]/Q[i]
      }
    }
    
    for(i in 2:nc-1){
      if(Q[i] == 0){
        Ct_k[,i] = 0}
      else{
        for(j in 1:k){Ct_k[j,i] = Qt_k[j,i]/Q[i]
        }
      }
    }
    
    # Density of water-sediment mixture (kg/m^3)
    p[-c(1,nc)] = pw*(1-Ct[-c(1,nc)])+ps*Ct[-c(1,nc)]
    
    # Fractional settling velocity of sediment particles in sediment-laden water (m/s)
    ws_k[1:k,2:nc-1] <- wso_k[1:k,2:nc-1]
    
    # Total and fractional entrainment rate per unit channel width (m/s)
    ws_star_k <- matrix(0,k,nc)
    Cb_star_k <- matrix(0,k,nc)
    E_k <- matrix(0,k,nc)
    for(j in 1:k){
      
      ## Dimensionless settling velocity of sediment particles
      ws_star_k[j,2:nc-1] = ws_k[j,2:nc-1]/ ((Rs*g*dk[j,2:nc-1])^0.5)
      
      ## Equilibrium near bed suspended sediment concentration
      for(i in 2:nc-1){
        if(Tau_star_k[j,i] >= Tau_star_c_k[j,i]){
          Cb_star_k[j,i] = pk[j,i]*(0.001*(Tau_star_k[j,i]/ws_star_k[j,i])^2)
        }
        else {Cb_star_k[j,i] = 0}
      }
      ## Fractional entrainment rate per unit channel width
      E_k[j,2:nc-1] = Cb_star_k[j,2:nc-1] * ws_k[j,2:nc-1]
    }
    E[2:nc-1] <- colSums(E_k[,2:nc-1])
    
    # Total and fractional deposition rate per unit channel width (m/s)
    Cb_k <- matrix(0,k,nc)
    D_k <- matrix(0,k,nc)
    kk <- (1-lambdao)/(C)
    kl <- rep(2,nc)
    alpha = pmin(kl,kk)
    for(i in 2:nc-1)
    {
      if(C[i] > 0.2)
      {
        Cb_k[1:k,i] = C_k[1:k,i]
      }
      else{
        Cb_k[1:k,i] = t(alpha[i] * t(C_k[1:k,i]))
      } 
      
    }
    D_k[,2:nc-1] <- Cb_k[,2:nc-1] * ws_k[,2:nc-1] 
    D[-c(1,nc)] <- colSums(D_k[,-c(1,nc)])
    
    # Total and fractional bedload transport rate gradient
    dQbdx[-c(1,nc)] <- (Qb[-c(1,nc)]-Qb[-c(1,2)])/dx
    dQb_kdx[,-c(1,nc)] <- (Qb_k[,-c(1,nc)]-Qb_k[,-c(1,2)])/dx
    
    # Concentration of total load transport gradient
    dCtdx[-c(1,nc)] <- (Ct[-c(1,nc)]-Ct[-c(1,2)])/dx
    
    ############# MORPHOLOGICAL EVOLUTION SUB-MODEL ###################
    ############### Exner Equation ##################
    # Fractional riverbed deformation (m^2)
    dAb_k[1:k,-c(1,nc)]=t((dt/(1-lambdao[-c(1,nc)]))*(t((-dQb_kndx[1:k,-c(1,nc)]+(t(b[-c(1,nc)]*(t(D_kn[1:k,-c(1,nc)]-E_kn[1:k,-c(1,nc)]))))))))
    
    # Total riverbed deformation (m^2)
    dAb[-c(1,nc)] = colSums(dAb_k[,-c(1,nc)])
    
    # Change in riverbed elevation (m)
    dzc[-c(1,nc)] = dAb[-c(1,nc)]/b[-c(1,nc)]
    zc[-c(1,nc)] = zcn[-c(1,nc)]+dzc[-c(1,nc)]
    
    # Adjust water surface elevation
    z[-c(1,nc)] = zc[-c(1,nc)]+h[-c(1,nc)]
    
    # Adjust riverbed elevation at ghost cells
    z[nc] = 2*z[nc-1]-z[nc-2]
    zc[nc] = z[nc]-h[nc-1]
    zc[1] = 2*zc[2]-zc[3]
    
    ######################### Riverbed Gradation ##########################
    # Riverbed mixing layer thickness (m)## removed 1000*d90
    delta_m = rep(1,nc)
    for(i in 2:nc-1){
      delta_m[i] = max(0.15*h[i], 2*d90_xc[i])
    }
    
    # Riverbed mixing layer cross-sectional area (m^2)
    A_m[-c(1,nc)] = b[-c(1,nc)]*delta_m[-c(1,nc)]
    
    # Fraction kth of sediment of size dk at time step n+1
    pk_star = matrix(0,k,nc)
    for(i in 2:nc-1){
      if(dAb[i] >= 0){
        pk_star[,i] = pk[,i]}
      else{pk_star[,i] = pk_sub[,i]
      }
    }
    pk[,-c(1,nc)] = t((1/A_m[-c(1,nc)])*t(dAb_k[,-c(1,nc)]-(t(dAb[-c(1,nc)]*(t(pk_star[,-c(1,nc)]))))))+pkn[,-c(1,nc)]
    for(j in 1:k){
      for(i in 2:nc-1){
        if(pk[j,i] < 0){
          pk[j,i] = 0}
      }
    }
    pk[,1] = pk[,2]
    pk[,nc] = pk[,nc-1]
    
    # Updated grain size distribution
    d50 <- rep(0,nc)
    d90 <- rep(0,nc)
    d65 <- rep(0,nc)
    d84 <- rep(0,nc)
    pk[pk < 0.0000001] <- 0.000001
    # Limits and representative grain sizes considered for each sediment class
    #dk_limits <- c(0.062, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16)
    dk_limits <- c(0.125,0.5,1,4,8,16,32,64)
    #dk_limits <- c(0.062, 1,2,3,4,8,16,32,64)
    pk_limits <- matrix(0,8,nc)
    pk_limits[1,] <- 0
    pk_limits[2,] <- pk[1,]+pk[2,]
    
    sum <- matrix(0,k,nc)
    sum[2,] <- pk_limits[2,]
    for(i in 3:k-1){
      pk_limits[i,] = sum[i-1,] + pk[i,]
      sum[i,] = pk_limits[i,]
    }
    pk_limits <- pk_limits*100
    
    for (i in 2:(nc - 1)) {
      d50[i] <- pchip(pk_limits[, i], dk_limits, 50) / 1000
      d65[i] <- pchip(pk_limits[, i], dk_limits, 65) / 1000
      d84[i] <- pchip(pk_limits[, i], dk_limits, 84) / 1000
      d90[i] <- pchip(pk_limits[, i], dk_limits, 90) / 1000
    }
    d50_xc <- d50
    d90_xc <- d90
    
    # Updated porosity of the riverbed surface layer
    lambdao[-c(1,nc)] = (0.13+0.21)/((d50_xc[-c(1,nc)]*1000+0.002)^0.21)
    lambdao[1] = lambdao[2]
    lambdao[nc] = lambdao[nc-1]
    
    ############# ADAPTIVE SUB-MODEL ###################
    # Froude number
    Fr[-c(1,nc)] = V[-c(1,nc)]/((g*h[-c(1,nc)])^0.5)
    
    # Dimensionless parameter Xi
    Xi[-c(1,nc)] = ((Ct[-c(1,nc)]*Rs)^0.5)/Fr[-c(1,nc)]
    
    # Adaptive numerical scheme
    for(i in 2:nc-1){
      if(Xi[i] >= Xi_cr[i]){
        Psi[i] = 1}
      else{Psi[i] = 0
      }
    }
    
    ## Hydraulic sub model
    # Average L2-norm of the water depth and flow rate for evaluating the change in the solution
    sum_h=0
    sum_Q=0;
    for(i in 2:nc-1){
      sum_h=sum_h+(h[i]-hn[i])^2
      sum_Q=sum_Q+(Q[i]-Qn[i])^2
    }
    L2_norm_h=((1/(nc-2))*sum_h)^(1/2)
    L2_norm_Q=((1/(nc-2))*sum_Q)^(1/2)
    
    ## saving the  solution 
    if (t > ts - dt){
      ts = ts + dts
      count_ts = count_ts + 1
      
      Cit = rbind(Cit, matrix(data=0,nrow = 1, ncol = ncol(Cit)))
      Cit[count_ts,] = C
      
    }
  }
  lt <- nrow(Cit)
  lt2 <- Cit[lt,dgo]
  observed <- SSC
  rmse <- sqrt(mean((lt2 - observed)^2))
  return(rmse) 
}

library(mcmc)

# MCMC-based optimization using Metropolis-Hastings algorithm with adaptive scaling
mcmc_optimization <- function(Num_model2, param_range, n_iter, initial_proposal_sd, scaling_factor) {
  # Function to calculate acceptance probability
  acceptance_prob <- function(curr_rmse, new_rmse) {
    return(min(1, exp(curr_rmse - new_rmse)))
  }
  
  # Initialize optimization variables
  best_param <- 0
  best_rmse <- Inf
  
  curr_param <- SSC[1]
  curr_rmse <- best_rmse
  
  proposal_sd <- initial_proposal_sd
  
  # Generate initial parameter value within the specified range
  curr_param <- runif(1, param_range[1], param_range[2])
  
  # Perform MCMC iterations
  for (i in 1:n_iter) {
    # Generate proposal parameter
    proposal_param <- rnorm(1, mean = curr_param, sd = proposal_sd)
    
    # Ensure proposal parameter is within the specified range
    proposal_param <- pmax(proposal_param, param_range[1])
    proposal_param <- pmin(proposal_param, param_range[2])
    
    # Calculate RMSE for the proposal parameter
    proposal_rmse <- Num_model2(proposal_param)
    
    # Accept or reject proposal based on acceptance probability
    if (runif(1) < acceptance_prob(curr_rmse, proposal_rmse)) {
      curr_param <- proposal_param
      curr_rmse <- proposal_rmse
    }
    
    # Update best parameter and RMSE
    if (curr_rmse < best_rmse) {
      best_param <- curr_param
      best_rmse <- curr_rmse
    }
    
    # Adaptively adjust proposal_sd based on acceptance rate
    acceptance_rate <- i / (i + 1)
    proposal_sd <- proposal_sd * (1 + scaling_factor * (acceptance_rate - 0.234))
  }
  
  # Return best parameter and RMSE
  return(list("best_param" = best_param, "best_rmse" = best_rmse))
}

# Define parameter range (min and max values)
param_range <- range(SSC)  # Range for the parameter

# Set number of iterations, initial proposal standard deviation, and scaling factor
n_iter <- 50          # Number of iterations
initial_proposal_sd <- 0.5   # Initial proposal standard deviation
scaling_factor <- 0.5   # Scaling factor for adaptive scaling

# Run MCMC-based optimization with adaptive scaling 11:18 to 2:53 (12506 sec)
tic()
result <- mcmc_optimization(Num_model2, param_range, n_iter, initial_proposal_sd, scaling_factor)
toc()
# Print the best parameter and RMSE
par1 <- result$best_param
rmse1 <- result$best_rmse

Num_1 <- data.frame(Css = par1, RMSE = rmse1)

write.csv(Num_1, file = "N_1.csv")
