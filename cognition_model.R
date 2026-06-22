require(deSolve)
require(tidyverse)
require(latex2exp)
require(cowplot)
require(scales)


odes <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    beta <- R0 * mu * (lambda + mu) * (alpha + gamma + mu) / (Gamma * lambda)
    
    dS_o <- Gamma - beta * ifelse(t > Tp, 1 - p1, 1) * S_o * I_o - mu * S_o
    dE_o <- beta * ifelse(t > Tp, 1 - p1, 1) * S_o * I_o - lambda * E_o - mu * E_o
    dI_o <- lambda * E_o - alpha * I_o - gamma * I_o - ifelse(t > Tp, p2 * I_o, 0) - mu * I_o
    dD_o <- alpha * I_o
    dR_o <- gamma * I_o + gammaQ * Q_o - mu * R_o
    dQ_o <- ifelse(t > Tp, p2 * I_o, 0) - gammaQ * Q_o
    dS <- Gamma - beta * (1 - eB * B - ifelse(t > Tp, p1, 0)) * S * I - mu * S
    dE <- beta * (1 - eB * B - ifelse(t > Tp, p1, 0)) * S * I - lambda * E - mu * E
    dI <- lambda * E - alpha * I - gamma * I - (kB * B + kC * C + ifelse(t > Tp, p2, 0)) * I - mu * I
    dD <- alpha * I + alphaQ * Q
    dR <- gamma * I + gammaQ * Q - mu * R
    dQ <- (kB * B + kC * C + ifelse(t > Tp, p2, 0)) * I - gammaQ * Q - alphaQ * Q - mu * Q
    dC <- ((cB1 - cB2 * B) * B + cD * D + cQ * Q + ifelse(t > Tp, p3, 0)) * (1 - C) - deltaC * C
    dB <- (bC * C + ifelse(t > Tp, p4, 0)) * (1 - B) - deltaB * B
    
    return(list(c(dS_o, dE_o, dI_o, dD_o, dR_o, dQ_o,
                  dS, dE, dI, dD, dR, dQ, dC, dB)))
  })
}

N <- 10^5 #population size
init_infections <- 10
y0 <- c(S_o=N-init_infections, E_o=init_infections, I_o=0, D_o=0, R_o=0, Q_o=0,
        S=N-init_infections, E=init_infections, I=0, D=0, R=0, Q=0, C=0, B=0)

parms_no_policy <- c(mu=1/(60*365), Gamma=10^5/(60*365), lambda=.3, gamma=.1, gammaQ=.1, alpha=.01, alphaQ=.01, delta=0, Tp=0,
                     R0=2, eB=.5, kB=0, kC=.05, deltaC=.2, cB1=.1, cB2=0,
                     cD=1/10^5, cQ=1/10^4, deltaB=.1, bC=.1,
                     p1=0, p2=0, p3=0, p4=0)
ts <- seq(0,3650,.1)
one_year <- which.min(abs(ts-365))


### NO POLICY

out_no_policy <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_no_policy)

out <- out_no_policy[1:one_year,]

I_max_t <- ts[which.max(out[,'I'])]

plot1 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('I','I_o','D','D_o','Q')) %>% 
  mutate(var=factor(var, levels=c('I','I_o','D','D_o','Q')),
         value=ifelse(var %in% c('I','I_o','D','D_o','Q'), value/12000, value)) %>% 
  ggplot() +
  ggtitle('lower risk perception (default)') +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value*12000, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(name='number of people', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('blue','blue','red','red','violet'), labels=c('infected population, I','I with C=B=0', 'cumulative deaths, D', 'D with C=B=0','quarantined population, Q')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12,'solid'), labels=c('infected population, I','I with C=B=0', 'cumulative deaths, D', 'D with C=B=0','quarantined population, Q')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust = 0.5))




## breakdown of effective reproduction number Re

plot2 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('C','B')) %>% 
  mutate(var=factor(var, levels=c('C','B'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(limits=c(0,1), expand=0) +
  scale_color_manual(labels=c('risk perception, C','NPI adoption, B'),
                     values=c('forestgreen','darkorange1')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


Re <-
  out[,'S']/N * parms_no_policy['R0'] * (1-parms_no_policy['eB'] * out[,'B']) *
  (parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu']) / (parms_no_policy['kC'] * out[,'C'] + parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu'])
impact_of_susceptibles <- out[,'S']/N
impact_of_behavior <- 1 - parms_no_policy['eB'] * out[,'B']
impact_of_testing <- (parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu']) / (parms_no_policy['kC'] * out[,'C'] + parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu'])


plot3 <- data.frame(value=c(Re, impact_of_susceptibles, impact_of_behavior, impact_of_testing),
                    time=ts[1:one_year],
                    var=rep(factor(c('Re','impact_of_susceptibles','impact_of_behavior','impact_of_testing'),levels=c('Re','impact_of_susceptibles','impact_of_behavior','impact_of_testing')), each=one_year)) %>% 
  ggplot() +
  geom_hline(yintercept=1, color='gray85', linewidth=.5) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(limits=c(0,2), expand=0) +
  scale_color_manual(labels=c(expression('effective reproduction,'~R[e]),expression('reduction in'~R[e]~' due to susceptible depletion'),expression('reduction in'~R[e]~' due to NPI adoption'),expression('reduction in'~R[e]~' due to testing')),
                     values=c('black','gray70','#ffc178','#a9e78f')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position="none",
        axis.title.y=element_blank())








parms_no_policy_high_risk_perception <- parms_no_policy
parms_no_policy_high_risk_perception['cD'] <- parms_no_policy['cD'] * 10
parms_no_policy_high_risk_perception['cQ'] <- parms_no_policy['cQ'] * 10

out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_no_policy_high_risk_perception)[1:one_year,]

I_max_t <- ts[which.max(out[,'I'])]

plot4 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('I','I_o','D','D_o','Q')) %>% 
  mutate(var=factor(var, levels=c('I','I_o','D','D_o','Q')),
         value=ifelse(var %in% c('I','I_o','D','D_o','Q'), value/12000, value)) %>% 
  ggplot() +
  ggtitle('higher risk perception') +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value*12000, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(name='number of people', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('blue','blue','red','red','violet'), labels=c('infected population, I','I with C=B=0', 'cumulative deaths, D', 'D with C=B=0','quarantined population, Q')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12,'solid'), labels=c('infected population, I','I with C=B=0', 'cumulative deaths, D', 'D with C=B=0','quarantined population, Q')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, .3, .3, 0), 'cm'),
        legend.text=element_text(size=10),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust = 0.5))




## breakdown of effective reproduction number Re

plot5 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('C','B')) %>% 
  mutate(var=factor(var, levels=c('C','B'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(limits=c(0,1), expand=0) +
  scale_color_manual(labels=c('risk perception, C','NPI adoption, B'),
                     values=c('forestgreen','darkorange1')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, .3, .3, 0), 'cm'),
        legend.text=element_text(size=10),
        legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())


Re <-
  out[,'S']/N * parms_no_policy['R0'] * (1-parms_no_policy['eB'] * out[,'B']) *
  (parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu']) / (parms_no_policy['kC'] * out[,'C'] + parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu'])
impact_of_susceptibles <- out[,'S']/N
impact_of_behavior <- 1 - parms_no_policy['eB'] * out[,'B']
impact_of_testing <- (parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu']) / (parms_no_policy['kC'] * out[,'C'] + parms_no_policy['alpha'] + parms_no_policy['gamma'] + parms_no_policy['mu'])


plot6 <- data.frame(value=c(Re, impact_of_susceptibles, impact_of_behavior, impact_of_testing),
                    time=ts[1:one_year],
                    var=rep(factor(c('Re','impact_of_susceptibles','impact_of_behavior','impact_of_testing'),levels=c('Re','impact_of_susceptibles','impact_of_behavior','impact_of_testing')), each=one_year)) %>% 
  ggplot() +
  geom_hline(yintercept=1, color='gray85', linewidth=.5) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(limits=c(0,2), expand=0) +
  scale_color_manual(labels=c(expression('effective reproduction,'~R[e]),expression('reduction in'~R[e]~' due to susceptible depletion'),expression('reduction in'~R[e]~' due to NPI adoption'),expression('reduction in'~R[e]~' due to testing')),
                     values=c('black','gray70','#ffc178','#a9e78f')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, .3, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        legend.title=element_blank(),
        axis.title.y=element_blank())




plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2, byrow=FALSE, align = 'vh', rel_heights=c(1,1,.85))
# 13x7




### no policy heatmaps



eB_kC <- data.frame()

parms <- parms_no_policy
for (eB in seq(0, parms_no_policy['eB']*2-parms_no_policy['eB']/50, parms_no_policy['eB']/50)) {
  print(eB)
  parms['eB'] <- eB
  for (kC in seq(0, parms_no_policy['kC']*2-parms_no_policy['kC']/50, parms_no_policy['kC']/50)) {
    parms['kC'] <- kC
    
    out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms)
    
    eB_kC <-
      eB_kC %>% 
      rbind(data.frame(cD = round(cD,7),
                       cQ = round(cQ,6),
                       deaths = unname(tail(out[,'D'],1)),
                       #I_max_t = ts[which.max(out[,'I'])],
                       #I_max = max(out[,'I']),
                       #B_eq = tail(out[,'B'],1),
                       B_max = max(out[,'B'])))
  }
}



cD_cQ <- data.frame()

parms <- parms_no_policy
for (cD in seq(0, parms_no_policy['cD']*2-parms_no_policy['cD']/50, parms_no_policy['cD']/50)) {
  print(cD)
  parms['cD'] <- cD
  for (cQ in seq(0, parms_no_policy['cQ']*2-parms_no_policy['cQ']/50, parms_no_policy['cQ']/50)) {
    parms['cQ'] <- cQ
    
    out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms)
    
    cD_cQ <-
      cD_cQ %>% 
      rbind(data.frame(cD = round(cD,7),
                       cQ = round(cQ,6),
                       deaths = unname(tail(out[,'D'],1)),
                       #I_max_t = ts[which.max(out[,'I'])],
                       #I_max = max(out[,'I']),
                       #B_eq = tail(out[,'B'],1),
                       B_max = max(out[,'B'])))
  }
}



plot1 <- eB_kC %>% 
  ggplot() +
  ggtitle("varying how \nC & B affect epi variables") +
  geom_raster(aes(x=eB,y=kC,fill=deaths), hjust=1, vjust=1) +
  scale_fill_gradientn(name='deaths after\nfirst wave',colors=rev(rainbow(7)[-7]),
                       limits=c(min(eB_kC$deaths, cD_cQ$deaths), max(eB_kC$deaths, cD_cQ$deaths)),
                       guide=guide_colorbar(frame.colour = 'black', ticks.colour = 'black',title.hjust=.5)) +
  scale_x_continuous(name=expression('effect of B on transmission,'~e[B]),expand=0) +
  scale_y_continuous(name=expression('effect of C on testing,'~k[C]),expand=0) +
  coord_fixed(10) +
  theme(axis.line=element_line(),
        plot.margin=unit(c(0,0,0,0),'cm'),
        plot.title=element_text(hjust=.5),
        legend.title=element_text(size=10))

plot2 <- eB_kC %>% 
  ggplot() +
  geom_raster(aes(x=eB,y=kC,fill=B_max), hjust=1, vjust=1) +
  scale_fill_gradientn(name='max value\nof B',colors=rev(rainbow(7)[-7]),
                       limits=c(min(eB_kC$B_max, cD_cQ$B_max), max(eB_kC$B_max, cD_cQ$B_max)),
                       guide=guide_colorbar(frame.colour = 'black', ticks.colour = 'black',title.hjust=.5)) +
  scale_x_continuous(name=expression('effect of B on transmission,'~e[B]),expand=0) +
  scale_y_continuous(name=expression('effect of C on testing,'~k[C]),expand=0) +
  coord_fixed(10) +
  theme(axis.line=element_line(),
        plot.margin=unit(c(0,0,0,0),'cm'),
        legend.title=element_text(size=10))



plot3 <- cD_cQ %>% 
  ggplot() +
  ggtitle("varying how\nC affects epi states") +
  geom_raster(aes(x=cD,y=cQ,fill=deaths), hjust=1, vjust=1) +
  scale_fill_gradientn(colors=rev(rainbow(7)[-7]),
                       limits=c(min(eB_kC$deaths, cD_cQ$deaths), max(eB_kC$deaths, cD_cQ$deaths))) +
  scale_x_continuous(name=expression('effect of D on C,'~c[D]),expand=0) +
  scale_y_continuous(name=expression('effect of Q on C,'~c[Q]),expand=0, limits=c(0,.0002)) +
  coord_fixed(.1) +
  theme(axis.line=element_line(),
        plot.margin=unit(c(0,0,0,0),'cm'),
        plot.title=element_text(hjust=.5),
        legend.position="none")

plot4 <- cD_cQ %>% 
  ggplot() +
  geom_raster(aes(x=cD,y=cQ,fill=B_max), hjust=1, vjust=1) +
  scale_fill_gradientn(colors=rev(rainbow(7)[-7]),
                       limits=c(min(eB_kC$B_max, cD_cQ$B_max), max(eB_kC$B_max, cD_cQ$B_max))) +
  scale_x_continuous(name=expression('effect of D on C,'~c[D]),expand=0) +
  scale_y_continuous(name=expression('effect of Q on C,'~c[Q]),expand=0, limits=c(0,.0002)) +
  coord_fixed(.1) +
  theme(axis.line=element_line(),
        plot.margin=unit(c(0,0,0,0),'cm'),
        legend.position="none")



plot_grid(plot1, plot3, plot2, plot4, nrow = 2, align = 'hv')
# save as 9.5 x 7.5









### POLICY

parms_policy <- function(p1=0, p2=0, p3=0, p4=0, Tp=0) {
  parms <- parms_no_policy
  parms['p1'] <- p1
  parms['p2'] <- p2
  parms['p3'] <- p3
  parms['p4'] <- p4
  parms['Tp'] <- Tp
  
  return(parms)
}

deaths_first_wave <- function(p1=0, p2=0, p3=0, p4=0, Tp=0) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1, p2=p2, p3=p3, p4=p4, Tp=Tp))[,'D'],1)))
}

deaths_no_policy <- deaths_first_wave()

deaths_red_p1 <- function(red, p1=0, p2=0, p3=0, p4=0, Tp=0) {
  if (p1+p2+p3+p4 == 0) {
    deaths_baseline <- deaths_no_policy
  } else {
    deaths_baseline <- deaths_first_wave(p1=p1, p2=p2, p3=p3, p4=p4, Tp=Tp)
  }
  
  return(uniroot(function(new_p1) deaths_first_wave(p1=new_p1, p2=p2, p3=p3, p4=p4, Tp=Tp) - (1-red) * deaths_baseline, c(0, 1-parms_no_policy['eB']))$root)
}

deaths_red_p2 <- function(red, p1=0, p2=0, p3=0, p4=0, Tp=0) {
  if (p1+p2+p3+p4 == 0) {
    deaths_baseline <- deaths_no_policy
  } else {
    deaths_baseline <- deaths_first_wave(p1=p1, p2=p2, p3=p3, p4=p4, Tp=Tp)
  }
  
  return(uniroot(function(new_p2) deaths_first_wave(p1=p1, p2=new_p2, p3=p3, p4=p4, Tp=Tp) - (1-red) * deaths_baseline, c(0, 10))$root)
}

deaths_red_p3 <- function(red, p1=0, p2=0, p3=0, p4=0, Tp=0) {
  if (p1+p2+p3+p4 == 0) {
    deaths_baseline <- deaths_no_policy
  } else {
    deaths_baseline <- deaths_first_wave(p1=p1, p2=p2, p3=p3, p4=p4, Tp=Tp)
  }
  
  return(uniroot(function(new_p3) deaths_first_wave(p1=p1, p2=p2, p3=new_p3, p4=p4, Tp=Tp) - (1-red) * deaths_baseline, c(0, 10))$root)
}

deaths_red_p4 <- function(red, p1=0, p2=0, p3=0, p4=0, Tp=0) {
  if (p1+p2+p3+p4 == 0) {
    deaths_baseline <- deaths_no_policy
  } else {
    deaths_baseline <- deaths_first_wave(p1=p1, p2=p2, p3=p3, p4=p4, Tp=Tp)
  }
  
  return(uniroot(function(new_p4) deaths_first_wave(p1=p1, p2=p2, p3=p3, p4=new_p4, Tp=Tp) - (1-red) * deaths_baseline, c(0, 10))$root)
}


I_max_t <- function(p1=0,p2=0,p3=0,p4=0,Tp=0) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1,p2=p2,p3=p3,p4=p4,Tp=Tp))[,'I'])])
}
I_max <- function(p1=0,p2=0,p3=0,p4=0,Tp=0) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1,p2=p2,p3=p3,p4=p4,Tp=Tp))[,'I']))
}

B_max <- function(p1=0,p2=0,p3=0,p4=0,Tp=0) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1,p2=p2,p3=p3,p4=p4,Tp=Tp))[,'B']))
}



### effects of each policy individually

I_max_t_p1s <- c()
I_max_t_p2s <- c()
I_max_t_p3s <- c()
I_max_t_p4s <- c()

I_max_p1s <- c()
I_max_p2s <- c()
I_max_p3s <- c()
I_max_p4s <- c()

B_max_p1s <- c()
B_max_p2s <- c()
B_max_p3s <- c()
B_max_p4s <- c()

deaths_p1s <- c()
deaths_p2s <- c()
deaths_p3s <- c()
deaths_p4s <- c()

p1s <- c()
p2s <- c()
p3s <- c()
p4s <- c()

reds <- seq(0,.6,.01)

for (red in reds) {
  print(red)
  
  p1 <- deaths_red_p1(red)
  p2 <- deaths_red_p2(red)
  p3 <- deaths_red_p3(red)
  p4 <- deaths_red_p4(red)
  
  p1s <- c(p1s, p1)
  p2s <- c(p2s, p2)
  p3s <- c(p3s, p3)
  p4s <- c(p4s, p4)
  
  I_max_t_p1s <- c(I_max_t_p1s, I_max_t(p1=p1))
  I_max_t_p2s <- c(I_max_t_p2s, I_max_t(p2=p2))
  I_max_t_p3s <- c(I_max_t_p3s, I_max_t(p3=p3))
  I_max_t_p4s <- c(I_max_t_p4s, I_max_t(p4=p4))
  
  I_max_p1s <- c(I_max_p1s, I_max(p1=p1))
  I_max_p2s <- c(I_max_p2s, I_max(p2=p2))
  I_max_p3s <- c(I_max_p3s, I_max(p3=p3))
  I_max_p4s <- c(I_max_p4s, I_max(p4=p4))

  B_max_p1s <- c(B_max_p1s, B_max(p1=p1))
  B_max_p2s <- c(B_max_p2s, B_max(p2=p2))
  B_max_p3s <- c(B_max_p3s, B_max(p3=p3))
  B_max_p4s <- c(B_max_p4s, B_max(p4=p4))
    
  deaths_p1s <- c(deaths_p1s, deaths_first_wave(p1=p1))
  deaths_p2s <- c(deaths_p2s, deaths_first_wave(p2=p2))
  deaths_p3s <- c(deaths_p3s, deaths_first_wave(p3=p3))
  deaths_p4s <- c(deaths_p4s, deaths_first_wave(p4=p4))
}

p1s_no_other_policies <- p1s
p2s_no_other_policies <- p2s
p3s_no_other_policies <- p3s
p4s_no_other_policies <- p4s

I_max_t_p1s_no_other_policies <- I_max_t_p1s
I_max_t_p2s_no_other_policies <- I_max_t_p2s
I_max_t_p3s_no_other_policies <- I_max_t_p3s
I_max_t_p4s_no_other_policies <- I_max_t_p4s

I_max_p1s_no_other_policies <- I_max_p1s
I_max_p2s_no_other_policies <- I_max_p2s
I_max_p3s_no_other_policies <- I_max_p3s
I_max_p4s_no_other_policies <- I_max_p4s

B_max_p1s_no_other_policies <- B_max_p1s
B_max_p2s_no_other_policies <- B_max_p2s
B_max_p3s_no_other_policies <- B_max_p3s
B_max_p4s_no_other_policies <- B_max_p4s

deaths_p1s_no_other_policies <- deaths_p1s
deaths_p2s_no_other_policies <- deaths_p2s
deaths_p3s_no_other_policies <- deaths_p3s
deaths_p4s_no_other_policies <- deaths_p4s

plot1 <- data.frame(red=reds,
                    I_max_t=c(I_max_t_p1s_no_other_policies, I_max_t_p2s_no_other_policies, I_max_t_p3s_no_other_policies, I_max_t_p4s_no_other_policies),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  ggtitle('no policy effects\non at baseline') +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='peak time of I (days)', expand=0, limits=c(0,550)) +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust=.5))


plot2 <- data.frame(red=reds,
                    I_max_t=c(I_max_p1s_no_other_policies, I_max_p2s_no_other_policies, I_max_p3s_no_other_policies, I_max_p4s_no_other_policies),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='max value of I', expand=0, limits=c(0,9000)) +
  theme_classic() +
  theme(legend.position='none')

plot3 <- data.frame(red=reds,
                    B_max=c(B_max_p1s_no_other_policies, B_max_p2s_no_other_policies, B_max_p3s_no_other_policies, B_max_p4s_no_other_policies),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=B_max,col=p_type)) +
  scale_x_continuous(name='reduction (%) in deaths due to focal policy effect', expand=0) +
  scale_y_continuous(name='max value of B', expand=0, limits=c(0,1)) +
  theme_classic() +
  theme(legend.position='none')



### effects of each policy with other policies each reducing 10% of deaths on their own

p1_def <- deaths_red_p1(.1)
p2_def <- deaths_red_p2(.1)
p3_def <- deaths_red_p3(.1)
p4_def <- deaths_red_p4(.1)

I_max_t_p1s <- c()
I_max_t_p2s <- c()
I_max_t_p3s <- c()
I_max_t_p4s <- c()

I_max_p1s <- c()
I_max_p2s <- c()
I_max_p3s <- c()
I_max_p4s <- c()

B_max_p1s <- c()
B_max_p2s <- c()
B_max_p3s <- c()
B_max_p4s <- c()

deaths_p1s <- c()
deaths_p2s <- c()
deaths_p3s <- c()
deaths_p4s <- c()

p1s <- c()
p2s <- c()
p3s <- c()
p4s <- c()

reds <- seq(0,.6,.01)

for (red in reds) {
  print(red)
  
  p1 <- deaths_red_p1(red, p1=p1_def, p2=p2_def, p3=p3_def, p4=p4_def)
  p2 <- deaths_red_p2(red, p1=p1_def, p2=p2_def, p3=p3_def, p4=p4_def)
  p3 <- deaths_red_p3(red, p1=p1_def, p2=p2_def, p3=p3_def, p4=p4_def)
  p4 <- deaths_red_p4(red, p1=p1_def, p2=p2_def, p3=p3_def, p4=p4_def)
  
  p1s <- c(p1s, p1)
  p2s <- c(p2s, p2)
  p3s <- c(p3s, p3)
  p4s <- c(p4s, p4)
  
  I_max_t_p1s <- c(I_max_t_p1s, I_max_t(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  I_max_t_p2s <- c(I_max_t_p2s, I_max_t(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  I_max_t_p3s <- c(I_max_t_p3s, I_max_t(p1=p1_def, p2=p2_def, p3=p3, p4=p4_def))
  I_max_t_p4s <- c(I_max_t_p4s, I_max_t(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4))
  
  I_max_p1s <- c(I_max_p1s, I_max(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  I_max_p2s <- c(I_max_p2s, I_max(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  I_max_p3s <- c(I_max_p3s, I_max(p1=p1_def, p2=p2_def, p3=p3, p4=p4_def))
  I_max_p4s <- c(I_max_p4s, I_max(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4))
  
  B_max_p1s <- c(B_max_p1s, B_max(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  B_max_p2s <- c(B_max_p2s, B_max(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  B_max_p3s <- c(B_max_p3s, B_max(p1=p1_def, p2=p2_def, p3=p3, p4=p4_def))
  B_max_p4s <- c(B_max_p4s, B_max(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4))
  
  deaths_p1s <- c(deaths_p1s, deaths_first_wave(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  deaths_p2s <- c(deaths_p2s, deaths_first_wave(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  deaths_p3s <- c(deaths_p3s, deaths_first_wave(p1=p1_def, p2=p2_def, p3=p3, p4=p4_def))
  deaths_p4s <- c(deaths_p4s, deaths_first_wave(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4))
}

p1s_low_other_policies <- p1s
p2s_low_other_policies <- p2s
p3s_low_other_policies <- p3s
p4s_low_other_policies <- p4s

I_max_t_p1s_low_other_policies <- I_max_t_p1s
I_max_t_p2s_low_other_policies <- I_max_t_p2s
I_max_t_p3s_low_other_policies <- I_max_t_p3s
I_max_t_p4s_low_other_policies <- I_max_t_p4s

I_max_p1s_low_other_policies <- I_max_p1s
I_max_p2s_low_other_policies <- I_max_p2s
I_max_p3s_low_other_policies <- I_max_p3s
I_max_p4s_low_other_policies <- I_max_p4s

B_max_p1s_low_other_policies <- B_max_p1s
B_max_p2s_low_other_policies <- B_max_p2s
B_max_p3s_low_other_policies <- B_max_p3s
B_max_p4s_low_other_policies <- B_max_p4s

deaths_p1s_low_other_policies <- deaths_p1s
deaths_p2s_low_other_policies <- deaths_p2s
deaths_p3s_low_other_policies <- deaths_p3s
deaths_p4s_low_other_policies <- deaths_p4s

plot4 <- data.frame(red=reds,
                    I_max_t=c(I_max_t_p1s_low_other_policies, I_max_t_p2s_low_other_policies, I_max_t_p3s_low_other_policies, I_max_t_p4s_low_other_policies),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  ggtitle('minor policy effects\non at baseline') +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='', expand=0, limits=c(0,550)) +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust=.5))


plot5 <- data.frame(red=reds,
                    I_max_t=c(I_max_p1s_low_other_policies, I_max_p2s_low_other_policies, I_max_p3s_low_other_policies, I_max_p4s_low_other_policies),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='', expand=0, limits=c(0,9000)) +
  scale_color_discrete(name='focal policy effect...', labels=c(TeX('reduces transmission ($p_1$)'),
                                                              TeX('increases testing ($p_2$)'),
                                                              TeX('increases C ($p_3$)'),
                                                              TeX('increases B ($p_4$)'))) +
  theme_classic()

plot6 <- data.frame(red=reds,
                    B_max=c(B_max_p1s_low_other_policies, B_max_p2s_low_other_policies, B_max_p3s_low_other_policies, B_max_p4s_low_other_policies),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=B_max,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='', expand=0, limits=c(0,1)) +
  theme_classic() +
  theme(legend.position='none')





plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, byrow=FALSE, nrow = 3, align = 'hv')

# saved as 9x8in





plot1 <- data.frame(red=reds,
           p=c(p1s_no_other_policies, p1s_low_other_policies),
           other_policies=rep(factor(c('no','low'),levels=c('no','low')),
                              each=length(reds))) %>%
  ggplot() +
  ggtitle(expression(p[1])) +
  geom_line(aes(x=red*100,y=p,linetype=other_policies), col=hue_pal()(4)[1]) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name=expression(atop(NA,atop(textstyle('policy effect strength'~p[i]~'needed to'),textstyle('get desired reduction in deaths')))), expand=0) +
  theme_classic() +
  theme(legend.position='none',
        plot.title=element_text(hjust=.5))

plot2 <- data.frame(red=reds,
           p=c(p2s_no_other_policies, p2s_low_other_policies),
           other_policies=rep(factor(c('no','low'),levels=c('no','low')),
                              each=length(reds))) %>%
  ggplot() +
  ggtitle(expression(p[2])) +
  geom_line(aes(x=red*100,y=p,linetype=other_policies), col=hue_pal()(4)[2]) +
  scale_x_continuous(name='reduction (%) in deaths due to focal intervention', expand=0) +
  scale_y_continuous(name='',
                     expand=0) +
  theme_classic() +
  theme(legend.position='none',
        plot.title=element_text(hjust=.5))

plot3 <- data.frame(red=reds,
           p=c(p3s_no_other_policies, p3s_low_other_policies,p4s_no_other_policies, p4s_low_other_policies),
           other_policies=rep(factor(c('no','low'),levels=c('no','low')),
                              each=length(reds)),
           p_type=rep(c('p3','p4'), each=2*length(reds))) %>%
  ggplot() +
  ggtitle(expression(p[3]~'and'~p[4])) +
  geom_line(aes(x=red*100,y=p,col=p_type, linetype=other_policies)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='',
                     expand=0, limits=c(0,.7)) +
  scale_color_manual(values=c(hue_pal()(4)[3],hue_pal()(4)[4]), guide='none') +
  scale_linetype_discrete(name='policy effects at baseline\neach reduce deaths by...', labels=c('0%','10%','20%')) +
  theme_classic() +
  theme(plot.title=element_text(hjust=.5),
        legend.position='bottom')




plot_grid(plot1, plot2, plot3, nrow = 1, align = 'hv')
#10x3.5









p1_big <- deaths_red_p1(.5)
p2_big <- deaths_red_p2(.5)
p3_big <- deaths_red_p3(.5)
p4_big <- deaths_red_p4(.5)
p3_small <- deaths_red_p3(.1)

deaths_red_p_vec <- function(red, p_vec, Tp=0) {
  if (p_vec[1] > 0) {
    mult_max <- (1-parms_no_policy['eB']) / p_vec[1]
  } else {
    mult_max <- 100
  }
  
  mult <- uniroot(function(mult) deaths_first_wave(p1=p_vec[1]*mult, p2=p_vec[2]*mult, p3=p_vec[3]*mult, p4=p_vec[4]*mult, Tp=Tp) - (1-red) * deaths_no_policy, c(0, mult_max))$root
  
  return(c(p1=mult*p_vec[1], p2=mult*p_vec[2], p3=mult*p_vec[3], p4=mult*p_vec[4]))
}





out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_big))
out_no_p <- out_no_policy
colnames(out_no_p) <- sapply(colnames(out_no_policy), function(var_name) paste0(var_name,"_no_p"))


I_max <- ts[which.max(out[,'I'])]

plot1 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('I','I_no_p','D','D_no_p','Q','Q_no_p')) %>% 
  mutate(var=factor(var, levels=c('I','I_no_p','D','D_no_p','Q','Q_no_p'))) %>% 
  ggplot() +
  ggtitle(expression(atop(NA,atop('improved ventilation',(p[1]*'>>0'))))) +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','blue','red','red','violet','violet'), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12,'solid',12), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))


plot2 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('C','C_no_p','B','B_no_p')) %>% 
  mutate(var=factor(var, levels=c('C','C_no_p','B','B_no_p'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','forestgreen','darkorange1','darkorange1'), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))



p_vec <- deaths_red_p_vec(.5,c(0,p2_big,p3_small,0))
out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p_vec[2],p3=p_vec[3]))


I_max <- ts[which.max(out[,'I'])]

plot3 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('I','I_no_p','D','D_no_p','Q','Q_no_p')) %>% 
  mutate(var=factor(var, levels=c('I','I_no_p','D','D_no_p','Q','Q_no_p'))) %>% 
  ggplot() +
  ggtitle(expression(atop(NA,atop('mandatory testing program',(p[2]*'>>0,'~p[3]>0))))) +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','blue','red','red','violet','violet'), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12,'solid',12), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))




plot4 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('C','C_no_p','B','B_no_p')) %>% 
  mutate(var=factor(var, levels=c('C','C_no_p','B','B_no_p'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','forestgreen','darkorange1','darkorange1'), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))




out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_big))


I_max <- ts[which.max(out[,'I'])]

plot5 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('I','I_no_p','D','D_no_p','Q','Q_no_p')) %>% 
  mutate(var=factor(var, levels=c('I','I_no_p','D','D_no_p','Q','Q_no_p'))) %>% 
  ggplot() +
  ggtitle(expression(atop(NA,atop('promoting risk awareness',(p[3]*'>>0'))))) +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','blue','red','red','violet','violet'), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12,'solid',12), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))


plot6 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('C','C_no_p','B','B_no_p')) %>% 
  mutate(var=factor(var, levels=c('C','C_no_p','B','B_no_p'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','forestgreen','darkorange1','darkorange1'), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))





p_vec <- deaths_red_p_vec(.5,c(0,0,p3_small,p4_big))
out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p_vec[3],p4=p_vec[4]))


I_max <- ts[which.max(out[,'I'])]

plot7 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('I','I_no_p','D','D_no_p','Q','Q_no_p')) %>% 
  mutate(var=factor(var, levels=c('I','I_no_p','D','D_no_p','Q','Q_no_p'))) %>% 
  ggplot() +
  ggtitle(expression(atop(NA,atop('promoting NPI use',(p[4]*'>>0,'~p[3]>0))))) +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','blue','red','red','violet','violet'), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12,'solid',12), labels=c('infected population, I','I with no policy', 'cumulative deaths, D', 'D with no policy','quarantined population, Q', 'Q with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size=18))


plot8 <- cbind(out[1:one_year,], out_no_p[1:one_year,]) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  filter(var %in% c('C','C_no_p','B','B_no_p')) %>% 
  mutate(var=factor(var, levels=c('C','C_no_p','B','B_no_p'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max, color='gray', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=var), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','forestgreen','darkorange1','darkorange1'), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  scale_linetype_manual(name='', values=c('solid',12,'solid',12), labels=c('perceived risk, C', 'C with no policy', 'NPI adoption, B', 'B with no policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size=18))









# export as 10x9
plot_grid(plot1, plot2, plot5, plot6, plot3, plot4, plot7, plot8, ncol = 2, byrow=FALSE, align='hv')





















p1_ex <- deaths_red_p1(.5)
p2_ex <- deaths_red_p2(.5)
p3_ex <- deaths_red_p3(.5)
p4_ex <- deaths_red_p4(.5)




Tps <- seq(0,200,1)


D365_Tp1 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[,'D'],1)))
}
D365_Tp2 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[,'D'],1)))
}
D365_Tp3 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[,'D'],1)))
}
D365_Tp4 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[,'D'],1)))
}

D365_Tp1s <- sapply(Tps, D365_Tp1)
D365_Tp2s <- sapply(Tps, D365_Tp2)
D365_Tp3s <- sapply(Tps, D365_Tp3)
D365_Tp4s <- sapply(Tps, D365_Tp4)

plot1 <- data.frame(Tp=Tps,
           D365=c(D365_Tp1s,D365_Tp2s,D365_Tp3s,D365_Tp4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(Tps))) %>%
  ggplot() +
  geom_line(aes(x=Tp,y=D365,col=p_type)) +
  scale_x_continuous(name=TeX('policy start time (days)'), expand=0) +
  scale_y_continuous(name=TeX('deaths after first wave'), expand=0, limits=c(2500,6000)) +
  theme_classic() +
  theme(legend.position="none")







I_max_time_Tp1 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[,'I'])])
}
I_max_time_Tp2 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[,'I'])])
}
I_max_time_Tp3 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[,'I'])])
}
I_max_time_Tp4 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[,'I'])])
}


I_max_time_Tp1s <- sapply(Tps, I_max_time_Tp1)
I_max_time_Tp2s <- sapply(Tps, I_max_time_Tp2)
I_max_time_Tp3s <- sapply(Tps, I_max_time_Tp3)
I_max_time_Tp4s <- sapply(Tps, I_max_time_Tp4)

plot2 <- data.frame(Tp=Tps,
           I_max_time=c(I_max_time_Tp1s,I_max_time_Tp2s,I_max_time_Tp3s,I_max_time_Tp4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(Tps))) %>%
  ggplot() +
  geom_line(aes(x=Tp,y=I_max_time,col=p_type)) +
  scale_x_continuous(name=TeX('policy start time (days)'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infected population (days)'), expand=0, limits=c(0,365)) +
  scale_color_discrete(name='policy effect that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes NPI use ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))


plot_grid(plot1, plot2, ncol = 2, align='hv')








I_max_Tp1 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[,'I']))
}
I_max_Tp2 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[,'I']))
}
I_max_Tp3 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[,'I']))
}
I_max_Tp4 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[,'I']))
}





I_max_Tp1s <- sapply(Tps, I_max_Tp1)
I_max_Tp2s <- sapply(Tps, I_max_Tp2)
I_max_Tp3s <- sapply(Tps, I_max_Tp3)
I_max_Tp4s <- sapply(Tps, I_max_Tp4)

data.frame(Tp=Tps,
           I_max=c(I_max_Tp1s,I_max_Tp2s,I_max_Tp3s,I_max_Tp4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(Tps))) %>%
  ggplot() +
  geom_line(aes(x=Tp,y=I_max,col=p_type)) +
  scale_x_continuous(name=TeX('policy start time (days)'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infected population (days)'), expand=0, limits=c(0,10000)) +
  scale_color_discrete(name='policy effect that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes NPI use ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))








B_max_Tp1 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[,'B']))
}
B_max_Tp2 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[,'B']))
}
B_max_Tp3 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[,'B']))
}
B_max_Tp4 <- function(Tp) {
  return(max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[,'B']))
}


B_max_Tp1s <- sapply(Tps, B_max_Tp1)
B_max_Tp2s <- sapply(Tps, B_max_Tp2)
B_max_Tp3s <- sapply(Tps, B_max_Tp3)
B_max_Tp4s <- sapply(Tps, B_max_Tp4)

data.frame(Tp=Tps,
           B_max=c(B_max_Tp1s,B_max_Tp2s,B_max_Tp3s,B_max_Tp4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(Tps))) %>%
  ggplot() +
  geom_line(aes(x=Tp,y=B_max,col=p_type)) +
  scale_x_continuous(name=TeX('policy start time (days)'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infected population (days)'), expand=0, limits=c(0,1)) +
  scale_color_discrete(name='policy effect that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes NPI use ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))

















# effects of freeriding and cB2





p1_ex <- deaths_red_p1(.3)
p2_ex <- deaths_red_p2(.3)
p3_ex <- deaths_red_p3(.3)
p4_ex <- deaths_red_p4(.3)






parms_freeriding <- function(cB2=0, p1=0, p2=0, p3=0, p4=0, Tp=0) {
  parms <- parms_policy(p1,p2,p3,p4,Tp)
  parms['cB2'] <- cB2
  
  return(parms)
}

out_p1_ex <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p1')
out_p2_ex <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>% 
  mutate(policy_type='p2')
out_p3_ex <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p3')
out_p4_ex <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p4')

out_p1_freeriding <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(cB2=.3, p1=p1_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p1')
out_p2_freeriding <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(cB2=.3, p2=p2_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p2')
out_p3_freeriding <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(cB2=.3, p3=p3_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p3')
out_p4_freeriding <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(cB2=.3, p4=p4_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p4')





policy_labels <- c('policy reducing\ntransmission directly (p1)',
                   'policy increasing\ntesting rate (p2)',
                   'policy increasing\nrisk perception (p3)',
                   'policy promoting\nNPI use (p4)')
names(policy_labels) <- c('p1','p2','p3','p4')


rbind(cbind(out_p1_ex, freeriding=FALSE),
      cbind(out_p2_ex, freeriding=FALSE),
      cbind(out_p3_ex, freeriding=FALSE),
      cbind(out_p4_ex, freeriding=FALSE),
      cbind(out_p1_freeriding, freeriding=TRUE),
      cbind(out_p2_freeriding, freeriding=TRUE),
      cbind(out_p3_freeriding, freeriding=TRUE),
      cbind(out_p4_freeriding, freeriding=TRUE)) %>%
  filter(var %in% c('I','D','C','B'), policy) %>% 
  mutate(var=factor(var, levels=c('I','D','C','B')),
         value=ifelse(var %in% c('I','D'), value/6000, value)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value*6000, color=interaction(freeriding,var), lty=interaction(freeriding,var), linewidth=interaction(freeriding,var))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0, limits=c(0,365)) +
  scale_y_continuous(name='number of people (I or D)', expand=0, limits=c(0,6000),
                     sec.axis = sec_axis( transform=~./6000, name='risk perception (C) or\nNPI adoption (B)')) +
  scale_color_manual(name='', values=c('blue','blue','red','red','forestgreen','forestgreen','darkorange1','darkorange1'), labels=c('infected population, I,\nwithout freeriding','infected population, I,\nwith freeriding','cumulative deaths, D,\nwithout freeriding','cumulative deaths, D,\nwith freeriding','risk perception, C,\nwithout freeriding','risk perception, C,\nwith freeriding','NPI adoption, B,\nwithout freeriding','NPI adoption, B,\nwith freeriding')) +
  scale_linetype_manual(name='', values=c('solid','24','solid','24','solid','24','solid','24'), labels=c('infected population, I,\nwithout freeriding','infected population, I,\nwith freeriding','cumulative deaths, D,\nwithout freeriding','cumulative deaths, D,\nwith freeriding','risk perception, C,\nwithout freeriding','risk perception, C,\nwith freeriding','NPI adoption, B,\nwithout freeriding','NPI adoption, B,\nwith freeriding')) +
  scale_linewidth_manual(name='', values=c(.5,.7,.5,.7,.5,.7,.5,.7), labels=c('infected population, I,\nwithout freeriding','infected population, I,\nwith freeriding','cumulative deaths, D,\nwithout freeriding','cumulative deaths, D,\nwith freeriding','risk perception, C,\nwithout freeriding','risk perception, C,\nwith freeriding','NPI adoption, B,\nwithout freeriding','NPI adoption, B,\nwith freeriding')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11),
        panel.spacing = unit(1, 'lines'),
        legend.key.spacing.y = unit(0.3, 'cm'),
        legend.text=element_text(size=10))




















### SUPPLEMENT

p1_ex <- deaths_red_p1(.5)
p2_ex <- deaths_red_p2(.5)
p3_ex <- deaths_red_p3(.5)
p4_ex <- deaths_red_p4(.5)


reds <- seq(.01,.6,.01)

p1s <- c()
p2s <- c()
p3s <- c()
p4s <- c()
for (red in reds) {
  print(red)
  
  vec12 <- deaths_red_p_vec(red, c(p1_ex,p2_ex,0,0))
  p1s <- c(p1s, vec12[1])
  p2s <- c(p2s, vec12[2])
  
  vec34 <- deaths_red_p_vec(red, c(0,0,p3_ex,p4_ex))
  p3s <- c(p3s, vec34[3])
  p4s <- c(p4s, vec34[4])
}


combined_ps <- data.frame()

for (red12 in 1:60) {
  print(reds[red12])
  for (red34 in 1:60) {
    for (Tp in c(0,50,100)) {
      out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1s[red12], p2=p2s[red12], p3=p3s[red34], p4=p4s[red34], Tp=Tp))
      
      combined_ps <-
        combined_ps %>% 
        rbind(data.frame(red12 = round(reds[red12],2),
                         red34 = round(reds[red34],2),
                         Tp = Tp,
                         D365 = unname(tail(out[, 'D'],1)),
                         I_max = max(out[,'I']),
                         I_max_time = ts[which.max(out[,'I'])],
                         B_max = max(out[,'B'])))
    }
  }
}



Tp_labels <- c('started at very beginning\nof epidemic (t=0)',
               'started early\nin epidemic (t=50)',
               'started late\nin epidemic (t=100)')
names(Tp_labels) <- c('0','50','100')

combined_ps %>% 
  ggplot() +
  geom_raster(aes(x=red12,y=red34,fill=B_max), hjust=1, vjust=1) +
  facet_wrap(~Tp, labeller = as_labeller(Tp_labels), scales='free') +
  scale_fill_gradientn(name='deaths after\n1 year',colors=rev(rainbow(7)[-7]),
                       guide=guide_colorbar(frame.colour = 'black', ticks.colour = 'black',title.hjust=.5)) +
  scale_x_continuous(name=expression(atop(NA,atop(textstyle('strength of top-down policies'),textstyle('(split equally between ' * p[1]*' and '*p[2]*')')))),expand=0) +
  scale_y_continuous(name=expression(atop(NA,atop(textstyle('strength of bottom-up policies'),textstyle('(split equally between ' * p[3]*' and '*p[4]*')')))),expand=0) +
  theme(axis.line=element_line(),
        legend.title=element_text(size=10),
        legend.position='right',
        panel.background = element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=11),
        aspect.ratio = 1)
