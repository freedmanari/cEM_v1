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
    dC <- ((cB - cB2 * B) * B + cD * D + cQ * Q + ifelse(t > Tp, p3, 0)) * (1 - C) - deltaC * C
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
                     R0=2, eB=.5, kB=0, kC=.05, deltaC=.2, cB=.1, cB2=0,
                     cD=1/10^5, cQ=1/10^4, deltaB=.1, bC=.1,
                     p1=0, p2=0, p3=0, p4=0)
ts <- seq(0,3650,.1)
one_year <- which.min(abs(ts-365))


##### NO POLICY


### Fig. 3 (outs_no_policy)

# low sensitivity to epi information

out_no_policy <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_no_policy)

out <- out_no_policy[1:one_year,]

I_max_t <- ts[which.max(out[,'I'])]

plot1 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='full_var',values_to='value') %>%
  mutate(var=substring(full_var,1,1),
         C_B=!grepl('o',full_var)) %>% 
  filter(var %in% c('I','D','Q')) %>% 
  mutate(var=factor(var, levels=c('I','D','Q')),
         C_B=factor(C_B, levels=c(FALSE,TRUE))) %>% 
  ggplot() +
  ggtitle('low sensitivity to\nepi information') +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var, linetype=C_B), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(name='number of people', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('model without C and B','model with C and B')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust = 0.5))


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
  scale_color_manual(labels=c('C (perceived risk)','B (NPI adoption)'),
                     values=c('forestgreen','darkorange1')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        legend.position='none',
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
        legend.position='none',
        axis.title.y=element_blank())




# high sensitivity to epi information

new_parms <- parms_no_policy
new_parms['cD'] <- parms_no_policy['cD'] * 10
new_parms['cQ'] <- parms_no_policy['cQ'] * 10

out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, new_parms)[1:one_year,]

I_max_t <- ts[which.max(out[,'I'])]

plot4 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='full_var',values_to='value') %>%
  mutate(var=substring(full_var,1,1),
         C_B=!grepl('o',full_var)) %>% 
  filter(var %in% c('I','D','Q')) %>% 
  mutate(var=factor(var, levels=c('I','D','Q'))) %>% 
  ggplot() +
  ggtitle('low sensitivity to\nepi information') +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var, linetype=C_B), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(name='number of people', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('model without C and B','model with C and B')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust = 0.5)) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

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
  scale_color_manual(labels=c('C (perceived risk)','B (NPI adoption)'),
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







### Fig. S1 (outs_no_policy_high_cB)


# low sensitivity to epi information, with high sensitivity to behavioral information

new_parms <- parms_no_policy
new_parms['cB'] <- parms_no_policy['cB'] * 10

out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, new_parms)[1:one_year,]

I_max_t <- ts[which.max(out[,'I'])]

plot1 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='full_var',values_to='value') %>%
  mutate(var=substring(full_var,1,1),
         C_B=!grepl('o',full_var)) %>% 
  filter(var %in% c('I','D','Q')) %>% 
  mutate(var=factor(var, levels=c('I','D','Q'))) %>% 
  ggplot() +
  ggtitle('low sensitivity to\nepi information') +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var, linetype=C_B), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(name='number of people', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('model without C and B','model with C and B')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust = 0.5))


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
  scale_color_manual(labels=c('C (perceived risk)','B (NPI adoption)'),
                     values=c('forestgreen','darkorange1')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        legend.position='none',
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
        legend.position='none',
        axis.title.y=element_blank())




# high sensitivity to epi information, with high sensitivity to behavioral information

new_parms <- parms_no_policy
new_parms['cD'] <- parms_no_policy['cD'] * 10
new_parms['cQ'] <- parms_no_policy['cQ'] * 10
new_parms['cB'] <- parms_no_policy['cB'] * 10

out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, new_parms)[1:one_year,]

I_max_t <- ts[which.max(out[,'I'])]

plot4 <- out %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='full_var',values_to='value') %>%
  mutate(var=substring(full_var,1,1),
         C_B=!grepl('o',full_var)) %>% 
  filter(var %in% c('I','D','Q')) %>% 
  mutate(var=factor(var, levels=c('I','D','Q'))) %>% 
  ggplot() +
  ggtitle('high sensitivity to\nepi information') +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.5) +
  geom_line(aes(x=time, y=value, color=var, linetype=C_B), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(name='number of people', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('model without C and B','model with C and B')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, .3, 0), 'cm'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title=element_text(hjust = 0.5)) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )

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
  scale_color_manual(labels=c('C (perceived risk)','B (NPI adoption)'),
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







### Fig. 4 (feedback_strengths)


feedback_strengths <- data.frame()
parms <- parms_no_policy

for (epi_to_CB_strength in c(.1,.2,.5,1,2,5,10)) {
  print(epi_to_CB_strength)
  
  parms['cD'] <- epi_to_CB_strength * parms_no_policy['cD']
  parms['cQ'] <- epi_to_CB_strength * parms_no_policy['cQ']
  
  for (CB_to_epi_strength in seq(0,2,.01)) {
    parms['eB'] <- CB_to_epi_strength * parms_no_policy['eB']
    parms['kC'] <- CB_to_epi_strength * parms_no_policy['kC']
    
    out <- ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms)
    
    feedback_strengths <-
      feedback_strengths %>% 
      rbind(data.frame(epi_to_CB_strength = epi_to_CB_strength,
                       CB_to_epi_strength = round(CB_to_epi_strength,2),
                       deaths = unname(tail(out[,'D'],1)),
                       I_max_t = ts[which.max(out[,'I'])],
                       I_max = max(out[,'I']),
                       B_eq = tail(out[,'B'],1),
                       B_max = max(out[,'B'])))
  }
}

plot1 <- feedback_strengths %>%
  ggplot() +
  geom_line(aes(x=CB_to_epi_strength, y=deaths, color=fct_rev(factor(epi_to_CB_strength)))) +
  scale_x_continuous(name='relative strength of cognitive-behavioral effects on transmission',expand=0,
                     labels=c('0','.5','1','1.5','2')) +
  scale_y_continuous(name='deaths after first wave', expand=0, limits=c(0,7500)) +
  scale_color_discrete(name='relative sensitivity of\nC to epi information', labels=c('10','5','2','1','.5','.2','.1')) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        legend.position='none')

plot2 <- feedback_strengths %>%
  ggplot() +
  geom_line(aes(x=CB_to_epi_strength, y=B_max, color=fct_rev(factor(epi_to_CB_strength)))) +
  scale_x_continuous(name='relative strength of cognitive-behavioral effects on transmission',expand=0,
                     labels=c('0','.5','1','1.5','2')) +
  scale_y_continuous(name='max value of B', expand=0, limits=c(0,.5), labels=c('0','.1','.2','.3','.4','.5')) +
  scale_color_discrete(name='relative sensitivity of\nC to epi information', labels=c('10','5','2','1','.5','.2','.1')) +
  theme_classic()

plot_grid(plot1, plot2, ncol=2, align='vh')













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






### Fig. 5 (example_policies)

p1_big <- deaths_red_p1(.5)
p2_big <- deaths_red_p2(.5)
p3_big <- deaths_red_p3(.5)
p4_big <- deaths_red_p4(.5)
p3_small <- deaths_red_p3(.1)


out_no_p <-
  out_no_policy %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(policy='no')
out <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_big)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(policy='yes')

I_max_t <- ts[which.max(pull(filter(out, var=='I'), value))]

plot1 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q')),
         policy=factor(policy, levels=c('no', 'yes'))) %>%
  ggplot() +
  ggtitle(expression(atop(NA,atop('improved ventilation',(p[1]))))) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('no policy','policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))


plot2 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B')),
         policy=factor(policy, levels=c('no','yes'))) %>%
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('no policy','policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))



out <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_big, p3=p3_small)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(policy='yes')


I_max_t <- ts[which.max(pull(filter(out, var=='I'), value))]

plot3 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q')),
         policy=factor(policy, levels=c('no','yes'))) %>%
  ggplot() +
  ggtitle(expression(atop(NA,atop('mandatory testing program',(p[2]*','~p[3]))))) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('no policy','policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))


plot4 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B')),
         policy=factor(policy, levels=c('no','yes'))) %>%
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('no policy','policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))




out <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_big)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(policy='yes')


I_max_t <- ts[which.max(pull(filter(out, var=='I'), value))]

plot5 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q')),
         policy=factor(policy, levels=c('no','yes'))) %>%
  ggplot() +
  ggtitle(expression(atop(NA,atop('promoting risk awareness',(p[3]))))) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('no policy','policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))


plot6 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B')),
         policy=factor(policy, levels=c('no','yes'))) %>%
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('no policy','policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank(),
        legend.position='none',
        plot.title = element_text(hjust = 0.5, size=18))




out <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_small, p4=p4_big)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(policy='yes')


I_max_t <- ts[which.max(pull(filter(out, var=='I'), value))]

plot7 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q')),
         policy=factor(policy, levels=c('no','yes'))) %>%
  ggplot() +
  ggtitle(expression(atop(NA,atop('promoting NPI use',(p[4]*','~p[3]))))) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(expand=0) +
  scale_y_continuous(expand=0, limits=c(0,9000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), guide='none') +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size=18))


plot8 <- 
  rbind(out_no_p,
        out) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B')),
         policy=factor(policy, levels=c('no','yes'))) %>%
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=policy), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('no policy','policy')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5, size=18)) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )



# export as 10x9
plot_grid(plot1, plot2, plot5, plot6, plot3, plot4, plot7, plot8, ncol = 2, byrow=FALSE, align='hv')













### Fig. 6 (deaths_red)


# no indirect bottom-up effects

reds <- seq(0,.6,.005)

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

p1s_without_bottom_up_effects <- p1s
p2s_without_bottom_up_effects <- p2s
p3s_without_bottom_up_effects <- p3s
p4s_without_bottom_up_effects <- p4s

I_max_t_p1s_without_bottom_up_effects <- I_max_t_p1s
I_max_t_p2s_without_bottom_up_effects <- I_max_t_p2s
I_max_t_p3s_without_bottom_up_effects <- I_max_t_p3s
I_max_t_p4s_without_bottom_up_effects <- I_max_t_p4s

I_max_p1s_without_bottom_up_effects <- I_max_p1s
I_max_p2s_without_bottom_up_effects <- I_max_p2s
I_max_p3s_without_bottom_up_effects <- I_max_p3s
I_max_p4s_without_bottom_up_effects <- I_max_p4s

B_max_p1s_without_bottom_up_effects <- B_max_p1s
B_max_p2s_without_bottom_up_effects <- B_max_p2s
B_max_p3s_without_bottom_up_effects <- B_max_p3s
B_max_p4s_without_bottom_up_effects <- B_max_p4s

deaths_p1s_without_bottom_up_effects <- deaths_p1s
deaths_p2s_without_bottom_up_effects <- deaths_p2s
deaths_p3s_without_bottom_up_effects <- deaths_p3s
deaths_p4s_without_bottom_up_effects <- deaths_p4s

plot1 <- data.frame(red=reds,
                    I_max_t=c(I_max_t_p1s_without_bottom_up_effects, I_max_t_p2s_without_bottom_up_effects, I_max_t_p3s_without_bottom_up_effects, I_max_t_p4s_without_bottom_up_effects),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  ggtitle('no indirect\nbottom-up effects') +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='peak time of I (days)', expand=0, limits=c(0,850)) +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust=.5))


plot2 <- data.frame(red=reds,
                    I_max_t=c(I_max_p1s_without_bottom_up_effects, I_max_p2s_without_bottom_up_effects, I_max_p3s_without_bottom_up_effects, I_max_p4s_without_bottom_up_effects),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='max value of I', expand=0, limits=c(0,9000)) +
  theme_classic() +
  theme(legend.position='none')

plot3 <- data.frame(red=reds,
                    B_max=c(B_max_p1s_without_bottom_up_effects, B_max_p2s_without_bottom_up_effects, B_max_p3s_without_bottom_up_effects, B_max_p4s_without_bottom_up_effects),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=B_max,col=p_type)) +
  scale_x_continuous(name='reduction (%) in deaths due to direct policy effect', expand=0) +
  scale_y_continuous(name='max value of B', expand=0, limits=c(0,.75)) +
  theme_classic() +
  theme(legend.position='none')




### effects of each policy with other policies each reducing 10% of deaths on their own

p1_def <- 0
p2_def <- 0
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
  
  I_max_t_p1s <- c(I_max_t_p1s, I_max_t(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  I_max_t_p2s <- c(I_max_t_p2s, I_max_t(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  I_max_t_p3s <- c(I_max_t_p3s, I_max_t(p1=p1_def, p2=p2_def, p3=p3+p3_def, p4=p4_def))
  I_max_t_p4s <- c(I_max_t_p4s, I_max_t(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4+p4_def))
  
  I_max_p1s <- c(I_max_p1s, I_max(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  I_max_p2s <- c(I_max_p2s, I_max(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  I_max_p3s <- c(I_max_p3s, I_max(p1=p1_def, p2=p2_def, p3=p3+p3_def, p4=p4_def))
  I_max_p4s <- c(I_max_p4s, I_max(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4+p4_def))
  
  B_max_p1s <- c(B_max_p1s, B_max(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  B_max_p2s <- c(B_max_p2s, B_max(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  B_max_p3s <- c(B_max_p3s, B_max(p1=p1_def, p2=p2_def, p3=p3+p3_def, p4=p4_def))
  B_max_p4s <- c(B_max_p4s, B_max(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4+p4_def))
  
  deaths_p1s <- c(deaths_p1s, deaths_first_wave(p1=p1, p2=p2_def, p3=p3_def, p4=p4_def))
  deaths_p2s <- c(deaths_p2s, deaths_first_wave(p1=p1_def, p2=p2, p3=p3_def, p4=p4_def))
  deaths_p3s <- c(deaths_p3s, deaths_first_wave(p1=p1_def, p2=p2_def, p3=p3+p3_def, p4=p4_def))
  deaths_p4s <- c(deaths_p4s, deaths_first_wave(p1=p1_def, p2=p2_def, p3=p3_def, p4=p4+p4_def))
}

p1s_with_bottom_up_effects <- p1s
p2s_with_bottom_up_effects <- p2s
p3s_with_bottom_up_effects <- p3s
p4s_with_bottom_up_effects <- p4s

I_max_t_p1s_with_bottom_up_effects <- I_max_t_p1s
I_max_t_p2s_with_bottom_up_effects <- I_max_t_p2s
I_max_t_p3s_with_bottom_up_effects <- I_max_t_p3s
I_max_t_p4s_with_bottom_up_effects <- I_max_t_p4s

I_max_p1s_with_bottom_up_effects <- I_max_p1s
I_max_p2s_with_bottom_up_effects <- I_max_p2s
I_max_p3s_with_bottom_up_effects <- I_max_p3s
I_max_p4s_with_bottom_up_effects <- I_max_p4s

B_max_p1s_with_bottom_up_effects <- B_max_p1s
B_max_p2s_with_bottom_up_effects <- B_max_p2s
B_max_p3s_with_bottom_up_effects <- B_max_p3s
B_max_p4s_with_bottom_up_effects <- B_max_p4s

deaths_p1s_with_bottom_up_effects <- deaths_p1s
deaths_p2s_with_bottom_up_effects <- deaths_p2s
deaths_p3s_with_bottom_up_effects <- deaths_p3s
deaths_p4s_with_bottom_up_effects <- deaths_p4s

plot4 <- data.frame(red=reds,
                    I_max_t=c(I_max_t_p1s_with_bottom_up_effects, I_max_t_p2s_with_bottom_up_effects, I_max_t_p3s_with_bottom_up_effects, I_max_t_p4s_with_bottom_up_effects),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  ggtitle(TeX('with indirect\nbottom-up effects ($p_3$, $p_4$)')) +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='', expand=0, limits=c(0,850)) +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust=.5))


plot5 <- data.frame(red=reds,
                    I_max_t=c(I_max_p1s_with_bottom_up_effects, I_max_p2s_with_bottom_up_effects, I_max_p3s_with_bottom_up_effects, I_max_p4s_with_bottom_up_effects),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=I_max_t,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='', expand=0, limits=c(0,9000)) +
  scale_color_discrete(name='direct policy effect...', labels=c(TeX('reduces transmission ($p_1$)'),
                                                               TeX('increases testing ($p_2$)'),
                                                               TeX('increases C ($p_3$)'),
                                                               TeX('increases B ($p_4$)'))) +
  theme_classic()

plot6 <- data.frame(red=reds,
                    B_max=c(B_max_p1s_with_bottom_up_effects, B_max_p2s_with_bottom_up_effects, B_max_p3s_with_bottom_up_effects, B_max_p4s_with_bottom_up_effects),
                    p_type=rep(c('p1','p2','p3','p4'),
                               each=length(reds))) %>%
  ggplot() +
  geom_line(aes(x=red*100,y=B_max,col=p_type)) +
  scale_x_continuous(name='', expand=0) +
  scale_y_continuous(name='', expand=0, limits=c(0,.75)) +
  theme_classic() +
  theme(legend.position='none')



plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, byrow=FALSE, nrow = 3, align = 'hv')

# saved as 9x8in




### Fig. S2 (total_deaths_red)

data.frame(red=reds,
           total_red=c(1-deaths_p1s_with_bottom_up_effects/deaths_no_policy,
                       1-deaths_p2s_with_bottom_up_effects/deaths_no_policy,
                       1-deaths_p3s_with_bottom_up_effects/deaths_no_policy,
                       1-deaths_p4s_with_bottom_up_effects/deaths_no_policy),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(reds))) %>% 
  ggplot() +
  geom_line(aes(x=red*100,y=total_red*100,col=p_type)) +
  geom_abline(slope=1,intercept=0,linetype='dashed') +
  scale_color_discrete(name='direct policy effect...', labels=c(TeX('reduces transmission ($p_1$)'),
                                                                TeX('increases testing ($p_2$)'),
                                                                TeX('increases C ($p_3$)'),
                                                                TeX('increases B ($p_4$)'))) +
  scale_x_continuous(name='reduction (%) in deaths due to\ndirect policy effect', expand=0) +
  scale_y_continuous(name='reduction (%) in deaths due to\ndirect and indirect policy effects', expand=0, limits=c(0,100)) +
  theme_classic()























### Fig. 7 (Tps)


p1_ex <- deaths_red_p1(.3)
p2_ex <- deaths_red_p2(.3)
p3_ex <- deaths_red_p3(.3)
p4_ex <- deaths_red_p4(.3)




Tps <- seq(0,200,1)


deaths_Tp1 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[,'D'],1)))
}
deaths_Tp2 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[,'D'],1)))
}
deaths_Tp3 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[,'D'],1)))
}
deaths_Tp4 <- function(Tp) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[,'D'],1)))
}

deaths_Tp1s <- sapply(Tps, deaths_Tp1)
deaths_Tp2s <- sapply(Tps, deaths_Tp2)
deaths_Tp3s <- sapply(Tps, deaths_Tp3)
deaths_Tp4s <- sapply(Tps, deaths_Tp4)

plot1 <- data.frame(Tp=Tps,
           deaths=c(deaths_Tp1s,deaths_Tp2s,deaths_Tp3s,deaths_Tp4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(Tps))) %>%
  ggplot() +
  geom_line(aes(x=Tp,y=deaths,col=p_type)) +
  scale_x_continuous(name=TeX('policy start time (days)'), expand=0) +
  scale_y_continuous(name=TeX('deaths after first wave'), expand=0, limits=c(2500,6000)) +
  theme_classic() +
  theme(legend.position='none')







I_max_t_Tp1 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[,'I'])])
}
I_max_t_Tp2 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[,'I'])])
}
I_max_t_Tp3 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[,'I'])])
}
I_max_t_Tp4 <- function(Tp) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[,'I'])])
}


I_max_t_Tp1s <- sapply(Tps, I_max_t_Tp1)
I_max_t_Tp2s <- sapply(Tps, I_max_t_Tp2)
I_max_t_Tp3s <- sapply(Tps, I_max_t_Tp3)
I_max_t_Tp4s <- sapply(Tps, I_max_t_Tp4)

plot2 <- data.frame(Tp=Tps,
           I_max_t=c(I_max_t_Tp1s,I_max_t_Tp2s,I_max_t_Tp3s,I_max_t_Tp4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(Tps))) %>%
  ggplot() +
  geom_line(aes(x=Tp,y=I_max_t,col=p_type)) +
  scale_x_continuous(name=TeX('policy start time (days)'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infected population (days)'), expand=0, limits=c(0,365)) +
  scale_color_discrete(name='policy effect that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes NPI use ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))


plot_grid(plot1, plot2, ncol = 2, align='hv')


















### Fig. 8 (freeriding)





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
  filter(var %in% c('I','D','C','B'), policy, time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','C','B')),
         value=ifelse(var %in% c('I','D'), value/6000, value)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value*6000, color=var, lty=freeriding)) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0, limits=c(0,365)) +
  scale_y_continuous(name='number of people (I or D)', expand=0, limits=c(0,6000),
                     sec.axis = sec_axis( transform=~./6000, name='risk perception (C) or\nNPI adoption (B)')) +
  scale_color_manual(name='', values=c('blue','red','forestgreen','darkorange1'), labels=c('I (infected population)','D (cumulative deaths)','C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c('solid','24'), labels=c('without freeriding','with freeriding')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11),
        panel.spacing = unit(1, 'lines'),
        legend.key.spacing.y = unit(0.3, 'cm'),
        legend.text=element_text(size=10))
# export as 9in x 6in



















### Fig. S3 (indirect_effects_examples)


p1_big <- deaths_red_p1(.5)
p2_big <- deaths_red_p2(.5)
p3_big <- deaths_red_p3(.5)
p4_big <- deaths_red_p4(.5)
p3_small <- deaths_red_p3(.1)
p4_small <- deaths_red_p4(.1)


out_no_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_big)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='no')
out_with_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p1=p1_big, p3=p3_small, p4=p4_small)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='yes')

I_max_t <- ts[which.max(pull(filter(out_with_indirect_effects, var=='I'), value))]

plot1 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q'))) %>% 
  ggplot() +
  ggtitle(expression('direct effect from'~p[1])) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,3000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('without indirect effects','with indirect effects')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

plot2 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('without indirect effects','with indirect effects')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position='none',
        axis.title.y=element_blank())




out_no_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_big)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='no')
out_with_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p2=p2_big, p3=p3_small, p4=p4_small)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='yes')

I_max_t <- ts[which.max(pull(filter(out_with_indirect_effects, var=='I'), value))]

plot3 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q'))) %>% 
  ggplot() +
  ggtitle(expression('direct effect from'~p[2])) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,3000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('without indirect effects','with indirect effects')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

plot4 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('without indirect effects','with indirect effects')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position='none',
        axis.title.y=element_blank())




out_no_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_big)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='no')
out_with_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_big+p3_small, p4=p4_small)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='yes')

I_max_t <- ts[which.max(pull(filter(out_with_indirect_effects, var=='I'), value))]

plot5 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q'))) %>% 
  ggplot() +
  ggtitle(expression('direct effect from'~p[3])) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,3000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('without indirect effects','with indirect effects')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position='none',
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

plot6 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('without indirect effects','with indirect effects')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.position='none',
        axis.title.y=element_blank())




out_no_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p4=p4_big)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='no')
out_with_indirect_effects <-
  ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_policy(p3=p3_small, p4=p4_big+p4_small)) %>% 
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>% 
  mutate(indirect_effects='yes')

I_max_t <- ts[which.max(pull(filter(out_with_indirect_effects, var=='I'), value))]

plot7 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('I','D','Q'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('I','D','Q'))) %>% 
  ggplot() +
  ggtitle(expression('direct effect from'~p[4])) +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,3000)) +
  scale_color_manual(name='', values=c('blue','red','violet'), labels=c('I (infected population)','D (cumulative deaths)','Q (quarantined population)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), guide='none') +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5))

plot8 <-
  rbind(out_no_indirect_effects,
        out_with_indirect_effects) %>% 
  filter(var %in% c('C','B'),
         time<=365) %>% 
  mutate(var=factor(var, levels=c('C','B'))) %>% 
  ggplot() +
  geom_vline(xintercept=I_max_t, color='gray85', linewidth=.4) +
  geom_line(aes(x=time, y=value, color=var, linetype=indirect_effects), linewidth=.7) +
  scale_x_continuous(name='time (days)', expand=0, limits=c(0,365)) +
  scale_y_continuous(expand=0, limits=c(0,1)) +
  scale_color_manual(name='', values=c('forestgreen','darkorange1'), labels=c('C (perceived risk)','B (NPI adoption)')) +
  scale_linetype_manual(name='', values=c(12,'solid'), labels=c('without indirect effects','with indirect effects')) +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), 'cm'),
        legend.text=element_text(size=10),
        axis.title.y=element_blank()) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  )


plot_grid(plot1, plot2, plot5, plot6, plot3, plot4, plot7, plot8, ncol = 2, byrow=FALSE, align='hv')
# export as 10x9








### Fig. S4 (cB2)


cB2s <- seq(0,.2,.002)
p1_ex <- deaths_red_p1(.3)
p2_ex <- deaths_red_p2(.3)
p3_ex <- deaths_red_p3(.3)
p4_ex <- deaths_red_p4(.3)

deaths_cB21 <- function(cB2) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p1=p1_ex, cB2=cB2))[,'D'],1)))
}
deaths_cB22 <- function(cB2) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p2=p2_ex, cB2=cB2))[,'D'],1)))
}
deaths_cB23 <- function(cB2) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p3=p3_ex, cB2=cB2))[,'D'],1)))
}
deaths_cB24 <- function(cB2) {
  return(unname(tail(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p4=p4_ex, cB2=cB2))[,'D'],1)))
}

deaths_cB21s <- sapply(cB2s, deaths_cB21)
deaths_cB22s <- sapply(cB2s, deaths_cB22)
deaths_cB23s <- sapply(cB2s, deaths_cB23)
deaths_cB24s <- sapply(cB2s, deaths_cB24)

data.frame(cB2=cB2s,
           deaths=c(deaths_cB21s,deaths_cB22s,deaths_cB23s,deaths_cB24s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(cB2s))) %>%
  ggplot() +
  geom_line(aes(x=cB2,y=deaths,col=p_type)) +
  scale_x_continuous(name=TeX('quadratic reduction term for effect of $B$ on $C$ ($c_{B2}$)'), expand=0) +
  scale_y_continuous(name=TeX('deaths after first wave'), expand=0, limits=c(4000,5000)) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes NPI use ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))





I_max_t_cB21 <- function(cB2) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p1=p1_ex, cB2=cB2))[,'I'])])
}
I_max_t_cB22 <- function(cB2) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p2=p2_ex, cB2=cB2))[,'I'])])
}
I_max_t_cB23 <- function(cB2) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p3=p3_ex, cB2=cB2))[,'I'])])
}
I_max_t_cB24 <- function(cB2) {
  return(ts[which.max(ode(rootfun=function(t, y, parms) max(y[8]+y[9]-init_infections,0)+max(365-t,0), y0, ts, odes, parms_freeriding(p4=p4_ex, cB2=cB2))[,'I'])])
}


I_max_t_cB21s <- sapply(cB2s, I_max_t_cB21)
I_max_t_cB22s <- sapply(cB2s, I_max_t_cB22)
I_max_t_cB23s <- sapply(cB2s, I_max_t_cB23)
I_max_t_cB24s <- sapply(cB2s, I_max_t_cB24)

data.frame(cB2=cB2s,
           I_max_t=c(I_max_t_cB21s,I_max_t_cB22s,I_max_t_cB23s,I_max_t_cB24s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(cB2s))) %>%
  ggplot() +
  geom_line(aes(x=cB2,y=I_max_t,col=p_type)) +
  scale_x_continuous(name=TeX('quadratic reduction term for effect of $B$ on $C$ ($c_{B2}$)'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infected population (days)'), expand=0, limits=c(0,350)) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes NPI use ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))



