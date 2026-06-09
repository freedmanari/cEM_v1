require(deSolve)
require(tidyverse)
require(latex2exp)

odes <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    beta <- R0 * mu * (lambda + mu) * (alpha + gamma + mu) / (Gamma * lambda)
      
    dS_o <- Gamma - beta * ifelse(t > Tp, 1 - eps1 * p1, 1) * S_o * I_o - mu * S_o
    dE_o <- beta * ifelse(t > Tp, 1 - eps1 * p1, 1) * S_o * I_o - lambda * E_o - mu * E_o
    dI_o <- lambda * E_o - alpha * I_o - gamma * I_o - ifelse(t > Tp, eps2 * p2 * I_o, 0) - mu * I_o
    dD_o <- alpha * I_o
    dR_o <- gamma * I_o + gammaQ * Q_o - mu * R_o
    dQ_o <- ifelse(t > Tp, eps2 * p2 * I_o, 0) - gammaQ * Q_o
    dS <- Gamma - beta * (1 - eB * B - ifelse(t > Tp, eps1 * p1, 0)) * S * I - mu * S
    dE <- beta * (1 - eB * B - ifelse(t > Tp, eps1 * p1, 0)) * S * I - lambda * E - mu * E
    dI <- lambda * E - alpha * I - gamma * I - (kB * B + kC * C + ifelse(t > Tp, eps2 * p2, 0)) * I - mu * I
    dD <- alpha * I + alphaQ * Q
    dR <- gamma * I + gammaQ * Q - mu * R
    dQ <- (kB * B + kC * C + ifelse(t > Tp, eps2 * p2, 0)) * I - gammaQ * Q - alphaQ * Q - mu * Q
    dC <- ((cB1 - cB2 * B) * B + aD * D + aQ * Q + ifelse(t > Tp, eps3 * p3, 0)) * (1 - C) - deltaC * C
    dB <- (bC * C + ifelse(t > Tp, eps4 * p4, 0)) * (1 - B) - deltaB * B
    
    return(list(c(dS_o, dE_o, dI_o, dD_o, dR_o, dQ_o,
                  dS, dE, dI, dD, dR, dQ, dC, dB)))
  })
}

N <- 10^5 #population size
init_infections <- 10
y0 <- c(S_o=N-init_infections, E_o=init_infections, I_o=0, D_o=0, R_o=0, Q_o=0,
        S=N-init_infections, E=init_infections, I=0, D=0, R=0, Q=0, C=0, B=0)

parms_no_policy <- c(mu=1/(60*365), Gamma=10^5/(60*365), lambda=.3, gamma=.1, gammaQ=.05, alpha=.01, alphaQ=.015, delta=0, Tp=0,
                     R0=2, eB=.5, eps1=0, kB=0, kC=.1, eps2=0, deltaC=.1, cB1=.1, cB2=0,
                     aD=1/10^5, aQ=1/10^4, eps3=0, deltaB=.05, bC=.1, eps4=0,
                     p1=0, p2=0, p3=0, p4=0)
ts <- seq(0,365,.1)

D365_full_p1 <- function(eps1) {
  parms <- parms_no_policy
  parms['eps1'] <- eps1
  parms['p1'] <- 1
  
  return(unname(ode(y0, ts, odes, parms)[length(ts), 'D']))
}
D365_full_p2 <- function(eps2) {
  parms <- parms_no_policy
  parms['eps2'] <- eps2
  parms['p2'] <- 1
  
  return(unname(ode(y0, ts, odes, parms)[length(ts), 'D']))
}
D365_full_p3 <- function(eps3) {
  parms <- parms_no_policy
  parms['eps3'] <- eps3
  parms['p3'] <- 1
  
  return(unname(ode(y0, ts, odes, parms)[length(ts), 'D']))
}
D365_full_p4 <- function(eps4) {
  parms <- parms_no_policy
  parms['eps4'] <- eps4
  parms['p4'] <- 1
  
  return(unname(ode(y0, ts, odes, parms)[length(ts), 'D']))
}

D365_no_policy <- D365_full_p1(0)

D365_red_eps1 <- function(red) {
  return(uniroot(function(eps1) D365_full_p1(eps1) - (1-red) * D365_no_policy, c(0, 1))$root)
}
D365_red_eps2 <- function(red) {
  return(uniroot(function(eps2) D365_full_p2(eps2) - (1-red) * D365_no_policy, c(0, 1))$root)
}
D365_red_eps3 <- function(red) {
  return(uniroot(function(eps3) D365_full_p3(eps3) - (1-red) * D365_no_policy, c(0, 1))$root)
}
D365_red_eps4 <- function(red) {
  return(uniroot(function(eps4) D365_full_p4(eps4) - (1-red) * D365_no_policy, c(0, .15))$root)
}


# default policy efficacy parameters, calibrated to reduce deaths by 75% after 1 year
red_def <- .75
eps1_def <- D365_red_eps1(red_def)
eps2_def <- D365_red_eps2(red_def)
eps3_def <- D365_red_eps3(red_def)
eps4_def <- D365_red_eps4(red_def)

parms_no_policy['eps1'] <- eps1_def
parms_no_policy['eps2'] <- eps2_def
parms_no_policy['eps3'] <- eps3_def
parms_no_policy['eps4'] <- eps4_def


parms_policy <- function(p1=0, p2=0, p3=0, p4=0, Tp=0) {
  parms <- parms_no_policy
  parms['p1'] <- p1
  parms['p2'] <- p2
  parms['p3'] <- p3
  parms['p4'] <- p4
  parms['Tp'] <- Tp
  
  return(parms)
}

D365_p1 <- function(p1) {
  return(unname(ode(y0, ts, odes, parms_policy(p1=p1))[length(ts), 'D']))
}
D365_p2 <- function(p2) {
  return(unname(ode(y0, ts, odes, parms_policy(p2=p2))[length(ts), 'D']))
}
D365_p3 <- function(p3) {
  return(unname(ode(y0, ts, odes, parms_policy(p3=p3))[length(ts), 'D']))
}
D365_p4 <- function(p4) {
  return(unname(ode(y0, ts, odes, parms_policy(p4=p4))[length(ts), 'D']))
}

ps <- seq(0,1,.01)
D365_p1s <- sapply(ps, D365_p1)
D365_p2s <- sapply(ps, D365_p2)
D365_p3s <- sapply(ps, D365_p3)
D365_p4s <- sapply(ps, D365_p4)

data.frame(p=ps,
           D365=c(D365_p1s,D365_p2s,D365_p3s,D365_p4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(ps))) %>%
  ggplot() +
  geom_line(aes(x=p,y=D365,col=p_type)) +
  scale_x_continuous(name=TeX('policy strength ($p_1$, $p_2$, $p_3$, or $p_4$)'), expand=0) +
  scale_y_continuous(name=TeX('deaths after 1 year'), expand=0, limits=c(0,6000 )) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes transmission-reducing behavior ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))







I_max_p1 <- function(p1) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p1=p1))[,'I'])])
}
I_max_p2 <- function(p2) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p2=p2))[,'I'])])
}
I_max_p3 <- function(p3) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p3=p3))[,'I'])])
}
I_max_p4 <- function(p4) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p4=p4))[,'I'])])
}


I_max_p1s <- sapply(ps, I_max_p1)
I_max_p2s <- sapply(ps, I_max_p2)
I_max_p3s <- sapply(ps, I_max_p3)
I_max_p4s <- sapply(ps, I_max_p4)

data.frame(p=ps,
           I_max=c(I_max_p1s,I_max_p2s,I_max_p3s,I_max_p4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(ps))) %>%
  ggplot() +
  geom_line(aes(x=p,y=I_max,col=p_type)) +
  scale_x_continuous(name=TeX('policy strength ($p_1$, $p_2$, $p_3$, or $p_4$)'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infection prevalence (days)'), expand=0, limits=c(0,365)) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes transmission-reducing behavior ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))








# policy strengths for the following example plots, calibrated to reduce deaths by 40% after 1 year
red_ex <- .5
p1_ex <- D365_red_eps1(red_ex) / eps1_def
p2_ex <- D365_red_eps2(red_ex) / eps2_def
p3_ex <- D365_red_eps3(red_ex) / eps3_def
p4_ex <- D365_red_eps4(red_ex) / eps4_def


## example policy starts at t=0

out_no_policy <-
  ode(y0, ts, odes, parms_no_policy) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=FALSE)

out_p1_ex <-
  ode(y0, ts, odes, parms_policy(p1=p1_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p1')
out_p2_ex <-
  ode(y0, ts, odes, parms_policy(p2=p2_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p2')
out_p3_ex <-
  ode(y0, ts, odes, parms_policy(p3=p3_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p3')
out_p4_ex <-
  ode(y0, ts, odes, parms_policy(p4=p4_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p4')



out_p_exs <-
  rbind(out_p1_ex,
        out_p2_ex,
        out_p3_ex,
        out_p4_ex)
policy_labels <- c('policy reducing\ntransmission directly',
            'policy increasing\ntesting rate',
            'policy increasing\nrisk perception',
            'policy promoting\ntransmission-reducing behavior')
names(policy_labels) <- c('p1','p2','p3','p4')

out_p_exs %>%
  filter(var %in% c('I','I_o')) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value, color=interaction(policy,behavior), lty=interaction(policy,behavior), linewidth=interaction(policy,behavior))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(name='infection prevalence', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('red','red','blue','blue'), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  scale_linetype_manual(name='', values=c('11','solid','11','solid'), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  scale_linewidth_manual(name='', values=c(.8,.4,.8,.4), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11))



out_p_exs %>%
  filter(var %in% c('I','C','B','Q')) %>% 
  mutate(var=factor(var, levels=c('I','C','B','Q')),
         value=ifelse(var %in% c('I','Q'), value/8000, value)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value*8000, color=interaction(policy,var), lty=interaction(policy,var), linewidth=interaction(policy,var))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(name='infection prevalence', expand=0, limits=c(0,8000),
                     sec.axis = sec_axis( trans=~./8000, name="risk perception (C) or\nbehavioral adoption (B)")) +
  scale_color_manual(name='', values=c('blue','blue','forestgreen','forestgreen','darkorange','darkorange','violet','violet'), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  scale_linetype_manual(name='', values=c('11','solid','11','solid','11','solid','11','solid'), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  scale_linewidth_manual(name='', values=c(.7,.5,.7,.5,.7,.5,.7,.5), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11),
        panel.spacing = unit(1, "lines"))















## example policy starts at t=50



out_late_p1_ex <-
  ode(y0, ts, odes, parms_policy(p1=p1_ex, Tp=50)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p1')
out_late_p2_ex <-
  ode(y0, ts, odes, parms_policy(p2=p2_ex, Tp=50)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p2')
out_late_p3_ex <-
  ode(y0, ts, odes, parms_policy(p3=p3_ex, Tp=50)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p3')
out_late_p4_ex <-
  ode(y0, ts, odes, parms_policy(p4=p4_ex, Tp=50)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p4')



out_late_p_exs <-
  rbind(out_late_p1_ex,
        out_late_p2_ex,
        out_late_p3_ex,
        out_late_p4_ex)

out_late_p_exs %>%
  filter(var %in% c('I','I_o')) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value, color=interaction(policy,behavior), lty=interaction(policy,behavior), linewidth=interaction(policy,behavior))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(name='infection prevalence', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('red','red','blue','blue'), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  scale_linetype_manual(name='', values=c('11','solid','11','solid'), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  scale_linewidth_manual(name='', values=c(.8,.4,.8,.4), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11))

out_late_p_exs %>%
  filter(var %in% c('I','C','B')) %>% 
  mutate(var=factor(var, levels=c('I','C','B')),
         value=ifelse(var=='I', value/8000, value)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value*8000, color=interaction(policy,var), lty=interaction(policy,var), linewidth=interaction(policy,var))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(name='infection prevalence', expand=0, limits=c(0,8000),
                     sec.axis = sec_axis( trans=~./8000, name="risk perception (C) or\nbehavioral adoption (B)")) +
  scale_color_manual(name='', values=c('blue','blue','forestgreen','forestgreen','darkorange','darkorange'), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  scale_linetype_manual(name='', values=c('11','solid','11','solid','11','solid'), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  scale_linewidth_manual(name='', values=c(.7,.5,.7,.5,.7,.5), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11),
        panel.spacing = unit(1, "lines"))






# example policy starts at t=100

out_very_late_p1_ex <-
  ode(y0, ts, odes, parms_policy(p1=p1_ex, Tp=100)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p1')
out_very_late_p2_ex <-
  ode(y0, ts, odes, parms_policy(p2=p2_ex, Tp=100)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p2')
out_very_late_p3_ex <-
  ode(y0, ts, odes, parms_policy(p3=p3_ex, Tp=100)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p3')
out_very_late_p4_ex <-
  ode(y0, ts, odes, parms_policy(p4=p4_ex, Tp=100)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  rbind(out_no_policy) %>% 
  mutate(policy_type='p4')




out_very_late_p_exs <-
  rbind(out_very_late_p1_ex,
        out_very_late_p2_ex,
        out_very_late_p3_ex,
        out_very_late_p4_ex)

out_very_late_p_exs %>%
  filter(var %in% c('I','I_o')) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value, color=interaction(policy,behavior), lty=interaction(policy,behavior), linewidth=interaction(policy,behavior))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(name='infection prevalence', expand=0, limits=c(0,12000)) +
  scale_color_manual(name='', values=c('red','red','blue','blue'), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  scale_linetype_manual(name='', values=c('11','solid','11','solid'), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  scale_linewidth_manual(name='', values=c(.8,.4,.8,.4), labels=c('without behavior, without policy','without behavior, with policy','with behavior, without policy','with behavior, with policy')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11))

out_very_late_p_exs %>%
  filter(var %in% c('I','C','B')) %>% 
  mutate(var=factor(var, levels=c('I','C','B')),
         value=ifelse(var=='I', value/8000, value)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value*8000, color=interaction(policy,var), lty=interaction(policy,var), linewidth=interaction(policy,var))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(name='infection prevalence', expand=0, limits=c(0,8000),
                     sec.axis = sec_axis( trans=~./8000, name="risk perception (C) or\nbehavioral adoption (B)")) +
  scale_color_manual(name='', values=c('blue','blue','forestgreen','forestgreen','darkorange','darkorange'), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  scale_linetype_manual(name='', values=c('11','solid','11','solid','11','solid'), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  scale_linewidth_manual(name='', values=c(.7,.5,.7,.5,.7,.5), labels=c('I, without policy','I, with policy','C, without policy','C, with policy','B, without policy','B, with policy')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11),
        panel.spacing = unit(1, "lines"))







Tps <- seq(0,200)


D365_Tp1 <- function(Tp) {
  return(unname(ode(y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[length(ts), 'D']))
}
D365_Tp2 <- function(Tp) {
  return(unname(ode(y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[length(ts), 'D']))
}
D365_Tp3 <- function(Tp) {
  return(unname(ode(y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[length(ts), 'D']))
}
D365_Tp4 <- function(Tp) {
  return(unname(ode(y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[length(ts), 'D']))
}

D365_Tp1s <- sapply(Tps, D365_Tp1)
D365_Tp2s <- sapply(Tps, D365_Tp2)
D365_Tp3s <- sapply(Tps, D365_Tp3)
D365_Tp4s <- sapply(Tps, D365_Tp4)

data.frame(Tp=Tps,
           D365=c(D365_Tp1s,D365_Tp2s,D365_Tp3s,D365_Tp4s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(Tps))) %>%
  ggplot() +
  geom_line(aes(x=Tp,y=D365,col=p_type)) +
  scale_x_continuous(name=TeX('policy start time'), expand=0) +
  scale_y_continuous(name=TeX('deaths after 1 year'), expand=0, limits=c(0,6000 )) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes transmission-reducing behavior ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))







I_max_Tp1 <- function(Tp) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p1=p1_ex, Tp=Tp))[,'I'])])
}
I_max_Tp2 <- function(Tp) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p2=p2_ex, Tp=Tp))[,'I'])])
}
I_max_Tp3 <- function(Tp) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p3=p3_ex, Tp=Tp))[,'I'])])
}
I_max_Tp4 <- function(Tp) {
  return(ts[which.max(ode(y0, ts, odes, parms_policy(p4=p4_ex, Tp=Tp))[,'I'])])
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
  scale_x_continuous(name=TeX('policy start time'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infection prevalence (days)'), expand=0, limits=c(0,365)) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes transmission-reducing behavior ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))











combined_ps <- data.frame()

for (p12 in seq(0,.99,.01)) {
  print(p12)
  for (p34 in seq(0,round(.99-p12,2),.01)) {
    for (Tp in c(0,50,100)) {
      out <- ode(y0, ts, odes, parms_policy(p1=p12/2, p2=p12/2, p3=p34/2, p4=p34/2, Tp=Tp))
      
      combined_ps <-
        combined_ps %>% 
        rbind(data.frame(p12 = round(p12,2),
                         p34 = round(p34,2),
                         Tp = Tp,
                         D365 = unname(out[length(ts), 'D']),
                         I_max = ts[which.max(out[,'I'])]))
    }
  }
}



Tp_labels <- c('started at very beginning\nof epidemic (t=0)',
               'started early\nin epidemic (t=50)',
               'started late\nin epidemic (t=100)')
names(Tp_labels) <- c('0','50','100')

combined_ps %>% 
  ggplot() +
  geom_raster(aes(x=p12,y=p34,fill=D365), hjust=1, vjust=1) +
  facet_wrap(~Tp, labeller = as_labeller(Tp_labels), scales='free') +
  scale_fill_gradientn(name='deaths after\n1 year',colors=rev(rainbow(7)[-7]),
                       guide=guide_colorbar(frame.colour = 'black', ticks.colour = 'black',title.hjust=.5)) +
  scale_x_continuous(name=expression(atop(NA,atop(textstyle('strength of transmission-controlling policies'),textstyle('(split equally between ' * p[1]*' and '*p[2]*')')))),expand=0) +
  scale_y_continuous(name=expression(atop(NA,atop(textstyle('strength of behavior-controlling policies'),textstyle('(split equally between ' * p[3]*' and '*p[4]*')')))),expand=0) +
  theme(axis.line=element_line(),
        legend.title=element_text(size=10),
        legend.position='right',
        panel.background = element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=11))


combined_ps %>% 
  ggplot() +
  geom_raster(aes(x=p12,y=p34,fill=I_max), hjust=1, vjust=1) +
  facet_wrap(~Tp, labeller = as_labeller(Tp_labels), scales='free') +
  scale_fill_gradientn(name='peak time of\ninfection prevalence\n(days)',colors=rainbow(7)[-7],
                       guide=guide_colorbar(frame.colour = 'black', ticks.colour = 'black',title.hjust=.5)) +
  scale_x_continuous(name=expression(atop(NA,atop(textstyle('strength of transmission-controlling policies'),textstyle('(split equally between ' * p[1]*' and '*p[2]*')')))),expand=0) +
  scale_y_continuous(name=expression(atop(NA,atop(textstyle('strength of behavior-controlling policies'),textstyle('(split equally between ' * p[3]*' and '*p[4]*')')))),expand=0) +
  theme(axis.line=element_line(),
        legend.title=element_text(size=10),
        legend.position='right',
        panel.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        strip.background=element_blank(),
        strip.text=element_text(size=11))















parms_freeriding <- function(cB2=0, p1=0, p2=0, p3=0, p4=0, Tp=0) {
  parms <- parms_policy(p1,p2,p3,p4,Tp)
  parms['cB2'] <- cB2
  
  return(parms)
}


out_p1_freeriding <-
  ode(y0, ts, odes, parms_freeriding(cB2=.15, p1=p1_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p1')
out_p2_freeriding <-
  ode(y0, ts, odes, parms_freeriding(cB2=.15, p2=p2_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p2')
out_p3_freeriding <-
  ode(y0, ts, odes, parms_freeriding(cB2=.15, p3=p3_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p3')
out_p4_freeriding <-
  ode(y0, ts, odes, parms_freeriding(cB2=.15, p4=p4_ex)) %>%
  as.data.frame %>% 
  pivot_longer(cols=-'time',names_to='var',values_to='value') %>%
  mutate(behavior=!str_detect(var, '_o'),
         policy=TRUE) %>%
  mutate(policy_type='p4')







rbind(cbind(out_p1_ex, freeriding=FALSE),
      cbind(out_p2_ex, freeriding=FALSE),
      cbind(out_p3_ex, freeriding=FALSE),
      cbind(out_p4_ex, freeriding=FALSE),
      cbind(out_p1_freeriding, freeriding=TRUE),
      cbind(out_p2_freeriding, freeriding=TRUE),
      cbind(out_p3_freeriding, freeriding=TRUE),
      cbind(out_p4_freeriding, freeriding=TRUE)) %>%
  filter(var %in% c('I','C','B'), policy) %>% 
  mutate(var=factor(var, levels=c('I','C','B')),
         value=ifelse(var=='I', value/8000, value)) %>% 
  ggplot() +
  geom_line(aes(x=time, y=value*8000, color=interaction(freeriding,var), lty=interaction(freeriding,var), linewidth=interaction(freeriding,var))) +
  facet_wrap(~ policy_type, labeller = as_labeller(policy_labels), scales='free') +
  scale_x_continuous(name='time (days)', expand=0) +
  scale_y_continuous(name='infection prevalence', expand=0, limits=c(0,8000),
                     sec.axis = sec_axis( trans=~./8000, name="risk perception (C) or\nbehavioral adoption (B)")) +
  scale_color_manual(name='', values=c('blue','blue','forestgreen','forestgreen','darkorange','darkorange'), labels=c('I, with policy,\nwithout freeriding','I, with policy,\nwith freeriding','C, without policy,\nwithout freeriding','C, with policy,\nwith freeriding','B, without policy,\nwithout freeriding','B, with policy,\nwith freeriding')) +
  scale_linetype_manual(name='', values=c('solid','24','solid','24','solid','24'), labels=c('I, with policy,\nwithout freeriding','I, with policy,\nwith freeriding','C, without policy,\nwithout freeriding','C, with policy,\nwith freeriding','B, without policy,\nwithout freeriding','B, with policy,\nwith freeriding')) +
  scale_linewidth_manual(name='', values=c(.5,.7,.5,.7,.5,.7), labels=c('I, with policy,\nwithout freeriding','I, with policy,\nwith freeriding','C, without policy,\nwithout freeriding','C, with policy,\nwith freeriding','B, without policy,\nwithout freeriding','B, with policy,\nwith freeriding')) +
  theme_classic() +
  theme(strip.background=element_blank(),
        strip.text=element_text(size=11),
        legend.key.spacing.y = unit(0.3, "cm"))







cB2s <- seq(0,.2,.002)


D365_cB21 <- function(cB2) {
  return(unname(ode(y0, ts, odes, parms_freeriding(p1=p1_ex, cB2=cB2))[length(ts), 'D']))
}
D365_cB22 <- function(cB2) {
  return(unname(ode(y0, ts, odes, parms_freeriding(p2=p2_ex, cB2=cB2))[length(ts), 'D']))
}
D365_cB23 <- function(cB2) {
  return(unname(ode(y0, ts, odes, parms_freeriding(p3=p3_ex, cB2=cB2))[length(ts), 'D']))
}
D365_cB24 <- function(cB2) {
  return(unname(ode(y0, ts, odes, parms_freeriding(p4=p4_ex, cB2=cB2))[length(ts), 'D']))
}

D365_cB21s <- sapply(cB2s, D365_cB21)
D365_cB22s <- sapply(cB2s, D365_cB22)
D365_cB23s <- sapply(cB2s, D365_cB23)
D365_cB24s <- sapply(cB2s, D365_cB24)

data.frame(cB2=cB2s,
           D365=c(D365_cB21s,D365_cB22s,D365_cB23s,D365_cB24s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(cB2s))) %>%
  ggplot() +
  geom_line(aes(x=cB2,y=D365,col=p_type)) +
  scale_x_continuous(name=TeX('quadratic reduction term for effect of $B$ on $C$ ($c_{B2}$)'), expand=0) +
  scale_y_continuous(name=TeX('deaths after 1 year'), expand=0, limits=c(0,6000 )) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes transmission-reducing behavior ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))







I_max_cB21 <- function(cB2) {
  return(ts[which.max(ode(y0, ts, odes, parms_freeriding(p1=p1_ex, cB2=cB2))[,'I'])])
}
I_max_cB22 <- function(cB2) {
  return(ts[which.max(ode(y0, ts, odes, parms_freeriding(p2=p2_ex, cB2=cB2))[,'I'])])
}
I_max_cB23 <- function(cB2) {
  return(ts[which.max(ode(y0, ts, odes, parms_freeriding(p3=p3_ex, cB2=cB2))[,'I'])])
}
I_max_cB24 <- function(cB2) {
  return(ts[which.max(ode(y0, ts, odes, parms_freeriding(p4=p4_ex, cB2=cB2))[,'I'])])
}


I_max_cB21s <- sapply(cB2s, I_max_cB21)
I_max_cB22s <- sapply(cB2s, I_max_cB22)
I_max_cB23s <- sapply(cB2s, I_max_cB23)
I_max_cB24s <- sapply(cB2s, I_max_cB24)

data.frame(cB2=cB2s,
           I_max=c(I_max_cB21s,I_max_cB22s,I_max_cB23s,I_max_cB24s),
           p_type=rep(c('p1','p2','p3','p4'),
                      each=length(cB2s))) %>%
  ggplot() +
  geom_line(aes(x=cB2,y=I_max,col=p_type)) +
  scale_x_continuous(name=TeX('quadratic reduction term for effect of $B$ on $dC/dt$ ($c_{B2}$)'), expand=0) +
  scale_y_continuous(name=TeX('peak time of infection prevalence (days)'), expand=0, limits=c(0,365)) +
  scale_color_discrete(name='policy that...', labels=c(TeX('reduces transmission rate directly ($p_1$)'),
                                                       TeX('increases testing rate ($p_2$)'),
                                                       TeX('increases risk perception ($p_3$)'),
                                                       TeX('promotes transmission-reducing behavior ($p_4$)'))) +
  theme_classic() +
  theme(legend.text=element_text(size=10))


Bs <- seq(0,1,.01)
cB1 <- parms_no_policy['cB1']
data.frame(B=Bs,
           cB2=factor(rep(c(0,.05,.1,.15,.2), each=length(Bs))),
           effect_of_B_on_dCdt=c(cB1 * Bs,
                                  (cB1-.05*Bs) * Bs,
                                  (cB1-.1*Bs) * Bs,
                                  (cB1-.15*Bs) * Bs,
                                  (cB1-.2*Bs) * Bs)) %>% 
  ggplot() +
  geom_line(aes(x=B,y=effect_of_B_on_dCdt,linetype=cB2)) +
  scale_x_continuous(name='behavioral adoption (B)', expand=0) +
  scale_y_continuous(name=expression(atop(NA,atop('effect of B on dC/dt,',(c[B1]-c[B2]*B)*B))), expand=0) +
  scale_linetype_discrete(name=TeX('quadratic reduction term for effect of $B$ on $dC/dt$ ($c_{B2}$)'),
                          guide=guide_legend(title.position="top")) +
  theme_classic() +
  theme(legend.position = 'bottom',
        legend.key.spacing.x = unit(0.5, "cm"),
        legend.key.width = unit(1.5, "line"),
        axis.title.y = element_text(size=15,angle=0,vjust=.5))








