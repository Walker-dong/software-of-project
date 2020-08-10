Sn = 1 * 10**5
Sq = 0
En = 109
Eq = 0
Im = 10
Io = 2
Is = 0
Iq = 21
A = 10
R = 0


C = 0
Release_quarantine = 0.00

# population
S = Sn + Sq
E = En + Eq
In = Im + Io + Is  
I = In + Iq 
N = S + E + I + R + A + C

# simuation Time / Day
T = 170
# susceptiable ratio
s = np.zeros([T])
# exposed ratio
e = np.zeros([T])
# infective ratio
i = np.zeros([T])
# remove ratio
r = np.zeros([T])


# contact rate between S and E'贝塔'
SE_contactRate = 9.08 * (10 ** -5)
#tracking probility of S'p'
S_trackProbility = 0.50
#Composition ratio of contacts‘司格马’
Infect_ratio = 0.50
#Relative to infected persons, ‘玉璞西戎’
# the transmission coefficient of the incubation group is between 0-1
E_infectRate = 0.50
#The rate at which people ‘阿拉法’
# in the incubation period progress to become infected
E_Transform = 1/5.2
#Symptomatic infections as a percentage of infected'倒钩'
I_symptomRatio = 0.70
#Ratio of minor infections'p1'
I_minorRatio = 0.50
#Ratio of ordinary infections'p2'
I_ordinaryRatio = 0.35
#Raio of serious infections 'p3'
I_seriousRatio = 0.15
#The rate of removal of infected persons is '伽马mos'
# the reciprocal of the course of response
I_removalRate = Ia_removalRate = Iq_removalRate = 0.10
#Probability of being diagnosed and isolated'sing他1'
Rm_quarantineRatio = 0.50
#Probability of being diagnosed and isolated'sing他2'
Ro_quarantineRatio = 0.80
#Probability of being diagnosed and isolated'sing他3'
Rs_quarantineRatio = 1.00
#communicator



# initial infective people
i[0] = 10.0 / N
s[0] = 1e7 / N
e[0] = 40.0 / N
for t in range(T-1):

    Sn[t+1] = -(1 - S_trackProbility) * SE_contactRate * Infect_ratio * (E_infectRate * En[t] + In[t] + A[t]) * Sn[t] \
        -S_trackProbility * SE_contactRate * (1 - Infect_ratio) * (E_infectRate * En[t] + In[t] + A[t]) * Sn[t] \
            -S_trackProbility * SE_contactRate * Infect_ratio * (E_infectRate * En[t] + In[t] + A[t]) * Sn[t] \
                  +Release_quarantine * Sq[t]
    En[t+1] = (1 - S_trackProbility) * SE_contactRate * Infect_ratio * (E_infectRate * En[t] + In[t] + A[t]) * Sn[t] \
        -E_Transform * En[t]
    Sq[t+1] = S_trackProbility * SE_contactRate * Infect_ratio * (E_infectRate * En[t] + In[t] + A[t]) * Sn[t] \
        -Release_quarantine * Sq[t]
    Eq[t+1] = S_trackProbility * SE_contactRate * (1 - Infect_ratio) * (E_infectRate * En[t] + In[t] + A[t]) * Sn[t] \
        -E_Transform * Eq[t]
    Io[t+1] = E_Transform * I_symptomRatio * I_minorRatio * En[t] - I_removalRate * Io[t]
    Im[t+1] = E_Transform * I_symptomRatio * I_ordinaryRatio * En[t] - I_removalRate * Im[t]
    Is[t+1] = E_Transform * I_symptomRatio * I_seriousRatio * En[t] - I_removalRate * Is[t]
    Iq[t+1] = E_Transform * Eq + Rm_quarantineRatio * I_removalRate * Io[t] \
        + Ro_quarantineRatio * I_removalRate * Im[t] \
            + Rs_quarantineRatio * I_removalRate * Is[t] \
                -Iq_removalRate * Iq[t]
    R[t+1] = Ia_removalRate * A[t] + (1 - Rm_quarantineRatio) * I_removalRate * Im[t] \
        + (1 - Ro_quarantineRatio) * I_removalRate * Io[t] \
            (1 - Rs_quarantineRatio) * I_removalRate * Is[t]
    A[t+1] = E_Transform * (1 - I_symptomRatio) * En[t] - Ia_removalRate * A[t]


    fig, ax = plt.subplots(figsize=(10,6))
ax.plot(s, c='b', lw=2, label='S')
ax.plot(e, c='orange', lw=2, label='E')
ax.plot(i, c='r', lw=2, label='I')
ax.plot(r, c='g', lw=2, label='R')
ax.set_xlabel('Day',fontsize=20)
ax.set_ylabel('Infective Ratio', fontsize=20)
ax.grid(1)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend();