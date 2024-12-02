# Aqueous species formulation, using Miron et al., 2017 database.


```math
G_{aq} = \sum_{i=1}^{n^{em}} n_i \mu_i
```

where $n_i$ the molar fraction and $\mu_i$ is the chemical potential of the species $i$.

```math
\mu_i  = G_i^0 + RT ln(\gamma_i m_i)
```

where $m_i$ is the molal fraction (!?)

And then develops into:



```math
\mu_i  = G_i^0 + RT( ln(\gamma_i) + ln(m_i))
```

```math
\mu_i  = G_i^0 + RT( ln(\gamma_i) + ln(\frac{n_i}{n_{H2O} M_{H2O}}))
```

```math
\mu_i  = G_i^0 + RT( ln(\gamma_i) + ln(\frac{n_i}{18.015 * n_{H2O}}))
```


---

For solute (aqueous species, excluding water), and assuming $\gamma_i$ = 1 (for now)

```math
\mu_i  = G_i^0 + RT ln(m_i)
```

where 

```math
m_i  = \frac{n_i}{m_{H2O}}
```

so here $n_i$ is again the molar fraction and $m_{H2O}$ is mass of water? Should it be water density?

I could find the following definition:

```math
m_i  = \frac{n_i}{n_{H2O} M_{H2O}}
```
and here $n_i$ is the molarity and $M_{H2O}$ is the Molar mass of water (18.015). However, here density of water like we discussed is not appearing, is that right?

Should it simply be:

```math
m_i  = \frac{n_i}{\rho_{H2O}}
```

---

For solvant (water):

```math
\mu_{H2O}  = G_{H2O }^0 + RT ln(\frac{1}{M_{H2O}})
```

!?
