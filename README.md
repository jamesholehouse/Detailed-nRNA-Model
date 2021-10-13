# Detailed-nRNA-Model

Here I detail the SSA for a detailed model of nascent RNA (nRNA), that incorporates:
- nRNA production rate `rho`,
- **delayed** nRNA degradation at time `tau` , 
- on/off gene switching at rates `\sigma_b` and `\sigma_u`,
- hopping rate between the cell cycle stages `k`,
- number of cell cycle stages `N`,
- the fraction denoting the gene dosage compensation after gene replication `kappa`.

The output of the SSA here is the number of nRNA at some time `t`. Note that the gene is assumed to replicate at stage `int(N/2)`, after which the transcription rate of each gene becomes `kappa`x`rho`. The different stages of the cell cycle pass with exponentially distributed waiting times at rate `k`. After cell cycle stage N (need a diagram for this) cell division occurs and the only one (random) copy of the gene remains, keeping only the nRNA that gene produced.

*Could later on incorporate an RNAp description, being produced at some rate and only removed via dilution at the binomial partitioning that occurs at cell division.*
