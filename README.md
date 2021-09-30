# Sandpiles

Code related to simulating different sandpile models. 


## Intro 

<!-- (introducing the theory of the sandpile models, copied from senior thesis?) -->

A broad class of systems exhibit scale-free phenomena in the form of avalanches of relaxation events. It is an open question whether these systems exhibit self-organized criticality (SOC) or if they are described by critical phenomena. In SOC systems the critical point is a dynamical attractor that is independent of external tuning or initial conditions. 

The code below simulates a simple cellular automata sandpile models which have been  proposed  as examples of SOC.This model was introduced as the Manna sandpile model and is generalized for a more broader study. 

## Model description 

<!-- (describe the general idea of the sandpile model)  -->

The Abelian sandpile model (AM) is cellular automata similar to the original BTW model. Typically defined on a $d$-dimensional hyper-cubic lattice. Each site, $\mathbf{n}$ is labelled by its height, $h_{\mathbf{n}}$. The dynamics can be described by a toppling matrix, $\Delta$ with matrix elements $\Delta_{\mathbf{n}, \mathbf{n}'}$ obeying the following rules: 
```math
\begin{equation}
	\Delta_{\mathbf{n}, \mathbf{n}} = h_\mathbf{n}^c > 0; \quad 
	\Delta_{\mathbf{n}, \mathbf{n}'} < 0; \quad
	\sum_{\mathbf{n}' = 1}^N \Delta_{\mathbf{n}, \mathbf{n}} \geq 0
\end{equation}
``` 
where $\mathbf{n}, \mathbf{n}' \in \{1,...,N\}$, $N$ is the number of sites on the lattice, and $h_\mathbf{n}^c > 0$ is the critical height threshold. 

![](Model_Images/lattice_init.png | width=100)
![](Model_Images/lattice_0.png | width=100)
![](Model_Images/lattice_1.png | width=100)
![](Model_Images/lattice_2.png | width=100)

### Describe differences between models

(how do the models differ) 

### Describe differences in code of thee same model
