# Gamma-Adaptive Initial-Projection Algorithm

## Overview

The Gamma-adaptive method automatically constructs initial Wannier projections from Gamma-point Bloch wavefunctions.

## Contents

- [1. Wannier initial guesses from Gamma-point wavefunctions](#1-wannier-initial-guesses-from-gamma-point-wavefunctions)
  - [1.1 Disentanglement](#11-disentanglement-selection-of-the-optimal-gamma-subspace)
  - [1.2 Jacobi rotation](#12-jacobi-rotation-foster-boys-localization-of-the-gamma-point-initial-guesses)
  - [1.3 Locate center](#13-locate-center-orbital-centers-of-the-gamma-point-initial-guesses)
  - [1.4 Phase fix](#14-phase-fix-global-phase-convention-for-the-gamma-point-initial-guesses)
- [2. Gamma-adaptive .amn matrix](#2-gamma-adaptive-amn-matrix)

## 1. Wannier initial guesses from Gamma-point wavefunctions

### 1.1 Disentanglement: selection of the optimal Gamma subspace

A $J$-dimensional subspace is selected from the active Bloch eigenstates at Gamma, where $J$ is the number of Wannier orbitals.

$$\ket{\phi_{n\Gamma}^{(S)}}=\sum_{m=1}^{N}\ket{\psi_{m\Gamma}} S_{mn}, n=1,\ldots,J.$$

Here,

$$S\in\mathbb C^{N\times J}, S^\dagger S=I_J.$$

The columns of $S$ provide an orthonormal basis for the selected subspace, with the corresponding projector

$$P=SS^\dagger.$$

The objective of disentanglement is to determine the projector $P$.

#### 1.1.1 Overlaps of Gamma-point wavefunctions

Although Gamma is the only independent crystal momentum, information about the periodic position operator can still be represented by a set of overlaps connecting Gamma to its reciprocal-lattice images. For each neighboring reciprocal-space vector $\mathbf b$, define

$$M^{(b)}_{mn}=\bra{\psi_{m\Gamma}}e^{-i\mathbf b\cdot\mathbf r}\ket{\psi_{n\Gamma}}.$$

Each neighboring vector has an associated geometric weight $w_b$. In the full active Bloch space,

$$M^{(b)}\in\mathbb C^{N\times N}.$$

For a selection matrix $S$, the overlap in the selected $J$-dimensional subspace is

$$M_S^{(b)}=S^\dagger M^{(b)}S\in\mathbb C^{J\times J}.$$

#### 1.1.2 Outer window and frozen window

The outer window defines the set of candidate Bloch states allowed to participate in the subspace selection:

$$\mathcal O=\lbrace m\mid E_{\mathrm{win}}^{\min}\leq\epsilon_m\leq E_{\mathrm{win}}^{\max}\rbrace.$$

The target subspace may be selected only from $\mathcal O$; therefore, for every $m\notin\mathcal O$,

$$S_{mn}=0.$$

The frozen window defines the set of Bloch eigenstates that must be retained in full in the target subspace:

$$\mathcal F=\lbrace m\in\mathcal O\mid E_{\mathrm{froz}}^{\min}\leq\epsilon_m\leq E_{\mathrm{froz}}^{\max}\rbrace.$$

Let:

- $N$: number of active Bloch states;
- $B=|\mathcal O|$: number of outer-window Bloch states;
- $J$: number of target Wannier orbitals;
- $M=|\mathcal F|$: number of frozen Bloch states;
- $\mathcal G=\mathcal O\setminus\mathcal F$: subspace formed by the non-frozen outer states;
- $L=J-M$: number of free directions to be selected from $\mathcal G$.

The dimensions satisfy

$$M\leq J\leq B\leq N.$$

Here, frozen denotes a subspace-inclusion constraint. For each $f\in\mathcal F$, let $\ket{e_f}$ be the unit coordinate vector of that Bloch eigenstate in the active-band basis. The constraint is

$$P\ket{e_f}=\ket{e_f}.$$

Let the frozen band indices be $f_1,\ldots,f_M$. Collecting the corresponding unit coordinate vectors as columns gives

$$F=(\ket{e_{f_1}},\ldots,\ket{e_{f_M}})\in\mathbb C^{N\times M}, F^\dagger F=I_M.$$

The projector onto the frozen subspace is

$$P_F=FF^\dagger.$$

The target space additionally requires the selection, from $\mathcal G$, of

$$L=J-M$$

orthogonal directions. Collect these directions in

$$C\in\mathbb C^{N\times L}, C^\dagger C=I_L, F^\dagger C=0.$$

All rows of $C$ outside $\mathcal G$ are required to vanish. The full subspace matrix is

$$S=(F,C).$$

The corresponding projector is

$$P=SS^\dagger=FF^\dagger+CC^\dagger.$$


The initial subspace contains all frozen states and, from $\mathcal G$, uses the $L$ lowest-energy Bloch states as the initial free directions.

#### 1.1.3 Gauge-invariant spread $\Omega_I$

For any candidate subspace $S$, define

$$\Omega_I[S]=\sum_b w_b\left[J-\sum_{m,n=1}^{J}|(M_S^{(b)})_{mn}|^2\right].$$

Using $P=SS^\dagger$, this can also be written as

$$\Omega_I[P]=\sum_b w_b\left\lbrace J-\mathrm{Tr}\left[PM^{(b)}PM^{(b)\dagger}\right]\right\rbrace.$$

$\Omega_I$ measures the leakage of the target subspace into its complement under the periodic-overlap mappings. Disentanglement selects the target subspace by minimizing $\Omega_I$.

#### 1.1.4 Self-consistent SMV iteration

Define

$$Z[P]=\sum_b w_bM^{(b)}PM^{(b)\dagger}.$$

The gauge-invariant spread can then be written as

$$\Omega_I[P]=J\sum_b w_b-\mathrm{Tr}\left[PZ[P]\right].$$

Thus, minimizing $\Omega_I[P]$ is equivalent to maximizing

$$\Phi[P]=\mathrm{Tr}\left[PZ[P]\right].$$

Because $Z[P]$ depends on the projector $P$ being sought, the problem is solved self-consistently. Given the projector $P^{(i-1)}$ from the preceding iteration, construct

$$Z^{(i)}=\sum_b w_bM^{(b)}P^{(i-1)}M^{(b)\dagger}.$$

Within the current iteration, $Z^{(i)}$ is treated as a fixed matrix. Denoting it simply by $Z$ below, the new target subspace is determined by maximizing

$$\mathrm{Tr}\left[PZ\right]$$

The projector onto the target subspace can be written as

$$P=P_F+CC^\dagger.$$

Substitution into the objective gives

$$\mathrm{Tr}\left[PZ\right]=\mathrm{Tr}\left[P_FZ\right]+\mathrm{Tr}\left[C^\dagger ZC\right].$$

The first term is fixed by the frozen subspace, so only the second term needs to be maximized within the non-frozen outer space $\mathcal G=\mathcal O\setminus\mathcal F$. Let the columns of $R_{\mathcal G}$ be the unit coordinate vectors in $\mathcal G$, and write

$$C=R_{\mathcal G}V, V^\dagger V=I_L.$$

Restricting $Z$ to the non-frozen outer space gives

$$Z_{\mathcal G}=R_{\mathcal G}^\dagger ZR_{\mathcal G}.$$

The objective correspondingly becomes

$$\mathrm{Tr}\left[C^\dagger ZC\right]=\mathrm{Tr}\left[V^\dagger Z_{\mathcal G}V\right].$$

By the Ky Fan maximum principle, the maximum of this objective is determined by the eigenvalues of $Z_{\mathcal G}$, specifically the $L$ largest. Let

$$Z_{\mathcal G}v_a=\lambda_av_a, \lambda_1\geq\lambda_2\geq\cdots,$$

then

$$\max_{V^\dagger V=I_L}\mathrm{Tr}\left[V^\dagger Z_{\mathcal G}V\right]=\sum_{a=1}^{L}\lambda_a.$$

Therefore,

$$V=(v_1,\ldots,v_L).$$

The free subspace, full subspace, and projector are then updated as

$$C=R_{\mathcal G}V, S=(F,C), P=SS^\dagger.$$

The updated $P$ is used to construct $Z$ for the next iteration, and the procedure is repeated. After each iteration, the change in the gauge-invariant spread,

$$\Delta\Omega_I^{(i)}=|\Omega_I^{(i)}-\Omega_I^{(i-1)}|,$$

and the change in the projector,

$$\Delta P^{(i)}=\|P^{(i)}-P^{(i-1)}\|_F.$$

are monitored. Once both $\Delta\Omega_I^{(i)}$ and $\Delta P^{(i)}$ have stabilized, the final subspace is obtained and the Gamma-periodic overlaps are projected into the selected $J$-dimensional subspace:

$$M_{\mathrm{S}}^{(b)}=S^\dagger M^{(b)}S.$$

### 1.2 Jacobi rotation: Foster-Boys localization of the Gamma-point initial guesses

Denote the $N_w=J$ orthonormal states obtained from disentanglement by $\ket{\psi_{m\Gamma}}$. The Jacobi rotation constructs localized Gamma-point orbitals within this fixed subspace:

$$\ket{\phi_{\alpha}^{\Gamma}}=\sum_{m=1}^{N_w}U_{m\alpha}^{(\Gamma)}\ket{\psi_{m\Gamma}}.$$

#### 1.2.1 Gamma-point Foster-Boys functional

The overlaps before the Jacobi rotation are

$$M_{mn}^{(0,\mathbf b)}=\bra{\psi_{m\Gamma}}e^{-i\mathbf b\cdot\mathbf r}\ket{\psi_{n\Gamma}}.$$

After rotation by $U^{(\Gamma)}$,

$$M^{(\mathbf b)}=U^{(\Gamma)\dagger}M^{(0,\mathbf b)}U^{(\Gamma)}.$$

The Gamma-point periodic Foster-Boys spread is

$$\Omega_{\Gamma}=\sum_{\mathbf b}w_{\mathbf b}\sum_{n=1}^{N_w}\left(1-|M_{nn}^{(\mathbf b)}|^2\right).$$

Thus, minimizing $\Omega_{\Gamma}$ is equivalent to maximizing

$$F=\sum_{\mathbf b}w_{\mathbf b}\sum_{n=1}^{N_w}|M_{nn}^{(\mathbf b)}|^2.$$

#### 1.2.2 Two-orbital $SU(2)$ rotation

The Jacobi method selects one pair of orbitals at a time and applies, within the corresponding two-dimensional subspace,

$$M'^{(\mathbf b)}=V^\dagger M^{(\mathbf b)}V, V\in SU(2).$$

Write the $2\times2$ overlap block corresponding to this orbital pair as

$$M^{(\mathbf b)}=a_0^{(\mathbf b)}I+\mathbf a^{(\mathbf b)}\cdot\boldsymbol{\sigma},$$

where

$$a_0^{(\mathbf b)}=\frac{1}{2}\left(M_{11}^{(\mathbf b)}+M_{22}^{(\mathbf b)}\right),$$

$$a_1^{(\mathbf b)}=\frac{1}{2}\left(M_{12}^{(\mathbf b)}+M_{21}^{(\mathbf b)}\right),$$

$$a_2^{(\mathbf b)}=\frac{i}{2}\left(M_{12}^{(\mathbf b)}-M_{21}^{(\mathbf b)}\right),$$

$$a_3^{(\mathbf b)}=\frac{1}{2}\left(M_{11}^{(\mathbf b)}-M_{22}^{(\mathbf b)}\right).$$

An $SU(2)$ rotation can be written as

$$V=e^{i\frac{\theta}{2}\mathbf n\cdot\boldsymbol{\sigma}}=\cos\frac{\theta}{2}I+i\sin\frac{\theta}{2}\mathbf n\cdot\boldsymbol{\sigma}.$$

It corresponds to an $SO(3)$ rotation $R$ in Pauli-vector space:

$$V^\dagger\left(\mathbf a\cdot\boldsymbol{\sigma}\right)V=(R\mathbf a)\cdot\boldsymbol{\sigma}.$$

The corresponding three-dimensional rotation is

$$R_{\mathbf n}(\theta)\mathbf a=\mathbf a\cos\theta+(\mathbf n\times\mathbf a)\sin\theta+\mathbf n(\mathbf n\cdot\mathbf a)(1-\cos\theta).$$

#### 1.2.3 Determination of the optimal rotation direction

The two diagonal elements after the rotation satisfy

$$|M_{11}'^{(\mathbf b)}|^2+|M_{22}'^{(\mathbf b)}|^2=2|a_0^{(\mathbf b)}|^2+2|a_3'^{(\mathbf b)}|^2.$$

Because $a_0^{(\mathbf b)}$ is invariant under the similarity transformation, the optimal rotation for the orbital pair is equivalent to maximizing

$$\sum_{\mathbf b}w_{\mathbf b}|a_3'^{(\mathbf b)}|^2.$$

Define

$$\mathbf u=R^T\mathbf e_3,$$

then

$$a_3'^{(\mathbf b)}=\mathbf u^T\mathbf a^{(\mathbf b)}.$$

Therefore,

$$\sum_{\mathbf b}w_{\mathbf b}|a_3'^{(\mathbf b)}|^2=\mathbf u^TW\mathbf u,$$

where $W$ is the $3\times3$ real symmetric matrix

$$W_{\alpha\beta}=\sum_{\mathbf b}w_{\mathbf b}\mathrm{Re}\left[a_{\alpha}^{(\mathbf b)}a_{\beta}^{(\mathbf b)\ast}\right].$$

By the Rayleigh-Ritz variational principle,

$$\max_{|\mathbf u|=1}\mathbf u^TW\mathbf u=\lambda_{\max},$$

and the optimal direction $\mathbf u_{\max}$ is the normalized eigenvector associated with the largest eigenvalue of $W$:

$$W\mathbf u_{\max}=\lambda_{\max}\mathbf u_{\max}.$$

#### 1.2.4 Accumulated unitary transformation

For $\mathbf u_{\max}=(u_x,u_y,u_z)$, define

$$c=\sqrt{\frac{1+u_z}{2}},z=\frac{u_x+iu_y}{2c},$$

then $V$ is

$$V=\begin{pmatrix}c & -z^\ast \\ 
z & c\end{pmatrix}.$$

For every $\mathbf b$, update the overlap for this orbital pair as

$$M^{(\mathbf b)}\rightarrow V^\dagger M^{(\mathbf b)}V.$$

At the same time, accumulate the Gamma-point Jacobi transformation:

$$U^{(\Gamma)}\rightarrow U^{(\Gamma)}V.$$

Rotate all orbital pairs in sequence and repeat the procedure until $\Omega_\Gamma$ has stabilized. Finally, combine the result with the disentanglement transformation $S$:

$$U_{\mathrm{total}}^{(\Gamma)}=SU^{(\Gamma)}.$$

### 1.3 Locate center: orbital centers of the Gamma-point initial guesses

For each localized Gamma-point orbital, its center is defined as the position at which a minimum-radius Fermi envelope contains a prescribed fraction of the charge.

#### 1.3.1 Real-space charge density

The charge density of orbital $\alpha$ is

$$\rho_\alpha(\mathbf r)=\sum_s|\phi_{\alpha s}^{\Gamma}(\mathbf r)|^2,$$

where a scalar wavefunction has one component, while a spinor wavefunction is summed over both spin components.

#### 1.3.2 Fermi envelope and minimum envelope radius

For each candidate envelope center $\mathbf r_c$, construct a Wigner-Seitz cell centered at $\mathbf r_c$. In what follows, $\mathbf r$ denotes the coordinate relative to the envelope center within this WS cell, so $|\mathbf r|$ is the distance from the envelope center.

The Fermi envelope is defined as

$$f_R(\mathbf r)=\left[1+\exp\left(\frac{|\mathbf r|-R}{\sigma}\right)\right]^{-1}.$$

The fraction of charge inside the envelope is

$$q_\alpha(\mathbf r_c,R)=\frac{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_R(\mathbf r)d\mathbf r}{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)d\mathbf r}.$$

For a prescribed target charge fraction $\eta$, the minimum envelope radius satisfying the condition is determined for each candidate center:

$$R_\alpha(\mathbf r_c)=\min\lbrace R\mid q_\alpha(\mathbf r_c,R)\geq\eta\rbrace.$$

The center search starts independently from multiple candidate positions, including the position of maximum orbital charge density and positions distributed uniformly throughout the periodic cell. The candidates are first ranked by the minimum envelope radius required to contain the target charge fraction. Full self-consistent center iterations are then performed, starting from candidates with smaller envelope radii.

#### 1.3.3 Self-consistent iteration of the envelope center

The displacement of the center of charge within the envelope relative to the envelope center is

$$\mathbf s_\alpha(\mathbf r_c)=\frac{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)\mathbf r d\mathbf r}{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)d\mathbf r}.$$

Move the envelope center toward the center of charge within the envelope:

$$\mathbf r_{\mathrm{trial}}=\mathbf r_c+\gamma\mathbf s_\alpha(\mathbf r_c).$$

At the trial center, reconstruct the WS cell and recompute the minimum envelope radius and envelope-center displacement. Accept the trial center if the envelope-center displacement decreases or the envelope radius becomes smaller. Repeat the procedure until

$$\left|\mathbf s_\alpha(\mathbf r_c)\right|<\epsilon_c.$$

For each converged center, compute the second moment of the charge within the Fermi envelope:

$$\mu_{2,\alpha}(\mathbf r_c)=\frac{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)|\mathbf r|^2d\mathbf r}{\int_{\mathrm{WS}}\rho_\alpha(\mathbf r)f_{R_\alpha}(\mathbf r)d\mathbf r}.$$

Finally, select the self-consistent center with the smallest envelope radius. When the envelope radii are close, select the center with the smaller second moment.

### 1.4 Phase fix: global phase convention for the Gamma-point initial guesses

#### 1.4.1 Phase convention for scalar Wannier orbitals

For a scalar orbital, take the position of maximum wavefunction amplitude,

$$\mathbf r_\alpha^{\max}=\mathrm\ast{argmax}_{\mathbf r}|\phi_\alpha^\Gamma(\mathbf r)|,$$

and apply the phase correction

$$p_\alpha=\frac{\phi_\alpha^\Gamma(\mathbf r_\alpha^{\max})^\ast}{|\phi_\alpha^\Gamma(\mathbf r_\alpha^{\max})|}.$$

The Gamma-point transformation matrix is correspondingly updated as

$$\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}=p_\alpha U_{\mathrm{total},m\alpha}^{(\Gamma)}.$$

#### 1.4.2 Phase convention for spinor Wannier orbitals

For a spinor orbital, write

$$\Phi_{\alpha,z}^\Gamma(\mathbf r)=\begin{pmatrix}a_{\alpha,z}(\mathbf r)\\
b_{\alpha,z}(\mathbf r)\end{pmatrix}.$$

All orbitals use a common spin basis $V$:

$$\Phi_{\alpha,V}^\Gamma(\mathbf r)=V^\dagger\Phi_{\alpha,z}^\Gamma(\mathbf r)=\begin{pmatrix}a_{\alpha,V}(\mathbf r)\\
b_{\alpha,V}(\mathbf r)\end{pmatrix}.$$

In this basis, choose the symmetric reference spinor

$$\chi_V(\theta,\phi)=e^{-i\pi/4}\begin{pmatrix}e^{-i\phi/2}\cos(\theta/2)\\
e^{i\phi/2}\sin(\theta/2)\end{pmatrix}.$$

The spinor function is represented as

$$\Phi_{\alpha,V}^\Gamma(\mathbf r)=\sqrt{\rho_\alpha(\mathbf r)}e^{i[\gamma_{\alpha,V}(\mathbf r)+\pi/4]}\chi_V(\theta_{\alpha,V}(\mathbf r),\phi_{\alpha,V}(\mathbf r)).$$

where

$$\rho_\alpha(\mathbf r)=|a_{\alpha,V}(\mathbf r)|^2+|b_{\alpha,V}(\mathbf r)|^2.$$

When both components are nonzero, the local phase is

$$\gamma_{\alpha,V}(\mathbf r)=\frac{\arg a_{\alpha,V}(\mathbf r)+\arg b_{\alpha,V}(\mathbf r)}{2}\pmod\pi.$$

When $a_{\alpha,V}$ or $b_{\alpha,V}$ is close to zero, the corresponding $\arg a_{\alpha,V}$ or $\arg b_{\alpha,V}$ becomes ill-defined. A common spin quantization axis must therefore be chosen such that the principal charge distribution remains as far as possible from the two poles of the new spin basis.

Define the local spin density

$$\mathbf s_\alpha(\mathbf r)=\Phi_{\alpha,z}^{\Gamma\dagger}(\mathbf r)\boldsymbol\sigma\Phi_{\alpha,z}^\Gamma(\mathbf r).$$

If the $z$-axis of the new basis points along the unit vector $\mathbf n$, then

$$\cos\theta_{\alpha,V}(\mathbf r)=\frac{\mathbf n\cdot\mathbf s_\alpha(\mathbf r)}{\rho_\alpha(\mathbf r)}.$$

The two poles correspond to $|\cos\theta_{\alpha,V}|=1$. The direction $\mathbf n$ is therefore selected by minimizing the mean polar weight over all orbitals:

$$\mathcal P(\mathbf n)=\frac{1}{N_w}\sum_{\alpha=1}^{N_w}\frac{\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)\cos^2\theta_{\alpha,V}(\mathbf r)d\mathbf r}{\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)d\mathbf r}.$$

This functional can be written as

$$\mathcal P(\mathbf n)=\mathbf n^TK\mathbf n,$$

where

$$K_{ij}=\frac{1}{N_w}\sum_{\alpha=1}^{N_w}\frac{\int_{\Omega_\alpha}s_{\alpha i}(\mathbf r)s_{\alpha j}(\mathbf r)/\rho_\alpha(\mathbf r)d\mathbf r}{\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)d\mathbf r}.$$

By the Rayleigh-Ritz variational principle, the optimal quantization axis is the normalized eigenvector associated with the smallest eigenvalue of $K$:

$$K\mathbf n=\kappa_{\min}\mathbf n.$$

The common spin basis is constructed from $\mathbf n=(n_x,n_y,n_z)$ as

$$c=\sqrt{\frac{1+n_z}{2}},$$

$$z=\frac{n_x+in_y}{2c},$$

$$V=\begin{pmatrix}c&-z^\ast\\
z&c\end{pmatrix}.$$

Once the common axis has been determined, use the doubled phase

$$D_\alpha(\mathbf r)=e^{2i[\gamma_{\alpha,V}(\mathbf r)+\pi/4]}.$$

In the north-pole region, use

$$D_\alpha^N(\mathbf r)=i\left[\frac{a_{\alpha,V}(\mathbf r)}{|a_{\alpha,V}(\mathbf r)|}\right]^2,$$

and in the south-pole region, use its time-reversal counterpart

$$D_\alpha^S(\mathbf r)=-i\left[\frac{b_{\alpha,V}(\mathbf r)}{|b_{\alpha,V}(\mathbf r)|}\right]^2.$$

Next, choose a common phase $e^{-i\vartheta_\alpha}$ for the entire orbital such that the integral over the WS cell of the squared real part relative to the reference spinor is maximized:

$$F_\alpha(\vartheta_\alpha)=\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)\left[\mathrm{Re}\left(e^{-i\vartheta_\alpha}e^{i[\gamma_{\alpha,V}(\mathbf r)+\pi/4]}\right)\right]^2d\mathbf r.$$

Define the doubled-phase integral

$$C_\alpha=\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)D_\alpha(\mathbf r)d\mathbf r.$$

Then

$$F_\alpha(\vartheta_\alpha)=\frac{1}{2}\int_{\Omega_\alpha}\rho_\alpha(\mathbf r)d\mathbf r+\frac{1}{2}\mathrm{Re}\left[e^{-2i\vartheta_\alpha}C_\alpha\right].$$

Therefore, $F_\alpha$ reaches its maximum at

$$\vartheta_\alpha=\frac{1}{2}\arg C_\alpha\pmod\pi$$

The final phase correction is

$$p_\alpha=e^{-i\vartheta_\alpha}=e^{-i\arg(C_\alpha)/2},$$

and the transformation is updated as

$$\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}=p_\alpha U_{\mathrm{total},m\alpha}^{(\Gamma)}.$$

The factor $e^{-i\pi/4}$ in the symmetric reference spinor ensures that

$$\Theta\chi_V(\theta,\phi)=\chi_V(\pi-\theta,\phi+\pi).$$

Therefore, for a Kramers pair, the phase-corrected orbitals automatically satisfy the following relation even without explicit pairing:

$$\widetilde\Phi_{\bar\alpha}^\Gamma=\pm\Theta\widetilde\Phi_\alpha^\Gamma.$$

## 2. Gamma-adaptive .amn matrix

### 2.1 Overlaps between Gamma-point and full-$\mathbf k$ wavefunctions

For Wannier orbital $\alpha$, use its center $\boldsymbol\tau_\alpha$ to construct the Wigner--Seitz cell $\Omega_{\mathrm{WS}}(\boldsymbol\tau_\alpha)$. Define the overlap between a full-$\mathbf k$ Bloch state and a Gamma-point Bloch state within this WS cell as

$$I_{nm}^{(\mathbf k;\alpha)}=\int_{\Omega_{\mathrm{WS}}(\boldsymbol\tau_\alpha)}d^3r\,\psi_{n\mathbf k}^\ast(\mathbf r)\psi_{m\Gamma}(\mathbf r).$$

For spinor wavefunctions, the inner product includes the sum over both spinor components:

$$I_{nm}^{(\mathbf k;\alpha)}=\sum_s\int_{\Omega_{\mathrm{WS}}(\boldsymbol\tau_\alpha)}d^3r\,\psi_{n\mathbf k s}^\ast(\mathbf r)\psi_{m\Gamma s}(\mathbf r).$$

Because the Bloch phases of the $\mathbf k$-point and Gamma-point wavefunctions do not cancel within the cell, the integrand
$$\psi_{n\mathbf k}^\ast(\mathbf r)\psi_{m\Gamma}(\mathbf r)=e^{-i\mathbf k\cdot\mathbf r}u_{n\mathbf k}^\ast(\mathbf r)u_{m\Gamma}(\mathbf r).$$

is not lattice-periodic. This overlap cannot be obtained directly from an inner product of periodic plane-wave coefficients and must instead be integrated explicitly within the selected real-space WS cell.

### 2.2 Construction of the .amn matrix

After disentanglement, Jacobi rotation, and phase fixing, the Gamma-adaptive Wannier initial guess is

$$\ket{g_{\alpha\mathbf 0}}=\sum_{m=1}^{N}\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}\ket{\psi_{m\Gamma}}.$$

The .amn matrix is defined as the overlap between a full-$\mathbf k$ Bloch state and this Gamma-adaptive initial guess:

$$A_{n\alpha}^{(\mathbf k)}=\braket{\psi_{n\mathbf k}|g_{\alpha\mathbf 0}}.$$

Substituting the Gamma-point initial guess gives

$$A_{n\alpha}^{(\mathbf k)}=\sum_{m=1}^{N}I_{nm}^{(\mathbf k;\alpha)}\widetilde U_{\mathrm{total},m\alpha}^{(\Gamma)}.$$
