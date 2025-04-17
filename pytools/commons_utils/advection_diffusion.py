#!/usr/bin/env python3

import numpy as np
import scipy.linalg.lapack as scilapack

#-----

def elm_soillayers(nlevgrnd=15, equal_thickness=False):
    #nlevgrnd = 15

    # For testing with equally distant soil layers
    if equal_thickness:
        dzsoi  = np.ones((nlevgrnd), dtype=float)*0.050
        zisoi  = np.insert(np.cumsum(dzsoi),0,0.0)
        zsoi   = np.cumsum(dzsoi) - 0.025

    else:
        
        # -- standard soil layers from ELM's 10- or 15-layer column--
        #  variable layer thickness
        """
        initVerticalMod.F90 --
        do j = 1, nlevgrnd
           zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
        enddo
        dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
        do j = 2,nlevgrnd-1
           dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
        enddo
        dzsoi(nlevgrnd) = zsoi(nlevgrnd)-zsoi(nlevgrnd-1)
        zisoi(0) = 0._r8
        do j = 1, nlevgrnd-1
           zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
        enddo
        zisoi(nlevgrnd) = zsoi(nlevgrnd) + 0.5_r8*dzsoi(nlevgrnd)
        """
        
        jidx = np.array(range(nlevgrnd))+1 
        zsoi = 0.025*(np.exp(0.5*(jidx-0.5))-1.0)       #ELM soil layer node depths - somewhere inside a layer but not centroid
        
        dzsoi= np.zeros_like(zsoi)
        dzsoi[0] = 0.5*(zsoi[0]+zsoi[1])                #thickness b/n two vertical interfaces (vertices)
        for j in range(1,nlevgrnd-1):
            dzsoi[j]= 0.5*(zsoi[j+1]-zsoi[j-1])
        dzsoi[nlevgrnd-1] = zsoi[nlevgrnd-1]-zsoi[nlevgrnd-2]
    
        zisoi= np.zeros((nlevgrnd+1), dtype=float)
        zisoi[0] = 0.0                                  # layer interface depth
        for j in range(1,nlevgrnd):
            zisoi[j] = 0.5*(zsoi[j-1]+zsoi[j])
        zisoi[nlevgrnd] = zsoi[nlevgrnd-1]+0.5*dzsoi[nlevgrnd-1]

    
    # fake zsoi, dzsoi from 1:nlevgrnd+1, so matching with Fortran style
    zsoi = np.insert(zsoi,0,np.nan)
    dzsoi= np.insert(dzsoi,0,np.nan)
    
    return zsoi, dzsoi, zisoi
 


#----- 
def a_step_explicit(nlevsoi, z, dz, c_prev, \
                    rain_conc, q_int, theta, watsat, sourcesink, d0, dt):
    """
    Author: Yaoping Wang, ESD, ORNL, 2025-04-15
    
    rain_conc - upper boundary condition (rain water chemistry) (mol/m3-water)
    c_prev - start concentration (mol/m3-water)
    q_int - water flux at grid boundaries (m/s)
    theta - soil moisture value (m3/m3)
    watsat - porosity
    sourcesink - source/sink strength (mol/m3-soil/s)
    d0 - diffusion coefficient in water (m2/s)
    dt - time step size (s)

    If we ignore spatial dependency, the equation becomes very simple

    Consider a generic cell C_i with adjacent cells
        C_{i-1} <--- Δ x_i ---> C_i, D_{eff,i} <--- Δ x_{i+1} ---> C_{i+1}

    Total mass change within C_i due to outflow
        Δ x  * θ * \frac{ dC_i }{ dt }

    Diffusion to the upper cell, if C_i > C_{i-1}
        dup = D_{eff,i} * I[C_i > C_{i-1}] * \frac{C_i - C_{i-1}}{Δ x_i}

    Diffusion to the lower cell, if C_i > C_{i+1}
        dlow = D_{eff,i} * I[C_i > C_{i+1}] * \frac{C_i - C_{i+1}}{Δ x_{i+1}}

    Advection to upper cell, if q_{in,i} < 0
        aup = I[q_{in,i} < 0] * abs( q_{in,i} ) * C_i

    Advection to lower cell, if q_{out,i} > 0
        alow = I[q_{out,i} > 0] * q_{out,i} * C_i

    Rain inflow
        rain = I[q_{in,i} > 0] * q_{in,0} * C_{rain}
    
    The relationship between self-change and the four outflow and internal source/sink, 
        after simplification, is

        Δ x * θ * \frac{ dC_i }{ dt } = 
          - ( \frac{ D_{eff,i} * I[C_i > C_{i-1}] }{Δ x_i} + 
              \frac{ D_{eff,i} * I[C_i > C_{i+1}] }{Δ x_{i+1}} + 
              I[q_{in,i} < 0] * abs( q_{in,i} ) + I[q_{out,i} > 0] * q_{out,i} ) * C_i
          + ( D_{eff,i} * I[C_i > C_{i-1}] \frac{C_{i-1}}{Δ x_i} )
          + ( D_{eff,i} * I[C_i > C_{i+1}] \frac{C_{i+1}}{Δ x_{i+1}} )
          + I[q_{in,i} > 0] * q_{in,0} * C_{rain}
          + R

        , which is of the form dC/dt = kC + r, and the analytical solution is just
             C = (C0 + r/k)e^{kt} - r/k
             on the interval [t, t + Δ t]

    The equations means when t → ∞, C → -r/k

    It is guaranteed that k < 0. Therefore, 

    if r > 0, the concentration will not be negative for large t (we can be sure )

    If r < 0, the concentration may be negative because the equilibrium value is
        negative. This can only happen when r < 0. In this case, the analytical 
        form may need to be assessed piecewise to insure non-negativity:

        Step 1. We evaluate the analytical form at (t + Δ t),
                if the resulting C_t' > C_{i-1} and C_t' > C_{i+1}, then C_t = C_t', 
                otherwise, we need to go to step 2
        Step 2. We assess piecewise: 
                first, find the dt' such that C(t + dt') = max(C_{i-1}, C_{i+1})
                then, on the interval [t+dt', t + Δ t], set the appropriate I[⋅]
                       term to 0, re-evaluate to t + Δ t
                third, if the result gives C(t + Δ t) > min(C_{i-1}, C_{i}), stop
                       otherwise, need a third piecewise, go to step 3
        Step 3. find the dt" such that C(t + dt") = min(C_{i-1}, C_{i+1}), on the
                interval [t+dt", t + Δ t], all the diffusion terms need to be 0.
                Re-integrate the pure-advection term to t + Δ t
    """
    def _analytical_c(C0, r, k, t1, t2):
        # analytical solution of dC/dt = kC + r, from t1 to t2
        c = (C0 + r/k)*np.exp(k*(t2-t1)) - r/k
        return c
    
    def _analytical_c_int(C0, r, k, t1, t2):
        # integral of the above analytical solution of C
        # over [t1, t2]
        c_int = (C0 + r/k)/k*(np.exp(k*(t2-t1))-1) - r/k*(t2-t1)
        return c_int

    def _analytical_dt(c1, c2, r, k):
        # analytical solution of the dt it takes for the above
        # analytical solution of C to evolve from c1 to c2
        dt = np.log((c2 + r/k) / (c1 + r/k)) / k
        return dt

    #
    zsoi = z[:nlevsoi]
    dz   = dz[:nlevsoi]
        
    #
    Deff = d0 * theta**(10/3) / watsat**2
    dx = np.append(np.insert(np.diff(zsoi), 0, 100), 21) # padding; doesn't matter, so long as >0

    scale = dz * theta # \Delta x * θ on the LHS

    c_next = np.zeros(nlevsoi)
    niter = np.zeros(nlevsoi) # keep track how many iteration was done

    dc_up = np.zeros(nlevsoi) # flow upward due to diffusion and advection
    dc_down = np.zeros(nlevsoi) # flow down due to diffusion and advection

    # convert sourcesink from mol/m3-soil/s to mol/m2/s
    sourcesink = sourcesink*dz

    if q_int[0] > 0:
        sourcesink[0] = sourcesink[0] + q_int[0]*rain_conc

    for i in range(nlevsoi):

        if i == 0:
            # top layer: no diffusion to above
            i1 = 0
            i2 = int(c_prev[i] > c_prev[i+1])

        elif i < (nlevsoi-1):
            # middle layers
            i1 = int(c_prev[i] > c_prev[i-1])
            i2 = int(c_prev[i] > c_prev[i+1])

        else:
            # bottom layer: no diffusion to below
            i1 = int(c_prev[i] > c_prev[i-1])
            i2 = 0

        i3 = int(q_int[i] < 0)
        i4 = int(q_int[i+1] > 0)

        k = - (i1*Deff[i]/dx[i] + i2*Deff[i]/dx[i+1] + i3*abs(q_int[i]) + i4*q_int[i+1]) / scale[i]

        if i == 0:
            r = (Deff[i]*i2*c_prev[i+1]/dx[i+1] + sourcesink[i]) / scale[i]
        elif i < (nlevsoi-1):
            r = (Deff[i]*i1*c_prev[i-1]/dx[i] + Deff[i]*i2*c_prev[i+1]/dx[i+1] \
                 + sourcesink[i]) / scale[i]
        else:
            r = (Deff[i]*i1*c_prev[i-1]/dx[i] + sourcesink[i]) / scale[i]
        
        # if there is no flow out of this cell, degrades to linear source
        if k == 0:
            c_next[i] = c_prev[i] + r*dt
            dc_up[i] = 0
            dc_down[i] = 0

            continue

        # otherwise, do actual calculations

        # step 1: find end-of-step solution
        c = _analytical_c(c_prev[i], r, k, 0, dt)

        if i == 0:
            cmax = c_prev[i+1]

            # end-of-step solution works, or there is no diffusion to begin with
            if (c > cmax) | (i2 == 0):
                c_next[i] = c
                niter[i] = 1

                # integration of c over [0, dt]
                c_int = _analytical_c_int(c_prev[i], r, k, 0, dt)

                # use c_int to calculate the fluxes
                dc_up[i] = i3*abs(q_int[i])*c_int
                dc_down[i] = i2*Deff[i]/dx[i+1] * (c_int - c_prev[i+1]*dt) \
                                + i4*q_int[i+1]*c_int

            else:
                # step 2: piecewise integrate to cmax, then update k, r and continue
                dt_p = _analytical_dt(c_prev[i], cmax, r, k)

                # integration of c over [0, dt']
                c_int_p = _analytical_c_int(c_prev[i], r, k, 0, dt_p)

                # use c_int to calculate the fluxes during [0, dt']
                dc_up[i] = i3*abs(q_int[i])*c_int_p
                dc_down[i] = i2*Deff[i]/dx[i+1] * (c_int_p - c_prev[i+1]*dt_p) \
                                + i4*q_int[i+1]*c_int_p

                # update k & r (no more diffusion)
                k = - (i3*abs(q_int[i]) + i4*q_int[i+1]) / scale[i]
                r = sourcesink[i] / scale[i]

                # itegration of c over [dt', dt]
                c_int_pp = _analytical_c_int(cmax, r, k, dt_p, dt)

                # add to the fluxes
                dc_up[i] = dc_up[i] + i3*abs(q_int[i])*c_int_pp
                dc_down[i] = dc_down[i] + i4*q_int[i+1]*c_int_pp

                # find final value
                c_next[i] = _analytical_c(cmax, r, k, dt_p, dt)
                niter[i] = 2

        elif i < (nlevsoi-1):

            if (i1 == 0) & (i2 == 0):
                # if there is no diffusion to begin with, end of solution works

                # final value found
                c_next[i] = c
                niter[i] = 1

                # integration of c over [0, dt]
                c_int = _analytical_c_int(c_prev[i], r, k, 0, dt)

                # calculate the fluxes during [0, dt] (all diffusion)
                dc_up[i] = i3*abs(q_int[i])*c_int
                dc_down[i] = i4*q_int[i+1]*c_int

            else:
                # if there is only one side diffusion
                if (i1 == 0):
                    cmax = c_prev[i+1] # diffuse down
                    cmin = 0
                    i1_keep = 0
                    i2_keep = 1

                elif (i2 == 0):
                    cmax = c_prev[i-1] # diffuse up
                    cmin = 0
                    i1_keep = 1
                    i2_keep = 0

                else:
                    # both sides have diffusion, decide which side to begin with

                    if c_prev[i-1] > c_prev[i+1]:
                        cmax = c_prev[i-1]
                        cmin = c_prev[i+1]
                        i1_keep = 0
                        i2_keep = 1

                    else:
                        cmax = c_prev[i+1]
                        cmin = c_prev[i-1]
                        i1_keep = 1
                        i2_keep = 0

                # if end of solution works
                if c > cmax:

                    # final value found
                    c_next[i] = c
                    niter[i] = 1

                    # integration of c over [0, dt]
                    c_int = _analytical_c_int(c_prev[i], r, k, 0, dt)

                    # calculate the fluxes during [0, dt] (all diffusion)
                    dc_up[i] = i1*Deff[i]/dx[i] * (c_int - c_prev[i-1]*dt) \
                            + i3*abs(q_int[i])*c_int
                    dc_down[i] = i2*Deff[i]/dx[i+1] * (c_int - c_prev[i+1]*dt) \
                            + i4*q_int[i+1]*c_int

                else:
                    # step 2: piecewise integrate to cmax, then update k, r and continue
                    dt_p = _analytical_dt(c_prev[i], cmax, r, k)

                    # integration of c over [0, dt']
                    c_int_p = _analytical_c_int(c_prev[i], r, k, 0, dt_p)

                    # calculate the fluxes during [0, dt'] (all diffusion)
                    dc_up[i] = i1*Deff[i]/dx[i] * (c_int_p - c_prev[i-1]*dt_p) \
                            + i3*abs(q_int[i])*c_int_p
                    dc_down[i] = i2*Deff[i]/dx[i+1] * (c_int_p - c_prev[i+1]*dt_p) \
                            + i4*q_int[i+1]*c_int_p

                    # update k & r (drop one side of diffusion)
                    k = - (i1_keep*i1*Deff[i]/dx[i] + i2_keep*i2*Deff[i]/dx[i+1] + \
                            i3*abs(q_int[i]) + i4*q_int[i+1]) / scale[i]
                    r = (Deff[i]*i1_keep*i1*c_prev[i-1]/dx[i] + \
                         Deff[i]*i2_keep*i2*c_prev[i+1]/dx[i+1] + sourcesink[i]) / scale[i]

                    # continue to end, but we need another check
                    c_p = _analytical_c(cmax, r, k, dt_p, dt)

                    if c_p > cmin:

                        # itegration of c over [dt', dt]
                        c_int_p = _analytical_c_int(cmax, r, k, dt_p, dt)

                        # add the fluxes during [dt', dt] (drop one side of diffusion)
                        dc_up[i] = dc_up[i] \
                            + i1_keep*i1*Deff[i]/dx[i] * (c_int_p - cmax*(dt-dt_p)) \
                            + i3*abs(q_int[i])*c_int_p
                        dc_down[i] = dc_down[i] \
                            + i2_keep*i2*Deff[i]/dx[i+1] * (c_int_p - cmax*(dt-dt_p)) \
                            + i4*q_int[i+1]*c_int_p

                        # final value at dt
                        c_next[i] = _analytical_c(cmax, r, k, dt_p, dt)
                        niter[i] = 2

                    else:
                        # step 3: do another piecewise integration

                        # second stop point
                        dt_pp = _analytical_dt(cmax, cmin, r, k)

                        # itegration of c over [dt', dt"]
                        c_int_pp = _analytical_c_int(c_p, r, k, dt_p, dt_pp)

                        # add the fluxes during [dt', dt"] (drop one side of diffusion)
                        dc_up[i] = dc_up[i] \
                            + i1_keep*i1*Deff[i]/dx[i] * (c_int_pp - cmin*(dt_pp-dt_p)) \
                            + i3*abs(q_int[i])*c_int_pp
                        dc_down[i] = dc_down[i] \
                            + i2_keep*i2*Deff[i]/dx[i+1] * (c_int_pp - cmin*(dt_pp-dt_p)) \
                            + i4*q_int[i+1]*c_int_pp

                        # update the k & r (drop all diffusion)
                        k = - (i3*abs(q_int[i]) + i4*q_int[i+1]) / scale[i]
                        r = sourcesink[i] / scale[i]

                        # itegration of c over [dt", dt]
                        c_int_pp = _analytical_c_int(cmin, r, k, dt_pp, dt)

                        # add the fluxes during [dt", dt]
                        dc_up[i] = dc_up[i] + i3*abs(q_int[i])*c_int_pp
                        dc_down[i] = dc_down[i] + i4*q_int[i+1]*c_int_pp

                        # final value
                        c_next[i] = _analytical_c(cmin, r, k, dt_pp, dt)
                        niter[i] = 3

        # last soil layer
        else:
            cmax = c_prev[i-1]

            # end-of-step solution works, or there is no diffusion to begin with 
            if (c > cmax) | (i1 == 0):
                c_next[i] = c
                niter[i] = 1
            
                # integration of c over [0, dt]
                c_int = _analytical_c_int(c_prev[i], r, k, 0, dt)

                # use c_int to calculate the fluxes
                dc_up[i] = i1*Deff[i]/dx[i] * (c_int - c_prev[i-1]*dt) \
                           + i3*abs(q_int[i])*c_int
                dc_down[i] = i4*q_int[i+1]*c_int

            else:
                # step 2: piecewise integrate to cmax, then update k, r and continue
                dt_p = _analytical_dt(c_prev[i-1], cmax, r, k)

                # integration of c over [0, dt']
                c_int_p = _analytical_c_int(c_prev[i], r, k, 0, dt_p)

                # calculate the fluxes during [0, dt'] (all diffusion)
                dc_up[i] = i1*Deff[i]/dx[i] * (c_int_p - c_prev[i-1]*dt_p) \
                           + i3*abs(q_int[i])*c_int_p
                dc_down[i] = i4*q_int[i+1]*c_int_p

                # update k & r (no more diffusion)
                k = - (i3*abs(q_int[i]) + i4*q_int[i+1]) / scale[i]
                r = sourcesink[i] / scale[i]

                # itegration of c over [dt', dt]
                c_int_pp = _analytical_c_int(cmax, r, k, dt_p, dt)

                # add to the fluxes
                dc_up[i] = dc_up[i] + i3*abs(q_int[i])*c_int_pp
                dc_down[i] = dc_down[i] + i4*q_int[i+1]*c_int_pp

                # find final value
                c_next[i] = _analytical_c(cmax, r, k, dt_p, dt)
                niter[i] = 2

    # we still need to catch the case when r is simply too negative
    # in that case, we really need to reduce r (secondary mineral
    # precipitation, etc.)
    # (TBD)

    # calculate the net between self-outflow and inflow fluxes
    # note the inflow fluxes need to be scaled by soil moisture
    # to get the correct concentration implications
    c_next[:-1] = c_next[:-1] + dc_up[1:]/theta[:-1]/dz[:-1]
    c_next[1:] = c_next[1:] + dc_down[:-1]/theta[1:]/dz[1:]

    for i in range(nlevsoi):
        if (c_next[i] < 0):
            print(i, c_next[i], dc_up[i+1]/theta[i]/dz[i], dc_down[i-1]/theta[i]/dz[i])

    return c_next, niter



#
def Patankar_solute_transport(nlevbed, dzsoi, zisoi, \
                        conc=np.empty((0)), conc_rate=np.empty((0)), \
                        adv=np.empty((0)), diffus=np.empty((0)), vwc=np.empty((0)), \
                        qsrc=np.empty((0)), conc_surf=0, dtime=1800.0):
    
    # soil bulk concentration, mol/m3-soil or g/m3-soil, is used, so vwc is needed to convert
    #  otherwise, i.e. in unit of per m3 water, vwc should be set to 1.0
    # 
    # Advection and diffusion for a single tracer in one column given diffusion coefficient, flow, and source-sink terms
    # Based on SoilLittVertTranspMod, which implements S. V. Patankar, Numerical Heat Transfer and Fluid Flow, Series in Computational Methods in Mechanics and Thermal Sciences, Hemisphere Publishing Corp., 1980. Chapter 5
    # 

    # naive solver of tridiagonal matrx
    def tridiag_naive(du, d, dl, r):
        '''
        Solves many tridiagonal matrix systems with diagonals du, d, dl and RHS vectors r
        '''
        assert du.shape == d.shape and du.shape == dl.shape and du.shape == r.shape
        n = du.shape[-1]
    
        d2 = d
        r2 = r
        for i in range(1, n):
            w = du[..., i] / d[..., i - 1]
            d2[..., i] += -w * dl[..., i - 1]
            r2[..., i] += -w * r[..., i - 1]
        out = np.empty_like(du)
        out[..., -1] = d2[..., -1] / d2[..., -1]
        for i in range(n - 2, -1, -1):
            out[..., i] = (r[..., i] - dl[..., i] * out[..., i + 1]) / d2[..., i]
        return out


    # Statement function of "A" functidln in Patankar Table 5.2, pg 95
    # pe: Peclet number (ratio of convection to diffusion)
    def aaa(pe):
        return max(0.0, (1.0 - 0.1 * abs(pe))**5)
     
    # mid-point of zisoil
    dzsoi[1:] = np.diff(zisoi)
    zsoi = zisoi-dzsoi/2.0 
          
    # Calculate the D and F terms in the Patankar algorithm
    # d: diffusivity,  m: layer above, p: layer below
    d_m1_zm1 = np.zeros((nlevbed+1),dtype=float)
    d_p1_zp1 = np.zeros((nlevbed+1),dtype=float)
    
    # f: flow for advection
    f_m1 = np.zeros((nlevbed+1),dtype=float)
    f_p1 = np.zeros((nlevbed+1),dtype=float)
    #
    pe_m1 = np.zeros((nlevbed+1),dtype=float)
    pe_p1 = np.zeros((nlevbed+1),dtype=float)
    # 
    for j in range(1, nlevbed+1):
        
        if j==nlevbed:
            print('checking!')
        
        if (j == 1):
            d_m1 = 0.
            d_m1_zm1[j] = 0. 
            
            # Weights for calculating harmonic mean of diffusivity
            #w_p1 = (zsoi[j+1]-zisoi[j])/dzsoi[j+1]
            w_p1 = (zsoi[j+1] - zisoi[j]) / (zsoi[j+1]-zsoi[j])
            if (diffus[j+1] > 0. and diffus[j] > 0.):
                #Harmonic mean of diffusivity
                d_p1 = 1./((1.- w_p1)/diffus[j] + w_p1/diffus[j+1])                
                #d_p1 = 1./((1.- w_p1)/diffus[j+1] + w_p1/diffus[j])                
            else:
                d_p1 = 0.         
            #d_p1_zp1[j] = d_p1 / dzsoi[j+1]
            d_p1_zp1[j] = d_p1 / (zsoi[j+1]-zsoi[j])

            #vwc_m1 = vwc[j]
            vwc_p1 = 1./((1.-w_p1)/vwc[j]+w_p1/vwc[j+1])
            #vwc_p1 = 1./((1.-w_p1)/vwc[j+1]+w_p1/vwc[j])   # appears this NOT makes large difference from above?
            
            f_m1[j] = 0
            f_p1[j] = adv[j+1]/vwc_p1
            pe_m1[j] = 0. 
            pe_p1[j] = f_p1[j]/d_p1_zp1[j]
        
        elif (j == nlevbed):
            # At the bottom, assume no gradient in d_z (i.e., they're the same)            
            #w_m1 = (zisoi[j-1]-zsoi[j-1])/dzsoi[j]
            w_m1 = (zisoi[j-1] - zsoi[j-1]) / (zsoi[j]-zsoi[j-1])            
            if (diffus[j]>0.0 and diffus[j-1]>0.0):
                d_m1 = 1./((1.-w_m1)/diffus[j]+w_m1/diffus[j-1])
                #d_m1 = 1./((1.-w_m1)/diffus[j-1]+w_m1/diffus[j])
            else:
                d_m1 = 0.           
            #d_m1_zm1[j] = d_m1 / dzsoi[j]
            d_m1_zm1[j] = d_m1 / (zsoi[j]-zsoi[j-1]) 

            d_p1 = 0.
            d_p1_zp1[j] = d_m1_zm1[j] # make sure 'pe_p1[j]' below not got nan

            #vwc_m1 = 1./((1.-w_m1)/vwc[j-1]+w_m1/vwc[j])
            vwc_m1 = 1./((1.-w_m1)/vwc[j]+w_m1/vwc[j-1])   # appears this NOT makes large difference from above?
            
            f_m1[j] = adv[j] / vwc_m1
            f_p1[j] = 0.            
            pe_m1[j] = f_m1[j]/d_m1_zm1[j]
            pe_p1[j] = f_p1[j]/d_p1_zp1[j]
      
        else:
            
            # Use distance from j-1 node to interface with j divided by distance between nodes
            #w_m1 = (zisoi[j-1] - zsoi[j-1]) / dzsoi[j]
            w_m1 = (zisoi[j-1] - zsoi[j-1]) / (zsoi[j]-zsoi[j-1])            
            if (diffus[j-1]>0. and diffus[j]> 0.):
                d_m1 = 1./((1.-w_m1)/diffus[j]+w_m1/diffus[j-1])
                #d_m1 = 1./((1.-w_m1)/diffus[j-1]+w_m1/diffus[j])
            else:
                d_m1 = 0. 
            #d_m1_zm1[j] = d_m1 / dzsoi[j]
            d_m1_zm1[j] = d_m1 / (zsoi[j]-zsoi[j-1]) 

          
            #w_p1 = (zsoi[j+1] - zisoi[j]) / dzsoi[j+1]
            w_p1 = (zsoi[j+1] - zisoi[j]) / (zsoi[j+1]-zsoi[j])
            if (diffus[j+1]>0. and diffus[j]>0.0):
                d_p1 = 1./((1.- w_p1)/diffus[j] + w_p1/diffus[j+1])
                #d_p1 = 1./((1.- w_p1)/diffus[j+1] + w_p1/diffus[j])
            else:
                #d_p1 = (1.-w_p1)*diffus[j]+w_p1*diffus[j+1]    # why different from d_m1 above?      
                d_p1 = 0.         
            #d_p1_zp1[j] = d_p1 / dzsoi[j+1]
            d_p1_zp1[j] = d_p1 / (zsoi[j+1]-zsoi[j])
            
            #vwc_m1 = 1./((1.-w_m1)/vwc[j-1]+w_m1/vwc[j])
            vwc_m1 = 1./((1.-w_m1)/vwc[j]+w_m1/vwc[j-1])   # appears this NOT makes large difference from above?
            vwc_p1 = 1./((1.-w_p1)/vwc[j]+w_p1/vwc[j+1])
            #vwc_p1 = 1./((1.-w_p1)/vwc[j+1]+w_p1/vwc[j])   # appears this NOT makes large difference from above?

            f_m1[j] = adv[j]   /vwc_m1         # adv[j] is from j-1 to j
            f_p1[j] = adv[j+1] /vwc_p1         # adv[j+1] is from j to j+1
            
            pe_m1[j] = f_m1[j]/d_m1_zm1[j]
            pe_p1[j] = f_p1[j]/d_p1_zp1[j]

    #


    # Calculate the tridiagonal coefficients
    # Coefficients of tridiagonal problem: a_i*x_(i-1) + b_i*(x_i) + c_i*x_(i+1) = r_i
    # Here, this is equivalent to Patankar equation 5.56 and 5.57 (but in one dimension):
    # a_P*phi_P = a_E*phi_E + a_W*phi_W + b [phi is concentration, = x in tridiagonal]. Converting East/West to above/below
    # -> -a_E*phi_E + a_P*phi_P - a_W+phi_W = b
    # -a_tri = a_above = D_above*A(Pe)+max(-F_above,0); D_above=diffus_above/dz
    # b_tri = a_above+a_below+rho*dz/dt
    # -c_tri = D_below*A(Pe)+max(F_below,0); D_below = diffus_below/dz
    # r_tri = b = source_const*dz + conc*rho*dz/dt
    
    a_tri = np.zeros((nlevbed+2),dtype=float)
    b_tri = np.zeros((nlevbed+2),dtype=float)
    c_tri = np.zeros((nlevbed+2),dtype=float)
    r_tri = np.zeros((nlevbed+2),dtype=float)
    
    for jx in range(nlevbed+2):

        if (jx > 0 and jx < nlevbed+1):
            a_p_0 =  dzsoi[jx]/dtime/vwc[jx]   # don't adjust using 'vwc', which already in 
            
      

        if (jx == 0): # top layer (atmosphere)
            a_tri[jx] = 0. 
            b_tri[jx] = 1. 
            c_tri[jx] = -1. 
            r_tri[jx] = 0. 
        
        elif (jx == 1):
            a_tri[jx] = -(d_m1_zm1[jx] * aaa(pe_m1[jx]) + max( f_m1[jx], 0.)) # Eqn 5.47 Patankar
            c_tri[jx] = -(d_p1_zp1[jx] * aaa(pe_p1[jx]) + max(-f_p1[jx], 0.))
            b_tri[jx] = -a_tri[jx] - c_tri[jx] + a_p_0
            
            # r_tri includes infiltration assuming same concentration as top layer.
            #  May want to change to either provide upper boundary condition or include in source term
            # r_tri[j] = source[j] * dzsoi[j] + (a_p_0 - adv[j]) * conc[j]
            r_tri[jx] = conc_rate[jx] * dzsoi[jx] + a_p_0 * conc[jx]
            if(qsrc[jx]>0): #infiltration
                r_tri[jx] = r_tri[jx] + qsrc[jx]*conc_surf
            elif(qsrc[jx]<0):         # upward flow to the surface or drainage
                r_tri[jx] = r_tri[jx] - qsrc[jx]*conc[jx]/vwc[jx]

        elif (jx < nlevbed+1):
            a_tri[jx] = -(d_m1_zm1[jx] * aaa(pe_m1[jx]) + max( f_m1[jx], 0.)) # Eqn 5.47 Patankar
            c_tri[jx] = -(d_p1_zp1[jx] * aaa(pe_p1[jx]) + max(-f_p1[jx], 0. ))
            b_tri[jx] = -a_tri[jx] - c_tri[jx] + a_p_0
            
            r_tri[jx] = conc_rate[jx] * dzsoi[jx] + a_p_0 * conc[jx] # Eq. 5.57
            # missing advection?
            
            
            # if lateral or bottom drainage occurs
            if(qsrc[jx]<0):        
                r_tri[jx] = r_tri[jx] - qsrc[jx]*conc[jx]
      
        else: # j==nlevbed+1; 0 concentration gradient at bottom
            a_tri[jx] = -1. 
            b_tri[jx] = 1. 
            c_tri[jx] = 0.  
            r_tri[jx] = 0. 

    #
    print('triag solver')
    
    ''''''
    du2, dl2, d2, x, info = \
        scilapack.dgtsv(c_tri[1:-1], b_tri[1:], a_tri[2:], r_tri[1:])

        #scilapack.dgtsv(c_tri[:-1], b_tri, a_tri[1:], r_tri)
        # (du, d, dl, b) - superdiagonal, main diagonal, subdiagonal, right-hand

    if(info < 0):
        print('dgtsv error in adv_diff line __LINE__: illegal argument')
        exit(info)
    if(info > 0):
        print('dgtsv error in adv_diff line __LINE__: singular matrix')
        exit(info)
    '''    
    x = tridiag_naive(a_tri, b_tri, c_tri, r_tri)
    '''
    
    conc_after = np.insert(x[:-1],0,np.nan)

    return conc_after

#

def mass_checking(nlevbed, dzsoi, vwc, \
                  conc, conc_after, conc_r, adv, conc_src, qsrc, dtime):
    
    mass_start = 0.  # mol/m2-soil
    mass_end   = 0.

    mass_in = 0.     # mol/m2-soil
    mass_out= 0.

    for j in range(1, nlevbed+1):        
        mass_start += conc[j]*dzsoi[j]                   
        mass_end += conc_after[j]*dzsoi[j]
        
        # reaction as src/sink
        if conc_r[j]>0:
            mass_in = mass_in + conc_r[j] * dzsoi[j] * dtime
        elif conc_r[j]<0:
            mass_out = mass_out + conc_r[j] * dzsoi[j] * dtime
            
        # flow src/sink, but out only assuming unknown external conc
        if qsrc[j]<0:
            mass_out = mass_out + conc[j]/vwc[j]*dtime

    # boundary        
    if adv[1]<0.0: # note < 0 is downwards
        mass_in -= adv[1]*conc_src*dtime        # mol/m2-soil: m/s * mol/m3 * s * m3/m3-soil
    elif adv[1]>0.0:
        mass_out += adv[1]*conc[1]*dtime        # mol/m2-soil: m/s * (mol/m3-soil) * m-soil/m * s
    
    if adv[nlevbed+1]<0.0:
        mass_out -= adv[nlevbed+1]*conc[nlevbed]*dtime    
        
    print('mass-change: ', mass_end - mass_start)
    print('mass-i/o: ', mass_in - mass_out)
    print('mass-balance: ', (mass_end - mass_start)-(mass_in - mass_out))
    
    return (mass_end - mass_start)-(mass_in - mass_out)
        
    
#
#--------------------------------------------------------------------
def test():

    """
    args = sys.argv[1:]
    input_path = args[0]
    output_path = args[1]
    """  
    input_path  = './'
    output_path = './'

    # default elm soil nodes and layers    
    zsoi, dzsoi, zisoi = elm_soillayers(nlevgrnd=15, equal_thickness=True)

    # ----- randomly generated inputs to test functions -------------------
    # insert np.nan into [0], so that indices matching with Fortran-style 
    
    # water flow caused transport (advection, leaching etc)
    OFF_QADV   = True
    NO_QSRC    = True
    UNI_VWC    = False
    
    # solute diffusion + reactions
    UNI_CONC   = False
    if UNI_CONC: UNI_VWC = True   # make sure conc in pore water uniform, because conc in unit of mol/m3-soil
    UNI_DIFUS = False
    ZERO_DIFUS = False

    
    NO_CONC_R  = True
    
    #-------------------
    nlevbed = 2
 
    porosity    = np.random.rand((nlevbed+1)) *0.10 + 0.40         # range: 0.4 ~ 0.5
    porosity[0] = np.nan
       
    # soil water and its flow
    sat         = np.random.rand((nlevbed+1)) *0.95 + 0.05         # saturation for Volumetric soil moisture in layer (m3/m3 or m/m)
    vwc         = sat * porosity                                   # range: 0.02~0.50
    if UNI_VWC or UNI_CONC: vwc[1:] = np.nanmean(vwc)
    sat[0]      = np.nan

    qadv        = (np.random.rand((nlevbed+2))-0.5)*1.e-6          # (mH2O/s), vertical into layer via interface (up: positive, down: negative)
    if OFF_QADV: qadv[...]   = 0.  # advection will be off
    qadv[0]     = np.nan

    qsrc        = np.random.rand((nlevbed+1))*1.e-8                # water Source/sink term (m/s: in +, out -)
    if NO_QSRC: qsrc[...]   = 0.       # close-water system
    #qsrc[1]     = 1.e-5    # infiltration
    #qsrc[1]     = -1.e-5  # exfiltration    
    #qsrc[-1]    = -1.e-5  # bottom drainage    
    qsrc[0]     = np.nan

    # soil pore water solution and its transport properties
    conc        = np.random.rand((nlevbed+1)) *1.e-7               # Bulk concentration (e.g. mol/m3-soil)
    if UNI_CONC: conc[...] = 1.e-7  # diffusion will be off 
    conc[0]     = np.nan
    
    diffus      = np.random.rand((nlevbed+1)) *1.e-9               # diffusivity (m2-soil/s)
    if ZERO_DIFUS: 
        diffus[...] = 0.  # another way of diffusion will be off 
    elif UNI_DIFUS:
        diffus[...] = 1.e-9
    diffus[0]   = np.nan
    
    conc_dt      = np.random.rand((nlevbed+1))*1.e-6             # Bulk concentration rate, e.g. reaction or other src/sink (mol/m3-soil/s)
    if NO_CONC_R: conc_dt[...] = 0.   # non-reaction
    conc_dt[0]   = np.nan

    conc_surf   = 0.e-6                             # Surface boundary layer concentration (for infiltration): mol/m3-water
    
    dt       = 1800.0                            # Time step (s)
    



    #---- solutions
    '''
    conc_perwater = conc/vwc # harmonize unit
    qdarcy = -(qadv[1:] * np.insert(vwc[1:], 0, 1)) # my x is defined to be positive downwards
    conc_next, niter = a_step_explicit(nlevbed, zsoi[1:], dzsoi[1:], conc_perwater[1:], \
        conc_surf, qdarcy, vwc[1:], porosity[1:], conc_dt[1:], diffus[1:], dt)
    conc_after = conc_next*vwc[1:]
    conc_after = np.insert(conc_after, 0, np.nan)
    
    
    '''
    conc_perwater = conc/vwc
    conc_r = conc_dt/vwc
    vwc_perwater = np.ones_like(conc_perwater)   # conc in unit of per vol. water, so no need to consider vwc
    conc_next = Patankar_solute_transport(nlevbed, dzsoi, zisoi, \
        conc_perwater, conc_r, qadv, diffus, \
        vwc_perwater, qsrc, conc_surf, \
        dt)
    conc_after = conc_next*vwc
    ''''''
      
    # ---- mass balance checking
    dmass_error = mass_checking(nlevbed, dzsoi, vwc, \
                                conc, conc_after, conc_dt, qadv, conc_surf, qsrc, dt)

    print('done transport solving with mass, mol/m2-soil, error of  : ', dmass_error)
    
if __name__ == '__main__':
    test()


