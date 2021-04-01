""" 
    Standar library containing all the basic types and functions used in modeling the solar home
"""
module BaseHEMS


#import modules that will be used inside the module
using JuMP,PyPlot,GLPK, SDDP;

using TypesHEMS # My types are defined in the file TypesHEMS.jl 

import Statistics;


export 
  build_scnr,
  mc_scrPrb,
  nMC_trans_tab,
  hp_hc,
  plot_solHomRslt,
  plot_OlfcVsMso,
  build_grid_snr,
  get_snr,
  checkGPO_state,
  ruleBasedControl;


#---------------------------------------------------------------------------------------------------
#'-------------------------------                 Functions            -----------------------------'
#---------------------------------------------------------------------------------------------------

"""
    initFixedVars(;<keyword arguments> )

    Initialize the fixed problem's variable 
`Return` a structure (composite type) FixedVars

# Keyword Arguments
- `dt::T=0.5`
- `Pv_max::T=3`: 
- `Eb_init::T=0`: [kWh] Initial storage for the simulation
- `Eb_min::T=0`: [kWh] Minimum reserved storage in nominal operation
- `Eb_max::T=8`: [kWh] Maximum Storage capacity
- `Pg_max::T=3`: [kW] Maximum power drawn from the grid when available
where {T <: AbstractFloat} 
"""
function initFixedVars(;dt::T =.5, Pv_max::T =3., Eb_init::T = 0., Eb_min::T = 0., Eb_max::T = 8.,
    Pg_max::T =3.) where {T<:AbstractFloat}
    
   return FixedVars(dt, Pv_max, Eb_init, Eb_min, Eb_max, Pg_max)
    
end


#***************************************************************************************************
"""
    initControledVars(Pg, Pb, Pw, Ps; <kwargs> )

    Initialize the controled variable and/or compute its sum.
`Return` either a structure (composite type) ControledVars or ControledVarsSum

# Arguments
- `Pg::Array{T, 1}`
- `Pb::Array{T, 1}` 
- `Pw::Array{T, 1}`
- `Ps::Array{T, 1}`
where {T <: AbstractFloat} 



# Keywords Arguments
- `opt::AbstractString`: Ignore or set to `"sum"` to compute the sum of the 
controled variable over the problème horizon
 
"""
function initControledVars(Pg::Array{T, 1}, Pb::Array{T, 1}, Pw::Array{T, 1}, Ps::Array{T, 1}; opt::AbstractString = " " ) where {T<:AbstractFloat}
  
  if isequal(opt," " )
    return ControledVars(Pg, Pb,Pw, Ps)
    
  elseif isequal(opt,"sum" )
    Pg_t, Ps_t, Pw_t,Cost_t = round.((Pg, Ps, Pw, Pg*0.2).*0.5 .|> sum, digits=2);
    return ControledVarsSum(Pg_t, Ps_t, Pw_t,Cost_t)
    
  else
    error("Wrong keywords argument. See the function description to set it right")
  end
    
end


#***************************************************************************************************
"""
    initControledVars(CtroledVars; <kwargs> )

    Compute the sum of the controlled variables
`Return` a structure (composite type) ControledVarsSum

# Arguments
- `ctroledVars::ControledVars`


# Keywords Arguments
- `opt::AbstractString="sum"`: ignore or set to `"sum"` to compute the 
sum of the controled variable over the problème horizon
 
"""
function initControledVars(ctroledVars::ControledVars; opt::AbstractString = "sum" )
  
  if isequal(opt,"sum" )
    ctrVr=ctroledVars;
    Pg_t, Ps_t, Pw_t,Cost_t = round.((ctrVr.Pg, ctrVr.Ps, ctrVr.Pw, ctrVr.Pg*0.2).*0.5 .|> sum, digits=2);
    return ControledVarsSum(Pg_t, Ps_t, Pw_t,Cost_t)
    
  else
    error("Wrong keywords argument. See the function description to set it right")
  
  end
    
end


#******************************************************************************************************
"""
    checkGPO_state(cur_ins, GpoPeriods=[0 0])

    Verify if a GPO is accuring at the current instant.
`Return`  boolean `true` if yes or `false` otherwise

# Arguments
- `cur_ins::T`:  Index of the current instant
- `GpoPeriods::Array{T,2}` =[0 0]: Array with each row defining a GPO periode
where {T <: Integer} 

# Example

julia> using BaseHEMS

julia> checkGPO_state( 100, [54 78; 120 153] )

true
"""
function checkGPO_state(cur_ins::T, GpoPeriods::Array{T,2}=[0 0]) where {T<:Integer}
    
    varr = false 
    
    for rowz in 1: size(GpoPeriods,1)
        
       if GpoPeriods[rowz,1] <= cur_ins <= GpoPeriods[rowz,2]
            varr = true 
            break
        else
            varr = false
        end
    end
    
    return varr

end


#***********************************************************************************************
"""
    ruleBasedControl(Pl_star, Pv_star, fixedVariables=initFixedVars(); <keyword arguments> )

Compute the energy dispatch using the Rulebased algorithm described in Figure 3.9 of the 
thesis.

`Return`  a subtype `AbstractControledVars` and  `E_b::Array{<:AbstractFloat,1+N}` with N the problem's Horizon.


# Arguments
- `Pl_star::Array{Real,1}`: load power request vector
- `Pv_star::Array{Real,1}`: solar power vector
- `fixedVariables::AbstractFixedVars=initFixedVars()`: Fied variables structure 

# Keywords arguments 
- `GpoPeriods::Array{T,2}=[0 0]`: Array with each row defining a GPO period `where {T<:Integer}`
- `powerSave::Real=0`: Power save option in percentage of the nominal load demande (0.5 for 
0.5Pl_star for instance). By default, no power save.
                

""" 
function ruleBasedControl(Pl_star::Array{U,1}, Pv_star::Array{U,1},
    fixedVariables::AbstractFixedVars= initFixedVars(); powerSave::Real=0, 
    GpoPeriods::Array{T,2}=[0 0])  where {T<:Integer, U<:Real}
    
  Pg_max = fixedVariables.Pg_max;
  dt = fixedVariables.dt;
  Eb_max = fixedVariables.Eb_max;
  
  N = length(Pl_star) #get the probleme horizon from the length of the vector Pl_star

  ## create output of the function ,i.e, control variables and storage evolution 
  Pg = zeros(N) ; Pw = zeros(N); Pb = zeros(N);  Eb = zeros(N+1); 
  Ps = zeros(N);  Pnl = zeros(N)
  
  Eb[1] = fixedVariables.Eb_init; #Initialize the first element (initial storage) of the storage vector 

    for k in 1:N

          if checkGPO_state(k, GpoPeriods)
              Pg_max = 0; 
              Pg[k] = 0;
              Eb_minF = 0;
              if powerSave != 0
                  Ps[k] = powerSave*Pl_star[k];
                  Pl_star[k] = Pl_star[k] - Ps[k];
              end
          else
              Pg_max = 3; 
              Eb_minF = fixedVariables.Eb_min;
              
              if k>=(N-20) Eb_minF=0 end  # Make sure to empty the battery to relaxe the 
              #minimum storage towards the ends of the simulation 
          end

          Pnl[k] = Pl_star[k] - Pv_star[k];
    
          if Pnl[k] > 0 
              Pw[k] = 0
              Diff = Eb[k] - Pnl[k]*dt;
              if Diff >= Eb_minF
                  Pb[k] = -Pnl[k]
                  Pg[k] = 0;
              elseif (Diff < Eb_minF) && (Pg_max ==3)
                  Pb[k] = 0;
                  Pg[k] = min(Pg_max,Pnl[k])
                  Ps[k] = Ps[k] + abs( min(0,Pg_max-Pnl[k]) )
              elseif Diff < Eb_minF && Pg_max ==0
                  Pg[k] = 0 ;
                  Pb[k] = -(Eb[k] - Eb_minF)/dt
                  Ps[k] = Ps[k] + (Pnl[k] + Pb[k])
              end
          else
              Pg[k] = 0 ;
              Diff = Eb[k] - Pnl[k]*dt;
              if Diff > Eb_max
                  Pb[k] = (Eb_max-Eb[k])/dt;
                  Pw[k] = (Diff - Eb_max)/dt;
              else
                  Pb[k] = - Pnl[k];
                  Pw[k] = 0 ;
              end 
          end
      Eb[k+1] = Eb[k] + Pb[k]*dt;
    end
  
  return initControledVars(Pg, Pb,Pw, Ps), Eb;
  
end     



#***********************************************************************************************
"""
    plot_solHomRslt(args...)

    Plot the resolved optimization problem on the the time frame considered 

# Arguments

- `args::`
- `args[1] ::Structure`: Output structure solarHome_data

"""
function plot_solHomRslt(args...; fig_size::Tuple{T,T}=(8,5)) where {T<: Real}
    
      #Unpack input
      input=args[1];
      t=input.t
      sim_tim = input.sim_tim
      nb_days = input.nb_days
  
      Pg_opt = input.ctrlVars.Pg
      Psh_opt = input.ctrlVars.Ps
      Ps_opt = input.ctrlVars.Pb
      Pcu_opt = input.ctrlVars.Pw
  
      P_nl = input.Pl_star - input.Pv_star
      P_nlSum = round((sum(P_nl)),digits=2)
      Es_opt = input.Eb_opt

      var = initControledVars(input.ctrlVars);
      P_gr_t, Psh_t, P_c_t, Cout_cons = var.Pg_t, var.Ps_t, var.Pw_t,var.Cost_t ;


  
  
      #Plot data on First subfigure------------------------------------------------ 
      fig, ax = subplots(2,1, figsize=fig_size,sharex=true)
      line_pgr, = ax[1,1].plot(t,Pg_opt', color="b", label=L"P_{g}");
      line_nl, = ax[1,1].plot(t,P_nl[1:sim_tim]',color="black",label=L"P_{nl}^*",alpha=0.6,lw=.7);
      line_sh, = ax[1,1].plot(t,Psh_opt', color="tab:red",label=L"P_{s}");
      line_cu, = ax[1,1].plot(t,Pcu_opt[1:sim_tim]',color="gold",label=L"P_{sp}",lw=2);
      line_sto, = ax[1,1].plot(t,Ps_opt', color="green",label=L"P_{b}",alpha=0.7,lw=2);

      
      ax[1,1].set(
          xlim=(0,nb_days+.025),
          ylim =(minimum(P_nl), 3.2),
          Yticks = ([-1,0,1,2,3]),
          ylabel="P(kW)")
      ax[1,1].grid();
      
  
      #Plot data on second subfigure ------------------------------------------------------ 
      t_sto = range(0,(sim_tim/2),step=.5)/24
  
      # If more thatn only one argument which is the structure SolarHome
      # compute at each instant Pnl_hor over the horizon given by the second argument.

      line_esto, = ax[2,1].plot(t_sto,Es_opt[1:end]',color="green",ls="-.",label=L"$E_{b}$",alpha=0.7,lw=2);

      ax[2,1].set(
          xlim=(0,nb_days+.025),
          xlabel="Time(Day)",
          ylabel=L"E_{b}(kWh)",
          Yticks=([0,2,4,6,8]),
          ylim=(-.3, 8.3)
      )
      ax[2,1].grid()
  
  
      #Highligh off peak periods
      ~,x_max = ax[1,1].get_xlim()
      for i in 0:1:x_max
          ax[1,1].axvspan(i, i+0.25,color="dodgerblue", alpha=0.1)
          ax[2,1].axvspan(i, i+0.25,color="dodgerblue", alpha=0.1)

      end
  
      # If more than one argument, color in dotted red the moment where the outage apears and 
      # in a a shade of red the period where the outage stays.  
      if length(args) >=3
              ax[1,1].axvspan(1, nb_days+.025,color="lightcoral",alpha=0.02)
              ax[2,1].axvspan(1, nb_days+.025,color="lightcoral",alpha=0.02)
              ax[1,1].axvline(1,color="r",ls=":")
              ax[2,1].axvline(1,color="r",ls=":")
      end

      
      Es_fin = round(Es_opt[end-1],sigdigits=4)
      fig.legend(ncol=6,loc="center",columnspacing=.8,handlelength=1,bbox_to_anchor=(0.48, 0.85))
       fig.suptitle("\$ P_{g} =$P_gr_t,  P_{nl}=$P_nlSum || kWh, Bill = $Cout_cons € \$ ");
      #fig.suptitle("ℙ fail: $p1fe3×10⁻³, shedding cost: $(C_shed[1]/1e3)×10³, 𝔼(J): $EJ")

end



  


  
end
