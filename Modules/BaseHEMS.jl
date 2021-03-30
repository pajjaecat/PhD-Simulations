""" 
    Standar library containing the definition of all the basic types and functions used in modeling the solar home
"""
module BaseHEMS


#import modules that will be used inside the module
using JuMP,PyPlot,GLPK, SDDP;
import Statistics;

#so that the function can be used via stoProg.functionName
export FixedVar, ControledVar, solar_HomData1, InitFixedVar, checkGPO_state,ruleBasedControl;


#---------------------------------------------------------------------------------------------------
# --------------------------                Abstrat types                 --------------------------
#---------------------------------------------------------------------------------------------------

abstract type AbstractProblemVars end #abstrat type to hold all the problem's variables

abstract type AbstractFixedVar<:AbstractProblemVars end  #abstrat type for fixed data of the problem

abstract type AbstractControledVar <:AbstractProblemVars end  #abstrat type for controlled variables 












#--------------------------------------------------------------------------------------------------
# --------------------                   Composite  types                 -------------------------
#--------------------------------------------------------------------------------------------------

struct FixedVar{T<:AbstractFloat} <: AbstractFixedVar 
  
    dt::T
    Pv_max::T
    Eb_init::T
    Eb_min::T
    Eb_max::T
    Pg_max::T
  
end


struct ControledVar{T<:AbstractFloat} <: AbstractControledVar 

    Pg::Array{T, 1}
    Pb::Array{T, 1}
    Pw::Array{T, 1}
    Ps::Array{T, 1}
  
end


struct solar_HomData1
        t
        sim_tim
        nb_days
        Pgr_opt
        Psh_opt
        Pst_opt
        Pnl
        Pcu_opt
        Est_opt
        Pgr_t 
        Psh_t
        Psu_t
        Pld_t
        Pcu_t
        Cout_cons
end 








#---------------------------------------------------------------------------------------------------
#'-------------------------------                 Functions            -----------------------------'
#---------------------------------------------------------------------------------------------------

"""
    InitFixedVar(;<keyword arguments> )

    Initialize the fixed problem's variable 
`Return` a structure (composite type) FixedVar

# Keyword Arguments
- `dt::T=0.5`
- `Pv_max::T=3`: 
- `Eb_init::T=0`: [kWh] Initial storage for the simulation
- `Eb_min::T=0`: [kWh] Minimum reserved storage in nominal operation
- `Eb_max::T=8`: [kWh] Maximum Storage capacity
- `Pg_max::T=3`: [kW] Maximum power drawn from the grid when available
where {T <: AbstractFloat} 
"""
function InitFixedVar(;dt::T =.5, Pv_max::T =3., Eb_init::T = 0., Eb_min::T = 0., Eb_max::T = 8.,
    Pg_max::T =3.) where {T<:AbstractFloat}
    
   return FixedVar(dt, Pv_max, Eb_init, Eb_min, Eb_max, Pg_max)
    
end


#***************************************************************************************************
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
    
    for rowzz in 1: size(GpoPeriods,1)
        
       if GpoPeriods[rowzz,1] <= cur_ins <= GpoPeriods[rowzz,2]
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
    ruleBasedControl(Pl_star, Pv_star, fixedVariables=InitFixedVar(); <keyword arguments> )

Compute the energy dispatch using the Rulebased algorithm described in Figure 3.9 of the 
thesis.


# Arguments
- `Pl_star::Array{Real,1}`: load power request vector
- `Pv_star::Array{Real,1}`: solar power vector
- `fixedVariables::AbstractFixedVar=InitFixedVar()`: Fied variables structure 

# Keywords arguments 
- `GpoPeriods::Array{T,2}=[0 0]`: Array with each row defining a GPO period `where {T<:Integer}`
- `powerSave::Real=0`: Power save option in percentage of the nominal load demande (0.5 for 0.5Pl_star for instance). By default, no power save.
                

""" 
function ruleBasedControl(Pl_star, Pv_star,fixedVariables::AbstractFixedVar= InitFixedVar(); powerSave::Real=0,GpoPeriods::Array{T,2}=[0 0])  where {T<:Integer}
    
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
  
  return ControledVar(Pg, Pb,Pw, Ps), Eb;
  
end     





  
end
