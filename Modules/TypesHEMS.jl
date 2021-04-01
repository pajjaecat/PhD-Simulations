### TypesHEMS.jl
###
### Standar library containing all the types used in modeling the solar home

""" 
    Standar library containing all the types used in modeling the solar home
"""
module TypesHEMS

  export
    # Types & aliases
    AbstractProblemVars,
    AbstractFixedVars,
    AbstractControledVars,
    FixedVars,
    ControledVars,
    ControledVarsSum,
    solvedProbData;





#---------------------------------------------------------------------------------------------------
# --------------------------                Abstrat types                 --------------------------
#---------------------------------------------------------------------------------------------------

abstract type AbstractProblemVars end #abstrat type to hold all the problem's variables

abstract type AbstractFixedVars<:AbstractProblemVars end  #abstrat type for fixed data of the problem

abstract type AbstractControledVars <:AbstractProblemVars end  #abstrat type for controlled variables 








#--------------------------------------------------------------------------------------------------
# --------------------                   Composite  types                 -------------------------
#--------------------------------------------------------------------------------------------------

  struct FixedVars{T<:AbstractFloat} <: AbstractFixedVars 

      dt::T
      Pv_max::T
      Eb_init::T
      Eb_min::T
      Eb_max::T
      Pg_max::T

  end
#---------------------------------------------------------------------------

  struct ControledVars{T<:AbstractFloat} <: AbstractControledVars 

      Pg::Array{T, 1}
      Pb::Array{T, 1}
      Pw::Array{T, 1}
      Ps::Array{T, 1}

  end
#---------------------------------------------------------------------------

  struct ControledVarsSum{T<:AbstractFloat} <: AbstractControledVars 

      Pg_t::T
      Ps_t::T
      Pw_t::T
      Cost_t::T

  end
#----------------------------------------------------------------------------

"""
    solvedProbData{T<:AbstractFloat} <: AbstractProblemVars

subtype of `AbstractProblemVars` to hold a solved problem data


# Arguments 
- `t`
- `sim_tim::Integer`
- `nb_days::Integer`
- `ctrlVars::ControledVars`
- `Eb_opt::Array{T,1}`    
- `Pl_star::Array{T,1}`  
- `Pv_star::Array{T,1}` 

where `{T <: AbstractFloat} `


"""
 struct solvedProbData{T<:AbstractFloat} <: AbstractProblemVars
      t
      sim_tim::Integer
      nb_days::Integer
      ctrlVars::ControledVars
      Eb_opt::Array{T,1}    
      Pl_star::Array{T,1}  
      Pv_star::Array{T,1}  
  end 



  
end
