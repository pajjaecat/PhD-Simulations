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
- `args[2] ::Any`: Any type of data is relevant

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








""" Build scenarios trees (Pgrid_max_scenario and Pgrid_max_scenario_prob)"""

function build_scnr(K, proba_fail)
 
# Uncertainty modeling: iid sequence of Bernoulli trials
    P_gm_prob = [1-proba_fail proba_fail]
    P_gm_n = length(P_gm_prob) # 2

    # Build scenarios for $P_{grid}^{max}$
    # nb scenarios
    S = P_gm_n^(K-1)

    P_gm_scenar_proba =  zeros(S)
    P_gm_scenar_state = zeros(Int, S,K)
    
    # create scenarios by enumerating all possibilities
    for k in 1:K
        step = P_gm_n^(K-k)
        state = 0
        for s in 1:step:S
            P_gm_scenar_state[s:s+step-1,k] .= state
            state = (state + 1) % 2
        end
    end
    
    # Translate [0,1] to [1,2]
    P_gm_scenar_state .+= 1

    P_gm_vals = [3. 0.0]
    # From state to real value
    P_gm_scenar = P_gm_vals[P_gm_scenar_state]

    return P_gm_scenar
    #return P_gm_scenar,P_gm_scenar_proba
end



#2---------------------------------------------------------------------------------
""" 
    Takes as INPUTS the initial state vector (init_st = [state1-ON state2-OFF]
                                                        [1 0]  i.e MC start at ON 
                                                    or  [0 1]) i.e MC start at OFF,
        the MC transition matrix (trans_mat) and the recursive horizon (rcsv_hr)
        RETURN a 2rows table which is the n step transition probability of being 
        in state1-ON (first row) and state2-OFF (second row)     
"""

function nMC_trans_tab(init_st, trans_mat,rcsv_hr)
    
    tabb = zeros(2,rcsv_hr);

    for k in 1:rcsv_hr
        tabb[:,k] = init_st*trans_mat^(k-1);
    end
    
    return tabb;   
end



#3-----------------------------------------------------------------------------------
"""Build the MC scenario probability tree"""
function mc_scrPrb(tran_mat,Pgm_scTree )
    # tran_mat : Transition probability matrice
    #  Pgm_scTree  : Pgrid_max scenario tree 

    λ = tran_mat[1,2];
    μ = tran_mat[2,1];
    row_s = length(Pgm_scTree[:,1])
    col_m = length(Pgm_scTree[1,:])
    mc_scnrProb = zeros(row_s,col_m-1) #build output table 
    
#I can clearly code this more efficiently.. TODO LATER   
    for i in 1:row_s # for each scenrio
        
        for j in 1: col_m-1 # for each colmn 
            
            # fill the mc_scnrProb according to the current and next state
            if (Pgm_scTree[i,j] ==3) && (Pgm_scTree[i,j+1]==0)     #ON-OFF
                mc_scnrProb[i,j] = λ;
                
            elseif (Pgm_scTree[i,j] ==3) && (Pgm_scTree[i,j+1]==3)  #ON-ON
                mc_scnrProb[i,j] = 1-λ;
                
            elseif (Pgm_scTree[i,j] ==0) && (Pgm_scTree[i,j+1]==3) #OFF-ON
                mc_scnrProb[i,j] = μ;
            
            else                                                   #OFF-OFF
                mc_scnrProb[i,j] = 1-μ;
            end
        end 
    end
    
    #Compute the product of transition prob at each instant of each scenario
    out_tab = prod(mc_scnrProb; dims=2);
    
    return out_tab;
end 


##4 -----------------------------------------------------------------

function hp_hc(args...)
    
#=This function creates a collumn vector that contains the price associated
to an off-peak and peak period of time over the number of days given as
input.
The day starts at 00h00 then, ends at 23h30 which corresponds to 48
slots of 30 minutes.  The first input 'nb_days' corresponds to the number
of days while the second input 'off_peak_itvl' is a row vector of two collums respectively associated to the starting and end time of the  off-peak period. The third input is the coef that multiply the base energy 
% prices i.e. to modify the price of the peak and off peak period =#
 
    nargin = length(args)
    off_peak_itvl = args[2];
    nb_days = args[1]
    
    if nargin < 3 # if no third argument create and instance it at 1; 
        coef = 1;
    end
    
    p_peak = .2*coef;
    p_offPeak = .1*coef;

    #Create a day of 48 slots(1/2 hour per slot) with a peak price associated
    #to each slot.
    day_one = p_peak*ones(1,48); 
    
    st_offP = off_peak_itvl[1];
    ed_offP = off_peak_itvl[2];

    #Set the off peak price for the off peak period.
    day_one[st_offP*2+1:ed_offP*2] .= p_offPeak;

    #Repeat day one for the number of days give as insput
    out_var = repeat(day_one,1,nb_days);
    
    return out_var
    
end




#5___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
### Set Plot Fonts

"""
   setPlotFonts(Font_Nrm_size, Font_Big_size)
Call this function when ploting. See 

# Arguments

- `Font_Nrm_size::String`: Normal Font size as a Sring. eg "14"
- `Font_Big_size::String`: Big Font size as a Sring. eg "16"


"""
function setPlotFonts(Font_Nrm_size, Font_Big_size,Font_BIG_size::Any=15)
  

    # Control figure front size------------------------------------------------
    #plt.rc("font", size=Font_Nrm_size)          # controls default text sizes
    plt.rc("axes", titlesize=Font_BIG_size)     # fontsize of the axes title
    plt.rc("axes", labelsize=Font_Big_size)     # fontsize of the x and y labels
    plt.rc("xtick", labelsize=Font_Nrm_size)    # fontsize of the tick labels
    plt.rc("ytick", labelsize=Font_Nrm_size)    # fontsize of the tick labels
    plt.rc("legend", fontsize=Font_BIG_size)    # legend fontsize
    plt.rc("figure", titlesize=Font_BIG_size)   # fontsize of the figure title
  
    #csfont = {"fontname":"Comic Sans MS"}
    #rc("font", **csfont)
  
end 






#7___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
### PLot Solar home results 

"""
    plot_OlfcVsMso(olfc,mso,mso1,nrows::real=2,ncol::real=3,titled)

    Plot data given as input 

# Arguments

- `olfc::`
- `mso ::Structure`: Output structure solarHome_data
- `mso1::Any`: Any type of data is relevant
- `nrows::real`: 
- `ncol::real`: 
- `titled::String`

"""
function plot_OlfcVsMso(def_rate,olfc,mso,mso1,nrows,ncol,titled,Hor_vect)
  
    #fig, ax = subplots(nrows,ncol,figsize=(8,5.5),sharex=true,sharey=true,)
      fig, ax = subplots(nrows,ncol,figsize=(10,6.5),sharex=true,sharey=true,)

    lin=1;
    colm = 1;
    
    for ii in 1:nrows*ncol
        if ii==ncol+1
            lin=2
            colm = 1;
        end 

        olfc_line, = ax[lin,colm].semilogx(def_rate,olfc[ii,:]',"x-",color="red")
        mso_line, = ax[lin,colm].semilogx(def_rate,mso[ii,:]',"x-",color="blue",lw=4,alpha=.3);
        mso_line2, = ax[lin,colm].semilogx(def_rate,mso1[ii,:]',"x:",color="green",lw=2,alpha=.8);

        ax[lin,colm].grid();
        kr = Hor_vect[ii]/2;
        ax[lin,colm].set_title("Hor(hour)=$kr");

        colm+=1;

        if ii== nrows*ncol
            fig.legend(
                [olfc_line, mso_line,mso_line2],
                ("OLFC","Mso_Recrs","MSO_Wait&See"), 
                loc="lower center",
        
                ncol=3,
                columnspacing=1,
                handletextpad=0.3);
        end
    end

    fig.tight_layout(rect=(0, 0.07, 1., .95));
    fig.suptitle(titled);
end 




#8___________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
### check_stateTrans 


""" 
check_stateTrans(curSta_in, Pg_max_prof)
    
Check the transition from the previous to the current state, returning  an 
integer value chanState_status as =>  1: On  - On 
                                      2: On  - off
                                      3: Off - off
                                      4: Off - On.

# Arguments
- `curSta_in::Int`: Current state indice: The indice of the current state of simulation`
- `Pg_max_prof::Array{Real,1}`:max grid power vector on the Mpc horizon
    
"""
function check_stateTrans(curSta_in, Pg_max_prof)
    # Check if the current Pgrid is the same as the previous one. Current
    # state being i+1 and previous i, because the state initial-1 has been 
    # added into the fisrt slot of Pg_max_prof_f which introduced the shift 
        
    if (Pg_max_prof[curSta_in] ==3) && (Pg_max_prof[curSta_in-1]==3) #ON-ON
        chanState_status = Int(1);
                
    elseif (Pg_max_prof[curSta_in] ==3) && (Pg_max_prof[curSta_in-1]==0) #ON-OFF
        chanState_status = Int(2);
                
    elseif (Pg_max_prof[curSta_in] ==0) && (Pg_max_prof[curSta_in-1]==0) #OFF-OFF
        chanState_status = Int(3)
            
    else    #OFF-ON
        chanState_status = Int(4) ;
    end
        
return chanState_status;

end



#10__________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
### build_grid_snr



# Function used in Sddp2VsOlfc located at Home/My_Jupiter/GitJupLab/JuliaMyLOve/DynProg/Sddp2VsOlfc.ipynb
function build_grid_snr(K,nb_snr,lambd,muu)
  
  λ = lambd; μ = muu;
  A = Array{Float64, 2}[[ 1.0 ]',
   [ 1-λ  λ]];
  B = kron(ones(K-2),[[1-λ  λ ; μ 1-μ]]);
  C = [A;B];
  
   modell = SDDP.MarkovianPolicyGraph(
      transition_matrices = C,
      sense = :Min,
      lower_bound = 0.0,
      optimizer = GLPK.Optimizer
  ) do modell, node
      # Unpack the stage and Markov index.
      t, markov_state = node

      pgmax = if markov_state == 1  # Grid available
          3.
      else  # Grid Unvailable
          0.
      end
      @variable(modell, 0 <= Eb <= 5, SDDP.State, initial_value = 0)
      # Define the control variables.
       @variables(modell, begin
          0 <= P_g <= pgmax
          P_b
      end)
      # Define the constraints
      @constraints(modell, begin
          Eb.out == Eb.in + P_b
      end)

      @stageobjective(modell, 10*P_b + 5*P_g)

      end;
  
  
    SDDP.train(modell; iteration_limit = 2);
    simulations = SDDP.simulate(modell,nb_snr);

    ab = [ [ elm[:node_index] for elm in sim] for sim in simulations]
    snr=zeros(length(ab),K)

    for j=1:length(ab)
        ll = [ell[2] for ell in ab[j]];
        ss = isequal.(ll,2)

        i=1;
        for elm in ss
            if elm==true
                snr[j,i] = 0
            else 
                snr[j,i] = 3
            end
            i+=1
        end
    end 
  return snr,simulations;
  
end
     


#11__________________________________________________________________________________________________________
#___________________________________________________________________________________________________________
### get_snr


function get_snr(simul,K)
  
    ab = [ [ elm[:node_index] for elm in sim] for sim in simul]
    snr=zeros(length(ab),K)

    for j=1:length(ab)
        ll = [ell[2] for ell in ab[j]];
        ss = isequal.(ll,2)

        i=1;
        for elm in ss
            if elm==true
                snr[j,i] = 0
            else 
                snr[j,i] = 3
            end
            i+=1
        end
    end 
  return snr;
  
end
            







#___________________________________________________________________________________________________________
### Forecasting


#----------------------------------------- for the thesis whriting ------------------------------------
## Forecasting algoritm 
"""
    forecasting(Pv_s_year,BaH,MpcH)

Computes the energy dispatch using Rulebased algorithm described in \\label{fig:RuleBasedConHEMS} of the 
thesis, returning the the control variables and storage evolution  Pg, Pw, Pb, Eb; 
NB: Failure of the grid is not acounted for here. Also Pg_max is choose such that it will always cover the 
maximum demand of the building, i.e., Pg_max > max(Pl_star) hence, the load shedding variable Ps is not needed.

# Arguments
- `Pv_s_year::Array{Real,1}`: Yearly load demand request vector
- `BaH::Real`: Backward Horizon
- `MpcH::Real`: Forward MPC horizon
- `k::Real`: current instant of the simulation
- `Eva::AbstractString`:="dist" Select the method to evaluate the distance. Must be "Dist" or "CorCoef" or "Both"

""" 

function forecasting(Pl_s_year, BaH, MpcH,k, Eva::AbstractString="Dist")
  LPv = length(Pl_s_year);
  Pls_frcst = zeros(1,MpcH);
  Pls_frcst2 = zeros(1,MpcH);
  PLs = 100*ones(366,MpcH) 
  PLs_cutErr = ones(365,1)
  PLs_cutErr2 = ones(365,1)

  
  Pl_s_yearR = [Pl_s_year;Pl_s_year[1:100]] 
  # Need to recode for MPC horizon different from 24h i.e. 48
  
  curr_day = Pl_s_year[k-BaH:k-BaH+MpcH-1];

  var1=1; # Count the current line in the PLs Table
    
  if k-BaH > MpcH # if there is enough data before current day to fill at least a whole MpcH 
  var2 = convert(Integer, floor((k-BaH-1)/MpcH)); # compute how many iterations is possible 
                                                  # backward from  the beginning of current day.
      for jj=1:var2
          PLs[jj,:] = Pl_s_year[k-BaH-jj*MpcH:k-BaH-jj*MpcH+MpcH-1];# I could have also use var1 to go
          # through the rows of PLs
          var1+=1;
      end
  end
  var3 = convert(Integer,floor(length(Pl_s_year[k-BaH+MpcH:end])/MpcH));# compute how many iterations
  #  is possible  forward  from  the end of the current day.
  var4=1;
  for kk=var1:var1+var3
      # begin the iteration where the last loop ended
      PLs[kk,:] = Pl_s_yearR[k-BaH+var4*MpcH:k-BaH+var4*MpcH+MpcH-1];
      var4+=1;
  end
  PLs[366,:] = PLs[365,:]; # fill the last row with the same data as the second to last row.
  
  if Eva =="Dist" #Using the distance
    PLs_cutErr= sum( (PLs[:,1:BaH] .- curr_day[1:BaH]').^2, dims=2)
    var7 = argmin(PLs_cutErr)[1] #find the index of the lowest element of PLs_cutErr which is the 
    # indice of the closest vector of curr_day[1:BaH] in PLs
    #println(var5[ii])
    var8 = convert(Integer,var7)
    Pls_frcst = vec(PLs')[MpcH*(var8-1)+BaH+1: MpcH*var8+BaH];
    return Pls_frcst,var7,PLs
  
  elseif Eva =="CorCoef" # Using the correlation coef
    PLs_cutErr2= Statistics.cor(curr_day[1:BaH]', PLs[:,1:BaH], dims=2);# Compute the corelation
    var7 = argmax(PLs_cutErr2)[2] #find the index of the highest element of PLs_cutErr which is the 
    # indice of the closest vector of curr_day[1:BaH] in PLs
    var8 = convert(Integer,var7)
    Pls_frcst2 = vec(PLs')[MpcH*(var8-1)+BaH+1: MpcH*var8+BaH];
    return Pls_frcst2,var7,PLs
    
  elseif Eva =="Both" #Using both previous
    #Using Distance------------------------------------------------
    PLs_cutErr= sum( (PLs[:,1:BaH] .- curr_day[1:BaH]').^2, dims=2)
    var7 = argmin(PLs_cutErr)[1]
    var8 = convert(Integer,var7)
    Pls_frcst = vec(PLs')[MpcH*(var8-1)+BaH+1: MpcH*var8+BaH];
    
    #Using Correlation coef-----------------------------------------
    PLs_cutErr2= Statistics.cor(curr_day[1:BaH]', PLs[:,1:BaH], dims=2)
    var7 = argmax(PLs_cutErr2)[2]
    var8 = convert(Integer,var7)
    Pls_frcst2 = vec(PLs')[MpcH*(var8-1)+BaH+1: MpcH*var8+BaH];
    return Pls_frcst, Pls_frcst2,var7,PLs
  else
      throw(UndefVarError(:Wrong_choice_of_input_Eva_MustBe_Mustbe_Dist_or_CorCoef_or_Both))
  end
end
    
  


  
end
