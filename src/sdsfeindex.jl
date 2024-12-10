
### get barely inefficiency  code 0 i.e. hsale * ut
### get  inefficiency  code 1 i.e. (I-wtau)^-1 * hsale * ut   if has
### get  inefficiency  code 2 i.e. (I-wrho)^-1 * (I-wtau)^-1 * hsale * ut,decompostion at (I-wrho)^-1 * (I-wtau)^-1 
#########################################
####        JLMS and BC index        ####
#########################################

#? --------------- Truncated Normal --------------




function  jlmsbckute0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, Wy::Matrix,
 PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})


 β  = rho[1:pos.endx]
τ  = rho[pos.begq:pos.endq]
phi = rho[pos.begphi:pos.endphi]

phi = reshape(phi,:,num.nofeta)
eps = EN-IV*phi
eta = rho[pos.begeta:pos.endeta]

δ2 = rho[pos.begw]  
γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
δ1 = rho[pos.begz]
gammap = rho[pos.beggamma]
gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));

hi  = exp.(Q*τ)

σᵤ²= exp(δ2) 
σᵤ= exp(0.5*δ2) 
σᵥ² = exp(γ)  
σᵥ = exp(0.5*γ)    
μ   = δ1
ϵ = PorC*(y - x * β )
T = size(rowIDT,1)

# sigs2 = zeros(eltype(y),T,1);
# mus = zeros(eltype(y),T,1);
# bc = zeros(eltype(y),size(hi,1),1);
# jlms = zeros(eltype(y),size(hi,1),1);

if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度

    @views N = rowIDT[1,2];
    @views Wyt = kron(I(T), Wy[1])


    @views invPi = 1.0 /σᵥ²;
    @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
    @views sigs2 = @. 1.0  / (hi^2 *invPi + 1.0 /σᵤ²) ;
    @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

    @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

elseif length(Wy)>1

    @views Wyt = BlockDiagonal([Wy...])
    @inbounds for ttt=1:T  
        @views N = rowIDT[ttt,2];

    end # for ttt=1:T

    @views invPi = 1.0/σᵥ²;
    @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
    @views sigs2 =@. 1.0 / (hi^2 *invPi + 1 /σᵤ²) ;
    @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

    @views jlms = @. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

end  #    if length(Wy)==1 

return jlms
end




function  jlmsbckut0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   Wy::Matrix,
    PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
   
    β  = rho[1:pos.endx]
   τ  = rho[pos.begq:pos.endq]
  
   δ2 = rho[pos.begw]  
   γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
   δ1 = rho[pos.begz]
   gammap = rho[pos.beggamma]
   gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
   
   hi  = exp.(Q*τ)

   σᵤ²= exp(δ2) 
   σᵤ= exp(0.5*δ2) 
   σᵥ² = exp(γ)  
   σᵥ = exp(0.5*γ)    
   μ   = δ1
   ϵ = PorC*(y - x * β )
   T = size(rowIDT,1)
   
   # sigs2 = zeros(eltype(y),T,1);
   # mus = zeros(eltype(y),T,1);
   # bc = zeros(eltype(y),size(hi,1),1);
   # jlms = zeros(eltype(y),size(hi,1),1);
   
   if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
   
       @views N = rowIDT[1,2];
       @views Wyt = kron(I(T), Wy[1])
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
 

       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 


   elseif length(Wy)>1
   
       @views Wyt = BlockDiagonal([Wy...])
       @inbounds for ttt=1:T  
        @views N = rowIDT[ttt,2];

       end # for ttt=1:T
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma.*Wyt*y  ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1 /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
  
   end  #    if length(Wy)==1 

   return jlms
   end
   
   




function jlmsbc0(::Type{SSFKUET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]
    jlms = jlmsbckute0(y, x, Q,  EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

    return jlms
 
 end

 function jlmsbc0(::Type{SSFKUT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]
    jlms = jlmsbckut0(y, x, Q,  Wy, PorC, pos, rho,  eigvalu, rowIDT )  
    
    return jlms
 
 end


 
function  jlmsbckuhe0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, Wy::Matrix,
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
   
    β  = rho[1:pos.endx]
   τ  = rho[pos.begq:pos.endq]
   phi = rho[pos.begphi:pos.endphi]
   
   phi = reshape(phi,:,num.nofeta)
   eps = EN-IV*phi
   eta = rho[pos.begeta:pos.endeta]
   
   δ2 = rho[pos.begw]  
   γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
#    δ1 = rho[pos.begz]
   gammap = rho[pos.beggamma]
   gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
   
   hi  = exp.(Q*τ)

   σᵤ²= exp(δ2) 
   σᵤ= exp(0.5*δ2) 
   σᵥ² = exp(γ)  
   σᵥ = exp(0.5*γ)    
   μ   = 0.0
   ϵ = PorC*(y - x * β )
   T = size(rowIDT,1)
   
   # sigs2 = zeros(eltype(y),T,1);
   # mus = zeros(eltype(y),T,1);
   # bc = zeros(eltype(y),size(hi,1),1);
   # jlms = zeros(eltype(y),size(hi,1),1);
   
   if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
   
       @views N = rowIDT[1,2];
       @views Wyt = kron(I(T), Wy[1])

       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

   
   elseif length(Wy)>1
   
       @views Wyt = BlockDiagonal([Wy...])
       @inbounds for ttt=1:T  
        @views N = rowIDT[ttt,2];

       end # for ttt=1:T
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

   end  #    if length(Wy)==1 

   return jlms
   end
   
   
   
   
   function  jlmsbckuh0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   Wy::Matrix,
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
      
       β  = rho[1:pos.endx]
      τ  = rho[pos.begq:pos.endq]
     
      δ2 = rho[pos.begw]  
      γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    #   δ1 = rho[pos.begz]
      gammap = rho[pos.beggamma]
      gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
      
      hi  = exp.(Q*τ)

      σᵤ²= exp(δ2) 
      σᵤ= exp(0.5*δ2) 
      σᵥ² = exp(γ)  
      σᵥ = exp(0.5*γ)    
      μ   = 0.0
      ϵ = PorC*(y - x * β )
      T = size(rowIDT,1)
      
      # sigs2 = zeros(eltype(y),T,1);
      # mus = zeros(eltype(y),T,1);
      # bc = zeros(eltype(y),size(hi,1),1);
      # jlms = zeros(eltype(y),size(hi,1),1);
      
      if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
      
          @views N = rowIDT[1,2];
          @views Wyt = kron(I(T), Wy[1])

      
          @views invPi = 1.0 /σᵥ²;
          @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
          @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
          @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

          @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

      
      elseif length(Wy)>1
      
        @views Wyt = BlockDiagonal([Wy...])
          @inbounds for ttt=1:T  
            @views N = rowIDT[ttt,2];

          end # for ttt=1:T
      
          @views invPi = 1.0 /σᵥ²;
          @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
          @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
          @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

          @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 


        end  #    if length(Wy)==1 

      return jlms   
      end
      
      
   
   
   
   
   function jlmsbc0(::Type{SSFKUEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
       PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
    
        Wy = _dicM[:wy]   
        jlms = jlmsbckuhe0(y, x, Q,  EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

       
        return jlms
    
    end



       
   function jlmsbc0(::Type{SSFKUH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]
    jlms = jlmsbckuh0(y, x, Q,  Wy, PorC, pos, rho,  eigvalu, rowIDT )  

     return jlms
 
 end



 
function  jlmsbckkte0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
    β  = rho[1:pos.endx]
    τ  = rho[pos.begq:pos.endq]
    phi = rho[pos.begphi:pos.endphi]
    
    phi = reshape(phi,:,num.nofeta)
    eps = EN-IV*phi
    eta = rho[pos.begeta:pos.endeta]
    
    δ2 = rho[pos.begw]  
    γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    δ1 = rho[pos.begz]
    
    hi  = exp.(Q*τ)
    σᵤ²= exp(δ2) 
    σᵤ= exp(0.5*δ2) 
    σᵥ² = exp(γ)  
    σᵥ = exp(0.5*γ)    
    μ   = δ1
    ϵ = PorC*(y - x * β )
    ID = size(rowIDT,1)
    
    sigs2 = zeros(eltype(y),ID,1);
    mus = zeros(eltype(y),ID,1);
    bc = zeros(eltype(y),size(hi,1),1);
    jlms = zeros(eltype(y),size(hi,1),1);
    
    @views invPi = 1.0 /σᵥ²  ;
    
    @floop begin
    @inbounds for idid=1:ID 

        @views ind = rowIDT[idid,1];
        @views ϵ[ind] = ϵ[ind] - PorC*(eps[ind,:]*eta) ;
        @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
        @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
        @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
    end # for idid=1:ID
    end # begin
   
    return jlms
    end
    
    
   
   
   
function  jlmsbckkt0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
       β  = rho[1:pos.endx]
       τ  = rho[pos.begq:pos.endq]
       
       δ2 = rho[pos.begw]  
       γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
       δ1 = rho[pos.begz]
       
       hi  = exp.(Q*τ)
       σᵤ²= exp(δ2) 
       σᵤ= exp(0.5*δ2) 
       σᵥ² = exp(γ)  
       σᵥ = exp(0.5*γ)    
       μ   = δ1
       ϵ = PorC*(y - x * β )
       ID = size(rowIDT,1)
       
       sigs2 = zeros(eltype(y),ID,1);
       mus = zeros(eltype(y),ID,1);
       bc = zeros(eltype(y),size(hi,1),1);
       jlms = zeros(eltype(y),size(hi,1),1);
       

       @views invPi = 1.0 /σᵥ²  ;
       @floop begin
       @inbounds for idid=1:ID 

           @views ind = rowIDT[idid,1];
           @views ϵ[ind] = ϵ[ind]  ;
           @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
           @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
           @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
       end # for idid=1:ID
       end # begin

       return jlms
    end    
      
      
    
function jlmsbc0(::Type{SSFKKET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms = jlmsbckkte0(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

    return jlms

end

function jlmsbc0(::Type{SSFKKT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms = jlmsbckkt0(y, x, Q, PorC, pos, rho,  eigvalu, rowIDT )  

    return jlms

end
   
   

 
function  jlmsbckkhe0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
    β  = rho[1:pos.endx]
    τ  = rho[pos.begq:pos.endq]
    phi = rho[pos.begphi:pos.endphi]
    
    phi = reshape(phi,:,num.nofeta)
    eps = EN-IV*phi
    eta = rho[pos.begeta:pos.endeta]
    
    δ2 = rho[pos.begw]  
    γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    # δ1 = rho[pos.begz]
    
    hi  = exp.(Q*τ)
    σᵤ²= exp(δ2) 
    σᵤ= exp(0.5*δ2) 
    σᵥ² = exp(γ)  
    σᵥ = exp(0.5*γ)    
    μ   = 0.0
    ϵ = PorC*(y - x * β )
    ID = size(rowIDT,1)
    
    sigs2 = zeros(eltype(y),ID,1);
    mus = zeros(eltype(y),ID,1);
    bc = zeros(eltype(y),size(hi,1),1);
    jlms = zeros(eltype(y),size(hi,1),1);
    
    @views invPi = 1.0 /σᵥ²  ;
    
    @floop begin
    @inbounds for idid=1:ID 

        @views ind = rowIDT[idid,1];
        @views ϵ[ind] = ϵ[ind] - PorC*(eps[ind,:]*eta) ;
        @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
        @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
        @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
    end # for idid=1:ID
    end # begin

    return jlms 
    end
    
    
   
   
   
function  jlmsbckkh0(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
       β  = rho[1:pos.endx]
       τ  = rho[pos.begq:pos.endq]
       
       δ2 = rho[pos.begw]  
       γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    #    δ1 = rho[pos.begz]
       
       hi  = exp.(Q*τ)
       σᵤ²= exp(δ2) 
       σᵤ= exp(0.5*δ2) 
       σᵥ² = exp(γ)  
       σᵥ = exp(0.5*γ)    
       μ   = 0.0
       ϵ = PorC*(y - x * β )
       ID = size(rowIDT,1)
       
       sigs2 = zeros(eltype(y),ID,1);
       mus = zeros(eltype(y),ID,1);
       bc = zeros(eltype(y),size(hi,1),1);
       jlms = zeros(eltype(y),size(hi,1),1);
       

       @views invPi = 1.0 /σᵥ²  ;

       @floop begin
       @inbounds for idid=1:ID 
           @views ind = rowIDT[idid,1];
           @views ϵ[ind] = ϵ[ind]  ;
           @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
           @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
           @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
       end # for idid=1:ID
       end # begin

       return jlms  
    end    
      
      
    
function jlmsbc0(::Type{SSFKKEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms = jlmsbckkhe0(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

    return jlms

end

function jlmsbc0(::Type{SSFKKH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms = jlmsbckkh0(y, x, Q, PorC, pos, rho,  eigvalu, rowIDT )  

    return jlms

end




### get barely inefficiency  code 0 i.e. hsale * ut
### get  inefficiency  code 1 i.e. (I-wtau)^-1 * hsale * ut   if has
### get  inefficiency  code 2 i.e. (I-wrho)^-1 * (I-wtau)^-1 * hsale * ut,decompostion at (I-wrho)^-1 * (I-wtau)^-1 
#########################################
####        JLMS and BC index        ####
#########################################

#? --------------- Truncated Normal --------------


function  jlmsbckute1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, Wy::Matrix,
 PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})


 β  = rho[1:pos.endx]
τ  = rho[pos.begq:pos.endq]
phi = rho[pos.begphi:pos.endphi]

phi = reshape(phi,:,num.nofeta)
eps = EN-IV*phi
eta = rho[pos.begeta:pos.endeta]

δ2 = rho[pos.begw]  
γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
δ1 = rho[pos.begz]
gammap = rho[pos.beggamma]
gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));

hi  = exp.(Q*τ)
σᵤ²= exp(δ2) 
σᵤ= exp(0.5*δ2) 
σᵥ² = exp(γ)  
σᵥ = exp(0.5*γ)    
μ   = δ1
ϵ = PorC*(y - x * β )
T = size(rowIDT,1)

# sigs2 = zeros(eltype(y),T,1);
# mus = zeros(eltype(y),T,1);
# bc = zeros(eltype(y),size(hi,1),1);
# jlms = zeros(eltype(y),size(hi,1),1);

if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度

    @views N = rowIDT[1,2];
    @views Wyt = kron(I(T), Wy[1])


    @views invPi = 1.0 /σᵥ²;
    @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
    @views sigs2 = @. 1.0  / (hi^2 *invPi + 1.0 /σᵤ²) ;
    @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

    @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

elseif length(Wy)>1

    @views Wyt = BlockDiagonal([Wy...])


    @views invPi = 1.0/σᵥ²;
    @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
    @views sigs2 =@. 1.0 / (hi^2 *invPi + 1 /σᵤ²) ;
    @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

    @views jlms = @. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

end  #    if length(Wy)==1 
@views TE_bc = exp.(-bc);
@views TE_jlms = exp.(-jlms);
return TE_jlms, TE_bc     
end




function  jlmsbckut1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   Wy::Matrix,
    PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
   
    β  = rho[1:pos.endx]
   τ  = rho[pos.begq:pos.endq]
  
   δ2 = rho[pos.begw]  
   γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
   δ1 = rho[pos.begz]
   gammap = rho[pos.beggamma]
   gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
   
   hi  = exp.(Q*τ)
   σᵤ²= exp(δ2) 
   σᵤ= exp(0.5*δ2) 
   σᵥ² = exp(γ)  
   σᵥ = exp(0.5*γ)    
   μ   = δ1
   ϵ = PorC*(y - x * β )
   T = size(rowIDT,1)
   
   # sigs2 = zeros(eltype(y),T,1);
   # mus = zeros(eltype(y),T,1);
   # bc = zeros(eltype(y),size(hi,1),1);
   # jlms = zeros(eltype(y),size(hi,1),1);
   
   if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
   
       @views N = rowIDT[1,2];
       @views Wyt = kron(I(T), Wy[1])
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 


   elseif length(Wy)>1
   
       @views Wyt = BlockDiagonal([Wy...])
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma.*Wyt*y  ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1 /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
   end  #    if length(Wy)==1 

   return jlms
   end
   
   




function jlmsbc1(::Type{SSFKUET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]
    jlms = jlmsbckute1(y, x, Q,  EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

     return jlms
 
 end

 function jlmsbc1(::Type{SSFKUT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]
    jlms = jlmsbckut1(y, x, Q,  Wy, PorC, pos, rho,  eigvalu, rowIDT )  

     return jlms
 
 end


 
function  jlmsbckuhe1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, Wy::Matrix,
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
   
    β  = rho[1:pos.endx]
   τ  = rho[pos.begq:pos.endq]
   phi = rho[pos.begphi:pos.endphi]
   
   phi = reshape(phi,:,num.nofeta)
   eps = EN-IV*phi
   eta = rho[pos.begeta:pos.endeta]
   
   δ2 = rho[pos.begw]  
   γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
#    δ1 = rho[pos.begz]
   gammap = rho[pos.beggamma]
   gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
   
   hi  = exp.(Q*τ)
   σᵤ²= exp(δ2) 
   σᵤ= exp(0.5*δ2) 
   σᵥ² = exp(γ)  
   σᵥ = exp(0.5*γ)    
   μ   = 0.0
   ϵ = PorC*(y - x * β )
   T = size(rowIDT,1)
   
   # sigs2 = zeros(eltype(y),T,1);
   # mus = zeros(eltype(y),T,1);
   # bc = zeros(eltype(y),size(hi,1),1);
   # jlms = zeros(eltype(y),size(hi,1),1);
   
   if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
   
       @views N = rowIDT[1,2];
       @views Wyt = kron(I(T), Wy[1])

       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

   
   elseif length(Wy)>1
   
       @views Wyt = BlockDiagonal([Wy...])
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

       @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
   end  #    if length(Wy)==1 

   return jlms
   end
   
   
   
   
   function  jlmsbckuh1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   Wy::Matrix,
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
      
       β  = rho[1:pos.endx]
      τ  = rho[pos.begq:pos.endq]
     
      δ2 = rho[pos.begw]  
      γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    #   δ1 = rho[pos.begz]
      gammap = rho[pos.beggamma]
      gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
      
      hi  = exp.(Q*τ)
      σᵤ²= exp(δ2) 
      σᵤ= exp(0.5*δ2) 
      σᵥ² = exp(γ)  
      σᵥ = exp(0.5*γ)    
      μ   = 0.0
      ϵ = PorC*(y - x * β )
      T = size(rowIDT,1)
      
      # sigs2 = zeros(eltype(y),T,1);
      # mus = zeros(eltype(y),T,1);
      # bc = zeros(eltype(y),size(hi,1),1);
      # jlms = zeros(eltype(y),size(hi,1),1);
      
      if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
      
          @views N = rowIDT[1,2];
          @views Wyt = kron(I(T), Wy[1])

      
          @views invPi = 1.0 /σᵥ²;
          @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
          @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
          @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

          @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
      
      elseif length(Wy)>1
      
        @views Wyt = BlockDiagonal([Wy...])

          @views invPi = 1.0 /σᵥ²;
          @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
          @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
          @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
          @views jlms =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 

        end  #    if length(Wy)==1 

      return jlms
      end
      
      
   
   
   
   
   function jlmsbc1(::Type{SSFKUEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
       PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
    
       Wy = _dicM[:wy]
        jlms = jlmsbckuhe1(y, x, Q,  EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

        return jlms
    
    end



       
   function jlmsbc1(::Type{SSFKUH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]
    jlms = jlmsbckuh1(y, x, Q,  Wy, PorC, pos, rho,  eigvalu, rowIDT )  

     return jlms
 
 end



 
function  jlmsbckkte1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
    β  = rho[1:pos.endx]
    τ  = rho[pos.begq:pos.endq]
    phi = rho[pos.begphi:pos.endphi]
    
    phi = reshape(phi,:,num.nofeta)
    eps = EN-IV*phi
    eta = rho[pos.begeta:pos.endeta]
    
    δ2 = rho[pos.begw]  
    γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    δ1 = rho[pos.begz]
    
    hi  = exp.(Q*τ)
    σᵤ²= exp(δ2) 
    σᵤ= exp(0.5*δ2) 
    σᵥ² = exp(γ)  
    σᵥ = exp(0.5*γ)    
    μ   = δ1
    ϵ = PorC*(y - x * β )
    ID = size(rowIDT,1)
    
    sigs2 = zeros(eltype(y),ID,1);
    mus = zeros(eltype(y),ID,1);
    bc = zeros(eltype(y),size(hi,1),1);
    jlms = zeros(eltype(y),size(hi,1),1);
    
    @views invPi = 1.0 /σᵥ²  ;
    
    @floop begin
    @inbounds for idid=1:ID 

        @views ind = rowIDT[idid,1];
        @views ϵ[ind] = ϵ[ind] - PorC*(eps[ind,:]*eta) ;
        @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
        @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
        @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
    end # for idid=1:ID
    end # begin
   

    return jlms
    end
    
    
   
   
   
function  jlmsbckkt1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
       β  = rho[1:pos.endx]
       τ  = rho[pos.begq:pos.endq]
       
       δ2 = rho[pos.begw]  
       γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
       δ1 = rho[pos.begz]
       
       hi  = exp.(Q*τ)
       σᵤ²= exp(δ2) 
       σᵤ= exp(0.5*δ2) 
       σᵥ² = exp(γ)  
       σᵥ = exp(0.5*γ)    
       μ   = δ1
       ϵ = PorC*(y - x * β )
       ID = size(rowIDT,1)
       
       sigs2 = zeros(eltype(y),ID,1);
       mus = zeros(eltype(y),ID,1);
       bc = zeros(eltype(y),size(hi,1),1);
       jlms = zeros(eltype(y),size(hi,1),1);
       

       @views invPi = 1.0 /σᵥ²  ;
       @floop begin
       @inbounds for idid=1:ID 

           @views ind = rowIDT[idid,1];
           @views ϵ[ind] = ϵ[ind]  ;
           @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
           @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
           @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
       end # for idid=1:ID
       end # begin

       return jlms
    end    
      
      
    
function jlmsbc1(::Type{SSFKKET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms = jlmsbckkte1(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

    return jlms

end

function jlmsbc1(::Type{SSFKKT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms = jlmsbckkt1(y, x, Q, PorC, pos, rho,  eigvalu, rowIDT )  

    return jlms

end
   
   

 
function  jlmsbckkhe1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
    β  = rho[1:pos.endx]
    τ  = rho[pos.begq:pos.endq]
    phi = rho[pos.begphi:pos.endphi]
    
    phi = reshape(phi,:,num.nofeta)
    eps = EN-IV*phi
    eta = rho[pos.begeta:pos.endeta]
    
    δ2 = rho[pos.begw]  
    γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    # δ1 = rho[pos.begz]
    
    hi  = exp.(Q*τ)
    σᵤ²= exp(δ2) 
    σᵤ= exp(0.5*δ2) 
    σᵥ² = exp(γ)  
    σᵥ = exp(0.5*γ)    
    μ   = 0.0
    ϵ = PorC*(y - x * β )
    ID = size(rowIDT,1)
    
    sigs2 = zeros(eltype(y),ID,1);
    mus = zeros(eltype(y),ID,1);
    bc = zeros(eltype(y),size(hi,1),1);
    jlms = zeros(eltype(y),size(hi,1),1);
    
    @views invPi = 1.0 /σᵥ²  ;
    
    @floop begin
    @inbounds for idid=1:ID 

        @views ind = rowIDT[idid,1];
        @views ϵ[ind] = ϵ[ind] - PorC*(eps[ind,:]*eta) ;
        @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
        @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
        @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
    end # for idid=1:ID
    end # begin
   
    return jlms  
    end
    
    
   
   
   
function  jlmsbckkh1(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
       β  = rho[1:pos.endx]
       τ  = rho[pos.begq:pos.endq]
       
       δ2 = rho[pos.begw]  
       γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    #    δ1 = rho[pos.begz]
       
       hi  = exp.(Q*τ)
       σᵤ²= exp(δ2) 
       σᵤ= exp(0.5*δ2) 
       σᵥ² = exp(γ)  
       σᵥ = exp(0.5*γ)    
       μ   = 0.0
       ϵ = PorC*(y - x * β )
       ID = size(rowIDT,1)
       
       sigs2 = zeros(eltype(y),ID,1);
       mus = zeros(eltype(y),ID,1);
       bc = zeros(eltype(y),size(hi,1),1);
       jlms = zeros(eltype(y),size(hi,1),1);
       

       @views invPi = 1.0 /σᵥ²  ;

       @floop begin
       @inbounds for idid=1:ID 
           @views ind = rowIDT[idid,1];
           @views ϵ[ind] = ϵ[ind]  ;
           @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
           @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
           @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
       end # for idid=1:ID
       end # begin

       return jlms
    end    
      
      
    
function jlmsbc1(::Type{SSFKKEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms_, bc_ = jlmsbckkhe1(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

    return jlms_, bc_

end

function jlmsbc1(::Type{SSFKKH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms = jlmsbckkh1(y, x, Q, PorC, pos, rho,  eigvalu, rowIDT )  

    return jlms

end







### get barely inefficiency  code 0 i.e. hsale * ut
### get  inefficiency  code 1 i.e. (I-wtau)^-1 * hsale * ut   if has
### get  inefficiency  code 2 i.e. (I-wrho)^-1 * (I-wtau)^-1 * hsale * ut,decompostion at (I-wrho)^-1 * (I-wtau)^-1 
#########################################
####        JLMS and BC index        ####
#########################################

#? --------------- Truncated Normal --------------




function  jlmsbckute(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, Wy::Matrix,
 PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})


 β  = rho[1:pos.endx]
τ  = rho[pos.begq:pos.endq]
phi = rho[pos.begphi:pos.endphi]

phi = reshape(phi,:,num.nofeta)
eps = EN-IV*phi
eta = rho[pos.begeta:pos.endeta]

δ2 = rho[pos.begw]  
γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
δ1 = rho[pos.begz]
gammap = rho[pos.beggamma]
gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));

hi  = exp.(Q*τ)

σᵤ²= exp(δ2) 
σᵤ= exp(0.5*δ2) 
σᵥ² = exp(γ)  
σᵥ = exp(0.5*γ)    
μ   = δ1
ϵ = PorC*(y - x * β )
T = size(rowIDT,1)

# sigs2 = zeros(eltype(y),T,1);
# mus = zeros(eltype(y),T,1);
# bc = zeros(eltype(y),size(hi,1),1);
# jlms = zeros(eltype(y),size(hi,1),1);

if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度

    @views N = rowIDT[1,2];
    @views Wyt = kron(I(T), Wy[1])
    @views Mgamma = (I(N)-gamma*Wy[1])\I(N)
    @views Mgammat = kron(I(T), Mgamma)

    @views invPi = 1.0 /σᵥ²;
    @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
    @views sigs2 = @. 1.0  / (hi^2 *invPi + 1.0 /σᵤ²) ;
    @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

    @views jlms1 =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
    @views jlms = Mgammat*jlms1;  
    @views jlms_direct = Diagonal(diag(Mgammat))*jlms1;  
    @views jlms_indirect = jlms - jlms_direct;

elseif length(Wy)>1

    @views Wyt = BlockDiagonal([Wy...])
    @views Mgammat_ = Array{Matrix}(undef, T);
    @inbounds for ttt=1:T  
        @views N = rowIDT[ttt,2];

        @views Mgammat_[ttt] = (I(N)-gamma*Wy[ttt])\I(N);
    end # for ttt=1:T
    @views Mgammat = BlockDiagonal([Mgammat_...])

    @views invPi = 1.0/σᵥ²;
    @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
    @views sigs2 =@. 1.0 / (hi^2 *invPi + 1 /σᵤ²) ;
    @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
    @views jlms1 = @. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
    @views jlms = Mgammat*jlms1;  
    @views jlms_direct = Diagonal(diag(Mgammat))*jlms1;  
    @views jlms_indirect = jlms - jlms_direct;
end  #    if length(Wy)==1 

return jlms,jlms_direct,jlms_indirect
end




function  jlmsbckut(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   Wy::Matrix,
    PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
   
    β  = rho[1:pos.endx]
   τ  = rho[pos.begq:pos.endq]
  
   δ2 = rho[pos.begw]  
   γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
   δ1 = rho[pos.begz]
   gammap = rho[pos.beggamma]
   gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
   
   hi  = exp.(Q*τ)
   σᵤ²= exp(δ2) 
   σᵤ= exp(0.5*δ2) 
   σᵥ² = exp(γ)  
   σᵥ = exp(0.5*γ)    
   μ   = δ1
   ϵ = PorC*(y - x * β )
   T = size(rowIDT,1)
   
   # sigs2 = zeros(eltype(y),T,1);
   # mus = zeros(eltype(y),T,1);
   # bc = zeros(eltype(y),size(hi,1),1);
   # jlms = zeros(eltype(y),size(hi,1),1);
   
   if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
   
       @views N = rowIDT[1,2];
       @views Wyt = kron(I(T), Wy[1])
       @views Mgamma = (I(N)-gamma*Wy[1])\I(N)
       @views Mgammat = kron(I(T), Mgamma)
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
 
       @views jlms1 =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
       @views jlms = Mgammat*jlms1;  
       @views jlms_direct = Diagonal(diag(Mgammat))*jlms1;  
       @views jlms_indirect = jlms - jlms_direct;
   elseif length(Wy)>1
   
       @views Wyt = BlockDiagonal([Wy...])
       @views Mgammat_ = Array{Matrix}(undef, T);
       @inbounds for ttt=1:T  
        @views N = rowIDT[ttt,2];

           @views Mgammat_[ttt] = (I(N)-gamma*Wy[ttt])\I(N);
       end # for ttt=1:T
       @views Mgammat = BlockDiagonal([Mgammat_...])
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma.*Wyt*y  ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1 /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

       @views jlms1 =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
       @views jlms =Mgammat* jlms1;  
       @views jlms_direct = Diagonal(diag(Mgammat))*jlms1;  
       @views jlms_indirect = jlms - jlms_direct;
   end  #    if length(Wy)==1 

   return jlms,jlms_direct,jlms_indirect
end
   
   




function jlmsbc(::Type{SSFKUET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]

    jlms,jlms_direct,jlms_indirect = jlmsbckute(y, x, Q,  EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

     return jlms,jlms_direct,jlms_indirect
 
 end

 function jlmsbc(::Type{SSFKUT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]

    jlms,jlms_direct,jlms_indirect = jlmsbckut(y, x, Q,  Wy, PorC, pos, rho,  eigvalu, rowIDT )  
    
     return jlms,jlms_direct,jlms_indirect
 
 end


 
function  jlmsbckuhe(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, Wy::Matrix,
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
   
    β  = rho[1:pos.endx]
   τ  = rho[pos.begq:pos.endq]
   phi = rho[pos.begphi:pos.endphi]
   
   phi = reshape(phi,:,num.nofeta)
   eps = EN-IV*phi
   eta = rho[pos.begeta:pos.endeta]
   
   δ2 = rho[pos.begw]  
   γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
#    δ1 = rho[pos.begz]
   gammap = rho[pos.beggamma]
   gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
   
   hi  = exp.(Q*τ)
   σᵤ²= exp(δ2) 
   σᵤ= exp(0.5*δ2) 
   σᵥ² = exp(γ)  
   σᵥ = exp(0.5*γ)    
   μ   = 0.0
   ϵ = PorC*(y - x * β )
   T = size(rowIDT,1)
   println("yy",y)
   # sigs2 = zeros(eltype(y),T,1);
   # mus = zeros(eltype(y),T,1);
   # bc = zeros(eltype(y),size(hi,1),1);
   # jlms = zeros(eltype(y),size(hi,1),1);
   
   if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
   
       @views N = rowIDT[1,2];
       @views Wyt = kron(I(T), Wy[1])
       @views Mgamma = (I(N)-gamma*Wy[1])\I(N)
       @views Mgammat = kron(I(T), Mgamma)
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

       @views jlms1 =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
    #    println(jlms1)
       @views jlms = Mgammat*jlms1;  
    #    println(jlms)

       @views jlms_direct = (diag(Mgammat)).*jlms1;  
       @views jlms_indirect = jlms - jlms_direct;
   elseif length(Wy)>1
   
       @views Wyt = BlockDiagonal([Wy...])
       @views Mgammat_ = Array{Matrix}(undef, T);
       @inbounds for ttt=1:T  
        @views N = rowIDT[ttt,2];

           @views Mgammat_[ttt] = (I(N)-gamma*Wy[ttt])\I(N);
       end # for ttt=1:T
       @views Mgammat = BlockDiagonal([Mgammat_...])
   
       @views invPi = 1.0 /σᵥ²;
       @views ϵ  = ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
       @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
       @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

       @views jlms1 =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
       @views jlms = Mgammat*jlms1;  
       @views jlms_direct = (diag(Mgammat)).*jlms1;  
       @views jlms_indirect = jlms - jlms_direct;
   end  #    if length(Wy)==1 

   return jlms,jlms_direct,jlms_indirect    
   end
   
   
   
   
   function  jlmsbckuh(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   Wy::Matrix,
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
      
       β  = rho[1:pos.endx]
      τ  = rho[pos.begq:pos.endq]
     
      δ2 = rho[pos.begw]  
      γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    #   δ1 = rho[pos.begz]
      gammap = rho[pos.beggamma]
      gamma  = eigvalu.rymin/(1.0 +exp(gammap))+eigvalu.rymax*exp(gammap)/(1.0 +exp(gammap));
      
      hi  = exp.(Q*τ)
      σᵤ²= exp(δ2) 
      σᵤ= exp(0.5*δ2) 
      σᵥ² = exp(γ)  
      σᵥ = exp(0.5*γ)    
      μ   = 0.0
      ϵ = PorC*(y - x * β )
      T = size(rowIDT,1)
      
      # sigs2 = zeros(eltype(y),T,1);
      # mus = zeros(eltype(y),T,1);
      # bc = zeros(eltype(y),size(hi,1),1);
      # jlms = zeros(eltype(y),size(hi,1),1);
      
      if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度
      
          @views N = rowIDT[1,2];
          @views Wyt = kron(I(T), Wy[1])
          @views Mgamma = (I(N)-gamma*Wy[1])\I(N)
          @views Mgammat = kron(I(T), Mgamma)
      
          @views invPi = 1.0 /σᵥ²;
          @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
          @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
          @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;

          @views jlms1 =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
          @views jlms = Mgammat*jlms1;  
          @views jlms_direct = (diag(Mgammat)).*jlms1;  
          @views jlms_indirect = jlms - jlms_direct;
      elseif length(Wy)>1
      
        @views Wyt = BlockDiagonal([Wy...])
          @views Mgammat_ = Array{Matrix}(undef, T);
          @inbounds for ttt=1:T  
            @views N = rowIDT[ttt,2];

              @views Mgammat_[ttt] = (I(N)-gamma*Wy[ttt])\I(N);
          end # for ttt=1:T
          @views Mgammat = BlockDiagonal([Mgammat_...])
      
          @views invPi = 1.0 /σᵥ²;
          @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
          @views sigs2 =@. 1.0  / (hi^2 *invPi + 1.0  /σᵤ²) ;
          @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
    
          @views jlms1 =@. hi *( mus + sqrt(sigs2) * normpdf(mus/sqrt(sigs2))/normcdf(mus/sqrt(sigs2))   ) 
          @views jlms = Mgammat*jlms1;  
          @views jlms_direct = (diag(Mgammat)).*jlms1;  
          @views jlms_indirect = jlms - jlms_direct;
        end  #    if length(Wy)==1 

      return jlms,jlms_direct,jlms_indirect
      end
      
      
   
   
   
   
   function jlmsbc(::Type{SSFKUEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
       PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
    
       Wy = _dicM[:wy]
       jlms,jlms_direct,jlms_indirect = jlmsbckuhe(y, x, Q,  EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

        return jlms,jlms_direct,jlms_indirect
    
    end



       
   function jlmsbc(::Type{SSFKUH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
 
    Wy = _dicM[:wy]
    jlms,jlms_direct,jlms_indirect = jlmsbckuh(y, x, Q,  Wy, PorC, pos, rho,  eigvalu, rowIDT )  

     return jlms,jlms_direct,jlms_indirect
 
 end



 
function  jlmsbckkte(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
    β  = rho[1:pos.endx]
    τ  = rho[pos.begq:pos.endq]
    phi = rho[pos.begphi:pos.endphi]
    
    phi = reshape(phi,:,num.nofeta)
    eps = EN-IV*phi
    eta = rho[pos.begeta:pos.endeta]
    
    δ2 = rho[pos.begw]  
    γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    δ1 = rho[pos.begz]
    
    hi  = exp.(Q*τ)
    σᵤ²= exp(δ2) 
    σᵤ= exp(0.5*δ2) 
    σᵥ² = exp(γ)  
    σᵥ = exp(0.5*γ)    
    μ   = δ1
    ϵ = PorC*(y - x * β )
    ID = size(rowIDT,1)
    
    sigs2 = zeros(eltype(y),ID,1);
    mus = zeros(eltype(y),ID,1);

    jlms = zeros(eltype(y),size(hi,1),1);
    jlms_direct = zeros(eltype(y),size(hi,1),1);
    jlms_indirect = zeros(eltype(y),size(hi,1),1);

    @views invPi = 1.0 /σᵥ²  ;
    
    @floop begin
    @inbounds for idid=1:ID 

        @views ind = rowIDT[idid,1];
        @views ϵ[ind] = ϵ[ind] - PorC*(eps[ind,:]*eta) ;
        @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
        @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
        @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
        @views jlms_direct = jlms
        @views jlms_indirect = jlms - jlms_direct;
    end # for idid=1:ID
    end # begin
   
    return jlms,jlms_direct,jlms_indirect
    end
    
    
   
   
   
function  jlmsbckkt(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
       β  = rho[1:pos.endx]
       τ  = rho[pos.begq:pos.endq]
       
       δ2 = rho[pos.begw]  
       γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
       δ1 = rho[pos.begz]
       
       hi  = exp.(Q*τ)
       σᵤ²= exp(δ2) 
       σᵤ= exp(0.5*δ2) 
       σᵥ² = exp(γ)  
       σᵥ = exp(0.5*γ)    
       μ   = δ1
       ϵ = PorC*(y - x * β )
       ID = size(rowIDT,1)
       
       sigs2 = zeros(eltype(y),ID,1);
       mus = zeros(eltype(y),ID,1);

       jlms = zeros(eltype(y),size(hi,1),1);
       jlms_direct = zeros(eltype(y),size(hi,1),1);
       jlms_indirect = zeros(eltype(y),size(hi,1),1);

       @views invPi = 1.0 /σᵥ²  ;
       @floop begin
       @inbounds for idid=1:ID 

           @views ind = rowIDT[idid,1];
           @views ϵ[ind] = ϵ[ind]  ;
           @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
           @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
           @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
           @views jlms_direct = jlms
           @views jlms_indirect = jlms - jlms_direct;
        end # for idid=1:ID
       end # begin
      

       return jlms,jlms_direct,jlms_indirect
    end    
      
      
    
function jlmsbc(::Type{SSFKKET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms,jlms_direct,jlms_indirect = jlmsbckkte(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

    return jlms,jlms_direct,jlms_indirect

end

function jlmsbc(::Type{SSFKKT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms,jlms_direct,jlms_indirect = jlmsbckkt(y, x, Q, PorC, pos, rho,  eigvalu, rowIDT )  

    return jlms,jlms_direct,jlms_indirect

end
   
   

 
function  jlmsbckkhe(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
    PorC::Int64, num::NamedTuple, pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
   
    β  = rho[1:pos.endx]
    τ  = rho[pos.begq:pos.endq]
    phi = rho[pos.begphi:pos.endphi]
    
    phi = reshape(phi,:,num.nofeta)
    eps = EN-IV*phi
    eta = rho[pos.begeta:pos.endeta]
    
    δ2 = rho[pos.begw]  
    γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    # δ1 = rho[pos.begz]
    
    hi  = exp.(Q*τ)
    σᵤ²= exp(δ2) 
    σᵤ= exp(0.5*δ2) 
    σᵥ² = exp(γ)  
    σᵥ = exp(0.5*γ)    
    μ   = 0.0
    ϵ = PorC*(y - x * β )
    ID = size(rowIDT,1)
    
    sigs2 = zeros(eltype(y),ID,1);
    mus = zeros(eltype(y),ID,1);

    jlms = zeros(eltype(y),size(hi,1),1);
    jlms_direct = zeros(eltype(y),size(hi,1),1);
    jlms_indirect = zeros(eltype(y),size(hi,1),1);

    @views invPi = 1.0 /σᵥ²  ;
    
    @floop begin
    @inbounds for idid=1:ID 

        @views ind = rowIDT[idid,1];
        @views ϵ[ind] = ϵ[ind] - PorC*(eps[ind,:]*eta) ;
        @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
        @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
        @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
        @views jlms_direct = jlms
        @views jlms_indirect = jlms - jlms_direct;
    end # for idid=1:ID
    end # begin
   

    return jlms,jlms_direct,jlms_indirect
    end
    
    
   
   
   
function  jlmsbckkh(y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,   
       PorC::Int64,  pos::NamedTuple, rho::Array{Float64, 1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})
      
       β  = rho[1:pos.endx]
       τ  = rho[pos.begq:pos.endq]
       
       δ2 = rho[pos.begw]  
       γ  = rho[pos.begv]  # May rho[po.begw : po.endw][1]
    #    δ1 = rho[pos.begz]
       
       hi  = exp.(Q*τ)
       σᵤ²= exp(δ2) 
       σᵤ= exp(0.5*δ2) 
       σᵥ² = exp(γ)  
       σᵥ = exp(0.5*γ)    
       μ   = 0.0
       ϵ = PorC*(y - x * β )
       ID = size(rowIDT,1)
       
       sigs2 = zeros(eltype(y),ID,1);
       mus = zeros(eltype(y),ID,1);

       jlms = zeros(eltype(y),size(hi,1),1);
       jlms_direct = zeros(eltype(y),size(hi,1),1);
       jlms_indirect = zeros(eltype(y),size(hi,1),1);


       @views invPi = 1.0 /σᵥ²  ;

       @floop begin
       @inbounds for idid=1:ID 
           @views ind = rowIDT[idid,1];
           @views ϵ[ind] = ϵ[ind]  ;
           @views sigs2[idid] = 1.0 /(hi[ind]'*hi[ind]*invPi  +1.0/σᵤ²);
           @views mus[idid] = (μ/σᵤ² - ϵ[ind]'*hi[ind] *invPi )*sigs2[idid] ;
           @views jlms[ind] = hi[ind] .*( mus[idid] + sqrt(sigs2[idid])* normpdf(mus[idid]/sqrt(sigs2[idid]))./normcdf(mus[idid]/sqrt(sigs2[idid]))   ) ;  
           @views jlms_direct = jlms
           @views jlms_indirect = jlms - jlms_direct;
        end # for idid=1:ID
       end # begin

       return jlms,jlms_direct,jlms_indirect
    end    
      
      
    
function jlmsbc(::Type{SSFKKEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms,jlms_direct,jlms_indirect = jlmsbckkhe(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

    return jlms,jlms_direct,jlms_indirect

end

function jlmsbc(::Type{SSFKKH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN, IV,
    PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

    jlms,jlms_direct,jlms_indirect = jlmsbckkh(y, x, Q, PorC, pos, rho,  eigvalu, rowIDT )  

    return jlms,jlms_direct,jlms_indirect

end