
#########################################
####        JLMS and BC index        ####
#########################################

#? --------------- Truncated Normal --------------






function prtlloglikekute( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix,
    Wy::Matrix, PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]
    println("tho",rho)
    phi = rho[po.begphi:po.endphi]
    phi = reshape(phi,:,num.nofeta)
    eps = zeros(eltype(EN),num.nofobs,num.nofeta);
    eps = EN-IV*phi
    eta = rho[po.begeta:po.endeta]
        ## calculate lky
        
        δ2 = rho[po.begw]  
        γ  = rho[po.begv]  # May rho[po.begw : po.endw][1]
        δ1 = rho[po.begz]
        gammap = rho[po.beggamma]
        gamma  = (eigvalu.rymin)/(1.0 +exp(gammap))+(eigvalu.rymax)*exp(gammap)/(1.0 +exp(gammap));
    
        hi  = exp.(Q*τ)
        σᵤ²= exp(δ2) 
        σᵤ= exp(0.5*δ2) 
        σᵥ² = exp(γ)            # todo: 重新换一下字母 
        σᵥ = exp(0.5*γ)  
        μ   = δ1
        ϵ = PorC*(y - x*β)
        T = size(rowIDT,1)
    

    if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度

        lik = zero(eltype(y));
        @views N = rowIDT[1,2];
        Wyt = kron(I(T), Wy[1])

        @views lndetIrhoW = log(det(I(N)-gamma*Wy[1]));  
        @views lndetIrhoWt = lndetIrhoW*T
        @views invPi = 1.0/σᵥ²;
        @views lndetPi = log(σᵥ²);
    
        @views ϵ = ϵ - PorC*gamma*Wyt*y - PorC*(eps*eta) ;
        @views sigs2 = @.  1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 =@. -0.5 * ϵ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
 
        # 检查 lik 是否为 NaN, 非实数, 或 Inf
        @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
        # 计算总和
        @views lik = sum(temp_1)
        @views lll =  lik+lndetIrhoWt

    elseif length(Wy)>1
        lik = zero(eltype(y));
        Wyt = kron(I(T), Wy[1])

        @floop begin
        @inbounds for ttt=1:T
          @views N = rowIDT[ttt,2];

            @views lndetIrhoW = log(det(I(N)-gamma*Wy[ttt]));    
            @views lndetIrhoWt += lndetIrhoW
        end # for ttt=1:T
        end # begin
        
        @views invPi = 1.0/σᵥ²;
        @views lndetPi = log(σᵥ²);
        @views ϵ  =  ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
        @views sigs2 = @.  1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 =@. -0.5 * ϵ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))

        # 检查 lik 是否为 NaN, 非实数, 或 Inf
        @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
        # 计算总和
        @views lik = sum(temp_1)
        @views lll =  lik+lndetIrhoWt

    end # if length(Wy)==1 

        return -lll
  
end

  

  
function prtlloglike(::Type{SSFKUET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN::Matrix, IV::Matrix,
  PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

  Wy = _dicM[:wy]
  println("tho1",rho)
  liky = prtlloglikekute(y, x, Q, EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

  return liky;       
end


function prtlloglikekuhe( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
  Wy::Matrix, PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
  β  = rho[1:po.endx]
  τ  = rho[po.begq:po.endq]
  phi = rho[po.begphi:po.endphi]
  phi = reshape(phi,:,num.nofeta)
  eps = EN-IV*phi
  eta = rho[po.begeta:po.endeta]
      ## calculate lky
      
      δ2 = rho[po.begw]  
      γ  = rho[po.begv]  # May rho[po.begw : po.endw][1]
      # δ1 = rho[po.begz]
      gammap = rho[po.beggamma]
      gamma  = (eigvalu.rymin)/(1.0 +exp(gammap))+(eigvalu.rymax)*exp(gammap)/(1.0 +exp(gammap));
  
      hi  = exp.(Q*τ)
      σᵤ²= exp(δ2) 
      σᵤ= exp(0.5*δ2) 
      σᵥ² = exp(γ)            # todo: 重新换一下字母 
      σᵥ = exp(0.5*γ)  
      μ   = 0.0
      ϵ = PorC*(y - x*β)
      T = size(rowIDT,1)
  

  if length(Wy)==1  # 可以传入单个cell的w，则默认cell的长度为时间的长度

      lik = zero(eltype(y));
      @views N = rowIDT[1,2];
      Wyt = kron(I(T), Wy[1])

      @views lndetIrhoW = log(det(I(N)-gamma*Wy[1]));  
      @views lndetIrhoWt = lndetIrhoW*T
      @views invPi = 1.0/σᵥ²;
      @views lndetPi = log(σᵥ²);
  
      @views ϵ = ϵ - PorC*gamma*Wyt*y - PorC*(eps*eta) ;
      @views sigs2 = @.  1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
      @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
      @views es2 =@. -0.5 * ϵ^2 *invPi;
      @views KK = -0.5*log(2 * π)-0.5*lndetPi;
      @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                      0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                      0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))

      # 检查 lik 是否为 NaN, 非实数, 或 Inf
      @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
      # 计算总和
      @views lik = sum(temp_1)
      @views lll =  lik+lndetIrhoWt

  elseif length(Wy)>1
      lik = zero(eltype(y));
      Wyt = kron(I(T), Wy[1])

      @floop begin
      @inbounds for ttt=1:T
        @views N = rowIDT[ttt,2];

          @views lndetIrhoW = log(det(I(N)-gamma*Wy[ttt]));    
          @views lndetIrhoWt += lndetIrhoW
      end # for ttt=1:T
      end # begin
      
      @views invPi = 1.0/σᵥ²;
      @views lndetPi = log(σᵥ²);
      @views ϵ  =  ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
      @views sigs2 = @.  1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
      @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
      @views es2 =@. -0.5 * ϵ^2 *invPi;
      @views KK = -0.5*log(2 * π)-0.5*lndetPi;
      @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                      0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                      0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))

      # 检查 lik 是否为 NaN, 非实数, 或 Inf
      @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
      # 计算总和
      @views lik = sum(temp_1)
      @views lll =  lik+lndetIrhoWt

  end # if length(Wy)==1 

      return -lll

end




function prtlloglike(::Type{SSFKUEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN::Matrix, IV::Matrix,
PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

  Wy = _dicM[:wy]

  liky = prtlloglikekuhe(y, x, Q, EN, IV, Wy, PorC, num, pos, rho,  eigvalu, rowIDT )  

  return liky;       
end



function prtlloglikekkhe( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
   PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
 
 
  β  = rho[1:po.endx]
  τ  = rho[po.begq:po.endq]
  phi = rho[po.begphi:po.endphi]
  @views phi = reshape(phi, :, num.nofeta)
  @views eps = EN- IV*phi

 
  eta = rho[po.begeta:po.endeta]
  δ2 = rho[po.begw]  
  γ  = rho[po.begv]  # May rho[po.begw : po.endw][1]
  # δ1 = rho[po.begz]


  hi  = exp.(Q*τ)
  σᵤ²= exp(δ2) 
  σᵤ= exp(0.5*δ2) 
  σᵥ² = exp(γ)            # todo: 重新换一下字母 
  σᵥ = exp(0.5*γ)  
  μ   = 0.0
  ϵ = PorC*(y - x*β)
  ID = size(rowIDT,1)


  lik = zero(eltype(y));
  @views invPi = 1.0/σᵥ²;

  @floop begin
  @inbounds  for iidd=1:ID  
    @views N = rowIDT[iidd,2];
    @views lndetPi = N*log(σᵥ²);
          @views ind = rowIDT[iidd,1];
          @views his = hi[ind];
          @views ϵs  = ϵ[ind]  - PorC*(eps[ind,:]*eta) ;
          @views sigs2 = 1.0 / (his'*his*invPi + 1.0 /σᵤ²) ;
          @views mus = (μ/σᵤ² - ϵs'*his*invPi)*sigs2 ;
          @views es2 = -0.5*ϵs'*ϵs*invPi ;
          @views KK = 0.5*N*log(2 * π)-0.5*lndetPi;

          @views temp = KK + es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                          0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                          0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                  if simple_check(temp)
                      lik += -1e19
                  else
                      lik += temp
                  end # simple_check(temp)
              end # for ttt=1:ID
      end # begin

  return -lik

end

    
  
    
function prtlloglike(::Type{SSFKKEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN::Matrix, IV::Matrix,
  PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

  liky = prtlloglikekkhe(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

  return liky;       
end
  
      

function prtlloglikekkte( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN::Matrix, IV::Matrix, 
  PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
 
 
  β  = rho[1:po.endx]
  τ  = rho[po.begq:po.endq]
  phi = rho[po.begphi:po.endphi]
  eps = zeros(eltype(EN),num.nofobs,num.nofeta);
  @views phi = reshape(phi, :, num.nofeta)
  @views eps = EN- IV*phi

 
  eta = rho[po.begeta:po.endeta]
  δ2 = rho[po.begw]  
  γ  = rho[po.begv]  # May rho[po.begw : po.endw][1]
  δ1 = rho[po.begz]


  hi  = exp.(Q*τ)
  σᵤ²= exp(δ2) 
  σᵤ= exp(0.5*δ2) 
  σᵥ² = exp(γ)            # todo: 重新换一下字母 
  σᵥ = exp(0.5*γ)  
  μ   = δ1
  ϵ = PorC*(y - x*β)
  ID = size(rowIDT,1)


  lik = zero(eltype(y));
  @views invPi = 1.0/σᵥ²;
  @floop begin
  @inbounds  for iidd=1:ID  
    @views N = rowIDT[iidd,2];
    @views lndetPi = N*log(σᵥ²);
          @views ind = rowIDT[iidd,1];
          @views his = hi[ind];
          @views ϵs  = ϵ[ind]  - PorC*(eps[ind,:]*eta) ;
          @views sigs2 = 1.0 / (his'*his*invPi + 1.0 /σᵤ²) ;
          @views mus = (μ/σᵤ² - ϵs'*his*invPi)*sigs2 ;
          @views es2 = -0.5*ϵs'*ϵs*invPi ;
          @views KK = 0.5*N*log(2 * π)-0.5*lndetPi;

          @views temp = KK + es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                          0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                          0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                  if simple_check(temp)
                      lik += -1e19
                  else
                      lik += temp
                  end # simple_check(temp)
              end # for ttt=1:ID
      end # begin

  return -lik

end

    
  
    
function prtlloglike(::Type{SSFKKET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w::Matrix, v::Matrix, z, EN::Matrix, IV::Matrix,
  PorC::Int64,  num::NamedTuple,  pos::NamedTuple, rho::Array{Float64,1}, eigvalu::NamedTuple, rowIDT::Matrix{Any})

  liky = prtlloglikekkte(y, x, Q, EN, IV, PorC, num, pos, rho,  eigvalu, rowIDT )  

  return liky;       
end
  
      



