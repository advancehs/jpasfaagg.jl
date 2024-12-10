# y x         σᵤ²  σᵥ² -- expo      
# y x         σᵤ²  σᵥ² -- half
# y x μ       σᵤ²  σᵥ² -- trun     
# y x μ  h    σᵤ²  σᵥ² -- scal  
# y x    h    σᵤ²  σᵥ² -- TFE_WH2010, half  
# y x μ  h    σᵤ²  σᵥ² -- TFE_WH2010, truncated  
# y x μ  g    σᵤ²  σᵥ² -- decay  
# y x         σᵤ²  σᵥ² -- panel half (2014 JoE)
# y x    σₐ²  σᵤ²  σᵥ² -- TRE  
# ------------------------------------------
# y x z  q    w    v   -- generic varname
#   β δ1 τ    δ2   γ   -- coeff 


function simple_check(xs)
    any(x -> isnan(x) || !isfinite(x) ,xs)
end






function ssdkuhe( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN,IV,
    Wy::Matrix, PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]
    phi = rho[po.begphi:po.endphi]
## calculate lkx
    nofiv = num.nofphi/num.nofeta
    eps = zeros(eltype(EN),num.nofobs,num.nofeta);

    # %%%%%%%%%%%%%%
    # @inbounds  for ii=1:num.nofeta
    #      eps[:,ii]=EN[:,ii]-IV*phi[((Int(ii)-1)*Int(nofiv)+1):(Int(ii)*Int(nofiv))];
    #  end
    @views phi = reshape(phi, :, num.nofeta)
    @views eps = EN- IV*phi

    @views LL = ((eps .- mean(eps, dims=1))' * (eps .- mean(eps, dims=1)) )/ num.nofobs
    @views logdetll = log(det(LL))
    @views invll = LL\I(num.nofeta)
    likx = zero(eltype(y));


   try
    @floop begin
    @inbounds for iitt =1:num.nofobs
                tempx=-0.5*num.nofeta*log(2*π)-0.5*logdetll-0.5*tr(invll*eps[iitt,:]'*eps[iitt,:]);
                if simple_check(tempx)
                    likx += -1e99
                else
                    likx += tempx
                end # simple_check(tempx)
    
        end # iitt =1:num.nofobs
    end # @floop begin
    ## calculate lky
    
        eta = rho[po.begeta:po.endeta]
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
    
        @views ϵ  =  ϵ-PorC*gamma*Wyt*y  - PorC*(eps*eta) ;
        @views sigs2 = @. 1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 =@. -0.5 * ϵ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                # print(size(temp_1))

        # 检查 lik 是否为 NaN, 非实数, 或 Inf
        @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
        # 计算总和
        @views lik = sum(temp_1)


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
        @views sigs2 = @. 1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 =@. -0.5 * ϵ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                # print(size(temp_1))

        # 检查 lik 是否为 NaN, 非实数, 或 Inf
        @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
        # 计算总和
        @views lik = sum(temp_1)

    end # if length(Wy)==1 
        return -lik-likx-lndetIrhoWt
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end
end


function ssdkuh( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,
    Wy::Matrix, PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]
   try
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
    
        @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
        @views sigs2 = @. 1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 =@.  -0.5 * ϵ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                # print(size(temp_1))

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
        @views ϵ  = ϵ-PorC*gamma*Wyt*y   ;
        @views sigs2 =@. 1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus = @. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 = @. -0.5 * ϵ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                # print(size(temp_1))

        # 检查 lik 是否为 NaN, 非实数, 或 Inf
        @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
        # 计算总和
        @views lik = sum(temp_1)
        @views lll =  lik+lndetIrhoWt

    end # if length(Wy)==1 
        return -lll
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end

end


function LL_T(::Type{SSFKUEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    Wy = _dicM[:wy]

    llt = ssdkuhe(y, x, Q, EN, IV, Wy, PorC, num, po, rho,  eigvalu, rowIDT )  

    return llt
end
    
function LL_T(::Type{SSFKUH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    Wy = _dicM[:wy]

    llt = ssdkuh(y, x, Q,  Wy, PorC, num, po, rho,  eigvalu, rowIDT ) 

    return llt
end



function ssdkute( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN,IV,
    Wy::Matrix, PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]
    phi = rho[po.begphi:po.endphi]
## calculate lkx
    nofiv = num.nofphi/num.nofeta
    eps = zeros(eltype(EN),num.nofobs,num.nofeta);

    # %%%%%%%%%%%%%%
    # @inbounds  for ii=1:num.nofeta
    #      eps[:,ii]=EN[:,ii]-IV*phi[((Int(ii)-1)*Int(nofiv)+1):(Int(ii)*Int(nofiv))];
    #  end
    @views phi = reshape(phi, :, num.nofeta)
    @views eps = EN- IV*phi

    @views LL = ((eps .- mean(eps, dims=1))' * (eps .- mean(eps, dims=1)) )/ num.nofobs
    @views logdetll = log(det(LL))
    @views invll = LL\I(num.nofeta)
    likx = zero(eltype(y));


   try
    @floop begin
        @inbounds for iitt =1:num.nofobs
                 tempx=-0.5*num.nofeta*log(2*π)-0.5*logdetll-0.5*tr(invll*eps[iitt,:]'*eps[iitt,:]);
                    if simple_check(tempx)
                        likx += -1e99
                    else
                        likx += tempx
                    end # simple_check(tempx)
        
            end # iitt =1:num.nofobs
        end # @floop begin
        ## calculate lky
    
        eta = rho[po.begeta:po.endeta]
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
        @views lll =  lik+likx+lndetIrhoWt

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
        @views lll =  lik+likx+lndetIrhoWt

    end # if length(Wy)==1 

        return -lll
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end
end


function ssdkut( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,
    Wy::Matrix, PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]
   try
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

        @views ϵ  =  ϵ-PorC*gamma*Wyt*y   ;
        @views sigs2 =@. 1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 =@. -0.5 * ϵ ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                # print(size(temp_1))

        # 检查 lik 是否为 NaN, 非实数, 或 Inf
        @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
        # 计算总和
        @views lik = sum(temp_1)
        @views lll = lik+lndetIrhoWt


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
        @views ϵ  =  ϵ-PorC*gamma*Wyt*y   ;
        @views sigs2 =@. 1.0 / (hi^2 *invPi + 1.0 /σᵤ²) ;
        @views mus =@. (μ/σᵤ² - ϵ * hi * invPi) * sigs2 ;
        @views es2 =@. -0.5 * ϵ ^2 *invPi;
        @views KK = -0.5*log(2 * π)-0.5*lndetPi;
        @views temp_1 =@. KK +  es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                        0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                        0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))

        # 检查 lik 是否为 NaN, 非实数, 或 Inf
        @views temp_1 = map(x -> x ≠ x ? -1e99  : isinf(x) ? -1e99  : x, temp_1)
        # 计算总和
        @views lik = sum(temp_1)

        @views lll = lik+lndetIrhoWt

    end # if length(Wy)==1 

        return -lll
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end

end


function LL_T(::Type{SSFKUET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    Wy = _dicM[:wy]
    llt = ssdkute(y, x, Q, EN, IV, Wy, PorC, num, po, rho,  eigvalu, rowIDT )  

    return llt
end
    
function LL_T(::Type{SSFKUT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    Wy = _dicM[:wy]

    llt = ssdkut(y, x, Q,  Wy, PorC, num, po, rho,  eigvalu, rowIDT )  
    return llt
end







    
function ssdkkhe( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN,IV,
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]
    phi = rho[po.begphi:po.endphi]
    ## calculate lkx
    nofiv = num.nofphi/num.nofeta
    eps = zeros(eltype(EN),num.nofobs,num.nofeta);

    # %%%%%%%%%%%%%%
    # @inbounds  for ii=1:num.nofeta
    #      eps[:,ii]=EN[:,ii]-IV*phi[((Int(ii)-1)*Int(nofiv)+1):(Int(ii)*Int(nofiv))];
    #  end
    @views phi = reshape(phi, :, num.nofeta)
    @views eps = EN- IV*phi

    @views LL = ((eps .- mean(eps, dims=1))' * (eps .- mean(eps, dims=1)) )/ num.nofobs
    @views logdetll = log(det(LL))
    @views invll = LL\I(num.nofeta)
    likx = zero(eltype(y));

   try
    @floop begin
    @inbounds for iitt =1:num.nofobs
                tempx=-0.5*num.nofeta*log(2*π)-0.5*logdetll-0.5*tr(invll*eps[iitt,:]'*eps[iitt,:]);
                if simple_check(tempx)
                    likx += -1e9
                else
                    likx += tempx
                end # simple_check(tempx)
    
        end # iitt =1:num.nofobs
    end # @floop begin
    ## calculate lky
    
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
                        lik += -1e9
                    else
                        lik += temp
                    end # simple_check(temp)
                end # for ttt=1:ID
        end # begin

    return -lik-likx
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end
end
    
    
    


    
function ssdkkh( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  
     PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]

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

    try
    @floop begin
    @inbounds  for iidd=1:ID  
        @views N = rowIDT[iidd,2];
        @views lndetPi = N*log(σᵥ²);

            @views ind = rowIDT[iidd,1];
            @views his = hi[ind];
            @views ϵs  = ϵ[ind]   ;
            @views sigs2 = 1.0 / (his'*his*invPi + 1.0 /σᵤ²) ;
            @views mus = (μ/σᵤ² - ϵs'*his*invPi)*sigs2 ;
            @views es2 = -0.5*ϵs'*ϵs*invPi ;
            @views KK = 0.5*N*log(2 * π)-0.5*lndetPi;

            @views temp = KK + es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                            0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                            0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                    if simple_check(temp)
                        lik += -1e9
                    else
                        lik += temp
                    end # simple_check(temp)
                end # for ttt=1:ID
        end # begin

    return -lik
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end
end
    
    


function LL_T(::Type{SSFKKEH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    llt = ssdkkhe(y, x, Q, EN, IV, PorC, num, po, rho,  eigvalu, rowIDT )  

    return llt
end
    
function LL_T(::Type{SSFKKH}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    llt = ssdkkh(y, x, Q, PorC, num, po, rho,  eigvalu, rowIDT ) 

    return llt
end



    
function ssdkkte( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  EN,IV,
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]
    phi = rho[po.begphi:po.endphi]
    ## calculate lkx
    nofiv = num.nofphi/num.nofeta
    eps = zeros(eltype(EN),num.nofobs,num.nofeta);

    # %%%%%%%%%%%%%%
    # @inbounds  for ii=1:num.nofeta
    #      eps[:,ii]=EN[:,ii]-IV*phi[((Int(ii)-1)*Int(nofiv)+1):(Int(ii)*Int(nofiv))];
    #  end
    @views phi = reshape(phi, :, num.nofeta)
    @views eps = EN- IV*phi

    @views LL = ((eps .- mean(eps, dims=1))' * (eps .- mean(eps, dims=1)) )/ num.nofobs
    @views logdetll = log(det(LL))
    @views invll = LL\I(num.nofeta)
    likx = zero(eltype(y));

   try
    @floop begin
    @inbounds for iitt =1:num.nofobs
                tempx=-0.5*num.nofeta*log(2*π)-0.5*logdetll-0.5*tr(invll*eps[iitt,:]'*eps[iitt,:]);
                if simple_check(tempx)
                    likx += -1e9
                else
                    likx += tempx
                end # simple_check(tempx)
    
        end # iitt =1:num.nofobs
    end # @floop begin
    ## calculate lky
    
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
        # println("T IS ",N)
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
                        lik += -1e9
                    else
                        lik += temp
                    end # simple_check(temp)
                end # for ttt=1:ID
        end # begin

    return -lik-likx
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end
end
    
    
    


    
function ssdkkt( y::Union{Vector,Matrix}, x::Matrix, Q::Matrix,  
     PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any} )
    β  = rho[1:po.endx]
    τ  = rho[po.begq:po.endq]

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

    try
    @floop begin
    @inbounds  for iidd=1:ID  
        @views N = rowIDT[iidd,2];
        @views lndetPi = N*log(σᵥ²);
            @views ind = rowIDT[iidd,1];
            @views his = hi[ind];
            @views ϵs  = ϵ[ind]   ;
            @views sigs2 = 1.0 / (his'*his*invPi + 1.0 /σᵤ²) ;
            @views mus = (μ/σᵤ² - ϵs'*his*invPi)*sigs2 ;
            @views es2 = -0.5*ϵs'*ϵs*invPi ;
            @views KK = 0.5*N*log(2 * π)-0.5*lndetPi;

            @views temp = KK + es2 + 0.5 * (((mus ^ 2) / sigs2) - (μ^2 / σᵤ²) ) +
                            0.5 * log(sigs2) + log(normcdf(mus / sqrt(sigs2))) -
                            0.5 * log(σᵤ²) - log(normcdf(μ / σᵤ))
                    if simple_check(temp)
                        lik += -1e9
                    else
                        lik += temp
                    end # simple_check(temp)
                end # for ttt=1:ID
        end # begin

    return -lik
    catch e
    # 处理异常的代码
    # println("操作失败，发生错误：$e")
        return 1e100
    end
end
    

function LL_T(::Type{SSFKKET}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    llt = ssdkkte(y, x, Q, EN, IV, PorC, num, po, rho,  eigvalu, rowIDT )  

    return llt
end
    
function LL_T(::Type{SSFKKT}, y::Union{Vector,Matrix}, x::Matrix, Q::Matrix, w, v, z, EN,IV, 
    PorC::Int64, num::NamedTuple, po::NamedTuple, rho,  eigvalu::NamedTuple, rowIDT::Matrix{Any}, ::Nothing) 

    llt = ssdkkt(y, x, Q, PorC, num, po, rho,  eigvalu, rowIDT )  
    return llt
end





