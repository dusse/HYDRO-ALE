#include "HeatManager.hpp"

using namespace std;
using namespace chrono;


HeatManager::HeatManager(std::shared_ptr<Loader> ldr, std::shared_ptr<GridManager> gridMnr):loader(move(ldr)), gridMgr(move(gridMnr)){
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[HeatManager] create...OK");
}

void HeatManager::initialize(){
    int xRes = loader->resolution[0], yRes = loader->resolution[1];
    int i,j;
    for ( i=0; i < xRes+2; i++){
        for ( j=0; j < yRes+2; j++){
            if(i < xRes/2+4 && i > xRes/2-4 && j < yRes/2+4 && j > yRes/2-4){
                currentTemperature.push_back(ScalarVar(TEMPERATURE, INITIAL_TEMPERATURE));
            }else{
                currentTemperature.push_back(ScalarVar(TEMPERATURE, INITIAL_TEMPERATURE));
            }
            
            currentFlux.push_back(VectorVar(FLUX, { 0.0, 0.0 }));
            currentDivFlux.push_back(ScalarVar(DIV_Q, 0.0));
        }
    }
    
}

/*
 *                                      i-1j+1 ----- ij+1 ---- i+1j+1
 *                                          |         |         |
 *    u' = -div(w)                          |    sηij |         |
 *                                          |         |_ ϕij    |
 *    w = -k∇T                              |         . |  sξij |
 *                                        i-1j ------ ij ---- i+1j
 *                                          |       |_._|       |
 *                                          |ϕi-1j-1  | ϕij-1   |
 *                                          |         |sη'ij    |
 *                                          |         |         |
 *                                      i-1j-1 ----- ij-1 --- i+1j-1
 *
 *              Un+1 - Un
 *   Ω(i,j) x  ----------  = - ( wξ(i+1,j) x sη(i+1,j) - wξ(i,j) x sη(i,j) + wη(i,j+1) x sξ(i,j+1) - wη(i,j) x sξ(i,j) )
 *                 Δt
 *
 *              wξ^2+ wη^2 + 2wξwηcos(ϕ)
 *      w^2 =  --------------------------
 *                     sin^2(ϕ)
 *
 *             w^2
 *     F =  ∫ ---- dΩ - 2∫u div(w) dΩ
 *          Ω   k        Ω
 *
 *                Ωij        (wξi+α,j)^2+ (wηi,j+β)^2 + (-1)^α+β 2(wξi+α,j)(wηi,j+β)cos(ϕi+α,j+β)
 *      F(w) =  Σ ---    Σ  --------------------------------------------------------------------- + 2uij  Σ (-1)^α( wξi+α,jsηi+α,j + wηi,j+αsξi,j+α )
 *             Ωij 4  α,β=0,1                        sin^2(ϕi+α,j+β)                                     α=0,1
 *
 *
 *                            Ωi-α,j                           α+β  Ω[i-α,j]cos(ϕ[i,j+β]{i-α,j})
 *              ( Σ   -------------------------) wξi,j +  Σ (-1)   ------------------------------  wηi-α,j+β) = - sηi,j( T(i,j) - T(i-1,j) )
 *             α,β=0,1 4sin^2(ϕ[i,j+β]{i-α,j})          α,β=0,1      4sin^2(ϕ[i,j+β]{i-α,j})
 *
 *
 *
 *                            Ωi,j-α                           α+β  Ω[i,j-α]cos(ϕ[i+β,j]{i,j-α})
 *              ( Σ   -------------------------) wηi,j +  Σ (-1)   ------------------------------  wξi+β,j-α) = - sξi,j( T(i,j) - T(i,j-1))
 *             α,β=0,1 4sin^2(ϕ[i+β,j]{i,j-α})          α,β=0,1      4sin^2(ϕ[i+β,j]{i,j-α})
 *      
 *   ________________________________________________________________________________________________________________________________________________________________
 *
 *              Aξ(i,j)wξ(i-1,j) + Cξ(i,j)wξ(i,j) + Bξ(i,j)wξ(i+1,j) = Fξ(i,j)[wη(i,j), wη(i,j+1), wη(i-1,j), wη(i-1,j+1)]
 *
 *                  Δt*sη(i-1,j)*sξ(i-1,j)
 *           A = - ----------------------
 *                        Ω(i-1,j)
 *
 *                  Δt*sη(i+1,j)*sξ(i+1,j)
 *           C = - ----------------------
 *                        Ω(i,j)
 *
 *                 1             Ω(i-α,j)                            1           1
 *           B =   Σ    -------------------------  + Δt*sη^2(i,j)( ------  +  ------- )
 *                α,β=0   4sin^2(ϕ(i,j+β){i-α,j})                  Ω(i,j)     Ω(i-1,j)
 *
 *
 *
 *                 1      α+β+1   Ω(i-α,j)cos(ϕ(i,j+β){i-α,j})    Δt*sη(i,j)*sξ(i-α,j+β)
 *           F =   Σ  (-1)     [ ------------------------------ + ---------------------- ] x wη(i-α,j+β) - sη(i,j) x ( T(i,j) - T(i-1,j) )
 *               α,β=0             4sin^2(ϕ(i,j+β){i-α,j})              Ω(i-α,j)
 *
 *          i = 1 : N-2
*           j = 0 : M-2
 *   ________________________________________________________________________________________________________________________________________________________________
 *
 *
 *              Aη(i,j)wη(i,j-1) + Cη(i,j)wη(i,j) + Bη(i,j)wη(i,j+1) = Fη(i,j)[wξ(i,j), wξ(i+1,j), wξ(i,j-1), wξ(i+1,j-1)]
 *
 *                  Δt*sξ(i,j-1)*sη(i,j-1)
 *           A = - ----------------------
 *                        Ω(i,j-1)
 *
 *                  Δt*sξ(i,j+1)*sη(i,j+1)
 *           C = - ----------------------
 *                        Ω(i,j)
 *
 *                 1             Ω(i,j-α)                            1           1
 *           B =   Σ    -------------------------  + Δt*sξ^2(i,j)( ------  +  ------- )
 *                α,β=0   4sin^2(ϕ(i+β,j){i,j-α})                  Ω(i,j)     Ω(i,j-1)
 *
 *
 *
 *                 1      α+β+1   Ω(i,j-α)cos(ϕ(i+β,j){i,j-α})    Δt*sξ(i,j)*sη(i+β,j-α)
 *           F =   Σ  (-1)     [ ------------------------------ + ---------------------- ] x wξ(i+β,j-α) - sξ(i,j) x ( T(i,j) - T(i,j-1) )
 *               α,β=0             4sin^2(ϕ(i+β,j){i,j-α})              Ω(i,j-α)
 *
 *
 *         i = 0 : N-2
 *         j = 1 : M-2
 *   ________________________________________________________________________________________________________________________________________________________________
 *
 *
 *
 */
void HeatManager::solve(){
    auto start_time = high_resolution_clock::now();
    

    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int i,j,k = 0, ijNode;
    
    initCurrentTemperature();
    

    int i1jNode,ij1Node;
    int i1j1Node;
    double Xij, Xi1j, Xij1, Xi1j1;
    double Yij, Yi1j, Yij1, Yi1j1;

    vector<VectorVar> coords = gridMgr->getExtendedGridCoordinates();
    vector<VectorVar> side;
    vector<VectorVar> sin;
    vector<ScalarVar> vol;
    
    for(int tmp = 0; tmp < (xRes+2)*(yRes+2); tmp++){
        vol.push_back(ScalarVar(VOLUME, 0.0));
        side.push_back(VectorVar(BASIS, {0.0, 0.0}));
        sin.push_back(VectorVar(SINUS, {0.0, 0.0, 0.0, 0.0}));
    }
    vector<VectorVar> ksiCoefs, ettaCoefs;
    double volij;
    double ksi_ij, etta_ij, ksi_ij1, etta_i1j;
    for ( i=0; i<xRes+1; i++){
        for ( j=0; j<yRes+1; j++){
            ijNode   = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            i1jNode  = IDX(i+1, j  , k, xRes+2, yRes+2, zRes);
            ij1Node  = IDX(i  , j+1, k, xRes+2, yRes+2, zRes);
            i1j1Node = IDX(i+1, j+1, k, xRes+2, yRes+2, zRes);
            
            Xij   = coords[ijNode].getValue()[0];
            Yij   = coords[ijNode].getValue()[1];
            Xi1j  = coords[i1jNode].getValue()[0];
            Yi1j  = coords[i1jNode].getValue()[1];
            Xij1  = coords[ij1Node].getValue()[0];
            Yij1  = coords[ij1Node].getValue()[1];
            Xi1j1 = coords[i1j1Node].getValue()[0];
            Yi1j1 = coords[i1j1Node].getValue()[1];
            
            ksi_ij  = sqrt(pow((Xi1j - Xij),2) + pow((Yi1j - Yij),2));// sξij = √{[x(i+1,j)-x(i,j)]^2+[y(i+1,j)-y(i,j)]^2}
            etta_ij = sqrt(pow((Xij1 - Xij),2) + pow((Yij1 - Yij),2));// sηij = √{[x(i,j+1)-x(i,j)]^2+[y(i,j+1)-y(i,j)]^2}
            side[ijNode]  = VectorVar(BASIS, { ksi_ij, etta_ij });
            
            ksi_ij1  = sqrt(pow((Xi1j1 - Xij1),2) + pow((Yi1j1 - Yij1),2));// sξij+1 = √{[x(i+1,j+1)-x(i,j+1)]^2+[y(i+1,j+1)-y(i,j+1)]^2}
            etta_i1j = sqrt(pow((Xi1j1 - Xi1j),2) + pow((Yi1j1 - Yi1j),2));// sηi+1j = √{[x(i+1,j+1)-x(i+1,j)]^2+[y(i+1,j+1)-y(i+1,j)]^2}
            
            volij = 0.5 * ((Xi1j1 - Xij)*(Yij1 - Yi1j) - (Xij1 - Xi1j)*(Yi1j1 - Yij));//ij
            vol[ijNode]  = ScalarVar(VOLUME, volij);
            if( volij < 0.0 ){
                string msg ="[HeatManager] Volume problem volij = "+to_string(volij);
                logger->writeMsg(msg.c_str(), CRITICAL);
            }

            vector<double> sinuses;
            double sinTmp = volij/ksi_ij /etta_ij;
            sinTmp = sinTmp > 1.0 ? 1.0 : sinTmp;
            sinuses.push_back(sinTmp) ; // i  ,j
            sinTmp = volij/ksi_ij /etta_i1j;
            sinTmp = sinTmp > 1.0 ? 1.0 : sinTmp;
            sinuses.push_back(sinTmp); // i+1,j
            sinTmp = volij/ksi_ij1/etta_ij;
            sinTmp = sinTmp > 1.0 ? 1.0 : sinTmp;
            sinuses.push_back(sinTmp) ; // i  ,j+1
            sinTmp = volij/ksi_ij1/etta_i1j;
            sinTmp = sinTmp > 1.0 ? 1.0 : sinTmp;
            sinuses.push_back(sinTmp); // i+1,j+1
            sin[ijNode] = VectorVar(SINUS, sinuses);
        }
    }
    

    
    ksiCoefs = getKsiCoeficients(side, sin, vol);
    
    //solve system
    vector<double> newFluxesKsi = solveTridiagonalMatrix( ksiCoefs );
    for ( i=1; i<xRes; i++){
        for ( j=0; j<yRes; j++){
            ijNode   = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            int idx = IDX(i-1  , j  , k, xRes-1, yRes, zRes);
            currentFlux[ijNode] = VectorVar(FLUX, {newFluxesKsi[idx],  currentFlux[ijNode].getValue()[1]});
        }
    }
    
    ettaCoefs = getEttaCoeficients(side, sin, vol);
    
    //solve system
    vector<double> newFluxesEtta = solveTridiagonalMatrix( ettaCoefs );
    for ( i=0; i<xRes; i++){
        for ( j=1; j<yRes; j++){
            ijNode   = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            int idx = IDX(i  , j-1  , k, xRes, yRes-1, zRes);
            currentFlux[ijNode] = VectorVar(FLUX, {currentFlux[ijNode].getValue()[0], newFluxesEtta[idx]});
        }
    }
    
    calculateDivQAndTemperature(side, vol);
    
    auto end_time = high_resolution_clock::now();
    string msg ="[HeatManager] solve() duration = "+to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str());
    
}

vector<VectorVar> HeatManager::getKsiCoeficients(vector<VectorVar> side, vector<VectorVar> sin, vector<ScalarVar> vol ){
    vector<VectorVar> ksiCoefs;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    double timeDelta = loader->getTimeStep();
    double ksiA, ksiB, ksiC, ksiF;
    int i,j,k = 0;
    int ijNode, i1jNode, ij1Node, i_1j1Node, i_1jNode;
    double ksiij, ettaij, tempij;
    double ksii_1j, ettai_1j, tempi_1j;
    double ksii_1j1, ksiij1, ksii1j, ettai1j;
    double volij, voli_1j;
    vector<double> sinij, sini_1j;
    double kappa = 0.0; //classical Spitzer-Harm plasma heat conductivity C*T^5/2/(ZlnΛ)
    
    for ( i=1; i<xRes+1; i++){
        for ( j=0; j<yRes+1; j++){
            ijNode   = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            i1jNode  = IDX(i+1, j  , k, xRes+2, yRes+2, zRes);
            ij1Node  = IDX(i  , j+1, k, xRes+2, yRes+2, zRes);
            
            i_1j1Node  = IDX(i-1, j+1, k, xRes+2, yRes+2, zRes);
            i_1jNode   = IDX(i-1, j  , k, xRes+2, yRes+2, zRes);
            
            ksii_1j  = side[i_1jNode].getValue()[0];
            ettai_1j = side[i_1jNode].getValue()[1];
            tempi_1j = currentTemperature[i_1jNode].getValue();
            
            ksiij  = side[ijNode].getValue()[0];
            ettaij = side[ijNode].getValue()[1];
            tempij = currentTemperature[ijNode].getValue();
            ksii1j  = side[i1jNode].getValue()[0];
            ettai1j = side[i1jNode].getValue()[1];
            ksii_1j1  = side[i_1j1Node].getValue()[0];
            ksiij1  = side[ij1Node].getValue()[0];
            volij  = vol[ijNode].getValue();
            voli_1j  =  vol[i_1jNode].getValue();
            sinij   = sin[ijNode].getValue();
            sini_1j = sin[i_1jNode].getValue();
            
            
            /*                  Δt*sη(i-1,j)*sξ(i-1,j)
             *           A = - ----------------------
             *                        Ω(i-1,j)
             */
            ksiA = -timeDelta*ettai_1j*ksii_1j/voli_1j;
            
            
            /*                 1             Ω(i-α,j)                            1           1
             *           B =    Σ    -------------------------  + Δt*sη^2(i,j)( ------  +  ------- )
             *                α,β=0   4sin^2(ϕ(i,j+β){i-α,j})                   Ω(i,j)     Ω(i-1,j)
             */
            kappa = getConductivity(tempij);
            ksiB  = 0.25*volij/(pow(sinij[0],2))/kappa;
            ksiB += 0.25*volij/(pow(sinij[2],2))/kappa;
            kappa = getConductivity(tempi_1j);
            ksiB += 0.25*voli_1j/(pow(sini_1j[0],2))/kappa;
            ksiB += 0.25*voli_1j/(pow(sini_1j[2],2))/kappa;
            ksiB += timeDelta*ettaij*ettaij*(1/volij + 1/voli_1j);
            
            
            /*                  Δt*sη(i+1,j)*sξ(i+1,j)
             *            C = - ----------------------
             *                        Ω(i,j)
             */
            ksiC = -timeDelta*ettai1j*ksii1j/volij;
            
            
            /*                 1      α+β+1   Ω(i-α,j)cos(ϕ(i,j+β){i-α,j})    Δt*sη(i,j)*sξ(i-α,j+β)
             *            F =   Σ  (-1)     [ ------------------------------ + ---------------------- ] x wη(i-α,j+β) - sη(i,j) x ( T(i,j) - T(i-1,j) )
             *                α,β=0             4sin^2(ϕ(i,j+β){i-α,j})              Ω(i-α,j)
             */
            kappa = getConductivity(tempij);
            ksiF  = -(0.25*volij  *sqrt((1 - pow(sinij  [0],2)))/(pow(sinij  [0],2))/kappa + timeDelta*ettaij*ksiij   )*currentFlux[ijNode].getValue()[1];   // α,β=0,0
            ksiF +=  (0.25*volij  *sqrt((1 - pow(sinij  [2],2)))/(pow(sinij  [2],2))/kappa + timeDelta*ettaij*ksiij1  )*currentFlux[ij1Node].getValue()[1];  // α,β=0,1
            kappa = getConductivity(tempi_1j);
            ksiF +=  (0.25*voli_1j*sqrt((1 - pow(sini_1j[0],2)))/(pow(sini_1j[0],2))/tempi_1j + timeDelta*ettaij*ksii_1j )*currentFlux[i_1jNode].getValue()[1]; // α,β=1,0
            ksiF += -(0.25*voli_1j*sqrt((1 - pow(sini_1j[2],2)))/(pow(sini_1j[2],2))/tempi_1j + timeDelta*ettaij*ksii_1j1)*currentFlux[i_1j1Node].getValue()[1];// α,β=1,1
            ksiF -= ettaij*(currentTemperature[ijNode].getValue() - currentTemperature[i_1jNode].getValue());
            ksiCoefs.push_back(VectorVar("KSI_COEF", { ksiA, ksiB, ksiC, ksiF }));
        }
    }
    return ksiCoefs;
}


/*
 *                       i-1j+1 ----- ij+1 ---- i+1j+1
 *                           |         |         |
 *                           |    sηij |         |
 *                           |         |_ ϕij    |
 *                           |         . |  sξij |
 *                         i-1j ------ ij ---- i+1j
 *                           |       |_._|       |
 *                           |ϕi-1j-1  | ϕij-1   |
 *                           |         |sη'ij    |
 *                           |         |         |
 *                        i-1j-1 ----- ij-1 --- i+1j-1
 */

vector<VectorVar> HeatManager::getEttaCoeficients(vector<VectorVar> side, vector<VectorVar> sin, vector<ScalarVar> vol ){
    vector<VectorVar> ettaCoefs;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int i,j,k = 0;
    double timeDelta = loader->getTimeStep();
    double ettaA, ettaB, ettaC, ettaF;
    int ijNode, i1jNode, ij1Node, i1j_1Node, ij_1Node;
    vector<double> sinij, sinij_1;

    double kappa = 0.0;
    
    double volij, volij_1;
    
    double ksiij_1, ettaij_1;
    double ksiij  , ettaij  ;
    double ksiij1 , ettaij1 ;
    double ksii1j , ettai1j ;
    double ettai1j_1;
    double tempij, tempij_1;
    
    for ( i=0; i<xRes+1; i++){
        for ( j=1; j<yRes+1; j++){
            ijNode   = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            i1jNode  = IDX(i+1, j  , k, xRes+2, yRes+2, zRes);
            ij1Node  = IDX(i  , j+1, k, xRes+2, yRes+2, zRes);
            
            i1j_1Node  = IDX(i+1, j-1, k, xRes+2, yRes+2, zRes);
            ij_1Node   = IDX(i  , j-1, k, xRes+2, yRes+2, zRes);
            
            ksiij  = side[ijNode].getValue()[0];
            ettaij = side[ijNode].getValue()[1];
            tempij = currentTemperature[ijNode].getValue();
            ksii1j  = side[i1jNode].getValue()[0];
            ettai1j = side[i1jNode].getValue()[1];
            
            
            ettai1j_1 = side[i1j_1Node].getValue()[1];
            
            ksiij1  = side[ij1Node].getValue()[0];
            ettaij1 = side[ij1Node].getValue()[1];
            ksiij_1  = side[ij_1Node].getValue()[0];
            ettaij_1 = side[ij_1Node].getValue()[1];
            tempij_1 = currentTemperature[ij_1Node].getValue();
            
            volij  = vol[ijNode].getValue();
            volij_1  =  vol[ij_1Node].getValue();
            sinij_1 = sin[ij_1Node].getValue();
            sinij = sin[ijNode].getValue();
            /*                  Δt*sξ(i,j-1)*sη(i,j-1)
             *           A = - ----------------------
             *                        Ω(i,j-1)
             */
            ettaA = -timeDelta*ksiij_1*ettaij_1/volij_1;
            
            
            /*                 1             Ω(i,j-α)                            1           1
             *            B =   Σ    -------------------------  + Δt*sξ^2(i,j)( ------  +  ------- )
             *                 α,β=0   4sin^2(ϕ(i+β,j){i,j-α})                  Ω(i,j)     Ω(i,j-1)
             */
            kappa = getConductivity(tempij);
            ettaB  = 0.25*volij/(pow(sinij[0],2))/kappa;
            ettaB += 0.25*volij/(pow(sinij[2],2))/kappa;
            kappa = getConductivity(tempij_1);
            ettaB += 0.25*volij_1/(pow(sinij_1[0],2))/kappa;
            ettaB += 0.25*volij_1/(pow(sinij_1[2],2))/kappa;
            ettaB += timeDelta*ettaij*ettaij*(1/volij + 1/volij_1);
            
            
            /*                  Δt*sξ(i,j+1)*sη(i,j+1)
             *            C = - ----------------------
             *                         Ω(i,j)
             */
            ettaC = -timeDelta*ksiij1*ettaij1/volij;
            
            /*                 1      α+β+1   Ω(i,j-α)cos(ϕ(i+β,j){i,j-α})    Δt*sξ(i,j)*sη(i+β,j-α)
             *            F =   Σ  (-1)     [ ------------------------------ + ---------------------- ] x wξ(i+β,j-α) - sξ(i,j) x ( T(i,j) - T(i,j-1) )
             *                α,β=0             4sin^2(ϕ(i+β,j){i,j-α})              Ω(i,j-α)
             */
            kappa = getConductivity(tempij);
            ettaF  = -(0.25*volij  *sqrt((1 - pow(sinij  [0],2)))/(pow(sinij  [0],2))/kappa + timeDelta*ksiij*ettaij   )*currentFlux[ijNode].getValue()[0];   // α,β=0,0
            ettaF +=  (0.25*volij  *sqrt((1 - pow(sinij  [1],2)))/(pow(sinij  [1],2))/kappa + timeDelta*ksiij*ettai1j  )*currentFlux[i1jNode].getValue()[0];  // α,β=0,1
            kappa = getConductivity(tempij_1);
            ettaF +=  (0.25*volij_1*sqrt((1 - pow(sinij_1[0],2)))/(pow(sinij_1[0],2))/kappa + timeDelta*ksiij*ettaij_1 )*currentFlux[ij_1Node].getValue()[0]; // α,β=1,0
            ettaF += -(0.25*volij_1*sqrt((1 - pow(sinij_1[1],2)))/(pow(sinij_1[1],2))/kappa + timeDelta*ksiij*ettai1j_1)*currentFlux[i1j_1Node].getValue()[0];// α,β=1,1
            ettaF -= ksiij*(currentTemperature[ijNode].getValue() - currentTemperature[ij_1Node].getValue());
            
            ettaCoefs.push_back(VectorVar("ETTA_COEF", { ettaA, ettaB, ettaC, ettaF }));
        }
    }
    return ettaCoefs;
}

double HeatManager::getConductivity(double temp){
    double kappa = 0.0;//classical Spitzer-Harm plasma heat conductivity C*T^5/2/(ZlnΛ)
    const double kBol = 1.3807e-16;// erg/K Boltzmann constant
    const double me = 9.1094e-28;// g
    const double e = 4.8032e-10;// statcoul
    const double Z =  3.0;//plasma mean ion charge
    const double deltaee = 0.095*(Z+0.24)/(1+0.24*Z);//electron-electron collision term
    const double C = 20*pow(2/PI, 1.5)*pow(kBol, 3.5)/sqrt(me)/pow(e, 4)*deltaee;
    const double lambdaLn = 10.0;
    kappa = C*pow(temp, 2.5)/Z/lambdaLn;
    return kappa;
}

void HeatManager::initCurrentTemperature(){
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int i,j,k = 0, ijNode, ijZone;
    vector<ScalarVar> presures  = gridMgr->getScalarVariableForAllZones("pres");
    vector<ScalarVar> densities = gridMgr->getScalarVariableForAllZones("dens");
    const double A = 27.0;//atomic mass
    const double Z =  3.0;//plasma mean ion charge
    const double Cpt = 0.9648e12;// [erg/eV/g] pressure/temperature conversion factor
    double mq = A/(Z+1)/Cpt;
    //                    A       p
    double Tcur;// T =  ----- x ------
    //                  Z + 1  C(p/T) ρ
    double pres, dens;
    for ( i=0; i < xRes+1; i++){
        for ( j=0; j < yRes+1; j++){
            ijNode = IDX(i, j, k, xRes+2, yRes+2, zRes);
            ijZone = IDX(i, j, k, xRes+1, yRes+1, zRes);
            dens = densities[ijZone].getValue();
            pres = presures [ijZone].getValue();
            Tcur = mq*pres/dens;
            currentTemperature[ijNode] = ScalarVar(TEMPERATURE, Tcur);
            
        }
    }
}

void HeatManager::calculateDivQAndTemperature(vector<VectorVar> side, vector<ScalarVar> vol){
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    double timeDelta = loader->getTimeStep();
    int i,j,k = 0, ijNode;
    int i1jNode, ij1Node;
    double ksiij, ettaij, ettaij1, ksii1j;
    double volij;
    double tempNext, tempPrev, fluxDiv;
    double fluxi1jKsi, fluxijKsi;
    double fluxij1Etta, fluxijEtta;
    for ( i=0; i<xRes+1; i++){
        for ( j=0; j<yRes+1; j++){
            ijNode   = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            i1jNode  = IDX(i+1, j  , k, xRes+2, yRes+2, zRes);
            ij1Node  = IDX(i  , j+1, k, xRes+2, yRes+2, zRes);
            
            fluxijKsi   = currentFlux[ijNode].getValue()[0];
            fluxijEtta  = currentFlux[ijNode].getValue()[1];
            ksiij  = side[ijNode].getValue()[0];
            ettaij = side[ijNode].getValue()[1];
            
            ksii1j  = side[i1jNode].getValue()[0];
            fluxi1jKsi  = currentFlux[i1jNode].getValue()[0];
            
            fluxij1Etta = currentFlux[ij1Node].getValue()[1];
            ettaij1 = side[ij1Node].getValue()[1];
           
            volij  = vol[ijNode].getValue();
            tempPrev = currentTemperature[ijNode].getValue();
            /*              Un+1 - Un
             *    Ω(i,j) x  ----------  = - ( wξ(i+1,j) x sξ(i+1,j) - wξ(i,j) x sξ(i,j) + wη(i,j+1) x sη(i,j+1) - wη(i,j) x sη(i,j) )
             *                  Δt
             */
            fluxDiv = fluxi1jKsi*ksii1j - fluxijKsi*ksiij + fluxij1Etta*ettaij1 - fluxijEtta*ettaij;
            currentDivFlux[ij1Node] = ScalarVar(DIV_Q, fluxDiv);
            tempNext = tempPrev - timeDelta*fluxDiv/volij;
            currentTemperature[ijNode] = ScalarVar(TEMPERATURE, tempNext);
        }
    }

}

/*
 // n is the number of unknowns
 
 |b0 c0 0 ||x0| |d0|
 |a1 b1 c1||x1|=|d1|
 |0  a2 b2||x2| |d2|
 
 1st iteration:
 b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->  
 x0 + g0x1 = r0
 g0 = c0/b0
 r0 = d0/b0
 
 2nd iteration:     
 a1x0 + b1x1 + c1x2 = d1
 from 1st it.: -| a1x0 + a1g0x1 = a1r0
 -----------------------------
 (b1 - a1g0)x1 + c1x2 = d1 - a1r0
 
 x1 + g1x2 = r1               
 g1=c1/(b1 - a1g0)
 r1 = (d1 - a1r0)/(b1 - a1g0)
 
 3rd iteration:      
 a2x1 + b2x2 = d2
 from 2st it. : -| a2x1 + a2g1x2 = a2r2
 -----------------------
 (b2 - a2g1)x2 = d2 - a2r2
 x2 = r2                      
 r2 = (d2 - a2r2)/(b2 - a2g1)
 
 Finally triangular matrix:
 |1  g0 0 ||x0| |r0|
 |0  1  g1||x1|=|r1|
 |0  0  1 ||x2| |r2|
 
 Condition: ||bi|| > ||ai|| + ||ci||????
 */
vector<double> HeatManager::solveTridiagonalMatrix( vector<VectorVar> coefs ) {
    
    vector<double> a, b, c, d;
    vector<VectorVar>::iterator iterator;
    for( iterator = coefs.begin(); iterator != coefs.end(); iterator++) {
        a.push_back(iterator->getValue()[0]);
        b.push_back(iterator->getValue()[1]);
        c.push_back(iterator->getValue()[2]);
        d.push_back(iterator->getValue()[3]);
        if(abs(iterator->getValue()[1]) < abs(iterator->getValue()[0]) + abs(iterator->getValue()[2])){
            string msg ="[HeatManager] CONDITION FAIL a ="+to_string(iterator->getValue()[0])+", b = "+to_string(iterator->getValue()[1])+", c = "+to_string(iterator->getValue()[2]);
            logger->writeMsg(msg.c_str(), CRITICAL);
        }
    }
    int n = (int) a.size();

    n--; // since we start from x0 (not x1)
    c[0] /= b[0];
    d[0] /= b[0];
    
    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }
    
    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);
    
    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
    return d;
}


vector<ScalarVar> HeatManager::getCurrentTemperature(){
    int idx;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int i,j,k = 0;
    vector<ScalarVar> temp;
    for ( i=1; i<xRes+2; i++){
        for ( j=1; j<yRes+2; j++){
            idx    = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            temp.push_back(ScalarVar(TEMPERATURE, currentTemperature[idx].getValue()));
        }
    }

    return temp;
}

vector<VectorVar> HeatManager::getCurrentFlux(){
    int idx;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int i,j,k = 0;
    vector<VectorVar> flux;
    for ( i=0; i<xRes+1; i++){
        for ( j=0; j<yRes+1; j++){
            idx    = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            flux.push_back(VectorVar(FLUX, currentFlux[idx].getValue()));
        }
    }
    return flux;
}


vector<ScalarVar> HeatManager::getFluxDivergence(){
    int idx;
    int xRes = loader->resolution[0], yRes = loader->resolution[1], zRes = 1;
    int i,j,k = 0;
    vector<ScalarVar> divQ;
    for ( i=0; i<xRes+1; i++){
        for ( j=0; j<yRes+1; j++){
            idx    = IDX(i  , j  , k, xRes+2, yRes+2, zRes);
            divQ.push_back(ScalarVar(DIV_Q, currentDivFlux[idx].getValue()));
        }
    }
    return divQ;
}


