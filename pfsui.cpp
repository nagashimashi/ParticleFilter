#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List PF(NumericVector y1, NumericVector y2,int n,int a,int b,int c,int lag){
    int T = y1.length();
    int veclen = 4;
    double sigma= pow(2,a); //観測モデルの分散
    double alpha1 = pow(10,b); //システムモデル（トレンド）の分散の比
    double alpha2 = pow(10,c); //システムモデル（係数）の分散の比
    
    List x(n); //各期の粒子を格納するリスト
    NumericVector w(n); //尤度から計算した粒子の重み
    NumericVector w_1(n); //和を1に正規化した重み
    double Log_Like = 0; //尤度
    NumericMatrix Pred(T,n); //予測分布
    NumericMatrix Filter(T,n); //トレンドのフィルタ分布
    NumericMatrix Filter_gamma(T,n); //係数のフィルタ分布
    NumericMatrix Filter_lag(T,n);
    List FixedLag(lag*n); //固定ラグベクトルを格納するリスト
    NumericMatrix Smooth(T,n); //トレンドの平滑化分布
    NumericMatrix Smooth_gamma(T,n); //係数の平滑化分布
    List z(n); //リサンプリングした値を格納
    NumericVector F(n+1); //リサンプリング時に用いる重みを足しあわせたベクトル
    
    for(int i=0; i<n; ++i){
        NumericVector initVec(veclen);
        initVec[0] = R::rnorm(0,alpha1*sigma);
        initVec[1] = R::rnorm(0,alpha1*sigma);
        initVec[2] = R::rnorm(1,alpha2*sigma);
        initVec[3] = R::rnorm(1,alpha2*sigma);
        x[i]= initVec;
        
	}
    
    for(int t=0; t<T; ++t){
        
        //観測モデルを動かして予測分布と尤度を計算
        for(int i=0; i<n; ++i){
            double v1; //システムモデル（トレンド）の誤差
            double v2; //システムモデル（係数）の誤差
            NumericVector x_t(2);
            double mu_t;
            double gamma_t;
            double predY;
            
            x_t = x[i];
            v1 = R::rnorm(0,alpha1*sigma); //システムモデルの誤差を発生
            mu_t = 2*x_t[0] - x_t[1] + v1; //予測分布を計算
            
            v2 = R::rnorm(0,alpha2*sigma);
            gamma_t = 2*x_t[2] - x_t[3] +v2;
            
            predY = mu_t + gamma_t*y2[t];
            
            Pred(t,i) = predY; //予測分布を格納
            w[i] = R::dnorm(y1[t],predY,sigma,false); //予測分布の尤度を計算
            
            NumericVector newVec(veclen);
            newVec[0] = mu_t;
            newVec[1] = x_t[0];
            newVec[2] = gamma_t;
            newVec[3] = x_t[2];
            x[i] = newVec;
            
            //ここで固定ラグを更新
            if(t<lag){
                    FixedLag[t+lag*i] = newVec;
            }
            else{
                    for(int j=0;j<lag-1;++j){
                        NumericVector lagVec(veclen);
                        lagVec = FixedLag[j+1+lag*i];
                        FixedLag[j+lag*i] = lagVec;
                    }
                    FixedLag[lag-1+lag*i] = newVec;
            }
        }
        //観測モデル動作完了
        
        //重みを計算する
        double sum_w = 0;
        for(int i=0;i<n;++i){
            sum_w = sum_w + w[i];
        }
        Log_Like = Log_Like +  log(sum_w);
        for(int i=0; i<n; ++i){
            w_1[i] = w[i]/sum_w; //重みを正規化
        }
        
        //これからリサンプリングを行う
        F[0] = 0;
        F[1] = w_1[0];
        for(int i=2; i<n+1; i++){
            F[i] = w_1[i-1]+F[i-1];
        }
        double seed;
        seed = R::runif(0,1);
        NumericVector ind(n);

        for(int k=0; k<n; k++){
            double u;
            u = (seed + k)/n; //櫛の値
            for(int j=0; j<n; ++j){
                if((u>F[j]) && (u<=F[j+1])){
                    NumericVector samplingVec(veclen);
                    samplingVec = x[j];
                    ind[k] = j;
                    z[k] = samplingVec; //櫛があたったところの予測分布の値を格納
                    break;
                }
            }
        }
        //リサンプリング完了
        
        //フィルタ分布を格納する
        for(int i=0;i<n;++i){
            NumericVector newVec(veclen);
            newVec = z[i];
            x[i] = newVec; //フィルタ分布を計算
            Filter(t,i) = newVec[0];
            Filter_gamma(t,i) = newVec[2];
            Filter_lag(t,i) = newVec[1];

        }
        
        //固定ラグの更新
        NumericMatrix newMat(lag*n,veclen);
        for(int i=0;i<n;++i){
            if(t<lag){
                for(int j=0;j<=t;++j){
                    NumericVector newVec(veclen);
                    newVec = FixedLag[j+lag*ind[i]];
                    newMat(j+lag*i,_) = newVec;
                }
            }
            else{
                for(int j=0;j<lag;++j){
                    NumericVector newVec(veclen);
                    newVec = FixedLag[j+lag*ind[i]];
                    newMat(j+lag*i,_) = newVec;
                }
            }
        }
        for(int i=0;i<n;++i){
            for(int j=0;j<lag;++j){
                FixedLag[j+lag*i] = newMat(j+lag*i,_);
            }
        }
        
        //平滑化分布の更新
        if(t>=lag-1){
            for(int i=0;i<n;++i){
                NumericVector smoothVec(veclen);
                smoothVec =FixedLag[lag*i];
                Smooth(t-lag+1,i) = smoothVec[0];
                Smooth_gamma(t-lag+1,i) = smoothVec[2];
            }
        }
        if(t==T-1){
            for(int i=0;i<n;++i){
                for(int j=0;j<lag;++j){
                    NumericVector smoothVec(veclen);
                    smoothVec = FixedLag[j+lag*i];
                    Smooth(T-lag+j,i) = smoothVec[0];
                    Smooth_gamma(T-lag+j,i) = smoothVec[2];
                }
            }
        }
        
    }
    
	return(List::create(Named("Predict")=Pred,Named("Filter")=Filter,Named("Filter_lag")=Filter_lag,Named("Filter_gamma")=Filter_gamma,Named("Smoother")=Smooth,Named("Smoother_gamma")=Smooth_gamma,Named("Log_Like")=Log_Like));
}

