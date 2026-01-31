




//BFL Linear//面對曲面使用Interpolation Bounce Method
    if( k == 3 || k == 4 ) {
        idx_xi = (k-3)*NYD6+j; //idx_xi哪裡宣告? evoplution.h84 
        if ( BFLReqF3_d[(k-3)*NYD6+j] == 1 ){//用1跟0來判斷是否需要壁面處理 
            //內插反彈技巧:為Half-way Bounce Back的變形 
            if(Q3_d[(k-3)*NYD6+j] > 0.5){ //Q3_d[(k-3)*NYD6+j]為正規化 邊界點Boundary node 到壁面上的距離 
                F3_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f4_old[index]+ ((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f3_old[index]; //(1/2q)*f4 + (2q-1)/2q*f3 , 3.4反方向 //線性外插
                F7_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f10_old[index]+((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f7_old[index];
                F8_in = (1/(2*Q3_d[(k-3)*NYD6+j]))*f9_old[index]+ ((2*Q3_d[(k-3)*NYD6+j]-1)/(2*Q3_d[(k-3)*NYD6+j]))*f8_old[index];
            }
            if(Q3_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(  f4_old,  F3_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6,  XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                //X_Y_XI_Intrpl7(f10_old, F7_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf37_0, XBFLf37_1, XBFLf37_2, XBFLf37_3, XBFLf37_4, XBFLf37_5, XBFLf37_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                //X_Y_XI_Intrpl7(f9_old,  F8_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf38_0, XBFLf38_1, XBFLf38_2, XBFLf38_3, XBFLf38_4, XBFLf38_5, XBFLf38_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                Y_XI_Intrpl7(  f4_old,  F3_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6,  XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                X_Y_XI_Intrpl7(f10_old, F7_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf37_0, XBFLf37_1, XBFLf37_2, XBFLf37_3, XBFLf37_4, XBFLf37_5, XBFLf37_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
                X_Y_XI_Intrpl7(f9_old,  F8_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf38_0, XBFLf38_1, XBFLf38_2, XBFLf38_3, XBFLf38_4, XBFLf38_5, XBFLf38_6, YBFLf3_0,  YBFLf3_1,  YBFLf3_2,  YBFLf3_3,  YBFLf3_4,  YBFLf3_5,  YBFLf3_6, XiBFLf3_0, XiBFLf3_1, XiBFLf3_2, XiBFLf3_3, XiBFLf3_4, XiBFLf3_5, XiBFLf3_6);
            }
            //F0_in = F0_in + ModifydRho_F378( F3_in, F7_in, F8_in, f4_old[index], f9_old[index], f10_old[index] ); 
        }
        if( BFLReqF4_d[(k-3)*NYD6+j] == 1 ){
            if(Q4_d[(k-3)*NYD6+j] > 0.5){
                F4_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f3_old[index]+ ((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f4_old[index];
                F9_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f8_old[index]+((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f9_old[index];
                F10_in = (1/(2*Q4_d[(k-3)*NYD6+j]))*f7_old[index]+ ((2*Q4_d[(k-3)*NYD6+j]-1)/(2*Q4_d[(k-3)*NYD6+j]))*f10_old[index]; 
            }
            if(Q4_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(  f3_old, F4_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                //X_Y_XI_Intrpl7(f8_old, F9_in,  i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf49_0,  XBFLf49_1,  XBFLf49_2,  XBFLf49_3,  XBFLf49_4,  XBFLf49_5,  XBFLf49_6,  YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                //X_Y_XI_Intrpl7(f7_old, F10_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, XBFLf410_0, XBFLf410_1, XBFLf410_2, XBFLf410_3, XBFLf410_4, XBFLf410_5, XBFLf410_6, YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                Y_XI_Intrpl7(  f3_old, F4_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf4_0,   YBFLf4_1,   YBFLf4_2,   YBFLf4_3,   YBFLf4_4,   YBFLf4_5,   YBFLf4_6,   XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                X_Y_XI_Intrpl7(f8_old, F9_in,  i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf49_0,  XBFLf49_1,  XBFLf49_2,  XBFLf49_3,  XBFLf49_4,  XBFLf49_5,  XBFLf49_6,  YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
                X_Y_XI_Intrpl7(f7_old, F10_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, XBFLf410_0, XBFLf410_1, XBFLf410_2, XBFLf410_3, XBFLf410_4, XBFLf410_5, XBFLf410_6, YBFLf4_0,  YBFLf4_1,  YBFLf4_2,  YBFLf4_3,  YBFLf4_4,  YBFLf4_5,  YBFLf4_6, XiBFLf4_0, XiBFLf4_1, XiBFLf4_2, XiBFLf4_3, XiBFLf4_4, XiBFLf4_5, XiBFLf4_6);
            }
            //F0_in = F0_in + ModifydRho_F4910( F4_in, F9_in, F10_in, f3_old[index], f7_old[index], f8_old[index] );
        }
        if ( BFLReqF15_d[(k-3)*NYD6+j] == 1 ){
            if(Q15_d[(k-3)*NYD6+j] > 0.5){
                F15_in = (1/(2*Q15_d[(k-3)*NYD6+j]))*f18_old[index]+ ((2*Q15_d[(k-3)*NYD6+j]-1)/(2*Q15_d[(k-3)*NYD6+j]))*f15_old[index];
            }
            if(Q15_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(f18_old, F15_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf15_0, YBFLf15_1, YBFLf15_2, YBFLf15_3, YBFLf15_4, YBFLf15_5, YBFLf15_6, XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6);
                Y_XI_Intrpl7(f18_old, F15_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf15_0, YBFLf15_1, YBFLf15_2, YBFLf15_3, YBFLf15_4, YBFLf15_5, YBFLf15_6, XiBFLf15_0, XiBFLf15_1, XiBFLf15_2, XiBFLf15_3, XiBFLf15_4, XiBFLf15_5, XiBFLf15_6);
            }
            //F0_in = F0_in + ModifydRho_F15( F15_in, f18_old[index] );
        }
        if ( BFLReqF16_d[(k-3)*NYD6+j] == 1 ){
            if(Q16_d[(k-3)*NYD6+j] > 0.5){
                F16_in = (1/(2*Q16_d[(k-3)*NYD6+j]))*f17_old[index]+ ((2*Q16_d[(k-3)*NYD6+j]-1)/(2*Q16_d[(k-3)*NYD6+j]))*f16_old[index];
            }
            if(Q16_d[(k-3)*NYD6+j] < 0.5) {
                //Y_XI_Intrpl7(f17_old, F16_in, i, j, k, (i-3), (j-3), 3, idx_xi, idx_xi, idx_xi, YBFLf16_0, YBFLf16_1, YBFLf16_2, YBFLf16_3, YBFLf16_4, YBFLf16_5, YBFLf16_6, XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6);
                Y_XI_Intrpl7(f17_old, F16_in, i, j, k, (i-3), (j-3), 3, i, j, idx_xi, YBFLf16_0, YBFLf16_1, YBFLf16_2, YBFLf16_3, YBFLf16_4, YBFLf16_5, YBFLf16_6, XiBFLf16_0, XiBFLf16_1, XiBFLf16_2, XiBFLf16_3, XiBFLf16_4, XiBFLf16_5, XiBFLf16_6);
            }
            //F0_in = F0_in + ModifydRho_F16( F16_in, f17_old[index] );
        }
    }





























double Lagrange_6th(//x被內插的點，x_i:對硬疊加指標，x1第一個連乘點
    double x,   double x_i, double x1,  double x2,  double x3,  double x4,  double x5,  double x6)//x1 : 從左邊數過來第一個被乘的
{
    double Para = (x - x1)/(x_i - x1)*(x - x2)/(x_i - x2)*(x - x3)/(x_i - x3)*(x - x4)/(x_i - x4)*(x - x5)/(x_i - x5)*(x - x6)/(x_i - x6);

    return Para;
}

void GetParameter_6th(
    double *Para_h[7],      double Position, //我一旦宣告了 Para_h[0]作為一個指標變數，那麼基於，Para_h[0]為Para_h[0][0]的初始存放位址，直接使用Para[0~6][?]的矩陣元素
    double *Pos,            int i,              int n  )
{   //Para_h[0+?][i] 對應到 Pos[n+?]
    Para_h[0][i] = Lagrange_6th(Position, Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);//Pos[n+0]對應Para_h[0][i]//Pos[n]=x_i//Pos[n+1]=x_1//....
    Para_h[1][i] = Lagrange_6th(Position, Pos[n+1], Pos[n],   Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[2][i] = Lagrange_6th(Position, Pos[n+2], Pos[n],   Pos[n+1], Pos[n+3], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[3][i] = Lagrange_6th(Position, Pos[n+3], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+4], Pos[n+5], Pos[n+6]);
    Para_h[4][i] = Lagrange_6th(Position, Pos[n+4], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+5], Pos[n+6]);
    Para_h[5][i] = Lagrange_6th(Position, Pos[n+5], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+6]);
    Para_h[6][i] = Lagrange_6th(Position, Pos[n+6], Pos[n],   Pos[n+1], Pos[n+2], Pos[n+3], Pos[n+4], Pos[n+5]);//但是他沒有宣告Para_h[0][?]的實際尺寸，即便Para_h[0]可以做為Para_h[0][0]的初始存放位址，我憑什麼也宣告了i，就可以使用Para_h[0][i]Para_h[0][i]這裡的I又是什麼意思，陣列尺寸還是第i個元素?A:信任機制ˋ
}//問題:i怎麼做使用 




void GetIntrplParameter_X() {
   
    for( int i = 3; i < NX6-3; i++ ){//包住整個X方向計算區域
        GetParameter_6th( XPara0_h, x_h[i]-minSize, x_h, i, i-3 );
        GetParameter_6th( XPara2_h, x_h[i]+minSize, x_h, i, i-3 );
    }

    for( int i = 0; i < 7; i++ ){
        CHECK_CUDA( cudaMemcpy(XPara0_d[i], XPara0_h[i], NX6*sizeof(double), cudaMemcpyHostToDevice) );
        CHECK_CUDA( cudaMemcpy(XPara2_d[i], XPara2_h[i], NX6*sizeof(double), cudaMemcpyHostToDevice) );
    }
    CHECK_CUDA( cudaDeviceSynchronize() );
}


XPara0_h為指標的指標  型別 double**
定義: &XPara0_h[0]
說明:為指標陣列(double* XPara0_h[7(尺寸)])的初始存放位址
特色:h代表 host 主機端(CPU端)




筆記起來

    nBytes = NX6 * sizeof(double);
    AllocateHostArray(  nBytes, 4,  &x_h, &Xdep_h[0], &Xdep_h[1], &Xdep_h[2]);//x_h的記憶體配置 作為一個表示格點座標的矩陣
    AllocateDeviceArray(nBytes, 4,  &x_d, &Xdep_d[0], &Xdep_d[1], &Xdep_d[2]);
    for( int i = 0; i < 7; i++ ){
        CHECK_CUDA( cudaMallocHost( (void**)&XPara0_h[i], nBytes ) );
        CHECK_CUDA( cudaMallocHost( (void**)&XPara2_h[i], nBytes ) );
        CHECK_CUDA( cudaMalloc( &XPara0_d[i], nBytes ) );
        CHECK_CUDA( cudaMalloc( &XPara2_d[i], nBytes ) );
    }




double GetDelta45Degree(
    const double z_pnt,      const double y_pnt,
    const double y_zdominate,const double y_ydominate  )//y_ydominate由y距離主導的y座標邊界 , y_zdominate由z距離主導的y座標邊界
{
    //初始化二分法上下界 
    double y_temp[2] = {y_zdominate, y_ydominate};//需告一個尺寸為2的陣列 
    double y_mid;
    int a = 0;  // 迭代计数器（当前未使用）

    // 二分法循环，寻找垂直距离与水平距离相等的 y_mid 点
    do{
        y_mid = (y_temp[0]+y_temp[1]) / 2.0; //先指定計算終點為上下界平均 

        // *** 错误警告 ***: 外层的 fabs() 作用于布尔值，这是不正确的
        // 应该改为: if( fabs(HillFunction(y_mid) - z_pnt) > fabs(y_mid - y_pnt) )
        // 当前逻辑意图: 如果垂直距离 > 水平距离，则更新下界
        if( fabs( fabs(HillFunction(y_mid) - z_pnt) > fabs(y_mid-y_pnt) ) ){
            y_temp[0] = y_mid;  // 向 y_ydominate 方向移动
        } else {
            y_temp[1] = y_mid;  // 向 y_zdominate 方向移动
        }
        a++;
    } while ( fabs( fabs(HillFunction(y_mid)-z_pnt) - fabs(y_mid-y_pnt) ) > 1e-13 );//你想要的相反作為判斷條件 

    // 计算调整后的网格间距
    // minSize 是允许的最小网格间距
    // 因子 2.0*fabs(y_pnt-y_mid) 表示从查询点到45度过渡点的距离
    double d = minSize - 2.0*fabs(y_pnt-y_mid);
    return d;
}

double GetDelta45Degree(
    const double z_pnt,      const double y_pnt,
    const double y_zdominate,const double y_ydominate  )
{
    double y_temp[2] = {y_zdominate, y_ydominate};
    double y_mid;
    int a = 0;
    do{
        y_mid = (y_temp[0]+y_temp[1]) / 2.0; //為(y_h[i]+y_h[i]+minSize)/2.0
        if( fabs(HillFunction(y_mid) - z_pnt) > fabs(y_mid-y_pnt) ){
            y_temp[0] = y_mid; //y_z = y_mid
        } else {
            y_temp[1] = y_mid;
        }
        a++;
    } while ( fabs( fabs(HillFunction(y_mid)-z_pnt) - fabs(y_mid-y_pnt) ) > 1e-13 ); //跳出迴圈判據地相反 
    double d = minSize - 2.0*fabs(y_pnt-y_mid); //若d<0代表，y目前離45度山丘遠，反之則近 
    return d;
}

//minSize為最小特徵尺度 



//超級重要， 
F0_Intrpl7( f0_old,  i, j, k);
    F1_Intrpl7( f1_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6 );
    F2_Intrpl7( f2_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6 );
    F3_Intrpl7( f3_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    F4_Intrpl7( f4_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    F5_Intrpl7( f5_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    F6_Intrpl7( f6_old,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    //7,8,9,10需要往三個方向內插
    X_Y_XI_Intrpl7( f7_old, F7_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    X_Y_XI_Intrpl7( f8_old, F8_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF3_0, XiF3_1, XiF3_2, XiF3_3, XiF3_4, XiF3_5, XiF3_6 );
    X_Y_XI_Intrpl7( f9_old, F9_in,  i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    X_Y_XI_Intrpl7(f10_old, F10_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF4_0, XiF4_1, XiF4_2, XiF4_3, XiF4_4, XiF4_5, XiF4_6 );
    //x_z平面方向之分量:先對X軸做內插，再對Z軸做內插
    X_XI_Intrpl7(  f11_old, F11_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    X_XI_Intrpl7(  f12_old, F12_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   XiF5_0, XiF5_1, XiF5_2, XiF5_3, XiF5_4, XiF5_5, XiF5_6 );
    X_XI_Intrpl7(  f13_old, F13_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0,   X0_1,   X0_2,   X0_3,   X0_4,   X0_5,   X0_6,   XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    X_XI_Intrpl7(  f14_old, F14_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X2_0,   X2_1,   X2_2,   X2_3,   X2_4,   X2_5,   X2_6,   XiF6_0, XiF6_1, XiF6_2, XiF6_3, XiF6_4, XiF6_5, XiF6_6 );
    //y_z平面方向之分量:先對Y軸做內插，再對Z軸做內插
    Y_XI_Intrpl7(  f15_old, F15_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF15_0, XiF15_1, XiF15_2, XiF15_3, XiF15_4, XiF15_5, XiF15_6 );
    Y_XI_Intrpl7(  f16_old, F16_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF16_0, XiF16_1, XiF16_2, XiF16_3, XiF16_4, XiF16_5, XiF16_6 );
    Y_XI_Intrpl7(  f17_old, F17_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y0_0,   Y0_1,   Y0_2,   Y0_3,   Y0_4,   Y0_5,   Y0_6,   XiF17_0, XiF17_1, XiF17_2, XiF17_3, XiF17_4, XiF17_5, XiF17_6 );
    Y_XI_Intrpl7(  f18_old, F18_in, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, Y2_0,   Y2_1,   Y2_2,   Y2_3,   Y2_4,   Y2_5,   Y2_6,   XiF18_0, XiF18_1, XiF18_2, XiF18_3, XiF18_4, XiF18_5, XiF18_6 );
    

F1_Intrpl7 函數參數說明
根據第 13-15 行的定義：

#define F1_Intrpl7(f, i, j, k, i_c, j_c, k_c, XPara_0, XPara_1, XPara_2, XPara_3, XPara_4, XPara_5, XPara_6)
參數列表（共 14 個參數）：
f - 要進行內插的分佈函數陣列（例如：f1_old）為非物理空間計算點的碰撞前一般態分佈函數
i - X方向的當前格點索引
j - Y方向的當前格點索引
k - Z方向的當前格點索引
i_c - X方向的計算起始索引（通常是 i-3）
j_c - Y方向的計算起始索引（通常是 j-3）
k_c - Z方向的單元格索引（例如：cell_z 或 3）
XPara_0 - 第0個Lagrange內插參數陣列（對應X方向）
XPara_1 - 第1個Lagrange內插參數陣列
XPara_2 - 第2個Lagrange內插參數陣列
XPara_3 - 第3個Lagrange內插參數陣列
XPara_4 - 第4個Lagrange內插參數陣列
XPara_5 - 第5個Lagrange內插參數陣列
XPara_6 - 第6個Lagrange內插參數陣列
使用範例（來自您的代碼第211行）：

F1_Intrpl7(f1_old, i, j, k, (i-3), (j-3), cell_z, i, j, idx_xi, X0_0, X0_1, X0_2, X0_3, X0_4, X0_5, X0_6);
填入說明：
f1_old - 舊的分佈函數 f1
i, j, k - 當前位置
(i-3), (j-3) - 內插起始位置（向前偏移3個格點）
cell_z - Z方向參數
i - X方向索引（用於取 XPara 陣列值）
j - Y方向索引
idx_xi - Xi方向索引
X0_0 到 X0_6 - 7個Lagrange 6階內插係數陣列
功能：
這個巨集執行 6階 Lagrange 多項式內插，在X方向上使用7個點（i-3 到 i+3）進行內插計算，結果存入 F1_in 變數中。


