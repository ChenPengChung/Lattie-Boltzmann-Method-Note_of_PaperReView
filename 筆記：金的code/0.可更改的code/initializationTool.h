#ifndef INITIALIZATIONTOOL_FILE
#define INITIALIZATIONTOOL_FILE

#define tanhFunction( L, LatticeSize, a, j, N )     \
(           \
    L/2.0 + LatticeSize/2.0 + ((L/2.0)/a)*tanh((-1.0+2.0*(double)(j)/(double)(N))/2.0*log((1.0+a)/(1.0-a)))     \
)

//上述為座標系統之映射，為甚麼要多加LatticeSize/2，因為要讓j所對應的物理空間計算點在網格中心點
//zi 的範圍 : [-L/2 ; L/2]
//\xi 的範圍 : [-1 ; 1]
//j的範圍 : [0 , N/2] 
//透過平移，讓真實Z座標從[0+LatticeSize/2]開始。

//找到最佳比例因子
double GetNonuniParameter() {
    double total = LZ - HillFunction( 0.0 ) - minSize;
    double a_temp[2] = {0.1, 1.0};
    double a_mid;

    double x_temp[2], dx;
    do{
        a_mid = (a_temp[0]+a_temp[1]) / 2.0;
        x_temp[0] = tanhFunction(total, minSize, a_mid, 0, (NZ6-7));//他是放NZ6-7不是放入NZ6-6，但是Z方向網格數是NZ6-6，他放NZ6-7的意義在於讓真正的網格數量j 一共有NZ6-6個。
        x_temp[1] = tanhFunction(total, minSize, a_mid, 1, (NZ6-7));
        dx = x_temp[1] - x_temp[0];
        if( dx - minSize >= 0.0 ){
            a_temp[0] = a_mid;
        } else {
            a_temp[1] = a_mid;
        }
    } while ( fabs( dx - minSize ) > 1e-14 );

    /*if( myid == 0 ){
        printf("a = %lf\n", a_mid);
    }*/

    return a_mid;
}

double Lagrange_2nd(
    double x,   double x_i,
    double x1,  double x2  )
{
    double Para = (x - x1)/(x_i - x1)*(x - x2)/(x_i - x2);

    return Para;
}

void GetParameter_2nd(
    double *Para_h[7],      double Position,
    double *Pos,            int i,              int n  )
{
    Para_h[0][i] = Lagrange_2nd(Position, Pos[n],   Pos[n+1], Pos[n+2]);
    Para_h[1][i] = Lagrange_2nd(Position, Pos[n+1], Pos[n]  , Pos[n+2]);
    Para_h[2][i] = Lagrange_2nd(Position, Pos[n+2], Pos[n]  , Pos[n+1]);
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

int IsBFLBCNeeded(const double Y, const double Z) {//真正的BFL條件判斷因子，
    const double hill = HillFunction( Y );
    if( hill > Z ){
        return 1;
    } else {
        return 0;
    }
}

double GetDeltaHorizontal(
    const double z_target,
    const double y_large,       const double y_small,       const double y_point )
{
    double y_temp[2] = {y_large, y_small};
    double y_mid;
    do{
        y_mid = (y_temp[0]+y_temp[1]) / 2.0;
        if( HillFunction(y_mid) >= z_target ){
            y_temp[0] = y_mid;
        } else {
            y_temp[1] = y_mid;
        }
    } while ( fabs( HillFunction(y_mid) - z_target ) > 1e-13 );

    double d = minSize - 2.0*fabs(y_point-y_mid);
    return d;
}

// 计算45度距离约束活跃点的网格间距
// 该函数在山丘表面搜索一个点，在该点处距离度量从垂直(z)距离主导
// 过渡到水平(y)距离主导
//
// 参数说明:
//   z_pnt, y_pnt        : 查询点的坐标 (y_pnt, z_pnt)
//   y_zdominate         : z距离主导的y边界（垂直约束活跃）
//   y_ydominate         : y距离主导的y边界（水平约束活跃）
//
// 算法原理:
//   使用二分法查找山丘表面上的 y_mid 点，使得:
//   |HillFunction(y_mid) - z_pnt| = |y_mid - y_pnt|
//   这是从点 (y_pnt, z_pnt) 出发的45度线与山丘表面的交点
//
// 返回值:
//   根据45度约束调整后的网格间距 delta
double GetDelta45Degree(
    const double z_pnt,      const double y_pnt,
    const double y_zdominate,const double y_ydominate  )
{
    // 初始化二分法的搜索区间
    double y_temp[2] = {y_zdominate, y_ydominate};
    double y_mid;
    int a = 0;  // 迭代计数器（当前未使用）

    // 二分法循环，寻找垂直距离与水平距离相等的 y_mid 点
    do{
        y_mid = (y_temp[0]+y_temp[1]) / 2.0;

        // *** 错误警告 ***: 外层的 fabs() 作用于布尔值，这是不正确的
        // 应该改为: if( fabs(HillFunction(y_mid) - z_pnt) > fabs(y_mid - y_pnt) )
        // 当前逻辑意图: 如果垂直距离 > 水平距离，则更新下界
        if( fabs( fabs(HillFunction(y_mid) - z_pnt) > fabs(y_mid-y_pnt) ) ){
            y_temp[0] = y_mid;  // 向 y_ydominate 方向移动
        } else {
            y_temp[1] = y_mid;  // 向 y_zdominate 方向移动
        }
        a++;
    } while ( fabs( fabs(HillFunction(y_mid)-z_pnt) - fabs(y_mid-y_pnt) ) > 1e-13 );

    // 计算调整后的网格间距
    // minSize 是允许的最小网格间距
    // 因子 2.0*fabs(y_pnt-y_mid) 表示从查询点到45度过渡点的距离
    double d = minSize - 2.0*fabs(y_pnt-y_mid);
    return d;
}

void GetBFLXiParameter(
    double *XiPara_h[7],    double pos_z,       double pos_y,
    double *Pos_xi,         int IdxToStore,     int k  )
{
    double L = LZ - HillFunction(pos_y) - minSize;
    double pos_xi = LXi * (pos_z - (HillFunction(pos_y)+minSize/2.0)) / L;
    double a = GetNonuniParameter();
    //double pos_xi = atanh((pos_z-HillFunction(pos_y)-minSize/2.0-L/2.0)/((L/2.0)/a))/log((1.0+a)/(1.0-a))*2.0;

    if( k >= 3 && k <= 4 ){
        GetParameter_6th( XiPara_h, pos_xi, Pos_xi, IdxToStore, 3 );
    }
}


#endif
