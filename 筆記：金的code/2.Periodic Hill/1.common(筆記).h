1.struct timeval结构体
//https://blog.csdn.net/lyc_daniel/article/details/11733715 
//该结构体是Linux系统中定义，struct timeval结构体在time.h中的定义为：
struct timeval {
    __time_t tv_sec ; //計算＿秒
    __suseconds_t tv_usec ; //計算_毫秒 
}

2.#ifndef #define #enif : 
#ifndef : "if not define" 的意思：當當前文件還沒被編譯過，就執行以下#define .....論述 ， 
防止定義被重複
#ifdef : "if define " 的意思：當前文件如果被編譯過。就執行以下 論述

#enif :結束定義：


3.inline : //想把 定義寫在header裡面就要加入inline 
cudaError_t: CUDA API 回傳的錯誤碼型別，描述每次 CUDA 呼叫的結果狀態。
cudaGetErrirString(errot) : 將 cudaError_t 轉成可讀的錯誤訊息字串，方便日誌或輸出。
checkCuda(call,__FILE__, __LINE__) : 常見包裝/巨集，執行 CUDA 呼叫並在失敗時列出檔名與行號。
__FILE__ : 編譯器內建巨集，展開為當前檔案名。
__LINE__ : 編譯器內建巨集，展開為當前行號。
默寫： inline void checkCuda(cudaError_t error, const char* file , int line)
{
    if(error!= cudaSSucess){
        cerr << "CUDA Error at : " file << ":" << line 
             << "-Code:" << error
             << ",Reason" << cudaGetErrorString(error) << emdl ;
            exit(1) ; //exit(1) 會立即停止整個程式，不管在哪個函數裡呼叫
    //- 回傳 1 給作業系統，表示程式因錯誤而結束
    //- 在你選取的那行，應該是在檢測到 CUDA 錯誤時使用
    }
}

#define CHECK_CUDA(call) checkCUDA(call,__FILE__,__LINE__) 

4.#define CHECH_MPI: 
(1). 為什麼 call 沒辦法宣告：
(2). mpi_error_string[MPI_MAX_ERROR_STRING];
(3). int mpi_error_string_length = 0; 是什麼？
(4). MPI_Error_string(mpi_status, mpi_error_string, &mpi_error_string_length);
(5). mpi_error_string : 
(6). #call : 
(7). __LINE__ : 編譯器內巨集，展開為當前行號 
(8). __FILE__ : 編譯器內建巨集，展開為當前檔案名。
(9). mpi_error_string : 
(10). mpi_status : vuhu3u3