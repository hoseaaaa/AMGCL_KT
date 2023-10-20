// #define BLOCK_MATRIX
// 
#define DEBUG

#ifdef DEBUG
#define SIM 1
#else
#define SIM 400
#endif
// #define SCALE_DIAGONAL

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

#include <iomanip>
#include<math.h>

#include <iostream>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai1.hpp>

#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/solver/idrs.hpp>

#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/relaxation/chebyshev.hpp>

#include <amgcl/solver/preonly.hpp>
#include <amgcl/coarsening/aggregation.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>

#include <amgcl/solver/idrs.hpp>
#include <amgcl/relaxation/damped_jacobi.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>

#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/relaxation/iluk.hpp>
#include <amgcl/relaxation/ilup.hpp>
#include <amgcl/relaxation/ilut.hpp>

#include <amgcl/solver/bicgstabl.hpp>
#include <amgcl/solver/lgmres.hpp>


#ifdef  BLOCK_MATRIX
#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#endif

#include <amgcl/value_type/static_matrix.hpp>
#include <amgcl/adapter/block_matrix.hpp>
#include <amgcl/make_block_solver.hpp>


#include <amgcl/io/mm.hpp>
#include <amgcl/profiler.hpp>

#include<tuple>


#include <amgcl/relaxation/as_preconditioner.hpp>


using namespace std ; 

void line2data(string s ,vector<ptrdiff_t> &a){
    int len_s = s.length();
	int i=0, j=0;
	while (i < len_s)
	{
		if (s[i] >= '0'&& s[i] <= '9')
		{
			j = i;
			int len = 0;
			while (s[i] >= '0'&& s[i] <= '9')
			{
				i++;
				len++;
			}
			string s0 = s.substr(j, len);//获取子串
            int num=0;//数字字符串转换为整型数字
			stringstream s1(s0);
			s1 >> num;
			a.push_back(num);
		}
		else
		{
			i++;
		}
	}
}
template <class I, class T>
void read_c_g_b( vector < I>& cout_point,
                   vector < I>& A_ptr,
                   vector < I>& A_col,
                   vector < T>& A_val,
                   string filename
                 ) {
	ifstream infile;
	infile.open(filename);   //将文件流对象与文件连接起来 
	if (!infile.is_open()){
        cout <<"error  open file " << filename <<endl ;   //若失败,则输出错误消息,并终止程序运行 
        exit(1) ;
    }
    string s ;
    //0 、
    if (getline(infile, s) ){
        line2data(s,cout_point);
    }

    //1、
    if (getline(infile, s) ){
        istringstream is(s);
        long int ptr ;
        while (is >> ptr) {
			A_ptr.push_back(ptr);
		}
    }
    //2、
    if (getline(infile, s) ){
        istringstream is(s);
        long int  col;
        while (is >> col) {
			A_col.push_back(col);
		}
    }
    //3、
    if (getline(infile, s) ){
        istringstream is(s);
        double  val;
        while (is >> val) {
			A_val.push_back(val);
		} 
    }
	infile.close();             //关闭文件输入流 
    // pan duan 
    if( ( A_ptr.size()-1 != cout_point[0] ) &&
         ( A_col.size()  != cout_point[2] ) && 
         ( A_val.size()  != cout_point[2] )  
        ){
            cout <<"error read file " << filename  <<" matrix size not right  "<<endl ;   //若失败,则输出错误消息,并终止程序运行 
            exit(1) ;
        }
}

template <class I, class T>
void csc_tocsr(const I n_row,
 	           const I n_col, 
               const vector< I>& Ap ,
               const vector< I>& Aj ,
	           const vector < T>& Ax,
	                 vector < I>& Bp,
	                 vector < I>& Bi,
	                 vector < T>& Bx)
{  
    const I nnz = Ap[n_row];

    //compute number of non-zero entries per column of A 
    // std::fill(Bp, Bp + n_col, 0);
    for (I n = 0; n < nnz; n++){            
        Bp[Aj[n]]++;
    }
    //cumsum the nnz per column to get Bp[]
    for(I col = 0, cumsum = 0; col < n_col; col++){     
        I temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz; 

    for(I row = 0; row < n_row; row++){
        for(I jj = Ap[row]; jj < Ap[row+1]; jj++){
            I col  = Aj[jj];
            I dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];
            Bp[col]++;
        }
    }  
    for(I col = 0, last = 0; col <= n_col; col++){
        I temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }
}   

template <class val>
void print_vector(const vector <val>& A){
    for ( auto i=A.begin() ; i< A.end() ; i++ ){
        cout  << * i << setw(20); 
    }
    cout << endl ;
}

template <class I, class T>
void read_u_t(     vector < I>& cout_point,
                   vector < T>& A_val,
                   I cols_0 ,
                   I cols_1 ,
                   string filename
                 ) {
	ifstream u_t_file;
	u_t_file.open(filename);   //将文件流对象与文件连接起来 

    if(!u_t_file.is_open()){
        cout << "error  open file " <<endl ;
        exit (1) ; 
    }

	string s;
    getline(u_t_file, s) ;  // matrix size
    line2data(s,cout_point) ;
    getline(u_t_file, s) ;  // drop  time step  

	while (getline(u_t_file, s)) {
		istringstream is(s); //将读出的一行转成数据流进行操作
		T d;
		while (!is.eof()) {
                int pam = 0 ; 
			    while (is >> d){
                    if ( pam == cols_0 ){
                        A_val.push_back(d);
                    }
                    if ( pam == cols_1){
                        A_val.at(A_val.size()-1) += d ;
                    }
                    pam++ ;
                }    
		}
		s.clear();
	}
	u_t_file.close();             //关闭文件输入流 
} 

int main(int argc ,char** argv) {    

    // if (argc < 4) {
    //     std::cerr << "Usage: " << argv[0] << " <matrix_C.mtx> <matrix_G.mtx> <matrix_B.mtx> <matrix_U.mtx> " << std::endl;
    //     return 1;
    // }
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <matrix_G.mtx> <matrix_U.mtx> " << std::endl;
        return 1;
    }
    // The profiler:
    
    double  H_step = 2.5e-11;

    amgcl::profiler<> prof("Ax=b");   
    cout.setf(std::ios::left);
    // Read the system matrix and the RHS:  csr
    // ptrdiff_t rows_C, cols_C ,nnz_C;
    // std::vector<ptrdiff_t> ptr_C, col_C ; 
    // std::vector<double> val_C ;

    ptrdiff_t rows_G, cols_G ,nnz_G;
    std::vector<ptrdiff_t> ptr_G,col_G ; 
    std::vector<double> val_G ;


    // ptrdiff_t rows_B, cols_B , nnz_B;
    // std::vector<ptrdiff_t> ptr_B, col_B ; 
    // std::vector<double> val_B ;

    ptrdiff_t cols_U,rows_U;
    std::vector<ptrdiff_t> ptr_U,col_U ; 
    std::vector<double> val_U  ;

    vector<ptrdiff_t>   cout_point ; 

    prof.tic("read");
    //read C_G
    // read_c_g_b<ptrdiff_t,double> (cout_point,ptr_C,col_C,val_C,argv[1]) ; 
    // rows_C = cout_point[0] ;    cols_C = cout_point[1] ;    nnz_C = cout_point[2] ;
    // cout_point.clear() ;

    read_c_g_b<ptrdiff_t ,double> (cout_point,ptr_G,col_G,val_G,argv[1]) ; 
    rows_G = cout_point[0] ;    cols_G = cout_point[1] ;    nnz_G = cout_point[2] ;
    cout_point.clear() ;

    //read b
    // read_c_g_b<ptrdiff_t ,double> (cout_point,ptr_B,col_B,val_B,argv[2]) ; 
    // rows_B = cout_point[0] ;    cols_B = cout_point[1] ;    nnz_B = cout_point[2] ;
    // cout_point.clear() ;

    // vector<ptrdiff_t >B_ptr_temp (rows_B+1); 
    // vector<ptrdiff_t >B_col_temp(nnz_B); 
    // vector<double> B_val_temp(nnz_B) ;

    // csc_tocsr<ptrdiff_t ,double>( cols_B  , rows_B  , 
    //                              ptr_B , col_B , val_B,
    //                              B_ptr_temp , B_col_temp , B_val_temp  ) ;
    // ptr_B.swap(B_ptr_temp)  ;
    // col_B.swap(B_col_temp)  ;
    // val_B.swap(B_val_temp)  ;

#ifdef DEBUG
    //print C 
    // std::cout << "Matrix C" << argv[1] << ": " << rows_C << "x" << cols_C << std::endl;
    //print G
    std::cout << "Matrix G" << argv[1] << ": " << rows_G << "x" << cols_G << std::endl;

    //print B
    // std::cout << "Matrix B" << argv[3] << ": " << rows_B << "x" << cols_B << std::endl;

    //print u
    // std::cout << " u  " << argv[4] << ": " << rows_U << "x" << cols_U << std::endl; 
#endif

    prof.toc("read");
    
    //calculate
    prof.tic("calculate");

    // auto C = std::tie(rows_C, ptr_C, col_C, val_C);
    auto G = std::tie(rows_C, ptr_G, col_G, val_G);

    // C sparse matrix +  C sparse matrix  C/h + G
    // auto sum_A = amgcl::backend::sum<double ,ptrdiff_t,ptrdiff_t>(1.0, C , H_step /2 , G ,true);  // may be change 
    auto sum_A  = G ;
    // auto sum_X_B = amgcl::backend::sum<double ,ptrdiff_t,ptrdiff_t>(1.0, C , -H_step /2 , G ,true);  // may be change 

    // C sparse matrix  * rhs dense   B * u 
    // auto B = std::tie(rows_B, ptr_B, col_B, val_B);

    prof.toc("calculate");

//////////////////////
//     //solver

    typedef amgcl::backend::builtin<double> SBackend;
    typedef amgcl::backend::builtin<double> PBackend;

    typedef amgcl::make_solver<
        amgcl::amg<
            PBackend,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai1
            >,
        amgcl::solver::bicgstabl<SBackend>
        > Solver;

    Solver::params prm;

    prm.solver.L = 4 ;
    prm.solver.tol = 2.5e-3 ;
    prm.solver.maxiter = 1000 ;

    // prm.solver.s = 2 ;
    // prm.solver.smoothing =true ;
#ifdef  DEBUG
    prm.solver.verbose = true ;
#endif
    prm.precond.direct_coarse  = true ; 
    prm.precond.coarsening.relax = 1.5f ;
    prm.precond.coarsening.power_iters = 1000 ; 
    prm.precond.coarsening.aggr.eps_strong = 0.0020 ;
    prm.precond.coarsening.estimate_spectral_radius = true ;

    // prm.precond.npre = 1 ;
    prm.precond.max_levels = 10;

    // may be need 
    // The -s option tells the solver to do that:
    
    vector<ptrdiff_t>  A_PTR (sum_A->ptr,sum_A->ptr+rows_C+1); 
    vector<ptrdiff_t>  A_COL (sum_A->col,sum_A->col + sum_A->ptr[rows_C] ) ; 
    vector<double>     A_VAL (sum_A->val,sum_A->val + sum_A->ptr[rows_C] ) ; 
    auto A = std::tie(rows_C,A_PTR ,A_COL ,A_VAL )  ;

    // Initialize the solver with the system matrix:
    prof.tic("setup");

    Solver solve(A,prm);

#ifdef DEBUG
        std::cout << solve  <<  std::endl;
#endif
    int iters;
    double error;
    // Show the mini-report on the constructed solver:
    // 设置阶段完成 。。尝试修改b


#ifndef DEBUG

#endif
    prof.toc("setup");

    std::vector<double> x(rows_C,0.0); 
    prof.tic("solve");
    ofstream outfile;
    outfile.open(argv[2]) ;
    for (auto i = 0 ;i < SIM ;i++ ){

        read_u_t<ptrdiff_t,double>(cout_point,val_U,i,i+1,argv[2]) ;  // read u0 
        rows_U = cout_point[0] ;    cols_U = cout_point[1];

        // sum_X_B * x(t)
        std::vector<double> new_val_C_x(rows_B);
        amgcl::backend::spmv( 1.0 , *sum_X_B , x, 0.0 , new_val_C_x);

        // B *u (t) + sum_X_B* x(t)

        amgcl::backend::spmv(  H_step/2, B , val_U , 1.0 , new_val_C_x);

        val_U.swap( new_val_C_x) ;

        // vector add vector 
        x.clear() ;
        std::tie(iters, error) = solve(A,val_U, x);

        // std::tuple<ptrdiff_t &, std::vector<ptrdiff_t> &, std::vector<ptrdiff_t> &, std::vector<double> &> C
        // Output the number of iterations, the relative error,
        // and the profiling data:

#ifndef DEBUG
       	if (!outfile.is_open()){
            cout <<"error  open file outx_data.txt"<<endl ;   //若失败,则输出错误消息,并终止程序运行 
            exit(1) ;
        }
        vector <ptrdiff_t> cout_sta ( cout_point.begin()+2,cout_point.end() ) ;
        if (i == 0 ) {
            outfile << "Time index:     \t"  <<setw(20); 
            for (auto cout_i = cout_sta.begin() ; cout_i<cout_sta.end() ; cout_i++){
                outfile << *cout_i <<setw(20);
            }
            outfile << endl ;
            for (auto i = 0 ;i <=cout_sta.size() ; i++ ){
                 outfile << "0.0 " <<setw(20) ; 
            }
            outfile << endl ;
        }
        outfile <<  2.5e-11 * (i+1) << setw(20); 
        for (auto i_cout = 0 ; i_cout< cout_sta.size() ; i_cout++ ){
            outfile <<   ( x [ cout_sta.at(i_cout) ]  ) <<setw(20);
        }
        outfile <<endl ;
        cout_sta.clear() ;
#endif
        val_U.clear() ;
        cout_point.clear() ;
    }
        outfile.close() ;
        prof.toc("solve");

        std::cout << "Iters: " << iters << std::endl
             << "Error: " << error << std::endl ;


    std::cout << prof << std::endl;
    return 0 ;
}
