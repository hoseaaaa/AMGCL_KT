#define DEBUG

#ifdef DEBUG
#define SIM 1
#else
#define SIM 400
#endif

#ifndef EPS_STRONG
#define EPS_STRONG 1.0
#endif

#ifndef RELAX
#define RELAX 0.08
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
#include <amgcl/relaxation/spai1.hpp>

#include <amgcl/relaxation/ilut.hpp>

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

#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/gauss_seidel.hpp>

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

#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/coarsening/smoothed_aggr_emin.hpp>
#include <amgcl/coarsening/pointwise_aggregates.hpp>
#include <amgcl/coarsening/ruge_stuben.hpp>
#include  <amgcl/coarsening/as_scalar.hpp>

#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/relaxation/chebyshev.hpp>
#include <amgcl/relaxation/detail/ilu_solve.hpp>


#include <amgcl/solver/richardson.hpp>
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
void read_c_g_b(   vector <I>& cout_point,
                   vector <I>& A_ptr,
                   vector <I>& A_col,
                   vector <T>& A_val,
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
    cout << "------s-----"  << s <<endl; 
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
            cout << " A_ptr.size : " <<  A_ptr.size()   << "cout_point[0]:  " << cout_point[0] << endl ;
            cout << " A_col.size : " <<  A_col.size()   << "cout_point[2]:  " << cout_point[2] << endl ;
            cout << " A_val.size : " <<  A_val.size()   << "cout_point[2]:  " << cout_point[2] << endl ;
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
    // getline(u_t_file, s) ;  // drop  time step  

	while (getline(u_t_file, s)) {
		istringstream is(s); //将读出的一行转成数据流进行操作
		T data;
		while (!is.eof()) {
                int pam = 0 ; 
			    while (is >> data){
                    if ( pam == cols_0 ){
                        A_val.push_back(data);
                    }
                    if ( pam == cols_1){
                        A_val.at(A_val.size()-1) += data ;
                    }
                    pam++ ;
                }    
		}
		s.clear();
	}
	u_t_file.close();             //关闭文件输入流 
} 

int main(int argc ,char** argv) {    
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <matrix_G.mtx> <matrix_Ut.mtx> " << std::endl;
        //argv[1]   <matrix_G.mtx>
        //argv[2]   <matrix_Ut.mtx>
        return 1;
    }
    
    double  H_step = 2.5e-11;
    amgcl::profiler<> prof("Ax=b");   
    cout.setf(std::ios::left);

    ptrdiff_t rows_G, cols_G ,nnz_G;
    std::vector<ptrdiff_t> ptr_G,col_G ; 
    std::vector<double> val_G ;

    ptrdiff_t cols_U,rows_U;
    std::vector<ptrdiff_t> ptr_U,col_U ; 
    std::vector<double> val_U  ;
    vector<ptrdiff_t>   cout_point ; 

    prof.tic("read");

    read_c_g_b<ptrdiff_t ,double> (cout_point,ptr_G,col_G,val_G,argv[1]) ; 
    rows_G = cout_point[0] ;    cols_G = cout_point[1] ;    nnz_G = cout_point[2] ;
    cout_point.clear() ;
#ifdef DEBUG
    std::cout << "Matrix G" << argv[1] << ": " << rows_G << "x" << cols_G << std::endl;
#endif
    prof.toc("read");
    prof.tic("calculate");
    auto G = std::tie(rows_G, ptr_G, col_G, val_G);
    auto sum_A  = G ;
    prof.toc("calculate");
//////////////////////
//     //solver

    typedef amgcl::backend::builtin<float> SBackend;
    typedef amgcl::backend::builtin<double> PBackend;

    typedef amgcl::make_solver<
        amgcl::amg<
            PBackend,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::gauss_seidel
            >,
        amgcl::solver::cg<SBackend>
        > Solver;
    // typedef amgcl::make_solver<
    //     amgcl::make_solver<
    //         amgcl::amg<
    //             PBackend,
    //             amgcl::coarsening::smoothed_aggregation,
    //             amgcl::relaxation::damped_jacobi
    //             >,
    //         amgcl::solver::cg<SBackend>
    //         >,
    //     amgcl::solver::cg<SBackend>
    //     > Solver;
    // typedef amgcl::make_solver<
    //     amgcl::amg<
    //         PBackend,
    //         amgcl::coarsening::smoothed_aggregation,
    //         amgcl::relaxation::spai0
    //         >,
    //     amgcl::solver::bicgstab<SBackend>
    //     > Solver;

    Solver::params prm;
    // prm.solver.tol  = 1e-6 ;
    // prm.solver.maxiter = 500 ;
    // prm.precond.coarsening.eps_strong = 0.485435426235199 ;

    // prm.solver.ns_search  = true ; 
    // prm.precond.relax.damping = 0.9 ;
    // prm.solver.verbose = true ;
    // prm.precond.direct_coarse  = true ; 
    prm.precond.coarsening.aggr.eps_strong = EPS_STRONG ;
    prm.precond.coarsening.relax = RELAX ;
    // prm.precond.coarsening.power_iters = 1 ; 
    // prm.precond.coarsening.estimate_spectral_radius = true ;
    // prm.precond.max_levels = 10;
    cout << "---------------flag2------" <<endl ; 
    auto A = sum_A;
    cout << "---------------flag3------" <<endl ; 

    prof.tic("setup");
    Solver solve(A,prm);
    cout << "---------------flag4------" <<endl ; 
    cout << "coarsening.relax : " <<prm.precond.coarsening.relax <<endl ;
    cout << "aggr.eps_strong  : " <<prm.precond.coarsening.aggr.eps_strong <<endl ;
    
#ifdef DEBUG
        std::cout << solve  <<  std::endl;
#endif
    int iters;
    double error;
    prof.toc("setup");

    std::vector<double> x(rows_G,0.0); 
    prof.tic("solve");
    ofstream outfile;

    read_u_t<ptrdiff_t,double>(cout_point, val_U, 0, 1, argv[2]) ;  // read u0 
    rows_U = cout_point[0] ;    cols_U = cout_point[1];

    // vector add vector 
    std::tie(iters, error) = solve(A,val_U, x);
    prof.toc("solve");
        std::cout << "Iters: " << iters << std::endl
            << "Error: " << error << std::endl ;
             
    std::cout << prof << std::endl;

    outfile.open(argv[3]) ;
    for (int i=0 ; i<rows_U;i++ ) {
        outfile << x[i] <<endl ; 
    }
    outfile.close() ;

    return 0 ;
}