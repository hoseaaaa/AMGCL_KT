#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>      
#include <numeric> // 包含了accumulate函数的头文件
#include <algorithm> // 包含了max_element函数的头文件
using namespace std ;
// 读取文件数据到动态数组



void read_data(  vector <double>& data,string filename) {
	ifstream infile;
	infile.open(filename);   //将文件流对象与文件连接起来 
	if (!infile.is_open()){
        cout <<"error  open file " << filename <<endl ;   //若失败,则输出错误消息,并终止程序运行 
        exit(1) ;
    }
    string s ;
    double  val;
    //0 、
    getline(infile, s) ;
    while (!infile.eof()){
        if (getline(infile, s) ){
            istringstream is(s);
            while (is >> val) {
                data.push_back(val);
            } 
        }
    }
	infile.close();             //关闭文件输入流 
}
void error_count (  vector <double>& dataA,vector <double>& dataB ) {
    vector <double> dataC; 
    double error ;
    for (int i = 0; i<dataA.size();i++){
        error = (dataA[i]-dataB[i])  ; 
        dataC.push_back ( error ) ;
    }
    
    std::cout << std::endl;
    double sum = std::accumulate(dataC.begin(), dataC.end(), 0.0);
    double average = sum / dataC.size();
    auto maxElement = std::max_element(dataC.begin(), dataC.end());
    cout << "avg_error_value " << average     <<endl ;
    cout << "max_error_value " << *maxElement  <<endl ;
    cout  <<endl ;

    // std::cout << "error  : ";
    // for (int i = 0; i < dataC.size(); ++i) {
    //     std::cout << dataC[i] << " "<<endl ;
    // }
}

int main(int argc ,char** argv ) {
    std::string filename1 = argv[1]; // 第一个文件名
    std::string filename2 = argv[2]; // 第二个文件名

    vector<double> array1 ; // 参考结果
    vector<double> array2 ; // 实际结果


    // 读取第一个文件数据到动态数组
    read_data(array1,filename1) ;
    read_data(array2,filename2) ;

/*
    // 打印第一个数组的内容
    std::cout << "Array 1: ";
    for (int i = 0; i < array1.size(); ++i) {
        std::cout << array1[i] << " ";
    }
    std::cout << std::endl;

    // 打印第二个数组的内容
    std::cout << "Array 2: ";
    for (int i = 0; i < array2.size(); ++i) {
        std::cout << array2[i] << " ";
    }
    std::cout << std::endl;
*/
    // 释放动态分配的内存
    error_count (array1,array2) ; 
    // cout << "array1.size() " << array1.size() << endl ; 
    // cout << "array2.size() " << array2.size() << endl ; 
    
    array1.clear();
    array2.clear();

    // delete[] array2;

    return 0;
}
