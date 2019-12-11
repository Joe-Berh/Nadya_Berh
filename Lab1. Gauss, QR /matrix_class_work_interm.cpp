
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <ctime>


using std:: cout;
using std:: cin;
using std:: string;

/* 1) ++i; 2) const!  3) constructor: Matrix b; - not exist; 4) отключить констр. коп. */
// вызов конструктора в цикле; (про приватность) я только ща понял, что внутри методов класса доступны приватные поля
    // оператор индекса и возврат из функции 
    // без move конструкторОВ инстр. a = Matrix(3,4) и a = zeros(3,4) аналогична


/* 1) конструктор с инит листом допускает разный размер у строк. Как обойти? или вообще мб лучше создать новую структуру вместо инит листа?
*/
template <typename T>
class Matrix {
    int rows, cols;
    T** ptr;

    public:

    friend void operator<< (std::ostream& out, const Matrix& x);
    friend Matrix operator+( const Matrix& lhs, const Matrix& rhs);
    friend Matrix operator+( const Matrix& lhs, Matrix&& rhs);
    friend Matrix operator+( Matrix&& lhs, const Matrix& rhs);
    friend Matrix operator-( const Matrix& lhs, const Matrix& rhs); // перегрузить на мув 
    friend Matrix operator+( const Matrix& lhs, const std:: initializer_list < std:: initializer_list <T> >& rhs);

    Matrix() : rows(0), cols(0), ptr(nullptr) {
        cout << "default ctor running for " << this << '\n';
    }
    Matrix(int _rows, int _cols) {
        cout << "Construtor (_rows, _cols) running for " << this << '\n';
        rows = _rows;
        cols = _cols;
        ptr = new T*[rows];
        //cout << "Memory at " << ptr << " allocated\n";

        for (int i = 0; i < rows; ++i){
            ptr[i] = new T[cols];
            //cout << "Memory at " << ptr[i] << " allocated\n";
            for (int j = 0; j < cols; ++j)
                ptr[i][j] = 0;
        }
    }
    void read_from_file(const string& file_name) {
        
        if (this->ptr != nullptr) {
            cout << "in read_from_file: deleting memory for " << this << '\n';
            for (int i = 0; i < rows; ++i) {
                cout << "Memory at " << ptr[i] << " deleted\n";
                delete[] ptr[i];
            }
            delete[] ptr;
            cout << "Memory at " << ptr << " deleted\n";
        }

        cout << "in read_from_file: allocating memory for " << this << '\n';
        std:: ifstream fin(file_name);
        if ( fin.is_open() ){
            string line, word;

            /* размеры матрицы */
            getline(fin, line, ' ');
            rows = std::stod(line);
            getline(fin, line, '\n');
            cols = std::stod(line);
            //cout << "Rows = " << rows << ", Cols = " << cols << '\n';
            
            ptr = new T*[rows];
            cout << "Memory at " << ptr << " allocated\n";
            for (int i = 0; i < rows; ++i){
                ptr[i] = new T[cols];
                cout << "Memory at " << ptr[i] << " allocated\n";
                for (int j = 0; j < cols; ++j)
                    ptr[i][j] = 0;
            }

            int i = 0;
            int j = 0;
            while( getline(fin, line, '\n') ) {
                std::istringstream stream(line);
                    //cout << "Reading line: " << line << '\n';  // разделитель входит или нет?
                    while ( getline(stream, word, ' ') ) {
                        //cout << "Reading word: " << word << '\n';
                        ptr[i][j] = std::stod(word);
                        j++;
                    }
                i++;
                j = 0;
            }
        }
        else
            cout<< "Error: file wasn't opened\n";
    }

    Matrix(const string& file_name){ // иницализировать матрицу файлом
        cout << "Constructor(file) running for " << this << '\n';
        std:: ifstream fin(file_name);
        if ( fin.is_open() ){
            string line, word;

            /* размеры матрицы */
            getline(fin, line, ' ');
            rows = std::stod(line);
            getline(fin, line, '\n');
            cols = std::stod(line);
            //cout << "Rows = " << rows << ", Cols = " << cols << '\n';
            
            ptr = new T*[rows];
            //cout << "Memory at " << ptr << " allocated\n";
            for (int i = 0; i < rows; ++i){
                ptr[i] = new T[cols];
                //cout << "Memory at " << ptr[i] << " allocated\n";
                for (int j = 0; j < cols; ++j)
                    ptr[i][j] = 0;
            }

            int i = 0;
            int j = 0;
            while( getline(fin, line, '\n') ) {
                std::istringstream stream(line);
                    //cout << "Reading line: " << line << '\n';  // разделитель входит или нет?
                    while ( getline(stream, word, ' ') ) {
                        //cout << "Reading word: " << word << '\n';
                        ptr[i][j] = std::stod(word);
                        j++;
                    }
                i++;
                j = 0;
            }
             
            // std::istringstream stream, _stream;
            //  string line, word;
            // getline(fin, line, '\n');
            // cout << "Is n: " << line.find('\n') << '\n';
            // cout << line << '\n';
            // stream.str(line);
            // cout << "Stream: " << stream.str() << '\n';
            // getline(stream, word, ';');
            // cout << word << '\n';
            // getline(stream, word, ';');
            // cout << word << '\n';
            // getline(fin, line, '\n');
            // cout << "Is n: " << line.find('\n') << '\n';
            // cout << line << '\n';
            // cout << "Trying to figure" << (EOF == line[line.size()]);
            // stream.str(line);
            // _stream.str("3;4");
            // cout << "_Stream: " << _stream.str() << '\n';
            // cout << "Stream: " << stream.str() << '\n';
            // cout << (bool)getline(stream, word, ';') << '\n';
            // cout << word << '\n';

        }
        else
            cout<< "Error: file wasn't opened\n";
    }
    
    Matrix(const std:: initializer_list < std:: initializer_list <T> >& lst) {
        cout << "Construtor with init list running for " << this << '\n';
        rows = lst.size();
        ptr = new T* [ rows ];
        //cout << "Memory at " << ptr << " allocated\n";
        int i = 0;
        for (auto str : lst) {
            ptr[i] = new T[ str.size()];
            //cout << "Memory at " << ptr[i] << " allocated\n";
            int j = 0;
            for (auto el : str) {
                cols = str.size();
                ptr[i][j] = el;
                j++;
            }
            i++;
        }
    }

    
    Matrix(const Matrix& rhs){
        cout << "Copy constructor running for " << this << ", rhs = " << &rhs << '\n';
        rows = rhs.rows;
        cols = rhs.cols;
        ptr = new T*[rows];
        //cout << "Memory at " << ptr << " allocated\n";

        for (int i = 0; i < rows; ++i){
            ptr[i] = new T[cols];
            //cout << "Memory at " << ptr[i] << " allocated\n";
            for (int j = 0; j < cols; j++)
                ptr[i][j] = rhs.ptr[i][j];
        }
    }

    Matrix( Matrix&& rhs ){
        cout << "Copy move constructor running for " << this << ", rhs = " << &rhs << '\n';
        ptr = rhs.ptr;
        rows = rhs.rows;
        cols = rhs.cols;
        rhs.ptr = nullptr;
    }
    
    void init(){
        for (int i = 0; i < rows; ++i){
            cout << "Enter " << i << " string:\n";
            for (int j = 0; j < cols; ++j)
                cin >> ptr[i][j];
        }
    }

    void operator=(const Matrix& rhs) {
        cout << "Im in op= for " << this << ", rhs = " << &rhs << '\n';
        /* !удаление старой памяти*/
        for (int i = 0; i < rows; ++i) {
            //cout << "Memory at " << ptr[i] << " deleted\n";
            delete[] ptr[i];
        }
        //cout << "Memory at " << ptr << " deleted\n";
        delete[] ptr;

        rows = rhs.rows;
        cols = rhs.cols;
        ptr = new T*[rows];
        //cout << "Memory at " << ptr << " allocated\n";
        for (int i = 0; i < rows; ++i){
            ptr[i] = new T[cols];
            //cout << "Memory at " << ptr[i] << " allocated\n";
            for (int j = 0; j < cols; j++)
                ptr[i][j] = rhs.ptr[i][j];
        }
    }

    void operator=(Matrix&& rhs){
        cout << "Im in move for " << this << ", rhs = " << &rhs << '\n';
        /* удалить старую память*/
        for (int i = 0; i < rows; ++i) {
            //cout << "Memory at " << ptr[i] << " deleted\n";
            delete[] ptr[i];
        }
        //cout << "Memory at " << ptr << " deleted\n";
        delete[] ptr;

        ptr = rhs.ptr;
        rows = rhs.rows; cols = rhs.cols;
        rhs.ptr = nullptr;
    }

    T operator() (int i, int j) const {
        return ptr[i][j];
    }
    T& operator() (int i, int j)  {
        return ptr[i][j];
    }

    T& operator() (int i) {
        if (rows == 1) 
            return ptr[0][i];
        if (cols == 1)
            return ptr[i][0];
        else
        {
            throw std:: invalid_argument("Think what ya doing when you do it!");
        }
    }
    T operator() (int i) const {
        if (rows == 1) 
            return ptr[0][i];
        if (cols == 1)
            return ptr[i][0];
        else
        {
            throw std:: invalid_argument("Think what ya doing when you do it!");
        }
    }

    Matrix operator-() const {
        Matrix res = *this;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                res.ptr[i][j] *= -1;
        return res;
    }

    void operator-= (const Matrix& rhs) {
        if ( (rows != rhs.rows) or (cols != rhs.cols))
            throw std::invalid_argument("Dimensions error!");
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++)
                ptr[i][j] -= rhs.ptr[i][j];
        }
    }
    void operator+= (const Matrix& rhs) {
        if ( (rows != rhs.rows) or (cols != rhs.cols))
            throw std::invalid_argument("Dimensions error!");
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++)
                ptr[i][j] += rhs.ptr[i][j];
        }
    }

    Matrix operator*(const Matrix& rhs) {
        if (cols != rhs.rows)
            throw std:: invalid_argument("Dimensions should correlate!");
        else {
            Matrix res(rows, rhs.cols);
            for (int i = 0; i < res.rows; i++) {
                for (int j = 0; j < res.cols; j++) {
                    for (int t = 0; t < cols; t++)
                        res.ptr[i][j] += this->ptr[i][t] * rhs.ptr[t][j];
                }
            }
            return res;
        }
    }

    /* геттеры */
    int get_rows() const {
        return rows;
    }
    int get_cols() const {
        return cols;
    }
    Matrix get_item(int i, const std::string& attr = "col") {
        if (attr == "col") {
            Matrix res(1, this -> rows);
            for (int j = 0; j < this -> rows; j++)
                res(j) = this-> ptr[j][i];
            return res;
        }
        else if (attr == "row") {
            Matrix res(1, this -> cols);
            for (int j = 0; j < this -> cols; j++)
                res(j) = this -> ptr[i][j];
            return res;
        }
        else 
            throw std::invalid_argument("No such attribute!");
    }

    /* преобразования строк*/  
    void sub(int i, int j){   // не забудь правильно проиндексировать!
        for (int t = 0; t < cols; ++t){
            ptr[j][t] = ptr[j][t] - ptr[i][t];
        }
    }
    void scale(int i, int _factor){
        for (int j = 0; j < cols; ++j)
            ptr[i][j] *= _factor;
    }
    void replace(int i, int j){
        T* tmp;
        tmp = ptr[i];
        ptr[i] = ptr[j];
        ptr[j] = tmp;
    }

    void transpose() {
        for (int i = 0; i < rows; i++) 
            for (int j = i + 1; j < cols; j++)
                std::swap(ptr[i][j], ptr[j][i]);
    }
    // Matrix transpose() {  // возможны случаи создания Matrix(n,1), т.е. не экономично 
    //     Matrix res(cols, rows);
    //     for (int i = 0; i < rows; i++) {
    //         for (int j = 0; j < cols; j++) {
    //             res.ptr[j][i] = ptr[i][j];
    //         }
    //     }
    //     return res;
    // }


    T norma() const { // наверное, else не нужен
        if (rows != 1)  throw "Norma is for vectors only!\n";
        else {
            T sum(0);
            for (int j = 0; j < cols; j++)
                sum += ptr[0][j]*ptr[0][j];
            return sqrt(sum);
        }
    }
    T norma_l2 () {
        T sum(0);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++)
                sum += ptr[i][j]*ptr[i][j];
        }
        return sqrt(sum);
    } 

    // T norma_l1 () {
    //     T sum(0);
    //     for (int i = 0; i < rows; i++) 
    //         for (int j = 0; j < cols; j++)
    //             sum += abs(ptr[i][j]);
    //     return sum;
    // }
    T norma_l1 () {
        T max = 0;
        T sum= 0;
        for (int i = 0; i < rows; i++) { // цикл по строчкам
            for (int j = 0; j < cols; j++) {
                sum += abs (this -> ptr[i][j] );
            }
            if (sum > max)
                max = sum;
            sum = 0;
        }
        return max;
    }
    T norma_inf () {
        T max = 0;
        T sum = 0;
        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < rows; i++) 
                sum += abs (this -> ptr[i][j]);
            if (sum > max)
                max = sum;
            sum = 0;
        }
        return max;
    }
    // T norma_inf () {
    //     T max = 0;
    //     for (int i = 0; i < rows; i++)
    //         for (int j = 0; j < cols; j++) {
    //             if (abs(ptr[i][j]) > max)
    //                 max = abs(ptr[i][j]); 
    //         }
    //     return max;
    // }


    ~Matrix(){
        cout << "Destr running for " << this << '\n';
        if (this->ptr != nullptr) {
            for (int i = 0; i < rows; ++i) {
                //cout << "Memory at " << ptr[i] << " deleted\n";
                delete[] ptr[i];
            }
            delete[] ptr;
            //cout << "Memory at " << ptr << " deleted\n";
        }
    }

    

};

template <typename T>
Matrix<T> operator+(const Matrix<T>& lhs, const Matrix<T>& rhs) {
    if ( (lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
            std::invalid_argument("Dimensions error!");
    Matrix<T> res = rhs;
    for (int i = 0; i < lhs.rows; i++) {
        for (int j = 0; j < lhs.cols; j++) 
            res.ptr[i][j] += lhs.ptr[i][j];
    }
    return res;
}
template <typename T>
Matrix<T> operator+( const Matrix<T>& lhs, Matrix<T>&& rhs) {
    if ( (lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
            throw std::invalid_argument("Dimensions error!");
    Matrix<T> res = std::move(rhs);
    for (int i = 0; i < lhs.rows; i++) {
        for (int j = 0; j < lhs.cols; j++) 
            res.ptr[i][j] += lhs.ptr[i][j];
    }
    return res;
}
template <typename T>
Matrix<T> operator+( Matrix<T>&& lhs, const Matrix<T>& rhs) {
    if ( (lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
            throw std::invalid_argument("Dimensions error!");
    Matrix<T> res = std::move(lhs);
    for (int i = 0; i < lhs.rows; i++) {
        for (int j = 0; j < lhs.cols; j++) 
            res.ptr[i][j] += rhs.ptr[i][j];
    }
    return res;
}

template <typename T>
Matrix<T> operator-( const Matrix<T>& lhs, const Matrix<T>& rhs) {
    return lhs + (-rhs);
}
// init list может иметь разный размер по строкам - нужна проверка
template <typename T>
Matrix<T> operator+( const Matrix<T>& lhs, const std:: initializer_list < std:: initializer_list <T> >& rhs) {
    Matrix<T> res = lhs;
    int i(0), j(0);
    for (auto str : rhs) {
        for (auto elem_of_str : str) {
            res.ptr[i][j] += elem_of_str;
            j++;
        }
        i++;
    }        
    return res;
}


template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& x){
    out << '\n';
    for (int i = 0; i < x.rows; ++i) {
        for (int j = 0; j < x.cols; ++j)
            out << x(i,j) << " ";
        out << '\n';
    }
    out << '\n';
    return out;
}

template <typename T>
Matrix<T> eye(int size){
    Matrix<T> res(size,size);
    for (int i = 0; i < size; ++i)
        res(i,i) = 1;
    return res;
}

template <typename T> 
Matrix<T> zeros(int _rows, int _cols){
    Matrix<T> res(_rows, _cols);
    return res;
}

/* МЕТОД ГАУССА И QR */
template <typename T>
int find_leading_element(const Matrix<T>& A, int k) {
    T max = abs(A(k,k));
        int row = k;
        for (int i = k + 1; i < A.get_rows(); i++) {
            if (abs(A(i,k)) > max) {
                max = abs(A(i,k));
                row = i;
            }
        }
    return row;
}
template <typename T>
void gauss_back_run(const Matrix<T>& A, const Matrix<T>& b,  Matrix<T>& X) {
    T sum(0);
    int rows_num = A.get_rows();
    for (int j = rows_num - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
        sum = 0;
        for (int k = j + 1; k < rows_num; k++)
            sum += X(k)*A(j, k);
        X(j) = (b(j) - sum)/A(j, j);
    }
}

//Gauss 3.0 c частичным выбором
template <typename T>
Matrix<T> Gauss( Matrix<T>& A, Matrix<T>& b, const T eps) {
    cout << "In Gauss\nArguments passed:\nA = " << A;
    cout << "b = " << b;

    int rows_num = A.get_rows();  // а может всегда вызывать get_...()? сколько занимает вызов метода (функции?)
    int cols_num = A.get_cols();
    T tmp1(0), tmp2(0);  // инициализировать нулями?

    /*  Прямой ход (c частичным выбором) */
    for (int k = 0; k < (rows_num - 1); k++) {
        //tmp1 = A(k,k);
        int row_of_leading_element = find_leading_element(A, k);
        /* проверк на вырожденность с заданной точностью eps */
        if ( abs(A(row_of_leading_element, k)) < eps )
            throw std::invalid_argument("Given Matrix is non-invertable");

        A.replace(k, row_of_leading_element);
        /* не забыть поменять строки вектора правой части */
        std::swap( b(k), b(row_of_leading_element));

        tmp1 = A(k,k);
        for (int i = k + 1; i < rows_num; i++) {
            tmp2 = A(i,k) / tmp1;

            for (int j = k; j < cols_num; j++) {
                A(i,j) = A(i,j) - A(k,j)*tmp2;
            } 

            /* соответствующим образом меняем вектор правой части */
            b(i) = b(i) - b(k)*tmp2;
        }
    }
    // cout << "После прямого хода:\nA= " << A;
    // cout << "b = " << b;

    /* Обратный ход */ 
    cout << "Create vector-matrix for Gauss-solution\n";
    Matrix<T> X(rows_num, 1);
    gauss_back_run(A, b, X);
    return X;
}
// template <typename T>
// void QR ( Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R) {
//     // cout << "In QR Factor\nArguments passed:\nA = " << _A;
//     // cout << "b = " << b;
    
//     //A.conctenate_with_column(b);
//     int rows_count = A.get_rows();
//     int cols_count = A.get_cols();

//     Matrix<T> str(1, cols_count);
//     Matrix<T> T = eye(rows_count);
//     Matrix<T> next = eye(rows_count);

//     for (int k = 0; k < rows_count - 1; k++) {
//         for (int i = k + 1; i < rows_count; i++) {
//             T norm = sqrt(A(k, k)*A(k, k) + A(i, k)*A(i, k));
//             T c = A(k, k) / norm;
//             T s  = A(i, k) / norm;
            
//             next(k, k) = c; next(k, i) = s;
//             next(i, k) = -s; next(i, i) = c;
//             // cout << "Next: " << next;
//             T = next * T;
//             next = eye(rows_count);

//             // Matrix str(1, cols_count);
//             for (int t = k; t < cols_count; t++)
//                 str(t) = A(k, t);
        
//             for (int j = k; j < cols_count; j++) 
//                 A(k,j) = c*str(j) + s*A(i,j);
        
//             for (int j = k; j < cols_count; j++) 
//                 A(i,j) = -s*str(j) + c*A(i,j);
//             //cout << A;
//         }
//     }
//     R = std:: move(A); // А больше не нужна
//     Q = T; Q.transpose(); // матрица Q - обратная к T, но T ортогональна

//     //Matrix X(1, rows_count);
//     //cout << "Confirm X initialized correctly: X = " << X; // what should it look like

//     // T sum(0);
//     // for (int j = rows_count - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
//     //     sum = 0;
//     //     for (int k = j + 1; k < rows_count; k++)
//     //         sum += X(k)*A(j, k);
//     //     X(j) = (A(j, cols_count - 1) - sum)/A(j,j);
//     // }
//     // return X;
// }




int main() {
    /* inputs */
    Matrix<double> A { {10, 6, 2, 0}, {5,1,-2,4}, {3,5,1,-1}, {0,6,-2,2}};
    Matrix<double>  b { {25},{14},{10},{8} }; 
    //ANSWER = {2, 1, -0.5, 0.5};

    // Matrix A { {1,1,1,1}, {0,1,1,1}, {0,0,1,1}, {0,0,0,1} };
    // Matrix b { {4,3,2,1} }; 
    

    // Matrix A { {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 1, 1}, {1, 1, 1, 1} };
    // Matrix b = { {1, 2, 3, 4} };

    // Matrix A { {1, 1, 1, 1}, {2, 3, 3, 3}, {2, 4, 4, 4}, {4, 5, 6, 7} };
    // Matrix b { {4, 11, 15, 22} };  // несовместная


    // Matrix A { {28.859, -0.008, 2.406, 19.240}, {14.436, -0.001, 1.203, 9.624}, {120.204, -0.032, 10.024, 80.144}, {-57.714, 0.016, -4.812, -38.478}};
    // Matrix b = { {30.459, 18.248, 128.156, -60.908} }; // огромный cond
    /**/
    const double eps = 1e-9; 
    const int n = 5;
    // std:: string A_input("A_input1.txt");
    // std:: string b_input("b_input1.txt");
    std:: string A_input("A_input2.txt");
    std:: string b_input("b_input2.txt");
    // std:: string A_input("A_triang.txt");
    // std:: string b_input("b_triang.txt");


    cout << A;
   //Matrix<double> Sol = Gauss(A,b, eps);
    //cout << "Sol SOLLO = " <<  Sol;
    // cout << "Creating source matrices: \n";
    // Matrix A, b;
    // cout << "Init source matrices from file\n";
    // A.read_from_file(A_input);
    // b.read_from_file(b_input);
    // cout << "Creating Gauss-solution-vector:\n";
    // Matrix Gauss_solution = Gauss(A, b, eps);
    // A.read_from_file(A_input);
    // b.read_from_file(b_input);



    // /* вычисление невязки */
    // cout << "Updating A,b\n"; A.read_from_file(A_input); b.read_from_file(b_input);
    // cout << "Computing residual:\n";
    // Gauss_solution = A*Gauss_solution;
    // Gauss_solution -= b;
    // cout << Gauss_solution;

    /* оценка числа обусловленности */
    // T x_delta, b_delta; // погрешности (относительные) 
    // T current; // текущие значение отношение погрешностей
    // T max = 0;
    // Gauss_solution = Gauss(A, b, eps);
    // T slv_norm = Gauss_solution.norma_l2();
    // A.read_from_file(A_input); b.read_from_file(b_input);
    // T b_norm = b.norma_l2();

    // cout << "Approx cond A:\n";
    // srand(time(NULL));
    // Matrix pertubation(A.get_rows(), 1);
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < A.get_rows(); j++){
    //         pertubation(j) = (T)(rand() % 21 + (-10)) / 1000;
    //     }

    //     b_delta = pertubation.norma_l2() / b_norm;
    //     pertubation += b;
    //     cout << "Pertub = " << pertubation;
    //     pertubation = Gauss(A, pertubation, eps);
    //     /* update A */ A.read_from_file(A_input);
    //     pertubation -= Gauss_solution;
    //     x_delta = pertubation.norma_l2() / slv_norm;
    
    //     current = x_delta / b_delta;
    //     if (current > max)
    //         max = current; 
    // }
    // cout << "СondA > " << max;



    // A.read_from_file("A_input1.txt");
    // b.read_from_file("b_input1.txt");
    // cout << " A = " << A;
    // Matrix Q, R;
    // QR(A, Q, R);
    // cout << Q*R;
    
    // A.read_from_file("A_input1.txt");
    // b.read_from_file("b_input1.txt");
    // Matrix solution = Gauss(A, b, eps);
    // cout << "Gauss Solution: " << solution;
    // A.read_from_file("A_input1.txt");
    // b.read_from_file("b_input1.txt");
    // Matrix Q,R;
    // QR(A, Q, R);
    // /* решение на основе полученного разложения */
    // Q.transpose();
    // b = Q*b;
    // Matrix QR_solution(1, R.get_rows());
    // gauss_back_run(R, b, QR_solution);
    // cout << QR_solution;

    // solution = A*solution;
    // solution -= b;
    // cout << "Residual: " << solution.norma_l2() << '\n';


    
    /* вычисление числа обусловленности 
        Найдем обратную матрицу, для этого решим n уравнений*/
    // Matrix E = eye(A.get_rows()); // единичная матрица
    // Matrix A_inversed;
    // for (int i = 0; i < A.get_rows(); i++) 
    //     A_inversed.conctenate_with_column( Gauss(A, E.get_item(i)) );
    // cout << "Make sure A_inversed*A == E: " << A_inversed*A;
    
    // cout << "Cond (l1) = " << A.norma_l1()*A_inversed.norma_l1() << '\n';
    // cout << "Cond (infinity) = " << A.norma_inf()*A_inversed.norma_inf() << '\n';

    return 0;
}

 /* 1) нужно ли думтаь об алгоритме через методы класса, если можно написать стандартный алгоритм, а класс использовать просто как некий контейнер?
    2) реализация векторов(?) (в гауссе - столбец, но на основании Matrix столбец - затратно (уточнить, почему, кстати) ) (метод concatenate ... использует специальную индексацию vec(i)) 
    3) Нужны ведь 2 версии оператора индекса - константная и нет?
    4) про init list in constrcts , мб см. скрин
    5) доступ к полям параметров в методах)
    6) деструктор не отработает (память не удалится), если выйти из проги по ошибке. Это вообще наскольок плохо? Ведь под экзещник операционка выделит отдельную область памяти и тд
    7) что все таки происходит здесь Matrix B = eye(2); и что должно было бы произойти? 
    8) зачем такой синтаксис: const Matrix some_method(...)
    */