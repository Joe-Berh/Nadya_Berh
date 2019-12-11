
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>


using std:: cout;
using std:: cin;
using std:: string;

/* 1) ++i; 2) const!  3) constructor: Matrix b; - not exist; 4) отключить констр. коп. */
// вызов конструктора в цикле; (про приватность) я только ща понял, что внутри методов класса доступны приватные поля
    // оператор индекса и возврат из функции 
    // без move конструкторОВ инстр. a = Matrix(3,4) и a = zeros(3,4) аналогична

class my_iter: public std:: iterator <std:: input_iterator_tag, int> {
    // friend class my_container
private:
    int* p;
public:
    my_iter(int* _p) : p(_p) {}
    my_iter(const my_iter& rhs) : p(rhs.p) {}
    bool operator!= (const my_iter& rhs) const{
        return p != rhs.p;
    }
    bool operator== (const my_iter& rhs) const{
        return p == rhs.p;
    }
    my_iter& operator++(){
        ++p;
        return *this;
    }
    int operator*() const {
        return *p;
    }

};

class my_container {
public:
    typedef my_iter iterator;
    my_container(int n) : size(n), data(new int[n]) {
        for (int i = 0; i < size; ++i)
            data[i] = i+1;
    }
    iterator begin(){
        return iterator(data);
    }
    iterator end(){
        return iterator(data + size);
    }
    ~my_container(){
        delete[] data;
    }
private:
    const size_t size;
    int* data;
};

// class mtrx_iter {
// private:
//     //friend class Matrix;
//     double* data;
//     double** to_rows_ptr;
//     int rows, cols;
//     int column_current = 0;
//     int row_current = 0;
// public:
//     mtrx_iter(double* _data, double** _to_rows, int _rows, int _cols) : data(_data),
//                 to_rows_ptr(_to_rows), rows(_rows), cols(_cols) {}
//     //mtrx_iter(const mtrx_iter& rhs) : data(rhs.data) {}  // what for?
//     bool operator!=(const mtrx_iter& rhs ) const {
//         return data != rhs.data;
//     }
//     double operator*(){
//         return *data;
//     }
//     mtrx_iter& operator++(){
//         if (column_current == (cols - 1)) {
//             column_current++;
//             data = to_rows_ptr[row_current] + column_current;
//         }
//         else {
//             row_current++;
//             column_current = 0;
//             data = to_rows_ptr[row_current];
//         }
//         return *this;
//     }
    

// };

/* 1) конструктор с инит листом допускает разный размер у строк. Как обойти? или вообще мб лучше создать новую структуру вместо инит листа?
*/
class Matrix {
    int rows, cols;
    double** ptr;

    public:

    friend void operator<< (std::ostream& out, const Matrix& x);
    friend Matrix operator+( const Matrix& lhs, const Matrix& rhs);
    friend Matrix operator+( const Matrix& lhs, Matrix&& rhs);
    friend Matrix operator+( Matrix&& lhs, const Matrix& rhs);
    friend Matrix operator-( const Matrix& lhs, const Matrix& rhs); // перегрузить на мув 
    friend Matrix operator+( const Matrix& lhs, const std:: initializer_list < std:: initializer_list <double> >& rhs);

    Matrix() : rows(0), cols(0), ptr(nullptr) {cout << "default ctor running for " << this << '\n';}
    Matrix(int _rows, int _cols) {
        cout << "Construtor running for " << this << '\n';
        rows = _rows;
        cols = _cols;
        ptr = new double*[rows];
        cout << "Memory at " << ptr << " allocated\n";

        for (int i = 0; i < rows; ++i){
            ptr[i] = new double[cols];
            cout << "Memory at " << ptr[i] << " allocated\n";
            for (int j = 0; j < cols; ++j)
                ptr[i][j] = 0;
        }
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
            
            ptr = new double*[rows];
            cout << "Memory at " << ptr << " allocated\n";
            for (int i = 0; i < rows; ++i){
                ptr[i] = new double[cols];
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
    
    Matrix(const std:: initializer_list < std:: initializer_list <double> >& lst) {
        cout << "Construtor with init list running for " << this << '\n';
        rows = lst.size();
        ptr = new double* [ rows ];
        cout << "Memory at " << ptr << " allocated\n";
        int i = 0;
        for (auto str : lst) {
            ptr[i] = new double[ str.size()];
            cout << "Memory at " << ptr[i] << " allocated\n";
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
        ptr = new double*[rows];
        cout << "Memory at " << ptr << " allocated\n";

        for (int i = 0; i < rows; ++i){
            ptr[i] = new double[cols];
            cout << "Memory at " << ptr[i] << " allocated\n";
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
            cout << "Memory at " << ptr[i] << " deleted\n";
            delete[] ptr[i];
        }
        cout << "Memory at " << ptr << " deleted\n";
        delete[] ptr;

        rows = rhs.rows;
        cols = rhs.cols;
        ptr = new double*[rows];
        cout << "Memory at " << ptr << " allocated\n";
        for (int i = 0; i < rows; ++i){
            ptr[i] = new double[cols];
            cout << "Memory at " << ptr[i] << " allocated\n";
            for (int j = 0; j < cols; j++)
                ptr[i][j] = rhs.ptr[i][j];
        }
    }

    void operator=(Matrix&& rhs){
        cout << "Im in move for " << this << ", rhs = " << &rhs << '\n';
        /* удалить старую память*/
        for (int i = 0; i < rows; ++i) {
            cout << "Memory at " << ptr[i] << " deleted\n";
            delete[] ptr[i];
        }
        cout << "Memory at " << ptr << " deleted\n";
        delete[] ptr;

        ptr = rhs.ptr;
        rows = rhs.rows; cols = rhs.cols;
        rhs.ptr = nullptr;
    }

    double& operator() (int i, int j) const {
        return ptr[i][j];
    }

    double& operator() (int i) const {
        if (this -> rows  > 1) { cout <<  "You should use this overload of op() for vectors (rows == 1) only!\n";}
        return ptr[0][i];
    }

    Matrix operator-() const {
        Matrix res = *this;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                res.ptr[i][j] *= -1;
        return res;
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
            Matrix res(1, this -> cols);
            for (int j = 0; j < this -> cols; j++)
                res(j) = this-> ptr[j][i];
            return res;
        }
        else if (attr == "row") {
            Matrix res(1, this -> rows);
            for (int j = 0; j < this -> rows; j++)
                res(j) = this -> ptr[i][j];
            return res;
        }
        else 
            throw "No such attribute!";
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
        double* tmp;
        tmp = ptr[i];
        ptr[i] = ptr[j];
        ptr[j] = tmp;
    }

    Matrix transpose() {  // возможны случаи создания Matrix(n,1), т.е. не экономично 
        Matrix res(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res.ptr[j][i] = ptr[i][j];
            }
        }
        return res;
    }

    template <typename T> // concats with column to the right; T should act like a vector;
    void conctenate_with_column(const T& vec) {  // внезапно решил писать все через this -> - меняется что-то?
        if (ptr != nullptr) {
            Matrix tmp = std::move(*this); 
            this -> rows = tmp.get_rows();
            this -> cols = tmp.cols + 1;  // "+1" - новый столбец под добавляемый вектор
            this -> ptr = new double* [ this -> rows];
            cout << "Memory at " << ptr << " allocated (in concatenate method)\n";
            for (int i = 0; i < this -> rows; i++) {
                this -> ptr[i] = new double[ this -> cols];
                cout << "Memory at " << ptr[i] << " allocated (in concatenate method)\n";
                for (int j = 0; j < this -> cols - 1; j++) {
                    this -> ptr[i][j] = tmp.ptr[i][j];
                }
                this -> ptr[i][cols - 1] = vec(i);   // ЗДЕСЬ НУЖНА ПОДДЕРЖКА ТАКОЙ ИНДЕКСАЦИИ ДЛЯ ВХОДА
            }
        }
        else {
            rows = vec.get_cols();  // Мда..... продумать реализацию ))
            cols = 1;
            ptr = new double* [ rows ];
            cout << "Memory at " << ptr << " allocated (in concatenate method)\n";
            for (int i = 0; i < rows; i++) {
                ptr[i] = new double;
                cout << "Memory at " << ptr[i] << " allocated (in concatenate method)\n";
                ptr[i][0] = vec(i);
            }
        }
    }

    double norma() const { // наверное, else не нужен
        if (rows != 1)  throw "Norma is for vectors only!\n";
        else {
            double sum(0);
            for (int j = 0; j < cols; j++)
                sum += ptr[0][j]*ptr[0][j];
            return sqrt(sum);
        }
    }
    double norma_l2 () {
        double sum(0);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++)
                sum += ptr[i][j]*ptr[i][j];
        }
        return sqrt(sum);
    } 

    double norma_l1 () {
        double sum(0);
        for (int i = 0; i < rows; i++) 
            for (int j = 0; j < cols; j++)
                sum += abs(ptr[i][j]);
        return sum;
    }
    double norma_inf () {
        double max = 0;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++) {
                if (abs(ptr[i][j]) > max)
                    max = abs(ptr[i][j]); 
            }
        return max;
    }

    ~Matrix(){
        cout << "Destr running for " << this << '\n';
        if (this->ptr != nullptr) {
            for (int i = 0; i < rows; ++i) {
                cout << "Memory at " << ptr[i] << " deleted\n";
                delete[] ptr[i];
            }
            delete[] ptr;
            cout << "Memory at " << ptr << " deleted\n";
        }
    }

    // typedef mtrx_iter iterator;
    // iterator begin(){ return iterator(ptr[0], ptr, rows, cols); }
    // iterator end(){ return iterator( ptr[rows - 1] + cols, ptr, rows, cols);}

};

Matrix operator+(const Matrix& lhs, const Matrix& rhs) {
    if ( (lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
            std::invalid_argument("Dimensions error!");
    Matrix res = rhs;
    for (int i = 0; i < lhs.rows; i++) {
        for (int j = 0; j < lhs.cols; j++) 
            res.ptr[i][j] += lhs.ptr[i][j];
    }
    return res;
}
Matrix operator+( const Matrix& lhs, Matrix&& rhs) {
    if ( (lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
            throw std::invalid_argument("Dimensions error!");
    Matrix res = std::move(rhs);
    for (int i = 0; i < lhs.rows; i++) {
        for (int j = 0; j < lhs.cols; j++) 
            res.ptr[i][j] += lhs.ptr[i][j];
    }
    return res;
}
Matrix operator+( Matrix&& lhs, const Matrix& rhs) {
    if ( (lhs.rows != rhs.rows) or (lhs.cols != rhs.cols))
            throw std::invalid_argument("Dimensions error!");
    Matrix res = std::move(lhs);
    for (int i = 0; i < lhs.rows; i++) {
        for (int j = 0; j < lhs.cols; j++) 
            res.ptr[i][j] += rhs.ptr[i][j];
    }
    return res;
}
Matrix operator-( const Matrix& lhs, const Matrix& rhs) {
    return lhs + (-rhs);
}
// init list может иметь разный размер по строкам - нужна проверка
Matrix operator+( const Matrix& lhs, const std:: initializer_list < std:: initializer_list <double> >& rhs) {
    Matrix res = lhs;
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

void operator<<(std::ostream& out, const Matrix& x){
    out << '\n';
    for (int i = 0; i < x.rows; ++i) {
        for (int j = 0; j < x.cols; ++j)
            out << x(i,j) << " ";
        out << '\n';
    }
    out << '\n';
}

Matrix eye(int size){
    Matrix res(size,size);
    for (int i = 0; i < size; ++i)
        res(i,i) = 1;
    return res;
}
Matrix zeros(int _rows, int _cols){
    Matrix res(_rows, _cols);
    return res;
}


// Matrix operator+ (const Matrix& lhs, const std::initializer_list <std::initializer_list<double>>& lst) {
//     Matrix res = lhs;
//         // нет проверки на соответствие размерностей
//     int i(0), j(0);
//         for (auto str : lst) {
//             j = 0; 
//             for (auto el : str) {
//                 res.ptr[i][j] += el;
//                 j++;
//             }
//             i++;
//         }
//     return res;
// }



// версия 1.0 работающая
// void Gauss( Matrix& A,  Matrix& b) {   
//     // как в итоге алгоритм связан с моим классом? - никак,
//         // он использ-ся просто как контейнер. А хотелось бы.
//     cout << "In Gauss\nArguments passed:\nA = " << A;
//     cout << "B = " << b;

//     int n = b.get_cols();
//     double tmp1, tmp2;  // инициализировать нулями?

//     /*  Прямой ход */
//     for (int k = 0; k < (n - 1); k++) {
//         tmp1 = A(k,k);
//         for (int i = k + 1; i < n; i++) {
//             tmp2 = A(i,k);
//             b(i) = b(i) - b(k)*A(i,k)/A(k,k);
//             for (int j = k; j < n; j++) {
//                 A(i,j) = A(i,j) - A(k,j)*tmp2/tmp1;
//             } 
//         }
//     }
//     cout << "После прямого хода:\nA= " << A;
//     cout << "b = " << b;
//     /* Обратный ход */
//     const Matrix X(1,n);
//     cout << "Confirm X initialized correctly: X = " << X; // what should it look like

//     double sum(0);
//     for (int j = n - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
//         sum = 0;
//         for (int k = j + 1; k < n; k++)
//             sum += X(k)*A(j, k);
//         X(j) = (b(j) - sum)/A(j,j);
//     }
//     cout << "X = " << X;
// }

// Gauss 2.0 with concat
// void Gauss( Matrix& A,  const Matrix& b) {
//     // как в итоге алгоритм связан с моим классом? - никак,
//         // он использ-ся просто как контейнер. А хотелось бы.

//     cout << "In Gauss\nArguments passed:\nA = " << A;
//     cout << "B = " << b;

//     // /* обьединим матрицы A и b в одну */
//     A.conctenate_with_column(b);
//     cout << "Confirm concatenation worked: A = " << A;


//     int rows_num = A.get_rows();  // а может всегда вызывать get_...()? сколько занимает вызов метода (функции?)
//     int cols_num = A.get_cols();
//     double tmp1, tmp2;  // инициализировать нулями?

//     /*  Прямой ход */
//     for (int k = 0; k < (rows_num - 1); k++) {
//         tmp1 = A(k,k);
//         for (int i = k + 1; i < rows_num; i++) {
//             tmp2 = A(i,k);
//             for (int j = k; j < cols_num; j++) {
//                 A(i,j) = A(i,j) - A(k,j)*tmp2/tmp1;
//             } 
//         }
//     }
//     cout << "После прямого хода:\nA= " << A;
//     cout << "b = " << b;
//     /* Обратный ход */
//     const Matrix X(1,rows_num);
//     cout << "Confirm X initialized correctly: X = " << X; // what should it look like

//     double sum(0);
//     for (int j = rows_num - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
//         sum = 0;
//         for (int k = j + 1; k < rows_num; k++)
//             sum += X(k)*A(j, k);
//         X(j) = (A(j, cols_num - 1) - sum)/A(j,j);
//     }
//     cout << "X = " << X;
// }

int find_leading_element(const Matrix& A, int k) {
    int max = abs(A(k,k));
        int row = k;
        for (int i = k + 1; i < A.get_rows(); i++) {
            if (abs(A(i,k)) > max) {
                max = abs(A(i,k));
                row = i;
            }
        }
    return row;
}

//Gauss 3.0 c частичным выбором
Matrix Gauss(const Matrix& _A,  const Matrix& b) {
    cout << "In Gauss\nArguments passed:\nA = " << _A;
    cout << "b = " << b;
    Matrix A = _A;
    // /* обьединим матрицы A и b в одну */
    A.conctenate_with_column(b);
    cout << "Confirm concatenation worked: A = " << A;


    int rows_num = A.get_rows();  // а может всегда вызывать get_...()? сколько занимает вызов метода (функции?)
    int cols_num = A.get_cols();
    double tmp1(0), tmp2(0);  // инициализировать нулями?

    /*  Прямой ход (c частичным выбором) */
    for (int k = 0; k < (rows_num - 1); k++) {
        //tmp1 = A(k,k);
        int row_of_leading_element = find_leading_element(A, k);
        // if (row_of_leading_element  == 666) 
        //     throw std::invalid_argument("Матрица вырождена!");
    
        A.replace(k, row_of_leading_element);
        tmp1 = A(k,k);

        for (int i = k + 1; i < rows_num; i++) {
            tmp2 = A(i,k);
            for (int j = k; j < cols_num; j++) {
                A(i,j) = A(i,j) - A(k,j)*tmp2/tmp1;
            } 
        }
    }
    cout << "После прямого хода:\nA= " << A;
    cout << "b = " << b;

    /* Обратный ход */ 
    const Matrix X(1, rows_num);
    cout << "Confirm X initialized correctly: X = " << X; // what should it look like

    double sum(0);
    for (int j = rows_num - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
        sum = 0;
        for (int k = j + 1; k < rows_num; k++)
            sum += X(k)*A(j, k);
        X(j) = (A(j, cols_num - 1) - sum)/A(j,j);
    }
    cout << "Solution = " << X;
    return X;
}

Matrix QR (const Matrix& _A, const Matrix& b) {
    cout << "In QR Factor\nArguments passed:\nA = " << _A;
    cout << "b = " << b;
    Matrix A = _A;
    A.conctenate_with_column(b);
    int rows_count = A.get_rows();
    int cols_count = A.get_cols();

    for (int k = 0; k < rows_count - 1; k++) {
        for (int i = k + 1; i < rows_count; i++) {
            double norm = sqrt(A(k, k)*A(k, k) + A(i, k)*A(i, k));
            double c = A(k, k) / norm;
            double s  = A(i, k) / norm;
            //
            Matrix str(1, cols_count);
            for (int t = k; t < cols_count; t++)
                str(t) = A(k, t);
            //
            for (int j = k; j < cols_count; j++) 
                A(k,j) = c*str(j) + s*A(i,j);
        
            for (int j = k; j < cols_count; j++) 
                A(i,j) = -s*str(j) + c*A(i,j);
            cout << A;
        }
    }
    const Matrix X(1, rows_count);
    cout << "Confirm X initialized correctly: X = " << X; // what should it look like

    double sum(0);
    for (int j = rows_count - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
        sum = 0;
        for (int k = j + 1; k < rows_count; k++)
            sum += X(k)*A(j, k);
        X(j) = (A(j, cols_count - 1) - sum)/A(j,j);
    }
    cout << "Solution = " << X;
    return X;
}


int main() {
    Matrix A { {10, 6, 2, 0}, {5,1,-2,4}, {3,5,1,-1}, {0,6,-2,2}};
    Matrix b { {25,14,10,8} }; 
    //ANSWER = {2, 1, -0.5, 0.5};

    // Matrix A { {1,1,1,1}, {0,1,1,1}, {0,0,1,1}, {0,0,0,1} };
    // Matrix b { {4,3,2,1} }; 
    

    // Matrix A { {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 1, 1}, {1, 1, 1, 1} };
    // Matrix vec = { {1, 2, 3, 4} };

    // Matrix A { {1, 1, 1, 1}, {2, 3, 3, 3}, {2, 4, 4, 4}, {4, 5, 6, 7} };
    // Matrix b { {4, 11, 15, 22} };  // несовместная
    


    // Matrix A { {28.859, -0.008, 2.406, 19.240}, {14.436, -0.001, 1.203, 9.624}, {120.204, -0.032, 10.024, 80.144}, {-57.714, 0.016, -4.812, -38.478}};
    // Matrix b = { {30.459, 18.248, 128.156, -60.908} }; // огромный cond

    Matrix slv = Gauss(A,b);
    
    // Matrix slv = Gauss(A, b);
    // cout << "Невязка (A*slv - b).norma_l2() = " << (A*(slv.transpose()) - b.transpose()).norma_l2();

    // /* оценка невязки */
    // const int num(3); // на основе числа num внесения погрешностей в правую часть
    // std:: vector <Matrix> disturbances;
    // std:: vector<Matrix> solves;

    // disturbances.push_back( {{0.01, 0.01, -0.01, 0.01 }});
    // disturbances.push_back( {{0.01, 0.1, 0.01, 0.01} });
    // disturbances.push_back( {{ -0.1, 0.1, 0.001, 0.01}});

    // double x_delta, b_delta; // погрешности (относительные) 

    // double slv_norm = slv.norma_l2();
    // double b_norm = b.norma_l2();
    // double current; // текущие значение отношение погрешностей
    // double max = 0;
    // for (int i = 0; i < num; i++) {
    //     solves.push_back ( Gauss(A, b + disturbances[i]));
    //     x_delta = (slv - solves[i]).norma_l2() / slv_norm;
    //     b_delta = disturbances[i].norma_l2() / b_norm;
    //     current = x_delta / b_delta;
    //     if (current > max)
    //         max = current; 
    // }
    // cout << "Оценка для cond: " << current << '\n';

    // /* вычисление числа обусловленности 
    //     Найдем обратную матрицу, для этого решим n уравнений*/
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