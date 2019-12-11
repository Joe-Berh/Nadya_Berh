
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


class Matrix;
Matrix Gauss( Matrix& A, Matrix& b, const float eps);
Matrix Gauss_by_value(Matrix A, Matrix b, const float eps);
 
class Matrix {
    int rows, cols;
    float** ptr;

    public:

    friend void operator<< (std::ostream& out, const Matrix& x);
    friend Matrix operator+( const Matrix& lhs, const Matrix& rhs);
    friend Matrix operator+( const Matrix& lhs, Matrix&& rhs);
    friend Matrix operator+( Matrix&& lhs, const Matrix& rhs);
    friend Matrix operator-( const Matrix& lhs, const Matrix& rhs); // перегрузить на мув 
    friend Matrix operator+( const Matrix& lhs, const std:: initializer_list < std:: initializer_list <float> >& rhs);

    Matrix() : rows(0), cols(0), ptr(nullptr) {
        cout << "default ctor running for " << this << '\n';
    }

    void allocate_memory(int _rows, int _cols) {
        ptr = new float*[_rows];
            // cout << "Memory at " << ptr << " allocated\n";
            for (int i = 0; i < _rows; ++i){
                ptr[i] = new float[_cols];
                // cout << "Memory at " << ptr[i] << " allocated\n";
                for (int j = 0; j < _cols; ++j)
                    ptr[i][j] = 0;
            }
        rows = _rows; cols = _cols;
    }
    void delete_memory(){
        // cout << "Destr running for " << this << '\n';
        if (this->ptr != nullptr) {
            for (int i = 0; i < rows; ++i) {
                // cout << "Memory at " << ptr[i] << " deleted\n";
                delete[] ptr[i];
            }
            delete[] ptr;
            // cout << "Memory at " << ptr << " deleted\n";
        }
    }

    Matrix(int _rows, int _cols) {
        // cout << "Construtor (_rows, _cols) running for " << this << '\n';
        rows = _rows;
        cols = _cols;
        ptr = new float*[rows];
        // cout << "Memory at " << ptr << " allocated\n";

        for (int i = 0; i < rows; ++i){
            ptr[i] = new float[cols];
            // cout << "Memory at " << ptr[i] << " allocated\n";
            for (int j = 0; j < cols; ++j)
                ptr[i][j] = 0;
        }
    }

    void read_from_file(const string& file_name) {
        
        std:: ifstream fin(file_name);
        if ( fin.is_open() ){
            string line, word;

            getline(fin, line, ' ');
            int _rows = std::stod(line);
            getline(fin, line, '\n');
            int _cols = std::stod(line);
            if ( !((rows == _rows) && (cols == _cols)))  // если с новыми размерами, то
                {
                    this->delete_memory();
                    this->allocate_memory(_rows, _cols);
                }

            int i = 0;
            int j = 0;
            while( getline(fin, line, '\n') ) {
                std::istringstream stream(line);
                    while ( getline(stream, word, ' ') ) {
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
            
            ptr = new float*[rows];
            //cout << "Memory at " << ptr << " allocated\n";
            for (int i = 0; i < rows; ++i){
                ptr[i] = new float[cols];
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
        }
        else
            cout<< "Error: file wasn't opened\n";
    }
    
    Matrix(const std:: initializer_list < std:: initializer_list <float> >& lst) {
        // cout << "Construtor with init list running for " << this << '\n';
        rows = lst.size();
        ptr = new float* [ rows ];
        // cout << "Memory at " << ptr << " allocated\n";
        int i = 0;
        for (auto str : lst) {
            ptr[i] = new float[ str.size()];
            // cout << "Memory at " << ptr[i] << " allocated\n";
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
        //cout << "Copy constructor running for " << this << ", rhs = " << &rhs << '\n';
        rows = rhs.rows;
        cols = rhs.cols;
        ptr = new float*[rows];
        //cout << "Memory at " << ptr << " allocated\n";

        for (int i = 0; i < rows; ++i){
            ptr[i] = new float[cols];
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
        // cout << "Im in op= for " << this << ", rhs = " << &rhs << '\n';
        /* !удаление старой памяти*/
        for (int i = 0; i < rows; ++i) {
            //cout << "Memory at " << ptr[i] << " deleted\n";
            delete[] ptr[i];
        }
        //cout << "Memory at " << ptr << " deleted\n";
        delete[] ptr;

        rows = rhs.rows;
        cols = rhs.cols;
        ptr = new float*[rows];
        //cout << "Memory at " << ptr << " allocated\n";
        for (int i = 0; i < rows; ++i){
            ptr[i] = new float[cols];
            //cout << "Memory at " << ptr[i] << " allocated\n";
            for (int j = 0; j < cols; j++)
                ptr[i][j] = rhs.ptr[i][j];
        }
    }

    void operator=(Matrix&& rhs){
        // cout << "Im in move for " << this << ", rhs = " << &rhs << '\n';
        /* удалить старую память*/
        for (int i = 0; i < rows; ++i) {
            // cout << "Memory at " << ptr[i] << " deleted\n";
            delete[] ptr[i];
        }
        // cout << "Memory at " << ptr << " deleted\n";
        delete[] ptr;

        ptr = rhs.ptr;
        rows = rhs.rows; cols = rhs.cols;
        rhs.ptr = nullptr;
    }

    float operator() (int i, int j) const {
        return ptr[i][j];
    }
    float& operator() (int i, int j)  {
        return ptr[i][j];
    }

    float& operator() (int i) {
        if (rows == 1) 
            return ptr[0][i];
        if (cols == 1)
            return ptr[i][0];
        else
        {
            throw std:: invalid_argument("Think what ya doing when you do it!");
        }
    }
    float operator() (int i) const {
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

    Matrix operator*(const Matrix& rhs) const{
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
        float* tmp;
        tmp = ptr[i];
        ptr[i] = ptr[j];
        ptr[j] = tmp;
    }

    void transpose() {
        if ( (rows != cols) or (cols == 0) or (rows == 0))
            throw std::invalid_argument("Transpose works for square matrices only");
        for (int i = 0; i < rows; i++) 
            for (int j = i + 1; j < cols; j++)
                std::swap(ptr[i][j], ptr[j][i]);
    }

    float norma_l2 () const {
        float sum(0);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++)
                sum += ptr[i][j]*ptr[i][j];
        }
        return sqrt(sum);
    } 

    float norma_l1 () {
        float max = 0;
        float sum= 0;
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
    float norma_inf () {
        float max = 0;
        float sum = 0;
        for (int j = 0; j < cols; j++) {
            for (int i = 0; i < rows; i++) 
                sum += abs (this -> ptr[i][j]);
            if (sum > max)
                max = sum;
            sum = 0;
        }
        return max;
    }
    


    Matrix inversed() {
        if ( (rows != cols) or (cols == 0) or (rows == 0))
            throw std::invalid_argument("Inverse works for square matrices only");
        int size = rows;
        Matrix res(size, size);
        Matrix rhs(1, size); rhs(0) = 1;

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++)  // rhs - соответствующий столбец единичной матрицы
                rhs(j) = 0; 
            rhs(i) = 1;
            rhs = Gauss_by_value(*this, rhs, 1e-9);
            for (int j = 0; j < size; j++)
            res(j, i) = rhs(j);
        }
        return res;
    }

    void write_to_file(std::ostream& out) {
        if (cols == 1) {
            out << ptr[0][0];
            for (int k = 1; k < rows; k++)
                out << ";" << ptr[k][0];
            out << '\n';
            return;
        }
        for (int i = 0; i < rows; i++) {
            out << ptr[i][0];
            for (int j = 1; j < cols; j++)
                out << ';' << ptr[i][j];
            out << '\n';
        }
    }

    ~Matrix(){
        // cout << "Destr running for " << this << '\n';
        if (this->ptr != nullptr) {
            for (int i = 0; i < rows; ++i) {
                // cout << "Memory at " << ptr[i] << " deleted\n";
                delete[] ptr[i];
            }
            delete[] ptr;
            // cout << "Memory at " << ptr << " deleted\n";
        }
    }

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
Matrix operator+( const Matrix& lhs, const std:: initializer_list < std:: initializer_list <float> >& rhs) {
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

/* МЕТОД ГАУССА И QR */

int find_leading_element(const Matrix& A, int k) {
    float max = abs(A(k,k));
        int row = k;
        for (int i = k + 1; i < A.get_rows(); i++) {
            if (abs(A(i,k)) > max) {
                max = abs(A(i,k));
                row = i;
            }
        }
    return row;
}

void gauss_back_run(const Matrix& A, const Matrix& b,  Matrix& X) {
    float sum(0);
    float Ajj(0), bj(0);
    int rows_num = A.get_rows();
    for (int j = rows_num - 1; j >= 0; j--) {    // идем с конца (с (n-1)й неизвестной)
        sum = 0;
        Ajj = A(j,j);
        bj = b(j);
        for (int k = j + 1; k < rows_num; k++)
            sum += X(k)*A(j, k);
        X(j) = (bj - sum)/Ajj;
    }
}

//Gauss 3.0 c частичным выбором
Matrix Gauss( Matrix& A, Matrix& b, const float eps) {
    //cout << "In Gauss\nArguments passed:\nA = " << A;
    //cout << "b = " << b;

    int rows_num = A.get_rows();  // а может всегда вызывать get_...()? сколько занимает вызов метода (функции?)
    int cols_num = A.get_cols();
    float tmp1(0), tmp2(0);  // инициализировать нулями?

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
    Matrix X(rows_num, 1);
    gauss_back_run(A, b, X);
    //cout << "Solution = " << X;
    return X;
}

Matrix Gauss_by_value(Matrix A, Matrix b, const float eps) {
    return Gauss(A, b, eps);
}

void QR ( Matrix& A, Matrix& Q) {
    // cout << "In QR Factor\nArguments passed:\nA = " << _A;
    // cout << "b = " << b;
    
    //A.conctenate_with_column(b);
    int rows_count = A.get_rows();
    int cols_count = A.get_cols();

    Matrix str(1, cols_count);
    Q = eye(rows_count);

    for (int k = 0; k < rows_count - 1; k++) {
        for (int i = k + 1; i < rows_count; i++) {
            double norm = sqrt(A(k, k)*A(k, k) + A(i, k)*A(i, k));
            double c = A(k, k) / norm;
            double s  = A(i, k) / norm;

            //cout << "Q WAS = " << Q;
            /* преобразование единич. матрицы E = Q(в нач. момент) в матрицу T, а затем, транспонированием, в Q*/
            for (int t = 0; t < cols_count; t++)
                str(t) = Q(k, t);

            for (int j = 0; j < cols_count; j++) 
                Q(k,j) = c*str(j) + s*Q(i,j);

            for (int j = 0; j < cols_count; j++) 
                Q(i,j) = -s*str(j) + c*Q(i,j);

            /* преобразование матрицы A в R */
            for (int t = k; t < cols_count; t++)
                str(t) = A(k, t);
        
            for (int j = k; j < cols_count; j++) 
                A(k,j) = c*str(j) + s*A(i,j);
        
            for (int j = k; j < cols_count; j++) 
                A(i,j) = -s*str(j) + c*A(i,j);
        }
    }
    Q.transpose();
    //Q = std:: move(T); Q.transpose(); // матрица Q - обратная к T, но T ортогональна
}

Matrix inv(const Matrix& A) {
    if ( (A.get_rows() != A.get_cols()) or (A.get_cols() == 0) or (A.get_rows() == 0))
            throw std::invalid_argument("Inverse works for square matrices only");
    int size = A.get_rows();
    Matrix Result(size, size), Q(size, size), rhs(size, 1);
    Matrix R(A);
    QR(R, Q);
    Q.transpose(); // Q -> T
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) // столбец единичной матрицы
            rhs(j) = 0; 
        rhs(i) = 1;

        gauss_back_run(R, Q*rhs, rhs);

        for (int j = 0; j < size; j++)
            Result(j, i) = rhs(j);
    }
    return Result;
}

// void compute_cond(Matrix& A, std::string A_input){
//     /* вычисление числа обусловленности */
//     const float eps=1e-9;
//     int size = A.get_rows();
//     A_inversed = A.reverse();
//     cout << A*A_inversed;
//     cout << "Cond (l1) = " << A.norma_l1()*A_inversed.norma_l1() << '\n';
//     cout << "Cond (infinity) = " << A.norma_inf()*A_inversed.norma_inf() << '\n';
// }

float compute_residual(const Matrix& A,  Matrix& b, const Matrix& solution, const float eps) {
    //Matrix _A = A; Matrix _b = b;
    cout << "In compute_residual:\n solution = " << solution;
    Matrix b_compare = A * solution;
    cout << "b_compare = " << b_compare;
    cout << "b = " << b;
    b_compare -= b;
    return b_compare.norma_l2();
    // return (A*Gauss(_A, _b, eps) - b).norma_l2();
}

float approx_cond(const Matrix& A, const Matrix& b, int n, const float eps) {
    float x_delta, b_delta; // погрешности (относительные) 
    float current; // текущие значение отношение погрешностей
    float max = 0;
    float pertubation(0); float pertubation_sum(0);
    
    Matrix _b = b; // можно заменить на чтение из файла
    Matrix Gauss_solution = Gauss_by_value(A, _b, eps);
    float slv_norm = Gauss_solution.norma_l2();

    srand(time(NULL));
    for (int i = 0; i < n; i++) {
        _b = b; 
        for (int j = 0; j < A.get_rows(); j++){
            pertubation = (float)(rand() % 21 + (-10)) / 1000;
            _b(j) += pertubation;
            pertubation_sum += pertubation * pertubation;
            // pertubation(j) = (double)(rand() % 21 + (-10)) / 1000;
        }
        b_delta = sqrt(pertubation_sum)/b.norma_l2();

        _b = Gauss_by_value(A, _b, eps);
        _b -= Gauss_solution;
        x_delta = _b.norma_l2() / slv_norm;
        current = x_delta / b_delta;
        //cout << "Current = " << current;
        if (current > max)
            max = current; 
    }
    //cout << "СondA > " << max;
    return max;
}

int main() {
    /* inputs */
    Matrix A { {10, 6, 2, 0}, {5,1,-2,4}, {3,5,1,-1}, {0,6,-2,2}};
    Matrix b { {25},{14},{10},{8} }; 
    //ANSWER = {2, 1, -0.5, 0.5};

    // Matrix A { {1,1,1,1}, {0,1,1,1}, {0,0,1,1}, {0,0,0,1} };
    // Matrix b { {4}, {3} , {2}, {1} }; 
    

    // Matrix A { {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 1, 1}, {1, 1, 1, 1} };
    // Matrix b = { {1, 2, 3, 4} };

    // Matrix A { {1, 1, 1, 1}, {2, 3, 3, 3}, {2, 4, 4, 4}, {4, 5, 6, 7} };
    // Matrix b { {4, 11, 15, 22} };  // несовместная


    // Matrix A { {28.859, -0.008, 2.406, 19.240}, {14.436, -0.001, 1.203, 9.624}, {120.204, -0.032, 10.024, 80.144}, {-57.714, 0.016, -4.812, -38.478}};
    // Matrix b = { {30.459}, {18.248}, {128.156}, {-60.908} }; // огромный cond
    /**/
    const float eps = 1e-9; 
    const int n = 5;
    std:: ofstream out("./Report/data_output_float.csv");
    std:: string A_input("A_input.txt");
    std:: string b_input("b_input.txt");
    // std:: string A_input("A_input2.txt");
    // std:: string b_input("b_input2.txt");
    // std:: string A_input("A_triang.txt");
    // std:: string b_input("b_triang.txt");

    //Matrix A, b; A.read_from_file(A_input); b.read_from_file(b_input);
    
    cout << A*inv(A);
    cout << A*A.inversed();

    Matrix Q;
    cout << "A = " << A;
    out << "System Matrix;;;\n"; A.write_to_file(out);
    Matrix A_copy = A; // для подсчета невязки
    QR(A,Q);
    cout << "QR = " << Q*A;
    out << "Q;\n"; Q.write_to_file(out);
    out << "R;\n"; A.write_to_file(out);
    out << "QR;\n"; (Q*A - A_copy).write_to_file(out);

    Matrix solution (A.get_rows(), 1);
    Q.transpose();
    gauss_back_run(A, Q*b, solution);
    out << "QR Solution;\n";
    solution.write_to_file(out);
    cout << "Solution = " << solution;
    out << "Residual;\n";
    out << compute_residual(A_copy, b, solution, eps) << ";\n";

    cout << "Approx CondA = " << approx_cond(A_copy, b, 5, eps) << '\n';
    cout << "CondA = " << A.norma_inf() * (A.inversed()).norma_inf();
    

    out.close();

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