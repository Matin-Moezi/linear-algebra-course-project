#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.141592654

typedef struct
{
    int col_no;
    int row_no;
    double **data;
} Matrix;

Matrix *matrix_init(int row_no, int col_no)
{
    Matrix *matrix = (Matrix *)malloc(sizeof(Matrix));
    matrix->col_no = col_no;
    matrix->row_no = row_no;
    matrix->data = (double **)malloc(row_no * sizeof(double *));
    for (int i = 0; i < row_no; i++)
        matrix->data[i] = (double *)malloc(col_no * sizeof(double));
    return matrix;
}

void matrix_set_all(Matrix *matrix, const double val)
{
    for (int i = 0; i < matrix->row_no; i++)
        for (int j = 0; j < matrix->col_no; j++)
            matrix->data[i][j] = val;
}

Matrix *tri_symm_matrix_init(const double a[], const double b[], int n)
{
    Matrix *matrix = matrix_init(n, n);
    matrix_set_all(matrix, 0);
    for (int i = 1; i < n - 1; i++)
    {
        matrix->data[i][i] = a[i];
        matrix->data[i][i - 1] = b[i - 1];
        matrix->data[i][i + 1] = b[i];
    }
    matrix->data[0][0] = a[0];
    matrix->data[0][1] = b[0];
    matrix->data[n - 1][n - 1] = a[n - 1];
    matrix->data[n - 1][n - 2] = b[n - 2];
    return matrix;
}

Matrix *trans_matrix(const Matrix *matrix)
{
    Matrix *trans = matrix_init(matrix->col_no, matrix->row_no);
    for (size_t i = 0; i < trans->row_no; i++)
    {
        for (size_t j = 0; j < trans->col_no; j++)
        {
            trans->data[i][j] = matrix->data[j][i];
        }
    }
    return trans;
}

double complex *tri_symm_matrix_eigen(const Matrix *matrix)
{
    if (matrix->row_no != matrix->col_no)
        return NULL;
    int n = matrix->row_no;
    double a = matrix->data[0][0];
    double b = matrix->data[1][0];
    double c = matrix->data[0][1];
    double complex *eigenvals = (double complex *)malloc(n * sizeof(double complex));
    double complex bc = b * c;
    for (int k = 0; k < n; k++)
        eigenvals[k] = a + 2 * csqrt(bc) * cos((k + 1) * PI / (n + 1));
    return eigenvals;
}

void print(const Matrix *matrix)
{
    for (int i = 0; i < matrix->row_no; i++)
    {
        for (int j = 0; j < matrix->col_no; j++)
            printf("%7.3g", matrix->data[i][j]);
        printf("\n");
    }
    printf("\n");
}

void vector_set_all(double *vec, double val, int n)
{
    for (size_t i = 0; i < n; i++)
        vec[i] = val;
}

Matrix *scale_matrix(Matrix *matrix, double val)
{
    Matrix *scaled = matrix_init(matrix->row_no, matrix->col_no);
    for (size_t i = 0; i < matrix->row_no; i++)
    {
        for (size_t j = 0; j < matrix->col_no; j++)
        {
            scaled->data[i][j] = matrix->data[i][j] * val;
        }
    }
    return scaled;
}

Matrix *multiple_matrix(Matrix *A, Matrix *B)
{
    if (A->col_no != B->row_no)
        return NULL;
    int m = A->col_no;
    Matrix *C = matrix_init(A->row_no, B->col_no);
    for (size_t i = 0; i < A->row_no; i++)
        for (size_t j = 0; j < B->col_no; j++)
            for (size_t k = 0; k < m; k++)
                C->data[i][j] += A->data[i][k] * B->data[k][j];
    return C;
}

Matrix *multiple_matrix_vec(Matrix *A, double *vec, int n)
{
    if (A->col_no != n)
        return NULL;
    int m = A->row_no;
    Matrix *result = matrix_init(m, 1);
    matrix_set_all(result, 0);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            result->data[i][0] += A->data[i][j] * vec[j];
        }
    }
    return result;
}

Matrix *multiple_vec_matrix(Matrix *A, double *vec, int n)
{
    if (A->row_no != n)
        return NULL;
    int m = A->col_no;
    Matrix *result = matrix_init(1, m);
    matrix_set_all(result, 0);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            result->data[0][i] += A->data[j][i] * vec[j];
        }
    }
    return result;
}

Matrix *identity_matrix(int n)
{
    Matrix *matrix = matrix_init(n, n);
    matrix_set_all(matrix, 0);
    for (size_t i = 0; i < n; i++)
    {
        matrix->data[i][i] = 1;
    }
    return matrix;
}

double *scale_vector(double *vector, double val, int n)
{
    double *scaled = (double *)malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++)
    {
        scaled[i] = vector[i] * val;
    }
    return scaled;
}

double norm2_vector(double *vec, int n)
{
    double result = 0;
    for (size_t i = 0; i < n; i++)
        result += vec[i] * vec[i];
    return sqrt(result);
}

double *vector_cpy(double *vector, int n)
{
    double *copy = (double *)malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++)
        copy[i] = vector[i];
    return copy;
}

void matrix_sub(Matrix *A, Matrix *B, Matrix *result)
{
    if (A->col_no != B->col_no || A->row_no != B->row_no)
        return;
    int n = A->row_no;
    int m = A->col_no;
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            result->data[i][j] = A->data[i][j] - B->data[i][j];
        }
    }
}

void matrix_sum(Matrix *A, Matrix *B, Matrix *result)
{
    if (A->col_no != B->col_no || A->row_no != B->row_no)
        return;
    int n = A->row_no;
    int m = A->col_no;
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < m; j++)
        {
            result->data[i][j] = A->data[i][j] + B->data[i][j];
        }
    }
}
void matrix_to_vec(Matrix *A, double *vec, int n)
{
    for (size_t i = 0; i < n; i++)
    {
        vec[i] = A->data[i][0];
    }
}

void vec_to_matrix(Matrix *A, double *vec, int n)
{
    for (size_t i = 0; i < n; i++)
    {
        A->data[i][0] = vec[i];
    }
}

void print_vec(double *vec, int n)
{
    printf("[");
    for (size_t i = 0; i < n; i++)
    {
        printf("%f ", vec[i]);
    }
    printf("]\n");
}

void symm_lanczos_method(Matrix *A, double *vec, Matrix *T, Matrix *V)
{
    if (A->col_no != A->row_no)
        return;
    int n = A->row_no;
    double a, b = 1;
    double **v = (double **)malloc((n + 1) * sizeof(double *)), *bv;
    Matrix *Av = matrix_init(n, 1), *bv_mat = matrix_init(n, 1), *vAv = matrix_init(1, 1), *iden = identity_matrix(n);
    v[0] = (double *)malloc(n * sizeof(double));
    vector_set_all(v[0], 0, n);
    double *r = vector_cpy(vec, n);
    for (size_t i = 0; i < n; i++)
    {
        v[i + 1] = (double *)malloc(n * sizeof(double));
        v[i + 1] = scale_vector(r, 1 / b, n);
        Av = multiple_matrix_vec(A, v[i + 1], n);
        vAv = multiple_vec_matrix(Av, v[i + 1], n);
        a = vAv->data[0][0];
        T->data[i][i] = a;
        Av = scale_matrix(iden, a);
        matrix_sub(A, Av, Av);
        Av = multiple_matrix_vec(Av, v[i + 1], n);
        bv = scale_vector(v[i], b, n);
        vec_to_matrix(bv_mat, bv, n);
        matrix_sub(Av, bv_mat, Av);
        matrix_to_vec(Av, r, n);
        b = norm2_vector(r, n);
        if (i < n - 1)
        {
            T->data[i + 1][i] = b;
            T->data[i][i + 1] = b;
        }
    }
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            V->data[i][j] = v[i + 1][j];
        }
    }
}

float determinant(const Matrix *A, float k)
{
    float s = 1, det = 0;
    Matrix *B = matrix_init(A->row_no, A->col_no);
    int i, j, m, n, c;
    if (k == 1)
        return A->data[0][0];
    else
    {
        det = 0;
        for (c = 0; c < k; c++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < k; i++)
            {
                for (j = 0; j < k; j++)
                {
                    B->data[i][j] = 0;
                    if (i != 0 && j != c)
                    {
                        B->data[m][n] = A->data[i][j];
                        if (n < (k - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (A->data[0][c] * determinant(B, k - 1));
            s = -1 * s;
        }
    }
    return (det);
}

/*Finding transpose of matrix*/
Matrix *transpose(Matrix *A, Matrix *fac, float r)
{
    int i, j;
    float d;
    Matrix *B = matrix_init(A->row_no, A->col_no);
    Matrix *inverse = matrix_init(A->row_no, A->col_no);
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            B->data[i][j] = fac->data[j][i];
        }
    }
    d = determinant(A, r);
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            inverse->data[i][j] = B->data[i][j] / d;
        }
    }
    return inverse;
}
Matrix *hilb(int n)
{
    Matrix *hilb = matrix_init(n, n);
    double tmp;
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            tmp = i + j + 1;
            hilb->data[i][j] = 1 / tmp;
        }
    }
    return hilb;
}
Matrix *inverse_matrix(Matrix *A, float f)
{
    Matrix *B = matrix_init(A->row_no, A->col_no);
    Matrix *fac = matrix_init(A->row_no, A->col_no);
    int p, q, m, n, i, j;
    for (q = 0; q < f; q++)
    {
        for (p = 0; p < f; p++)
        {
            m = 0;
            n = 0;
            for (i = 0; i < f; i++)
            {
                for (j = 0; j < f; j++)
                {
                    if (i != q && j != p)
                    {
                        B->data[m][n] = A->data[i][j];
                        if (n < (f - 2))
                            n++;
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            fac->data[q][p] = pow(-1, q + p) * determinant(B, f - 1);
        }
    }
    return transpose(A, fac, f);
}

Matrix *diag(Matrix *A)
{
    int n = A->row_no;
    Matrix *diag = matrix_init(n, n);
    matrix_set_all(diag, 0);
    for (size_t i = 0; i < n; i++)
    {
        diag->data[i][i] = A->data[i][i];
    }
    return diag;
}

Matrix *upper(Matrix *A)
{
    int n = A->row_no;
    int m = A->col_no;
    Matrix *upper = matrix_init(n, m);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = i + 1; j < m; j++)
        {
            upper->data[i][j] = A->data[i][j];
        }
    }
    return upper;
}

Matrix *lower(Matrix *A)
{
    int n = A->row_no;
    int m = A->col_no;
    Matrix *lower = matrix_init(n, m);
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            lower->data[i][j] = A->data[i][j];
        }
    }
    return lower;
}

double *jacobi(Matrix *A, double *b, double *x0, const double tol)
{
    int n = A->row_no;
    double *r_vec = (double *)malloc(n * sizeof(double));
    double *x = vector_cpy(x0, n);
    Matrix *b_mat = matrix_init(n, 1), *r = matrix_init(n, 1);
    Matrix *ul = matrix_init(n, n);
    Matrix *x_mat = matrix_init(n, 1);
    Matrix *D = diag(A);
    Matrix *L = lower(A);
    Matrix *U = upper(A);
    matrix_sum(U, L, ul);
    Matrix *d_inv = inverse_matrix(D, n);
    Matrix *B_j = scale_matrix(d_inv, -1);
    B_j = multiple_matrix(B_j, ul);
    Matrix *b_j = multiple_matrix_vec(d_inv, b, n);
    Matrix *Ax0 = multiple_matrix_vec(A, x0, n);
    vec_to_matrix(b_mat, b, n);
    matrix_sub(Ax0, b_mat, r);
    matrix_to_vec(r, r_vec, n);

    while (norm2_vector(r_vec, n) >= tol)
    {
        Matrix *Bx = multiple_matrix_vec(B_j, x, n);
        matrix_sum(Bx, b_j, x_mat);
        matrix_to_vec(x_mat, x, n);
        Ax0 = multiple_matrix_vec(A, x, n);
        matrix_sub(Ax0, b_mat, r);
        matrix_to_vec(r, r_vec, n);
    }
    return x;
}
void vector_sum(double *vec1, double *vec2, double *result, int n)
{
    for (size_t i = 0; i < n; i++)
    {
        result[i] = vec1[i] + vec2[i];
    }
}
double *conj_grad(Matrix *A, double *b, double *x0, const double tol)
{
    int n = A->row_no;
    double beta, a;
    double *r_vec = (double *)malloc(n * sizeof(double));
    double *r_tmp = (double *)malloc(n * sizeof(double));
    double *ap;
    double *x = vector_cpy(x0, n);
    Matrix *b_mat = matrix_init(n, 1), *r = matrix_init(n, 1), *p = matrix_init(n, 1);
    Matrix *r_tmp_mat, *w, *pt, *ptw, *aw, *bp;
    Matrix *Ax0 = multiple_matrix_vec(A, x0, n);
    vec_to_matrix(b_mat, b, n);
    matrix_sub(b_mat, Ax0, r);
    matrix_to_vec(r, r_vec, n);
    double *p_vec = vector_cpy(r_vec, n);
    while (norm2_vector(r_vec, n) * norm2_vector(r_vec, n) >= tol)
    {
        w = multiple_matrix_vec(A, p_vec, n);
        vec_to_matrix(p, p_vec, n);
        pt = trans_matrix(p);
        ptw = multiple_matrix(pt, w);
        a = norm2_vector(r_vec, n) * norm2_vector(r_vec, n) / ptw->data[0][0];
        ap = scale_vector(p_vec, a, n);
        vector_sum(x, ap, x, n);
        aw = scale_matrix(w, a);
        r_tmp = vector_cpy(r_vec, n);
        matrix_sub(r, aw, r);
        matrix_to_vec(r, r_vec, n);
        beta = norm2_vector(r_vec, n) * norm2_vector(r_vec, n) / (norm2_vector(r_tmp, n) * norm2_vector(r_tmp, n));
        bp = scale_matrix(p, beta);
        matrix_sum(r, bp, r);
        matrix_to_vec(p, p_vec, n);
        matrix_to_vec(r, r_vec, n);
    }
    return x;
}
double *gauss_seidel(Matrix *A, double *b, double *x0, const double tol)
{
    int n = A->row_no;
    double *r_vec = (double *)malloc(n * sizeof(double));
    double *x = vector_cpy(x0, n);
    Matrix *b_mat = matrix_init(n, 1), *r = matrix_init(n, 1);
    Matrix *dl = matrix_init(n, n);
    Matrix *x_mat = matrix_init(n, 1);
    Matrix *D = diag(A);
    Matrix *L = lower(A);
    Matrix *U = upper(A);
    matrix_sum(D, L, dl);
    Matrix *dl_inv = inverse_matrix(dl, n);
    Matrix *B_gs = scale_matrix(dl_inv, -1);
    B_gs = multiple_matrix(B_gs, U);
    Matrix *b_gs = multiple_matrix_vec(dl_inv, b, n);
    Matrix *Ax0 = multiple_matrix_vec(A, x0, n);
    vec_to_matrix(b_mat, b, n);
    matrix_sub(Ax0, b_mat, r);
    matrix_to_vec(r, r_vec, n);

    while (norm2_vector(r_vec, n) >= tol)
    {
        Matrix *Bx = multiple_matrix_vec(B_gs, x, n);
        matrix_sum(Bx, b_gs, x_mat);
        matrix_to_vec(x_mat, x, n);
        Ax0 = multiple_matrix_vec(A, x, n);
        matrix_sub(Ax0, b_mat, r);
        matrix_to_vec(r, r_vec, n);
    }
    return x;
}

double *sor_method(Matrix *A, double *b, double *x0, const double w, const double tol)
{
    int n = A->row_no;
    double *r_vec = (double *)malloc(n * sizeof(double));
    double *x = vector_cpy(x0, n);
    Matrix *b_mat = matrix_init(n, 1), *r = matrix_init(n, 1);
    Matrix *dl = matrix_init(n, n), *du = matrix_init(n, n);
    Matrix *x_mat = matrix_init(n, 1);
    Matrix *D = diag(A);
    Matrix *L = lower(A);
    Matrix *U = upper(A);
    Matrix *wL = scale_matrix(L, w);
    Matrix *wU = scale_matrix(U, w);
    Matrix *wD = scale_matrix(D, 1 - w);
    matrix_sum(D, wL, dl);
    Matrix *dl_inv = inverse_matrix(dl, n);
    matrix_sub(wD, wU, du);
    Matrix *B_sor = multiple_matrix(dl_inv, du);
    Matrix *wdl_inv = scale_matrix(dl_inv, w);
    Matrix *b_sor = multiple_matrix_vec(wdl_inv, b, n);

    Matrix *Ax0 = multiple_matrix_vec(A, x0, n);
    vec_to_matrix(b_mat, b, n);
    matrix_sub(Ax0, b_mat, r);
    matrix_to_vec(r, r_vec, n);

    while (norm2_vector(r_vec, n) >= tol)
    {
        Matrix *Bx = multiple_matrix_vec(B_sor, x, n);
        matrix_sum(Bx, b_sor, x_mat);
        matrix_to_vec(x_mat, x, n);
        Ax0 = multiple_matrix_vec(A, x, n);
        matrix_sub(Ax0, b_mat, r);
        matrix_to_vec(r, r_vec, n);
    }
    return x;
}

int main()
{

    // Matrix *iden = identity_matrix(2);
    // Matrix *one = matrix_init(2, 2);
    // double vec[] = {2, -1};
    // matrix_set_all(one, 1);
    // Matrix *mat = matrix_init(3, 3);
    // mat->data[0][0] = 1;
    // mat->data[0][1] = 2;
    // mat->data[0][2] = 3;
    // mat->data[1][0] = 2;
    // mat->data[1][1] = 3;
    // mat->data[1][2] = 4;
    // mat->data[2][0] = 3;
    // mat->data[2][1] = 4;
    // mat->data[2][2] = 5;
    // Matrix *c = upper(mat);
    // print(c);

    /************************ Tridiagonal Symmetric Matrix Eigenvalues *************************************/
    // double a[] = {2, 2, 2, 2}, b[] = {9, 4, 7};
    // Matrix *mat = tri_symm_matrix_init(a, b, 4);
    // printf("The matrix is:\n");
    // print(mat);
    // double complex *eigens = tri_symm_matrix_eigen(mat);
    // if(eigens == NULL)
    //     return -1;
    // printf("The eigenvalues is:\n");
    // for(int i = 0; i < mat->col_no; i++)
    //     printf("%4.4f%+4.4fi\n", creal(eigens[i]), cimag(eigens[i]));
    // return 0;

    /************************* Lanczos Method *****************************************************************/
    // Matrix *mat = matrix_init(3, 3);
    // mat->data[0][0] = 1;
    // mat->data[0][1] = 2;
    // mat->data[0][2] = 3;
    // mat->data[1][0] = 2;
    // mat->data[1][1] = 3;
    // mat->data[1][2] = 4;
    // mat->data[2][0] = 3;
    // mat->data[2][1] = 4;
    // mat->data[2][2] = 5;
    // double vec[] = {1, 0, 0};
    // Matrix *T = matrix_init(3, 3);
    // Matrix *V = matrix_init(3, 3);
    // matrix_set_all(T, 0);
    // printf("Matrix A is:\n");
    // print(mat);
    // symm_lanczos_method(mat, vec, T, V);
    // printf("Matrix T is:\n");
    // print(T);
    // printf("Matrix V is:\n");
    // print(V);

    /************************************ SOR Method ******************************/
    Matrix *mat = matrix_init(3, 3);
    mat->data[0][0] = 5;
    mat->data[0][1] = 1;
    mat->data[0][2] = 1;
    mat->data[1][0] = 1;
    mat->data[1][1] = 5;
    mat->data[1][2] = 1;
    mat->data[2][0] = 1;
    mat->data[2][1] = 1;
    mat->data[2][2] = 5;
    Matrix *A = hilb(5);
    double b[] = {1, 1, 1, 1, 1};
    double x0[] = {1, 1, 1, 1, 1};
    double *res = conj_grad(A, b, x0, 0.0001);


    print_vec(res, 5);
}