/*-*-c++-*-*************************************************************************************************************
* Copyright 2016-2022 Inesonic, LLC.
* 
* This file is licensed under two licenses.
*
* Inesonic Commercial License, Version 1:
*   All rights reserved.  Inesonic, LLC retains all rights to this software, including the right to relicense the
*   software in source or binary formats under different terms.  Unauthorized use under the terms of this license is
*   strictly prohibited.
*
* GNU Public License, Version 2:
*   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
*   License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
*   version.
*   
*   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
*   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
*   details.
*   
*   You should have received a copy of the GNU General Public License along with this program; if not, write to the Free
*   Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
********************************************************************************************************************//**
* \file
*
* This header defines the \ref Mat::Api pure-virtual class.
***********************************************************************************************************************/

/* .. sphinx-project inemat */

#ifndef MAT_API_H
#define MAT_API_H

#include "mat_common.h"

#if (defined(__cplusplus))

    #include <cstdint>
    #include <cstdlib>

    using std::uint8_t;
    using std::uint32_t;
    using std::uint64_t;

    using std::int8_t;
    using std::int32_t;
    using std::int64_t;

    extern "C" {

#else

    #include <stdint.h>
    #include <stdlib.h>

#endif

/**
 * Type used to represent a memory allocation size, in bytes.
 */
typedef size_t MatSizeInBytes;

/**
 * Type used to represent a number of elements, rows, or columns.
 */
typedef unsigned long MatNumberElements;

/**
 * Type used to represent an integer value.
 */
typedef long long MatInteger;

/**
 * Type used to hold logical index values.
 */
typedef long long MatLogical;

/**
 * Type used to represent a complex value.
 */
typedef struct _MatComplex {
    /**
     * Real portion
     */
    double r;

    /**
     * Imaginary portion
     */
    double i;
} MatComplex;

/**
 * Enumeration used to indicate column major or row major.
 */
typedef enum _MatMatrixMode {
    /**
     * Indicates row major.
     */
    ROW_MAJOR = 0,

    /**
     * Indicates column major.
     */
    COLUMN_MAJOR = 1
} MatMatrixMode;

/**
 * Enumeration used to indicate the Cholesky decomposition type.
 */
typedef enum _MatCholeskyType {
    /**
     * Indicates lower Cholesky decomposition.
     */
    LOWER = 'L',

    /**
     * Indicates upper Cholesky decomposition.
     */
    UPPER = 'U'
} MatCholeskyType;

/**
 * Enumeration used to indicate a desired operation.
 */
typedef enum _MatOperation {
    /**
     * No operation
     */
    NO_OPERATION = 'N',

    /**
     * Transpose
     */
    TRANSPOSE = 'T',

    /**
     * Conjugate (no transpose)
     */
    CONJUGATE = 'R',

    /**
     * Conjugate transpose
     */
    CONJUGATE_TRANSPOSE = 'C'
} MatOperation;

/**
 * Enumeration of eigenvalue/eigenvector balancing operations.
 */
typedef enum _MatEigenBalanceJob {
    /**
     * Matrix will be neither permuted or scaled.
     */
    NONE = 'N',

    /**
     * Matrix will be permuted.
     */
    PERMUTED = 'P',

    /**
     * Matrix will be scaled.
     */
    SCALED = 'S',

    /**
     * Matrix will be permuted and scaled.
     */
    BALANCED = 'B'
} MatEigenBalanceJob;

/**
 * Enumeration of Schur computation modes.
 */
typedef enum _MatSchurMode {
    /**
     * Indicates that the routine should only calculate eignevalues of the matrix.
     */
    EIGENVALUES_ONLY = 0,

    /**
     * Indicates that the routine should calculate the Schur form without a supplied Q matrix.
     */
    SCHUR_WITHOUT_Q_MATRIX = 1,

    /**
     * Indicates that the routine should calculate the Schur form using a supplied Q matrix.
     */
    SCHUR_WITH_Q_MATRIX = 2
} MatSchurMode;

/**
 * Enumeration of eigenvector sides.
 */
typedef enum _MatSide {
    /**
     * Value indicating left eigenvectors.
     */
    LEFT = 'L',

    /**
     * Value indicating right eigenvectors.
     */
    RIGHT = 'R'
} MatSide;

/**
 * Enumeration of eigenvector computation modes.
 */
typedef enum _MatEigenMode {
    /**
     * Indicates all eigenvectors for the specified side should be computed.
     */
    ALL_FOR_SIDE = 'A',

    /**
     * Indicates all eigenvectors for the specified side should be computed and back-transformed.
     */
    ALL_FOR_SIDE_WITH_BACKTRANSFORM = 'B',

    /**
     * Indicates only select eigenvectors to be computed.
     */
    SELECT = 'S'
} MatEigenMode;

/**
 * Enumeration of machine parameters.
 */
typedef enum _MatMachineParameter {
    /**
     * Machine epsilon.
     */
    EPSILON = 'E',

    /**
     * The safe minimum value, safe_min.  Value such that 1/safe_min does not overflow.
     */
    SAFE_MINIMUM = 'S',

    /**
     * The machine base.
     */
    BASE = 'B',

    /**
     * Machine precision.
     */
    PRECISION = 'P',

    /**
     * Mantissa size - number of mantissa digits supported.
     */
    MANTISSA_SIZE = 'n',

    /**
     * Indicates if rounding occurs in addition. 1.0 indicates rounding, 0.0 indicates no rounding.
     */
    ROUNDING = 'R',

    /**
     * Minimum exponent before the start of gradual underflow
     */
    MINIMUM_EXPONENT = 'M',

    /**
     * Underflow threshold
     */
    UNDERFLOW_THRESHOLD = 'U',

    /**
     * Largest exponent before values will overflow.
     */
    LARGEST_EXPONENT = 'L',

    /**
     * Overflow threshold
     */
    OVERFLOW_THRESHOLD = 'O'
} MatMachineParameter;

/**
 * Type used to represent a matrix allocation function.  This function must be reentrant.
 *
 * \param[in] sizeInBytes      The required allocation size, in bytes.
 *
 * \param[in] alignmentInBytes The required matrix alignment, in bytes.  Value must be a power of 2.
 *
 * \return Returns a pointer to the newly allocated space.  A null pointer is returned on error.
 */
typedef void* (*MatAllocatorFunction)(MatSizeInBytes sizeInBytes, unsigned alignmentInBytes);

/**
 * Type used to represent a matrix deallocation function.  This function must be reentrant.
 *
 * \param[in] p The pointer to the matrix to be deallocated.
 */
typedef void (*MatDeallocatorFunction)(void* p);

/**
 * Function that copies one double precision vector to another -- See BLAS dcopy
 *
 * \param[in]     numberElements    The number of vector elements.
 *
 * \param[in]     source            The source vector.
 *
 * \param[in]     strideSource      The increment or stride to apply to the source elements.
 *
 * \param[in,out] destination       The destination vector.  Space for the destination vector will be allocated before
 *                                  this function is called.
 *
 * \param[in]     strideDestination The stride to impose on the destination vector.
 */
typedef void (*MatBlasDoubleCopy)(
    MatNumberElements numberElements,
    const double*     source,
    MatNumberElements strideSource,
    double*           destination,
    MatNumberElements strideDestination
);

/**
 * Function that computes a matrix-matrix product with double precision dense matrices.  See BLAS dgemm.
 *
 * This function computes \f$ C = \alpha op _ A \left ( A \right ) op _ B \left ( B \right ) + \beta C \f$
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     operationA     The desired operation to perform on A.
 *
 * \param[in]     operationB     The desired operation to perform on B.
 *
 * \param[in]     numberRowsC    The number of rows in matrix A after operation A is applied.  Also represents the
 *                               number of rows in matrix C.
 *
 * \param[in]     numberColumnsC The number of columns in matrix B after operation B is applied.  Also represents the
 *                               number of columns in matrix C.
 *
 * \param[in]     numberColumnsA The number of columns in matrix A after operation A is applied.  Also represents the
 *                               number of rows in matrix B after operation B is applied.
 *
 * \param[in]     alpha          Scalar alpha term.
 *
 * \param[in]     a              Matrix A.
 *
 * \param[in]     strideA        Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     b              Matrix B.
 *
 * \param[in]     strideB        Number of major elements between minor elements for matrix B.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     beta           Scalar beta term.
 *
 * \param[in,out] c              Matrix C term, also contains the result of the operation.
 *
 * \param[in]     strideC        Number of major elements between minor elements for matrix C.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 */
typedef void (*MatBlasDoubleMultplyAdd)(
    MatMatrixMode     isColumnMajor,
    MatOperation      operationA,
    MatOperation      operationB,
    MatNumberElements numberRowsC,
    MatNumberElements numberColumnsC,
    MatNumberElements numberColumnsA,
    double            alpha,
    const double*     a,
    MatNumberElements strideA,
    const double*     b,
    MatNumberElements strideB,
    double            beta,
    double*           c,
    MatNumberElements strideC
);

/**
 * Function that copies one complex vector to another -- See BLAS zcopy
 *
 * \param[in]     numberElements    The number of vector elements.
 *
 * \param[in]     source            The source vector.
 *
 * \param[in]     strideSource      The increment or stride to apply to the source elements.
 *
 * \param[in,out] destination       The destination vector.  Space for the destination vector will be allocated before
 *                                  this function is called.
 *
 * \param[in]     strideDestination The stride to impose on the destination vector.
 */
typedef void (*MatBlasComplexCopy)(
    MatNumberElements numberElements,
    const MatComplex* source,
    MatNumberElements strideSource,
    MatComplex*       destination,
    MatNumberElements strideDestination
);

/**
 * Function that computes a matrix-matrix product with double precision dense matrices.  See BLAS zgemm.
 *
 * This function computes \f$ C = \alpha op _ A \left ( A \right ) op _ B \left ( B \right ) + \beta C \f$
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     operationA     The desired operation to perform on A.
 *
 * \param[in]     operationB     The desired operation to perform on B.
 *
 * \param[in]     numberRowsC    The number of rows in matrix A after operation A is applied.  Also represents the
 *                               number of rows in matrix C.
 *
 * \param[in]     numberColumnsC The number of columns in matrix B after operation B is applied.  Also represents the
 *                               number of columns in matrix C.
 *
 * \param[in]     numberColumnsA The number of columns in matrix A after operation A is applied.  Also represents the
 *                               number of rows in matrix B after operation B is applied.
 *
 * \param[in]     alpha          Scalar alpha term.
 *
 * \param[in]     a              Matrix A.
 *
 * \param[in]     strideA        Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     b              Matrix B.
 *
 * \param[in]     strideB        Number of major elements between minor elements for matrix B.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     beta           Scalar beta term.
 *
 * \param[in,out] c              Matrix C term, also contains the result of the operation.
 *
 * \param[in]     strideC        Number of major elements between minor elements for matrix C.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 */
typedef void (*MatBlasComplexMultplyAdd)(
    MatMatrixMode     isColumnMajor,
    MatOperation      operationA,
    MatOperation      operationB,
    MatNumberElements numberRowsC,
    MatNumberElements numberColumnsC,
    MatNumberElements numberColumnsA,
    const MatComplex* alpha,
    const MatComplex* a,
    MatNumberElements strideA,
    const MatComplex* b,
    MatNumberElements strideB,
    const MatComplex* beta,
    MatComplex*       c,
    MatNumberElements strideC
);

/**
 * Function that performs a PLU decomposition of a dense (general) matrix.  See LAPACK dgetrf.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     numberRows    Number of input matrix rows.
 *
 * \param[in]     numberColumns Number of input matrix columns.
 *
 * \param[in,out] a             The input matrix.  Holds the L and U matrix components on exit.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] pivotData     Returned data used to construct the P matrix. The array should be pre-allocated with at
 *                              least min(numberRows, numberColumns) elements.
 *
 * \return Returns 0 on success.  Returns a negative value if an invalid parameter was supplied.  Returns a positive
 *         value if the matrix is singular.  The value at \f$ a _ { i,i } \f$ is zero.
 */
typedef MatInteger (*MatLapackDoublePlu)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    double*           a,
    MatNumberElements lda,
    MatInteger*       pivotData
);

/**
 * Function that calculates the inverse of an LU factorized square matrix.  See LAPACK getri
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     matrixOrder   The number of rows/columns of the matrix.
 *
 * \param[in,out] a             The input matrix.  Holds the L and U matrix components calculated using
 *                              MatLapackDoublePlu (getrf).  On exit, holds the calculated inverse matrix.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     pivotData     Pivot data calculated by MatLapackDoublePlu (getrf).
 *
 * \return Returns 0 on success.  Returns a negative value if an invalid parameter was supplied.  Returns a positive
 *         value if the matrix is singular.  The value at \f$ a _ { i,i } \f$ is zero.
 */
typedef MatInteger (*MatLapackDoubleLuInverse)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    double*           a,
    MatNumberElements lda,
    const MatInteger* pivotData
);

/**
 * Function that calculates QR factorization of a matrix.  See LAPACK deqrf.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     numberRows    The number of rows in the input matrix.
 *
 * \param[in]     numberColumns The number of columns in the input matrix.
 *
 * \param[in,out] a             The input matrix.  Holds the L and U matrix components calculated using
 *                              MatLapackDoublePlu (getrf).  On exit, holds the calculated inverse matrix.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] tau           Array holding at least min(numberRows, numberColumns) that can be used to generate the
 *                              elementary reflectors for the Q matrix.
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackDoubleQrFactorization)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    double*           a,
    MatNumberElements lda,
    double*           tau
);

/**
 * Function that calculates the Q matrix from the QR factorization generated by MatLapackDoubleQrFactorization (LAPACK
 * deqrf).  See LAPACK dorgqr.
 *
 * \param[in]     isColumnMajor    If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                 operation should assume a row major matrix.
 *
 * \param[in]     numberRows       The number of rows in the input matrix.
 *
 * \param[in]     numberColumns    The number of columns in the input matrix.
 *
 * \param[in]     numberReflectors The number of Q matrix reflectors.
 *
 * \param[in,out] a                The matrix data returned by MatLapackDoubleQrFactorization (LAPACK deqrf).  The
 *                                 matrix will hold the Q orthogonal matrix on exit.
 *
 * \param[in]     lda              Number of major elements between minor elements for matrix A.  The value is used to
 *                                 determine the stride.  In column major mode, this will represent the number of rows
 *                                 of space between columns which should always be greater than or equal to the number
 *                                 of rows.  In row major mode, this will represent the number of columns of space
 *                                 between rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     tau              Array holding the elementary reflector terms returned by
 *                                 MatLapackDoubleQrFactorization (LAPACK deqrf).
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackDoubleGenerateQFromQrMatrix)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    MatNumberElements numberReflectors,
    double*           a,
    MatNumberElements lda,
    const double*     tau
);

/**
 * Function that calculates the Cholesky factorization of a Hermitian positive-definite matrix. See LAPACK dpotrf.
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     choleskyType   Indicates if lower or upper Cholesky should be computed.
 *
 * \param[in]     matrixOrder    The matrix order.
 *
 * \param[in,out] a              The matrix data to compute the Cholesky decomposition of.  Returns the upper or lower
 *                               triangular matrix on exit.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackDoubleCholesky)(
    MatMatrixMode     isColumnMajor,
    MatCholeskyType   choleskyType,
    MatNumberElements matrixOrder,
    double*           a,
    MatNumberElements lda
);

/**
 * Function that calculates the upper Hessenberg form of a matrix.  See LAPACK dgehrd.
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     matrixOrder    The matrix order.
 *
 * \param[in,out] a              The matrix data to compute the upper Hessenberg form from.  On exit, this array will
 *                               hold the upper Hessenberg matrix.  The elements below the subdiagonal elements will
 *                               contain the Q matrix as the product of n elementary reflectors.  Use the LAPACK dorghr
 *                               function to calculate the Q matrix.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] tau            An array that will contain the elementry reflectors for the Q matrix.  The array
 *                               will be allocated on entry and must have at least matrixOrder - 1 entries (minimum of
 *                               1 entry).
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackDoubleUpperHessenberg)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    double*           a,
    MatNumberElements lda,
    double*           tau
);

/**
 * Function that calculates the Q orthogonal matrix determined by MatLapackDoubleHessenberg (LAPACK dgehrd).  See
 * LAPACK dorghr.
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     matrixOrder    The matrix order.
 *
 * \param[in,out] a              The matrix data returned by the MatLapackDoubleHessenberg (LAPACK dgehrd) function.
 *                               On exit, this array will contain the Q matrix.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     tau            An array holding the elementry reflectors returned by MatLapackDoubleHessenberg
 *                               (LAPACK dgehrd).
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackDoubleUpperHessenbergQMatrix)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    double*           a,
    MatNumberElements lda,
    const double*     tau
);

/**
 * Function that calculates the singular value decomposition of a matrix.  See LAPACK dgesvd.
 *
 * Note that this function will return all coumns of U and V (jobu, jobvt 'A' option).
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     numberRows     The number of rows in the input matrix.
 *
 * \param[in]     numberColumns  The number of columns in the input matrix.
 *
 * \param[in,out] a              The input matrix A.  The array will be destroyed on exit.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] s              Array to hold the S matrix on exit.  The array will be allocated on entry.
 *
 * \param[in,out] u              Array to hold the U matrix on exit.  The array will be allocated on entry.
 *
 * \param[in]     ldu            Number of major elements between minor elements for matrix U.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     vt             Array to hold the V-transpose matrix on exit.  The array will be allocated on entry.
 *
 * \param[in]     ldvt           Number of major elements between minor elements for matrix V-transpose.  The value is
 *                               used to determine the stride.  In column major mode, this will represent the number of
 *                               rows of space between columns which should always be greater than or equal to the
 *                               number of rows.  In row major mode, this will represent the number of columns of space
 *                               between rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] suberb         If the result does not converge, this array will hold the unconverged
 *                               super-diagnonals such that the diagnonal matrix, B, defined by these values satisfies
 *                               \f$ A = U B V^T \f$.  This array will be allocated on entry.
 *
 * \return Returns 0 on success.  A negative value represents an illegal parameter value.  A positive value indicates
 *         that the decomposition did not converge and the value represents the number of super-diagnonals in suberb
 *         that did not converge to zero.
 */
typedef MatInteger (*MatLapackDoubleSvd)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    double*           a,
    MatNumberElements lda,
    double*           s,
    double*           u,
    MatNumberElements ldu,
    double*           vt,
    MatNumberElements ldvt,
    double*           superb
);

/**
 * Function that calculates matrix row/column scale factors to equilibrate a matrix.   See LAPACK dgeequ.
 *
 * \param[in]     isColumnMajor      If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                   operation should assume a row major matrix.
 *
 * \param[in]     numberRows         The number of rows in the input matrix.
 *
 * \param[in]     numberColumns      The number of columns in the input matrix.
 *
 * \param[in,out] a                  The input matrix A.  The array will be destroyed on exit.
 *
 * \param[in]     lda                Number of major elements between minor elements for matrix A.  The value is used
 *                                   to determine the stride.  In column major mode, this will represent the number of
 *                                   rows of space between columns which should always be greater than or equal to the
 *                                   number of rows.  In row major mode, this will represent the number of columns of
 *                                   space between rows which will always be greater than or equal to the number of
 *                                   columns.
 *
 * \param[in,out] rowScaleFactors    On exit, this array holding scale factors for each row.  The array will be
 *                                   allocated on entry.
 *
 * \param[in,out] columnScaleFactors On exit, this array holding scale factors for each column.  The array will be
 *                                   allocated on entry.
 *
 * \param[out]    rowRatio           On exit, this value will hold the ratio of the largest row scale factor to the
 *                                   smallest row scale factor.
 *
 * \param[out]    columnRatio        On exit, this value will hold the ratio of the lartest column scale factor to the
 *                                   smallest column scale factor.
 *
 * \param[out]    aMax               The absolute value of the largest element of A.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  Returns a positive value less than the number of
 *         rows indicating the row number of the first row that is exactly 0.
 */
typedef MatInteger (*MatLapackDoubleEquilibrate)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const double*     a,
    MatNumberElements lda,
    double*           rowScaleFactors,
    double*           columnScaleFactors,
    double*           rowRatio,
    double*           columnRatio,
    double*           aMax
);

/**
 * Function that calculates matrix row/column scale factors to equilibrate a matrix to a power of 2.   See LAPACK
 * dgeequb.
 *
 * \param[in]     isColumnMajor      If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                   operation should assume a row major matrix.
 *
 * \param[in]     numberRows         The number of rows in the input matrix.
 *
 * \param[in]     numberColumns      The number of columns in the input matrix.
 *
 * \param[in,out] a                  The input matrix A.  The array will be destroyed on exit.
 *
 * \param[in]     lda                Number of major elements between minor elements for matrix A.  The value is used
 *                                   to determine the stride.  In column major mode, this will represent the number of
 *                                   rows of space between columns which should always be greater than or equal to the
 *                                   number of rows.  In row major mode, this will represent the number of columns of
 *                                   space between rows which will always be greater than or equal to the number of
 *                                   columns.
 *
 * \param[in,out] rowScaleFactors    On exit, this array holding scale factors for each row.  The array will be
 *                                   allocated on entry.
 *
 * \param[in,out] columnScaleFactors On exit, this array holding scale factors for each column.  The array will be
 *                                   allocated on entry.
 *
 * \param[out]    rowRatio           On exit, this value will hold the ratio of the largest row scale factor to the
 *                                   smallest row scale factor.
 *
 * \param[out]    columnRatio        On exit, this value will hold the ratio of the lartest column scale factor to the
 *                                   smallest column scale factor.
 *
 * \param[out]    aMax               The absolute value of the largest element of A.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  Returns a positive value less than the number of
 *         rows indicating the row number of the first row that is exactly 0.
 */
typedef MatInteger (*MatLapackDoubleEquilibratePowerOf2)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const double*     a,
    MatNumberElements lda,
    double*           rowScaleFactors,
    double*           columnScaleFactors,
    double*           rowRatio,
    double*           columnRatio,
    double*           aMax
);

/**
 * Function that calculates the solution to a system of linear equations.  Solves for \f$ X \f$ given \f$ A X = B \f$.
 * See LAPACK dsgesv.
 *
 * \param[in]     isColumnMajor        If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                     operation should assume a row major matrix.
 *
 * \param[in]     matrixAOrder         The order of matrix A.
 *
 * \param[in]     numberRightHandSides The number of terms in the right hand side of the equation.  This value is also
 *                                     the same as the number of columns in matrix B.
 *
 * \param[in,out] a                    The square coefficient matrix.  The matrix will be destroyed on exit.
 *
 * \param[in]     lda                  Number of major elements between minor elements for matrix A.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[in,out] pivotData            An array holding matrixOrder pivot data terms.  The array will be allocated on
 *                                     entry.
 *
 * \param[in]     b                    An array of size ldb * numberRightHandSides for column major or
 *                                     lbd * matrixOrder for row major holding the B matrix terms.
 *
 * \param[in]     ldb                  Number of major elements between minor elements for matrix B.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[in,out] x                    Array to hold the resulting X matrix values.  The array will be allocated on
 *                                     entry.
 *
 * \param[in]     ldx                  Number of major elements between minor elements for matrix X.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[out]    iterationDetails     A value indicating the result of the solver.  The value is currently ignored.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  A positive value indicates the upper triangular
 *         matrix in the PLU decomposition was singular so a solution could not be computed.
 */
typedef MatInteger (*MatLapackDoubleSolve)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    MatNumberElements numberRightHandSides,
    double*           a,
    MatNumberElements lda,
    MatInteger*       pivotData,
    const double*     b,
    MatNumberElements ldb,
    double*           x,
    MatNumberElements ldx,
    MatInteger*       iterationDetails
);

/**
 * Function that calculates the solution to an over-determined or under-determined system of linear equations.  Solves
 * for \f$ X \f$ given \f$ A X = B \f$.  See LAPACK dgels.
 *
 * \param[in]     isColumnMajor        If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                     operation should assume a row major matrix.
 *
 * \param[in]     transposed           Flag indicating if the A matrix is transposed.  A non-zero value indicates
 *                                     transposed.
 *
 * \param[in]     numberRows           The number of rows in matrix A.
 *
 * \param[in]     numberColumns        The number of columns in matrix A.
 *
 * \param[in]     numberRightHandSides The number of terms in the right hand side of the equation.  This value is also
 *                                     the same as the number of columns in matrix B.
 *
 * \param[in,out] a                    The coefficient matrix.  The matrix will be destroyed on exit.
 *
 * \param[in]     lda                  Number of major elements between minor elements for matrix A.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[in,out] b                    An array of size ldb * numberRightHandSides for column major or
 *                                     lbd * matrixOrder * max(numberRows, numberColumns) for row major holding the B
 *                                     matrix terms.  On exit, this array will hold either the least squares fit
 *                                     solution or a minimum norm solution (if the a matrix is not full rank).
 *
 * \param[in]     ldb                  Number of major elements between minor elements for matrix B.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  A positive value indicates a matrix is not full
 *         rank so a least squares solution could not be determined.  In this scenario, b will contain the minimum
 *         norm solution vectors.
 */
typedef MatInteger (*MatLapackDoubleLeastSquaresSolve)(
    MatMatrixMode     isColumnMajor,
    MatOperation      transposed,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    MatNumberElements numberRightHandSides,
    double*           a,
    MatNumberElements lda,
    double*           b,
    MatNumberElements ldb
);

/**
 * Function that performs a PLU decomposition of a dense (general) matrix.  See LAPACK zgetrf.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     numberRows    Number of input matrix rows.
 *
 * \param[in]     numberColumns Number of input matrix columns.
 *
 * \param[in,out] a             The input matrix.  Holds the L and U matrix components on exit.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] pivotData     Returned data used to construct the P matrix. The array should be pre-allocated with at
 *                              least min(numberRows, numberColumns) elements.
 *
 * \return Returns 0 on success.  Returns a negative value if an invalid parameter was supplied.  Returns a positive
 *         value if the matrix is singular.  The value at \f$ a _ { i,i } \f$ is zero.
 */
typedef MatInteger (*MatLapackComplexPlu)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    MatComplex*       a,
    MatNumberElements lda,
    MatInteger*       pivotData
);

/**
 * Function that calculates the inverse of an LU factorized square matrix.  See LAPACK getri
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     matrixOrder   The number of rows/columns of the matrix.
 *
 * \param[in,out] a             The input matrix.  Holds the L and U matrix components calculated using
 *                              MatLapackComplexPlu (getrf).  On exit, holds the calculated inverse matrix.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     pivotData     Pivot data calculated by MatLapackComplexPlu (getrf).
 *
 * \return Returns 0 on success.  Returns a negative value if an invalid parameter was supplied.  Returns a positive
 *         value if the matrix is singular.  The value at \f$ a _ { i,i } \f$ is zero.
 */
typedef MatInteger (*MatLapackComplexLuInverse)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    MatComplex*       a,
    MatNumberElements lda,
    const MatInteger* pivotData
);

/**
 * Function that calculates QR factorization of a matrix.  See LAPACK zgeqrf.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     numberRows    The number of rows in the input matrix.
 *
 * \param[in]     numberColumns The number of columns in the input matrix.
 *
 * \param[in,out] a             The input matrix.  Holds the L and U matrix components calculated using
 *                              MatLapackComplexPlu (getrf).  On exit, holds the calculated inverse matrix.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] tau           Array holding at least min(numberRows, numberColumns) that can be used to generate the
 *                              elementary reflectors for the Q matrix.
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackComplexQrFactorization)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    MatComplex*       a,
    MatNumberElements lda,
    MatComplex*       tau
);

/**
 * Function that calculates the Q matrix from the QR factorization generated by MatLapackComplexQrFactorization (LAPACK
 * deqrf).  See LAPACK zungqr.
 *
 * \param[in]     isColumnMajor    If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                 operation should assume a row major matrix.
 *
 * \param[in]     numberRows       The number of rows in the input matrix.
 *
 * \param[in]     numberColumns    The number of columns in the input matrix.
 *
 * \param[in]     numberReflectors The number of Q matrix reflectors.
 *
 * \param[in,out] a                The matrix data returned by MatLapackComplexQrFactorization (LAPACK zeqrf).  The
 *                                 matrix will hold the Q orthogonal matrix on exit.
 *
 * \param[in]     lda              Number of major elements between minor elements for matrix A.  The value is used to
 *                                 determine the stride.  In column major mode, this will represent the number of rows
 *                                 of space between columns which should always be greater than or equal to the number
 *                                 of rows.  In row major mode, this will represent the number of columns of space
 *                                 between rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     tau              Array holding the elementary reflector terms returned by
 *                                 MatLapackComplexQrFactorization (LAPACK zeqrf).
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackComplexGenerateQFromQrMatrix)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    MatNumberElements numberReflectors,
    MatComplex*       a,
    MatNumberElements lda,
    const MatComplex* tau
);

/**
 * Function that calculates the Cholesky factorization of a Hermitian positive-definite matrix. See LAPACK zpotrf.
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     choleskyType   Indicates if lower or upper Cholesky should be computed.
 *
 * \param[in]     matrixOrder    The matrix order.
 *
 * \param[in,out] a              The matrix data to compute the Cholesky decomposition of.  Returns the upper or lower
 *                               triangular matrix on exit.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackComplexCholesky)(
    MatMatrixMode     isColumnMajor,
    MatCholeskyType   choleskyType,
    MatNumberElements matrixOrder,
    MatComplex*       a,
    MatNumberElements lda
);

/**
 * Function that calculates the upper Hessenberg form of a matrix.  See LAPACK zgehrd.
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     matrixOrder    The matrix order.
 *
 * \param[in]     ilo            The calculated LAPACK zgebal ilo value.  Set to 1 for normal use.
 *
 * \param[in]     ihi            The calculated LAPACK zgebal ihi value.  Set to matrixOrder for normal use.
 *
 * \param[in,out] a              The matrix data to compute the upper Hessenberg form from.  On exit, this array will
 *                               hold the upper Hessenberg matrix.  The elements below the subdiagonal elements will
 *                               contain the Q matrix as the product of n elementary reflectors.  Use the LAPACK zorghr
 *                               function to calculate the Q matrix.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] tau            An array that will contain the elementry reflectors for the Q matrix.  The array
 *                               will be allocated on entry and must have at least matrixOrder - 1 entries (minimum of
 *                               1 entry).
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackComplexUpperHessenberg)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    MatInteger        ilo,
    MatInteger        ihi,
    MatComplex*       a,
    MatNumberElements lda,
    MatComplex*       tau
);

/**
 * Function that calculates the Q orthogonal matrix determined by MatLapackComplexHessenberg (LAPACK zgehrd).  See
 * LAPACK zunghr.
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     matrixOrder    The matrix order.
 *
 * \param[in]     ilo            The calculated LAPACK zgebal ilo value.  Set to 1 for normal use.
 *
 * \param[in]     ihi            The calculated LAPACK zgebal ihi value.  Set to matrixOrder for normal use.
 *
 * \param[in,out] a              The matrix data returned by the MatLapackComplexHessenberg (LAPACK zgehrd) function.
 *                               On exit, this array will contain the Q matrix.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     tau            An array holding the elementry reflectors returned by MatLapackComplexHessenberg
 *                               (LAPACK zgehrd).
 *
 * \return Returns 0 on success.  Returns a non-zero value on error.
 */
typedef MatInteger (*MatLapackComplexUpperHessenbergQMatrix)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    MatInteger        ilo,
    MatInteger        ihi,
    MatComplex*       a,
    MatNumberElements lda,
    const MatComplex* tau
);

/**
 * Function that calculates the singular value decomposition of a matrix.  See LAPACK zgesvd.
 *
 * Note that this function will return all coumns of U and V (jobu, jobvt 'A' option).
 *
 * \param[in]     isColumnMajor  If non-zero, the operation should assume a column major matrix.  If zero, the
 *                               operation should assume a row major matrix.
 *
 * \param[in]     numberRows     The number of rows in the input matrix.
 *
 * \param[in]     numberColumns  The number of columns in the input matrix.
 *
 * \param[in,out] a              The input matrix A.  The array will be destroyed on exit.
 *
 * \param[in]     lda            Number of major elements between minor elements for matrix A.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] s              Array to hold the S matrix on exit.  The array will be allocated on entry.
 *
 * \param[in,out] u              Array to hold the U matrix on exit.  The array will be allocated on entry.
 *
 * \param[in]     ldu            Number of major elements between minor elements for matrix U.  The value is used to
 *                               determine the stride.  In column major mode, this will represent the number of rows of
 *                               space between columns which should always be greater than or equal to the number of
 *                               rows.  In row major mode, this will represent the number of columns of space between
 *                               rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     vh             Array to hold the V-conjugate transpose matrix on exit.  The array will be allocated
 *                               on entry.
 *
 * \param[in]     ldvh           Number of major elements between minor elements for matrix V-conjugate transpose.  The
 *                               value is used to determine the stride.  In column major mode, this will represent the
 *                               number of rows of space between columns which should always be greater than or equal
 *                               to the number of rows.  In row major mode, this will represent the number of columns
 *                               of space between rows which will always be greater than or equal to the number of
 *                               columns.
 *
 * \param[in,out] suberb         If the result does not converge, this array will hold the unconverged
 *                               super-diagnonals such that the diagnonal matrix, B, defined by these values satisfies
 *                               \f$ A = U B V^T \f$.  This array will be allocated on entry.
 *
 * \return Returns 0 on success.  A negative value represents an illegal parameter value.  A positive value indicates
 *         that the decomposition did not converge and the value represents the number of super-diagnonals in suberb
 *         that did not converge to zero.
 */
typedef MatInteger (*MatLapackComplexSvd)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    MatComplex*       a,
    MatNumberElements lda,
    double*           s,
    MatComplex*       u,
    MatNumberElements ldu,
    MatComplex*       vt,
    MatNumberElements ldvt,
    double*           superb
);

/**
 * Function that calculates matrix row/column scale factors to equilibrate a matrix.   See LAPACK zgeequ.
 *
 * \param[in]     isColumnMajor      If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                   operation should assume a row major matrix.
 *
 * \param[in]     numberRows         The number of rows in the input matrix.
 *
 * \param[in]     numberColumns      The number of columns in the input matrix.
 *
 * \param[in,out] a                  The input matrix A.  The array will be destroyed on exit.
 *
 * \param[in]     lda                Number of major elements between minor elements for matrix A.  The value is used
 *                                   to determine the stride.  In column major mode, this will represent the number of
 *                                   rows of space between columns which should always be greater than or equal to the
 *                                   number of rows.  In row major mode, this will represent the number of columns of
 *                                   space between rows which will always be greater than or equal to the number of
 *                                   columns.
 *
 * \param[in,out] rowScaleFactors    On exit, this array holding scale factors for each row.  The array will be
 *                                   allocated on entry.
 *
 * \param[in,out] columnScaleFactors On exit, this array holding scale factors for each column.  The array will be
 *                                   allocated on entry.
 *
 * \param[out]    rowRatio           On exit, this value will hold the ratio of the largest row scale factor to the
 *                                   smallest row scale factor.
 *
 * \param[out]    columnRatio        On exit, this value will hold the ratio of the lartest column scale factor to the
 *                                   smallest column scale factor.
 *
 * \param[out]    aMax               The absolute value of the largest element of A.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  Returns a positive value less than the number of
 *         rows indicating the row number of the first row that is exactly 0.
 */
typedef MatInteger (*MatLapackComplexEquilibrate)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const MatComplex* a,
    MatNumberElements lda,
    double*           rowScaleFactors,
    double*           columnScaleFactors,
    double*           rowRatio,
    double*           columnRatio,
    double*           aMax
);

/**
 * Function that calculates matrix row/column scale factors to equilibrate a matrix to a power of 2.   See LAPACK
 * dgeequb.
 *
 * \param[in]     isColumnMajor      If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                   operation should assume a row major matrix.
 *
 * \param[in]     numberRows         The number of rows in the input matrix.
 *
 * \param[in]     numberColumns      The number of columns in the input matrix.
 *
 * \param[in,out] a                  The input matrix A.  The array will be destroyed on exit.
 *
 * \param[in]     lda                Number of major elements between minor elements for matrix A.  The value is used
 *                                   to determine the stride.  In column major mode, this will represent the number of
 *                                   rows of space between columns which should always be greater than or equal to the
 *                                   number of rows.  In row major mode, this will represent the number of columns of
 *                                   space between rows which will always be greater than or equal to the number of
 *                                   columns.
 *
 * \param[in,out] rowScaleFactors    On exit, this array holding scale factors for each row.  The array will be
 *                                   allocated on entry.
 *
 * \param[in,out] columnScaleFactors On exit, this array holding scale factors for each column.  The array will be
 *                                   allocated on entry.
 *
 * \param[out]    rowRatio           On exit, this value will hold the ratio of the largest row scale factor to the
 *                                   smallest row scale factor.
 *
 * \param[out]    columnRatio        On exit, this value will hold the ratio of the lartest column scale factor to the
 *                                   smallest column scale factor.
 *
 * \param[out]    aMax               The absolute value of the largest element of A.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  Returns a positive value less than the number of
 *         rows indicating the row number of the first row that is exactly 0.
 */
typedef MatInteger (*MatLapackComplexEquilibratePowerOf2)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const MatComplex* a,
    MatNumberElements lda,
    double*           rowScaleFactors,
    double*           columnScaleFactors,
    double*           rowRatio,
    double*           columnRatio,
    double*           aMax
);

/**
 * Function that calculates the solution to a system of linear equations.  Solves for \f$ X \f$ given \f$ A X = B \f$.
 * See LAPACK zcgesv.
 *
 * \param[in]     isColumnMajor        If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                     operation should assume a row major matrix.
 *
 * \param[in]     matrixAOrder         The order of matrix A.
 *
 * \param[in]     numberRightHandSides The number of terms in the right hand side of the equation.  This value is also
 *                                     the same as the number of columns in matrix B.
 *
 * \param[in,out] a                    The square coefficient matrix.  The matrix will be destroyed on exit.
 *
 * \param[in]     lda                  Number of major elements between minor elements for matrix A.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[in,out] pivotData            An array holding matrixOrder pivot data terms.  The array will be allocated on
 *                                     entry.
 *
 * \param[in]     b                    An array of size ldb * numberRightHandSides for column major or
 *                                     lbd * matrixOrder for row major holding the B matrix terms.
 *
 * \param[in]     ldb                  Number of major elements between minor elements for matrix B.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[in,out] x                    Array to hold the resulting X matrix values.  The array will be allocated on
 *                                     entry.
 *
 * \param[in]     ldx                  Number of major elements between minor elements for matrix X.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[out]    iterationDetails     A value indicating the result of the solver.  The value is currently ignored.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  A positive value indicates the upper triangular
 *         matrix in the PLU decomposition was singular so a solution could not be computed.
 */
typedef MatInteger (*MatLapackComplexSolve)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements matrixOrder,
    MatNumberElements numberRightHandSides,
    MatComplex*       a,
    MatNumberElements lda,
    MatInteger*       pivotData,
    const MatComplex* b,
    MatNumberElements ldb,
    MatComplex*       x,
    MatNumberElements ldx,
    MatInteger*       iterationDetails
);

/**
 * Function that calculates the solution to an over-determined or under-determined system of linear equations.  Solves
 * for \f$ X \f$ given \f$ A X = B \f$.  See LAPACK zgels.
 *
 * \param[in]     isColumnMajor        If non-zero, the operation should assume a column major matrix.  If zero, the
 *                                     operation should assume a row major matrix.
 *
 * \param[in]     conjugateTransposed  Flag indicating if the A matrix is conjugate-transposed.  A non-zero value
 *                                     indicates conjugate-transposed.
 *
 * \param[in]     numberRows           The number of rows in matrix A.
 *
 * \param[in]     numberColumns        The number of columns in matrix A.
 *
 * \param[in]     numberRightHandSides The number of terms in the right hand side of the equation.  This value is also
 *                                     the same as the number of columns in matrix B.
 *
 * \param[in,out] a                    The coefficient matrix.  The matrix will be destroyed on exit.
 *
 * \param[in]     lda                  Number of major elements between minor elements for matrix A.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \param[in,out] b                    An array of size ldb * numberRightHandSides for column major or
 *                                     lbd * matrixOrder * max(numberRows, numberColumns) for row major holding the B
 *                                     matrix terms.  On exit, this array will hold either the least squares fit
 *                                     solution or a minimum norm solution (if the a matrix is not full rank).
 *
 * \param[in]     ldb                  Number of major elements between minor elements for matrix B.  The value is used
 *                                     to determine the stride.  In column major mode, this will represent the number
 *                                     of rows of space between columns which should always be greater than or equal to
 *                                     the number of rows.  In row major mode, this will represent the number of
 *                                     columns of space between rows which will always be greater than or equal to the
 *                                     number of columns.
 *
 * \return Returns 0 on success.  Returns a negative value on error.  A positive value indicates a matrix is not full
 *         rank so a least squares solution could not be determined.  In this scenario, b will contain the minimum
 *         norm solution vectors.
 */
typedef MatInteger (*MatLapackComplexLeastSquaresSolve)(
    MatMatrixMode     isColumnMajor,
    MatOperation      conjugateTransposed,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    MatNumberElements numberRightHandSides,
    MatComplex*       a,
    MatNumberElements lda,
    MatComplex*       b,
    MatNumberElements ldb
);

/**
 * Function that (re)balances a matrix to improve accuracy of eigenvalue/eigenvector calculations.  See LAPACK zgebal.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     job           The desired balancing job.
 *
 * \param[in]     matrixOrder   The matrix order (rows/columns) of A.
 *
 * \param[in,out] a             On entry, this array will hold the matrix to be balanced.  On exit, this array will
 *                              hold the newly balanced matrix.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[out]    ilo           The calculated LAPACK zgebal ilo value.
 *
 * \param[out]    ihi           The calculated LAPACK zgebal ihi value.
 *
 * \param[in,out] scaleTerms    Array returning the calculated permutation and scaling factors.
 *
 * \return Return 0 on success.  Returns a negative value on error.
 */
typedef MatInteger (*MatLapackComplexEigenBalance)(
    MatMatrixMode      isColumnMajor,
    MatEigenBalanceJob job,
    MatNumberElements  matrixOrder,
    MatComplex*        a,
    MatNumberElements  lda,
    MatInteger*        ilo,
    MatInteger*        ihi,
    double*            scaleTerms
);

/**
 * Function that calculates the eigenvalues and Schur factorization of a matrix.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     schurMode     The Schur computation mode (combines job and compz parameter).
 *
 * \param[in]     matrixOrder   The matrix order (rows/columns) of H.
 *
 * \param[out]    ilo           The calculated LAPACK zgebal ilo value.
 *
 * \param[out]    ihi           The calculated LAPACK zgebal ihi value.
 *
 * \param[in,out] h             On entry, this array will hold the Hessenberg matrix to use to calculate the
 *                              eigenvalues and Schur form from.
 *
 * \param[in]     ldh           Number of major elements between minor elements for matrix H.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] w             Array that will hold the computed eigenvalues on exit.  The array will be initialized
 *                              to hold at least matrixOrder terms on entry.
 *
 * \param[in,out] z             If schurMode == MatSchurMode::SCHUR_WITH_Q_MATRIX then this matrix should hold the
 *                              Q matrix from the reduction to Hessenberg form.
 *
 *                              For schurMode == MatSchurMode::SHUR_WITHOUT_Q_MATRIX, this array must be pre-allocated
 *                              to hold the Z matrix.
 *
 *                              On exit, this array will hold either the Z matrix or the Q * Z matrix depending on the
 *                              value of schurMode.  The array is not referenced if
 *                              schurMode == MatSchurMode::EIGENVALUES_ONLY.
 *
 * \param[in]     ldz           Number of major elements between minor elements for matrix Z.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \return Return 0 on success.  Returns a negative value on error.  A positive value indicates that this function was
 *         unable to compute all the eigenvalues of H.
 */
typedef MatInteger (*MatLapackComplexSchur)(
    MatMatrixMode      isColumnMajor,
    MatSchurMode       schurMode,
    MatNumberElements  matrixOrder,
    MatInteger         ilo,
    MatInteger         ihi,
    MatComplex*        h,
    MatNumberElements  ldh,
    MatComplex*        w,
    MatComplex*        z,
    MatNumberElements  ldz
);

/**
 * Function that restores a matrix (re)balanced by MatLapackComplexEigenBalance (LAPACK zgebal) back to an original,
 * unbalanced, non-symmetric matrix.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     job           The desired balancing job.
 *
 * \param[in]     side          The desired side (left eigenvectors or right eigenvectors).
 *
 * \param[in]     numberRows    The number of rows of the matrix of eigenvectors.
 *
 * \param[in]     numberColumns The number of columns of the matrix of eigenvectors.
 *
 * \param[out]    ilo           The calculated LAPACK zgebal ilo value.
 *
 * \param[out]    ihi           The calculated LAPACK zgebal ihi value.
 *
 * \param[in]     scaleTerms    Array holding the calculated permutation and scaling factors.
 *
 * \param[in,out] v             Matrix holding the eigenvectors to be transformed.  Returns the transformed
 *                              eigenvectors on exit.
 *
 * \param[in]     ldv           Number of major elements between minor elements for matrix V.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \return Return 0 on success.  Returns a negative value on error.
 */
typedef MatInteger (*MatLapackComplexEigenUnbalance)(
    MatMatrixMode      isColumnMajor,
    MatEigenBalanceJob job,
    MatSide            side,
    MatNumberElements  numberRows,
    MatNumberElements  numberColumns,
    MatInteger         ilo,
    MatInteger         ihi,
    double*            scaleTerms,
    MatComplex*        v,
    MatNumberElements  ldv
);

/**
 * Function that calculates eigenvectors of an upper triangular matrix computed via Schur decomposition.  See LAPACK
 * ztrevc.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     side          The desired side (left eigenvectors or right eigenvectors).
 *
 * \param[in]     howMany       Indicates the mode used to calculate eigenvectors.
 *
 * \param[in,out] select        Array of eigenvectors to be computed.  A null pointer can be used when this parameter
 *                              is not needed.  This array is modified on exit.
 *
 * \param[in]     matrixOrder   The order ot the supplied T matrix.
 *
 * \param[in,out] t             On entry, this value holds the T upper triangular matrix.Schur matrix. This method
 *                              may modify the contents of T.
 *
 * \param[in]     ldt           Number of major elements between minor elements for matrix T.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     vl            On entry, this matrix holds the Q matrix.  If the side is set to MatSide::LEFT, then
 *                              this array will hold the eigenvectors on exit.
 *
 * \param[in]     ldvl          Number of major elements between minor elements for matrix Vl.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     vr            On entry, this matrix holds the Q matrix.  This array is not reference if side is set
 *                              to MatSide::LEFT.  If MatSide is set to RIGHT, then this array will how the right
 *                              eigenvectors on exit.
 *
 * \param[in]     ldvr          Number of major elements between minor elements for matrix Vr.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     mm            The number of eigenvectors to be returned.  Set to either matrixOrder or the number of
 *                              eintries in select.
 *
 * \param[out]    m             The returned number of selected eigenvectors.
 *
 * \return Returns 0 on success.  Returns a negative value on error.
 */
typedef MatInteger (*MatLapackComplexEigenvectors)(
    MatMatrixMode      isColumnMajor,
    MatSide            side,
    MatEigenMode       howMany,
    const MatLogical*  select,
    MatNumberElements  matrixOrder,
    MatComplex*        t,
    MatNumberElements  ldt,
    MatComplex*        vl,
    MatNumberElements  ldvl,
    MatComplex*        vr,
    MatNumberElements  ldvr,
    MatNumberElements  mm,
    MatInteger*        m
);

/**
 * Function that returns LAPACK machine parameters.
 *
 * \param[in] parameter The desired machine parameter.  Currently only MatMachineParameter::SAFE_MINIMUM is used.
 *
 * \return Returns the requested parameter.
 */
typedef double (*MatLapackDoubleMachineParameter)(MatMachineParameter parameter);

/**
 * Function that performs scaling and transposition or copying of a real matrix.  Function is similar to the MKL
 * mkl_domatcopy function.
 *
 * This function performs \f$ B = \alpha op \left ( A \right ) \f$.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     operation     The desired operation.
 *
 * \param[in]     numberRows    The number of input matrix rows.
 *
 * \param[in]     numberColumns The number of input matrix columns.
 *
 * \param[in]     alpha         Scale factor to apply to A or A transpose.
 *
 * \param[in]     a             Pointer to the matrix, A, data.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] result        Pointer to the result matrix, B.  Space will be allocated prior to calling this
 *                              function.
 *
 * \param[in]     ldResult      Number of major elements between minor elements for matrix B.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 */
typedef void (*MatDoubleScaleCopy)(
    MatMatrixMode     isColumnMajor,
    MatOperation      operation,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    double            alpha,
    const double*     a,
    MatNumberElements lda,
    double*           result,
    MatNumberElements ldResult
);

/**
 * Function that performs scaling, transposition on two real matrices, adding the results.  Results are placed in a
 * third matrix.  Function is similar to the MKL mkl_domatadd function.
 *
 * This function performs \f$ C = \alpha op _ A \left ( A \right ) + \beta op _ B \left ( B \right ) \f$.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     operationA    The desired operation on matrix A.
 *
 * \param[in]     operationB    The desired operation on matrix B.
 *
 * \param[in]     numberRows    The number of input matrix rows.
 *
 * \param[in]     numberColumns The number of input matrix columns.
 *
 * \param[in]     alpha         Scale factor to apply to A or A transpose.
 *
 * \param[in]     a             Pointer to the matrix A data.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     beta          Scale factor to apply to B or B transpose.
 *
 * \param[in]     b             Pointer to the matrix B data.
 *
 * \param[in]     ldb           Number of major elements between minor elements for matrix B.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] result        Pointer to the result matrix C.  Space will be allocated prior to calling this
 *                              function.
 *
 * \param[in]     ldResult      Number of major elements between minor elements for matrix C.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 */
typedef void (*MatDoubleScaleAdd)(
    MatMatrixMode     isColumnMajor,
    MatOperation      operationA,
    MatOperation      operationB,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    double            alpha,
    const double*     a,
    MatNumberElements lda,
    double            beta,
    const double*     b,
    MatNumberElements ldb,
    double*           result,
    MatNumberElements ldResult
);

/**
 * Function that calculates the per-element linear fraction of vectors.  This function calculates
 * \f$ z _ i = \frac{m _ x x _ i + b _ x}{m _ y y _ i + b _ y} \f$
 *
 * \param[in]     numberElements The number of elements to calculate over.
 *
 * \param[in]     x              Input vector x.
 *
 * \param[in]     y              Input vector y.
 *
 * \param[in]     mx             Scale multiplier for x
 *
 * \param[in]     bx             Offset value for x
 *
 * \param[in]     my             Scale multiplier for y
 *
 * \param[in]     by             Offset value for y
 *
 * \param[in,out] z              Result values.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as either input vector.
 */
typedef void (*MatDoubleLinearFraction)(
    MatNumberElements numberElements,
    const double*     x,
    const double*     y,
    double            mx,
    double            bx,
    double            my,
    double            by,
    double*           z
);

/**
 * Function that fills a vector with a value.
 *
 * \param[in]     numberElements The number of elements to be populated.
 *
 * \param[in]     value          The value to fill each element with.
 *
 * \param[in,out] destination    The buffer to hold the destination.  The buffer must be pre-allocated.
 */
typedef void(*MatIntegerFill)(MatNumberElements numberElements, MatInteger value, MatInteger* destination);

/**
 * Function that fills a vector with a value.
 *
 * \param[in]     numberElements The number of elements to be populated.
 *
 * \param[in]     value          The value to fill each element with.
 *
 * \param[in,out] destination    The buffer to hold the destination.  The buffer must be pre-allocated.
 */
typedef void(*MatDoubleFill)(MatNumberElements numberElements, double value, double* destination);

/**
 * Function that fills a vector with a value.
 *
 * \param[in]     numberElements The number of elements to be populated.
 *
 * \param[in]     value          The value to fill each element with.
 *
 * \param[in,out] destination    The buffer to hold the destination.
 */
typedef void(*MatComplexFill)(MatNumberElements numberElements, const MatComplex& value, MatComplex* destination);

/**
 * Function that performs per-element multiply and add against scalar quantities.
 * \f$ z _ i = m x _ i + b \f$
 *
 * \param[in]     numberElements The number of elements to calculate over.
 *
 * \param[in]     x              Input vector x.
 *
 * \param[in]     mx             Scale multiplier for x
 *
 * \param[in]     bx             Offset value for x
 *
 * \param[in,out] z              Result values.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as either input vector.
 */
typedef void (*MatDoubleVectorScalarMultiplyAdd)(
    MatNumberElements numberElements,
    const double*     x,
    double            m,
    double            b,
    double*           z
);

/**
 * Function that calculates the floor of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the floor of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function. The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleFloor)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the Ceiling of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the ceiling of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function. The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleCeiling)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the closest integer for each element of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the nearest integers of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleNearestInteger)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the per-element natural log of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the natural log of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleLog)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the per-element exponential of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the exponential of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleExponential)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the per-element square root of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the per-element square root of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleSquareRoot)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the per-element sine of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the per-element sine of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleSine)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the per-element cosine of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the per-element cosine of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleCosine)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that calculates the per-element tangent of a double precision vector.
 *
 * \param[in]     numberElements The number of elements to calculate the per-element cosine of.
 *
 * \param[in]     source         Array of input vectors.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleTangent)(MatNumberElements numberElements, const double* source, double* destination);

/**
 * Function that packs values from one vector into another vector.
 *
 * \param[in]     numberElements The number of elements to be packed.  This determines the length of the result.
 *
 * \param[in]     source         Input vector.
 *
 * \param[in]     increment      The source increment, in elements.  A value of 1 will copy elements.
 *
 * \param[in,out] destination    Result vector.  The array will be allocated before calling this function.
 */
typedef void (*MatDoublePack)(
    MatNumberElements numberElements,
    const double*     source,
    MatNumberElements increment,
    double*           destination
);

/**
 * Function that calculates the dot product of two double precision vectors.
 *
 * \param[in]     numberElements The number of elements to calculate the natural log of.
 *
 * \param[in]     a              The first vector to multiply.
 *
 * \param[in]     b              The second vector to multiply.
 *
 * \param[in,out] result         Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatDoubleDotProduct)(
    MatNumberElements numberElements,
    const double*     a,
    const double*     b,
    double*           result
);

/**
 * Function that performs scaling and transposition or copying of a complex matrix.  Function is similar to the MKL
 * mkl_zomatcopy function.
 *
 * This function performs \f$ B = \alpha op \left ( A \right ) \f$.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     operation     The desired operation.
 *
 * \param[in]     numberRows    The number of input matrix rows.
 *
 * \param[in]     numberColumns The number of input matrix columns.
 *
 * \param[in]     alpha         Scale factor to apply to A or A transpose.
 *
 * \param[in]     a             Pointer to the matrix, A, data.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] result        Pointer to the result matrix, B.  Space will be allocated prior to calling this
 *                              function.
 *
 * \param[in]     ldResult      Number of major elements between minor elements for matrix B.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 */
typedef void (*MatComplexScaleCopy)(
    MatMatrixMode     isColumnMajor,
    MatOperation      operation,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const MatComplex* alpha,
    const MatComplex* a,
    MatNumberElements lda,
    MatComplex*       result,
    MatNumberElements ldResult
);

/**
 * Function that performs scaling, transposition on two complex matrices, adding the results.  Results are placed in a
 * third matrix.  Function is similar to the MKL mkl_zomatadd function.
 *
 * This function performs \f$ C = \alpha op _ A \left ( A \right ) + \beta op _ B \left ( B \right ) \f$.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     operationA    The desired operation on matrix A.
 *
 * \param[in]     operationB    The desired operation on matrix B.
 *
 * \param[in]     numberRows    The number of input matrix rows.
 *
 * \param[in]     numberColumns The number of input matrix columns.
 *
 * \param[in]     alpha         Scale factor to apply to A or A transpose.
 *
 * \param[in]     a             Pointer to the matrix A data.
 *
 * \param[in]     lda           Number of major elements between minor elements for matrix A.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in]     beta          Scale factor to apply to B or B transpose.
 *
 * \param[in]     b             Pointer to the matrix B data.
 *
 * \param[in]     ldb           Number of major elements between minor elements for matrix B.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 * \param[in,out] result        Pointer to the result matrix C.  Space will be allocated prior to calling this
 *                              function.
 *
 * \param[in]     ldResult      Number of major elements between minor elements for matrix C.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 */
typedef void (*MatComplexScaleAdd)(
    MatMatrixMode     isColumnMajor,
    MatOperation      operationA,
    MatOperation      operationB,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const MatComplex* alpha,
    const MatComplex* a,
    MatNumberElements lda,
    const MatComplex* beta,
    const MatComplex* b,
    MatNumberElements ldb,
    MatComplex*       result,
    MatNumberElements ldResult
);

/**
 * Function that calculates the dot product of two complex vectors.
 *
 * \param[in]     numberElements The number of elements to calculate the natural log of.
 *
 * \param[in]     a              The first vector to multiply.
 *
 * \param[in]     b              The second vector to multiply.
 *
 * \param[in,out] result         Result vector.  The array will be allocated before calling this function.  The
 *                               destination array can be the same as the source array.
 */
typedef void (*MatComplexDotProduct)(
    MatNumberElements numberElements,
    const MatComplex* a,
    const MatComplex* b,
    MatComplex*       result
);

/**
 * Function that converts a vector of 64-bit integer to a floating point value between [0, 1].
 *
 * \param[in]     numberElements Number of elements to be converted.
 *
 * \param[in,out] p              Pointer to a vector of 64-bit unsigned integer values to be converted.  Values should
 *                               range from 0 to 0xFFFFFFFFFFFFFFFF.  This vector will contain floating point values
 *                               on exit.
 */
typedef void (*MatIntegerToFloatInclusive)(MatNumberElements numberElements, void* p);

/**
 * Function that converts a vector of 64-bit integer to a floating point value between (0, 1].
 *
 * \param[in]     numberElements Number of elements to be converted.
 *
 * \param[in,out] p              Pointer to a vector of 64-bit unsigned integer values to be converted.  Values should
 *                               range from 0 to 0xFFFFFFFFFFFFFFFF.  This vector will contain floating point values
 *                               on exit.
 */
typedef void (*MatIntegerToFloatExclusiveInclusive)(MatNumberElements numberElements, void* p);

/**
 * Function that converts a vector of 64-bit integer to a floating point value between [0, 1).
 *
 * \param[in]     numberElements Number of elements to be converted.
 *
 * \param[in,out] p              Pointer to a vector of 64-bit unsigned integer values to be converted.  Values should
 *                               range from 0 to 0xFFFFFFFFFFFFFFFF.  This vector will contain floating point values
 *                               on exit.
 */
typedef void (*MatIntegerToFloatInclusiveExclusive)(MatNumberElements numberElements, void* p);

/**
 * Function that converts a vector of 64-bit integer to a floating point value between (0, 1).
 *
 * \param[in]     numberElements Number of elements to be converted.
 *
 * \param[in,out] p              Pointer to a vector of 64-bit unsigned integer values to be converted.  Values should
 *                               range from 0 to 0xFFFFFFFFFFFFFFFF.  This vector will contain floating point values
 *                               on exit.
 */
typedef void (*MatIntegerToFloatExclusive)(MatNumberElements numberElements, void* p);

/**
 * Function that returns a true random number.  This function must be reentrant.
 *
 * \return Returns a true random value.
 */
typedef uint32_t (*MatTrueRandomValue)();

/**
 * Function that returns an array of true random values.  This function must be reentrant.
 *
 * \param[in,out] array The array to contain the data.  The array is expected to be allocated on entry.
 *
 * \param[in]     numberTerms The number of values to be read.
 */
typedef void (*MatTrueRandomArray)(uint32_t* array, unsigned long numberTerms);

/**
 * Function that specifies the space required for the MT216091 PRNG value array, in 64-bit values.
 *
 * \return Returns the space needed for the MT216091 value array, in 64-bit values.
 */
typedef unsigned (*MatMt216091ValueArraySize)();

/**
 * Function that seeds a new MT216091 value array.
 *
 * \param[in,out] valueArray  The MT216091 value array to be seeded.  The size of the array in 64-bit values must match
 *                            the value returned by MatMt216091ValueArraySize.  The array should also be allocated
 *                            aligned for SIMD instructions using the MatAllocatorFunction.
 *
 * \param[in]     seeds       Array of seeds.
 *
 * \param[in]     numberSeeds The number of seed entries.
 */
typedef void (*MatMt216091Seed)(uint64_t* valueArray, const uint64_t* seeds, unsigned numberSeeds);

/**
 * Function that updates an MT216091 value array with new values.   The old array should be used as an input.
 *
 * \param[in,out] valueArray  The MT216091 value array to be seeded.  The size of the array in 64-bit values must match
 *                            the value returned by MatMt216091ValueArraySize.  The array should also be allocated
 *                            aligned for SIMD instructions using the MatAllocatorFunction.
 */
typedef void (*MatMt216091Update)(uint64_t* valueArray);

/**
 * Function that calculates a forward DFT of a one or two dimensional matrix.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     numberRows    The number of input rows.
 *
 * \param[in]     numberColumns The number of input columns.
 *
 * \param[in]     i             The input data.  Data can be assumed to be memory consecutive for row and column
 *                              matrices.
 *
 * \param[in]     ldi           Number of major elements between minor elements for matrix I.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 *                              This value is ignored for row/column matrices.
 *
 * \param[in,out] o             The output data. The array will be allocated prior to calling this method.
 *
 * \param[in]     ldo           Number of major elements between minor elements for matrix O.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 *                              This value is ignored for row/column matrices.
 *
 * \return Returns 0 on success.  Returns non-zero on failure.
 */
typedef MatInteger (*MatComplexDft)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const MatComplex* i,
    MatNumberElements ldi,
    MatComplex*       o,
    MatNumberElements ldo
);

/**
 * Function that calculates a backward (inverse) DFT of a one or two dimensional matrix.
 *
 * \param[in]     isColumnMajor If non-zero, the operation should assume a column major matrix.  If zero, the operation
 *                              should assume a row major matrix.
 *
 * \param[in]     numberRows    The number of input rows.
 *
 * \param[in]     numberColumns The number of input columns.
 *
 * \param[in]     i             The input data.  Data can be assumed to be memory consecutive for row and column
 *                              matrices.
 *
 * \param[in]     ldi           Number of major elements between minor elements for matrix I.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 *                              This value is ignored for row/column matrices.
 *
 * \param[in,out] o             The output data. The array will be allocated prior to calling this method.
 *
 * \param[in]     ldo           Number of major elements between minor elements for matrix O.  The value is used to
 *                              determine the stride.  In column major mode, this will represent the number of rows of
 *                              space between columns which should always be greater than or equal to the number of
 *                              rows.  In row major mode, this will represent the number of columns of space between
 *                              rows which will always be greater than or equal to the number of columns.
 *
 *                              This value is ignored for row/column matrices.
 *
 * \return Returns 0 on success.  Returns non-zero on failure.
 */
typedef MatInteger (*MatComplexInverseDft)(
    MatMatrixMode     isColumnMajor,
    MatNumberElements numberRows,
    MatNumberElements numberColumns,
    const MatComplex* i,
    MatNumberElements ldi,
    MatComplex*       o,
    MatNumberElements ldo
);

/**
 * Function that calculates a forward (type 2) or the reverse (type 3) DCT of a vector.
 *
 * \param[in]     numberElements The number of input rows.
 *
 * \param[in]     i              The input data.  Data can be assumed to be memory consecutive for row and column
 *
 * \param[in,out] o              The output data. The array will be allocated prior to calling this method.
 *
 * \return Returns 0 on success.  Returns non-zero on failure.
 */
typedef MatInteger (*MatRealDct)(MatNumberElements numberElements, const double* i, double* o);

/**
 * Function that calculates a Hilbert transform of a vector.
 *
 * \param[in]     numberElements The number of input rows.
 *
 * \param[in]     i              The input data.  Data can be assumed to be memory consecutive for row and column
 *
 * \param[in,out] o              The output data. The array will be allocated prior to calling this method.
 *
 * \return Returns 0 on success.  Returns non-zero on failure.
 */
typedef MatInteger (*MatComplexHilbertTransform)(MatNumberElements numberElements, const double* i, MatComplex* o);

/**
 * Structure used to define entry points to useful BLAS, LAPACK, and similar matrix functions.  Most functions should
 * map directly to standard BLAS/LAPACK calls.
 */
typedef struct _MatApi {
    /*************************************************
     * Constants.
     */

    /**
     * Memory alignment constraint for this matrix library.  Matrices will be aligned to multiples of this value,
     * in bytes.  Value is used to force memory alignment of matrices for improved performance with SIMD instructions.
     * Set to 0 to disable forced memory alignment.
     */
    uint8_t memoryAlignmentRequirementBytes;

    /*************************************************
     * Memory management
     */

    /**
     * Pointer to the matrix allocator function.  The pointer should typically be mapped to the system malloc function.
     */
    MatAllocatorFunction allocateMemory;

    /**
     * Pointer to the matrix deallocator function.  The pointer should typically be mapped to the system free function.
     */
    MatDeallocatorFunction releaseMemory;

    /*************************************************
     * BLAS functions
     */

    /**
     * Vector copy (double precision)
     */
    MatBlasDoubleCopy blasDoubleCopy;

    /**
     * BLAS dgemm
     */
    MatBlasDoubleMultplyAdd blasDoubleMultiplyAdd;

    /**
     * Vector copy (complex)
     */
    MatBlasComplexCopy blasComplexCopy;

    /**
     * BLAS zgemm
     */
    MatBlasComplexMultplyAdd blasComplexMultiplyAdd;

    /*************************************************
     * LAPACK functions
     */

    /**
     * LAPACK dgetrf - Compute PLU factorization
     */
    MatLapackDoublePlu lapackDoublePlu;

    /**
     * LAPACK dgetri - Compute inverse from LU factorization.
     */
    MatLapackDoubleLuInverse lapackDoubleLuInverse;

    /**
     * LAPACK dgeqrf - Compute QR factorization of general matrix.
     */
    MatLapackDoubleQrFactorization lapackDoubleQrFactorization;

    /**
     * LAPACK dorgqr - Compute orthogonal Q matrix from QR factorization.
     */
    MatLapackDoubleGenerateQFromQrMatrix lapackDoubleGenerateQFromQrMatrix;

    /**
     * LAPACK dpotrf - Compute upper or lower Cholesky decomposition of a semi-definite Hermitian matrix.
     */
    MatLapackDoubleCholesky lapackDoubleCholesky;

    /**
     * LAPACK dgehrd - Compute upper Hessenberg of a dense matrix.
     */
    MatLapackDoubleUpperHessenberg lapackDoubleUpperHessenberg;

    /**
     * LAPACK dorghr - Compute orthogonal Q matrix from Hessenberg.
     */
    MatLapackDoubleUpperHessenbergQMatrix lapackDoubleUpperHessenbergQMatrix;

    /**
     * LAPACK degsvd - Compute singular value decomposition of a general matrix.
     */
    MatLapackDoubleSvd lapackDoubleSvd;

    /**
     * LAPACK dgeequ - Compute row/column scale factors to equilibrate a matrix.
     */
    MatLapackDoubleEquilibrate lapackDoubleEquilibrate;

    /**
     * LAPACK dgeequ - Compute row/column scale factors to equilibrate a matrix to a power of 2.
     */
    MatLapackDoubleEquilibratePowerOf2 lapackDoubleEquilibratePowerOf2;

    /**
     * LAPACK dsgesv - Solve systems of simultaneous equations.
     */
    MatLapackDoubleSolve lapackDoubleSolve;

    /**
     * LAPACK dgels - Solve least squares solution to an over-determined/under-determined system.
     */
    MatLapackDoubleLeastSquaresSolve lapackDoubleLeastSquaresSolve;

    /**
     * LAPACK zgetrf - Compute PLU factorization (complex)
     */
    MatLapackComplexPlu lapackComplexPlu;

    /**
     * LAPACK zgetri - Compute inverse from LU factorization (complex)
     */
    MatLapackComplexLuInverse lapackComplexLuInverse;

    /**
     * LAPACK dgeqrf - Compute QR factorization of general matrix (complex)
     */
    MatLapackComplexQrFactorization lapackComplexQrFactorization;

    /**
     * LAPACK dorgqr - Compute orthogonal Q matrix from QR factorization (complex)
     */
    MatLapackComplexGenerateQFromQrMatrix lapackComplexGenerateQFromQrMatrix;

    /**
     * LAPACK dpotrf - Compute upper or lower Cholesky decomposition of a semi-definite Hermitian matrix (complex).
     */
    MatLapackComplexCholesky lapackComplexCholesky;

    /**
     * LAPACK dgehrd - Compute upper Hessenberg of a dense matrix (complex).
     */
    MatLapackComplexUpperHessenberg lapackComplexUpperHessenberg;

    /**
     * LAPACK dorghr - Compute orthogonal Q matrix from Hessenberg (complex).
     */
    MatLapackComplexUpperHessenbergQMatrix lapackComplexUpperHessenbergQMatrix;

    /**
     * LAPACK degsvd - Compute singular value decomposition of a general matrix (complex).
     */
    MatLapackComplexSvd lapackComplexSvd;

    /**
     * LAPACK dgeequ - Compute row/column scale factors to equilibrate a matrix (complex).
     */
    MatLapackComplexEquilibrate lapackComplexEquilibrate;

    /**
     * LAPACK dgeequ - Compute row/column scale factors to equilibrate a matrix to a power of 2 (complex).
     */
    MatLapackComplexEquilibratePowerOf2 lapackComplexEquilibratePowerOf2;

    /**
     * LAPACK dsgesv - Solve systems of simultaneous equations (complex).
     */
    MatLapackComplexSolve lapackComplexSolve;

    /**
     * LAPACK dgels - Solve least squares solution to an over-determined/under-determined system (complex).
     */
    MatLapackComplexLeastSquaresSolve lapackComplexLeastSquaresSolve;

    /**
     * LAPACK zgebal - Balances a matrix prior to performing eigenvector/eigenvalue computation.
     */
    MatLapackComplexEigenBalance lapackComplexEigenBalance;

    /**
     * LAPACK zhseqr - Calculate Eigenvalues and Schur factorization of a Hessenberg matrix.
     */
    MatLapackComplexSchur lapackComplexSchur;

    /**
     * LAPACK zgebak - Restores a balanced matrix of eigenvectors back to its original, unbalanced form.
     */
    MatLapackComplexEigenUnbalance lapackComplexEigenUnbalance;

    /**
     * LAPACK ztrevc - Calculate eigenvectors from a Schur triangular matrix.
     */
    MatLapackComplexEigenvectors lapackComplexEigenvectors;

    /**
     * LAPACK dlamch - Returns double precision machine parameters.
     */
    MatLapackDoubleMachineParameter lapackDoubleMachineParameter;

    /*************************************************
     * Specialized functions
     */

    /**
     * Matrix scale & copy method, double precision.
     */
    MatDoubleScaleCopy doubleScaleCopy;

    /**
     * Matrix scale & add method, double precision.
     */
    MatDoubleScaleAdd doubleScaleAdd;

    /**
     * Vector linear fraction function, double precision.
     */
    MatDoubleLinearFraction doubleLinearFraction;

    /**
     * Vector fill, 64-bit integer.
     */
    MatIntegerFill integerFill;

    /**
     * Vector fill, double precision.
     */
    MatDoubleFill doubleFill;

    /**
     * Vector fill, complex.
     */
    MatComplexFill complexFill;

    /**
     * Vector per-element multiply & add against scaler, double precision.
     */
    MatDoubleVectorScalarMultiplyAdd doubleVectorScalarMultiplyAdd;

    /**
     * Vector dot product, double precision.
     */
    MatDoubleDotProduct doubleDotProduct;

    /**
     * Vector floor function, double precision.
     */
    MatDoubleFloor doubleFloor;

    /**
     * Vector floor function, double precision.
     */
    MatDoubleCeiling doubleCeiling;

    /**
     * Vector nearest integer function, double precision.
     */
    MatDoubleNearestInteger doubleNearestInteger;

    /**
     * Vector natural log function (per-element), double precision.
     */
    MatDoubleLog doubleLog;

    /**
     * Vector exponential function (per-element), double precision.
     */
    MatDoubleExponential doubleExponential;

    /**
     * Vector square root function (per-element), double precision
     */
    MatDoubleSquareRoot doubleSquareRoot;

    /**
     * Vector per-element sine, double precision.
     */
    MatDoubleSine doubleSine;

    /**
     * Vector per-element cosine, double precision.
     */
    MatDoubleCosine doubleCosine;

    /**
     * Vector per-element tangent, double precision.
     */
    MatDoubleTangent doubleTangent;

    /**
     * Vector packer, double precision.
     */
    MatDoublePack doublePack;

    /**
     * Matrix scale & copy method, complex.
     */
    MatComplexScaleCopy complexScaleCopy;

    /**
     * Matrix scale & add method, complex.
     */
    MatComplexScaleAdd complexScaleAdd;

    /**
     * Vector dot product, complex.
     */
    MatComplexDotProduct complexDotProduct;

    /**
     * Vector integer to float conversion over the range [0, 1]
     */
    MatIntegerToFloatInclusive integerToFloatInclusive;

    /**
     * Vector integer to float conversion over the range (0, 1].
     */
    MatIntegerToFloatExclusiveInclusive integerToFloatExclusiveInclusive;

    /**
     * Vector integer to float conversion over the range [0, 1).
     */
    MatIntegerToFloatInclusiveExclusive integerToFloatInclusiveExclusive;

    /**
     * Vector integer to float conversio over the range (0, 1).
     */
    MatIntegerToFloatExclusive integerToFloatExclusive;

    /*************************************************
     * RNG functions.
     */

    /**
     * Returns a true random number.
     *
     * \return Returns a true random value.
     */
    MatTrueRandomValue trueRandomValue;

    /**
     * Populates an array of true-random values.
     */
    MatTrueRandomArray trueRandomArray;

    /**
     * Determines the size of the MT216091 value array, in 64-bit entries.
     */
    MatMt216091ValueArraySize mt216091ValueArraySize;

    /**
     * Function that seeds a new MT216091 value array.
     */
    MatMt216091Seed mt216091Seed;

    /**
     * Function that updates an MT216091 value array.
     */
    MatMt216091Update mt216091Update;

    /*************************************************
     * Transform functions.
     */

    /**
     * One and two dimensional DFT, complex.
     */
    MatComplexDft complexDft;

    /**
     * One and two dimensional inverse DFT, complex.
     */
    MatComplexInverseDft complexInverseDft;

    /**
     * Forward (type 2) DCT, real.
     */
    MatRealDct realDctType2;

    /**
     * Inverse (type 3) DCT, real.
     */
    MatRealDct realDctType3;

    /**
     * Complex Hilbert transform, complex.
     */
    MatComplexHilbertTransform complexHilbertTransform;
} MatApi;

#if (defined(__cplusplus))

    }

#endif

#endif
