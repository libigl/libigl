#define NO_TAUCS

#ifdef NO_TAUCS
  typedef double taucsType;
#  include <vector>
#  include <map>
#else
#  include "taucsaddon.h"
#endif

class SparseSolver {
protected:
	// is set to true iff m_A is symmetric positive definite
	bool  m_SPD;
        
#ifndef NO_TAUCS
	taucs_ccs_matrix * m_A;
	taucs_ccs_matrix * m_AT;
	taucs_ccs_matrix * m_ATA;
#endif

#ifndef NO_TAUCS
	void * m_factorATA;
	void * m_factorA;
	void * m_etree; // for the factor update
#endif
	
	// placeholder, so that we don't need to allocate space every time
	// the space is allocated when a factor for ATA is created
	std::vector<taucsType> m_ATb;

	// helper map to create A
	std::vector< std::map<int,taucsType> > m_colsA;

	// the dimensions of the A matrix
	int m_numRows;
	int m_numCols; 

public:

	SparseSolver(const int numRows = 0, const int numCols = 0, const bool isSPD = false) 
		: 
		m_SPD(isSPD)
                ,
#ifndef NO_TAUCS
                  m_A(NULL)
		, m_AT(NULL)
		, m_ATA(NULL)
		, m_factorATA(NULL)
		, m_factorA(NULL)
		, m_etree(NULL)
		, m_ATb(NULL)
		, 
#endif
		  m_colsA(numCols)
                , m_numRows(numRows)
		, m_numCols(numCols)
	{};

	~SparseSolver() { ClearObject();}

//Operators
public:
	/** Copies the matrix.

		However, factorizations are not copied.
		This is because the factorizations are actually owned by taucs
		and we have to call taucs to release them.
		If we would release them from two different objects,
		we would run into trouble.
	*/
	SparseSolver& operator=(const SparseSolver& other)
	{
		//Reset ourselves
		Reset(other.GetNumRows(), other.GetNumColumns());

		//Copy the sparse matrix
		m_colsA = other.m_colsA;
		m_SPD = other.m_SPD;
		m_numRows = other.m_numRows;
		m_numCols = other.m_numCols;

		return *this;
	}

//Methods
public:

	///Resets all internal data structures and sets the new size of this matrix
	void Reset(const int numRows, const int numCols);

	///Returns the number of rows of the \c A matrix.
	inline int GetNumRows() const {return m_numRows;}

	///Returns the number of columns of the \c A matrix.
	inline int GetNumColumns() const {return m_numCols;}

	///Returns the dimensions of the \c A matrix.
	inline void GetAMatrixSize(int & numRows, int & numCols) const
	{
		numRows = m_numRows;
		numCols = m_numCols;
	}

	///Returns true if no rows or columns are in the matrix.
	inline bool IsEmpty() const {return m_numRows == 0 || m_numCols == 0;}

	///Clears all data structures.
	void ClearObject();

	// adds a new entry to the sparse matrix A
	// i and j are 0-based
	// if the A^T A matrix had been factored, this will destroy the factor
	void AddIJV(const int i, const int j, const taucsType v);

	// updates the entry in A by adding the value v to the current value
	// i and j are 0-based
	// if the A^T A matrix had been factored, this will destroy the factor
	void AddToIJV(const int i, const int j, const taucsType v);

        // updates the entries in A by adding the values v in B to the current
        // Value
	// if the A^T A matrix had been factored, this will destroy the factor
	void AddMatrix(SparseSolver & B);
        // updates the entries in A by adding the values v in B to the current
        // Value in A offset by (i,j)
        // like matlab's A[i:m,j:n] = B, where [m,n] = size(B);
        // i and j are 0-based
	// if the A^T A matrix had been factored, this will destroy the factor
	void AddMatrix(const int i, const int j, SparseSolver & B);

	///Removes all entries from the given row.
	///If the A^T A matrix had been factored, this will destroy the factor.
	inline bool ClearRow(const int iRow)
	{
		if (iRow >= 0 && iRow < m_numRows)
		{
#ifndef NO_TAUCS
			ClearFactorATA();
			ClearFactorA();
			ClearMatricesATA();
			ClearMatricesA();
#endif

			//Go over all columns, try to find the respective row and delete it
			for(std::vector< std::map<int,taucsType> >::iterator itCol=m_colsA.begin();itCol!=m_colsA.end();itCol++)
			{
				//Find the row
				std::map<int,taucsType>::iterator itRow = itCol->find(iRow);
				if (itRow != itCol->end())
				{
					itCol->erase(itRow);
				}
			}

			return true;
		}

		return false;
	}

	///Removes all entries from the given column.
	///If the A^T A matrix had been factored, this will destroy the factor.
	inline bool ClearColumn(const int iColumn)
	{
		if (iColumn >= 0 && iColumn < m_numCols)
		{
#ifndef NO_TAUCS
			ClearFactorATA();
			ClearFactorA();
			ClearMatricesATA();
			ClearMatricesA();
#endif

			m_colsA[iColumn].clear();
			return true;
		}

		return false;
	}

	// Sets the status of a matrix - whether SPD or not
	void SetSPD(bool isSPD) {m_SPD = isSPD;}
        // Return true only if symmetric, automatically true if isSPD flag is
        // set to true, otherwise actually determine if matrix is symmetric
        bool IsSymmetric();

#ifndef NO_TAUCS
	// allows to add an anchored vertex without destroying the factor, if there was one
	// i is the anchor's number (i.e. the index of the mesh vertex that is anchored is i)
	// w is the weight of the anchor in the original Ax=b system. HAS TO BE POSITIVE!!
	// if there was no factor for A^T A, this will simply add an anchor row to the Ax=b system
	// if there was a factor, the factor will be updated, approximately at the cost of
	// a single solve.
	void AddAnchor(const int i, const taucsType w);

	// Will create A^T A and its factorization, as well as store A^T
	// returns true on success, false otherwise
	bool FactorATA();

	// Will create A and its factorization. If A is spd then Cholesky factor, otherwise LU
	// returns true on success, false otherwise
	// Assumes A is square and invertible
	bool FactorA();
	// Will create *symbolic* factorization of ATA and return it; returns NULL on failure
	void * SymbolicFactorATA();

	taucs_ccs_matrix * GetATA_Copy();

	// Will create the actual numerical factorization of the present ATA matrix
	// with the same non-zero pattern as created by SymbolicFactorATA()
	bool FactorATA_UseSymbolic(void * symbFactor,
	  taucs_ccs_matrix ** PAP,
	  int * perm,
	  int * invperm);

	// Solves the normal equations for Ax = b, namely A^T A x = A^T b
	// it is assumed that enough space is allocated in x
	// Uses the factorization provided in factor
	// numRhs should be the the number of right-hand sides stored in b (and respectively,
	// there should be enough space allocated in the solution vector x)
	// returns true on success, false otherwise
	bool SolveATA_UseFactor(void * factor, 
							taucs_ccs_matrix * PAP,
							int * perm,
							int * invperm,
							const taucsType * b, taucsType * x, const int numRhs);

	// Solves the normal equations for Ax = b, namely A^T A x = A^T b
	// it is assumed that enough space is allocated in x
	// If a factorization of A^T A exists, it will be used, otherwise it will be created
	// numRhs should be the the number of right-hand sides stored in b (and respectively,
	// there should be enough space allocated in the solution vector x)
	// returns true on success, false otherwise
	bool SolveATA(const taucsType * b, taucsType * x, const int numRhs); 

	// Same as SolveATA only that this solves the actual system Ax = b
	// It is assumed that A is rectangular and invertible
	bool SolveA(const taucsType * b, taucsType * x, const int numRhs); 

	// The version of SolveATA with 3 right-hand sides
	// returns true on success, false otherwise
	bool SolveATA3(const taucsType * bx, const taucsType * by, const taucsType * bz,
						 taucsType * x,        taucsType * y,        taucsType * z);

	// The version of SolveATA with 2 right-hand sides
	// returns true on success, false otherwise
	bool SolveATA2(const taucsType * bx, const taucsType * by,
						 taucsType * x,        taucsType * y);
#endif

// Matrix utility routines

	// Stores the result of Matrix*v in result.
	// v can have numCols columns; then the result is stored column-wise as well
	// It is assumed space is allocaed in result!
	// v cannot be the same pointer as result
	void MultiplyMatrixVector(const taucsType * v, taucsType * result, const int numCols) const;

	// Multiplies the matrix by diagonal matrix D from the left, stores the result
	// in the columns of res
	void MultiplyDiagMatrixMatrix(const taucsType * D, SparseSolver & res) const;

	// Multiplies the matrix by the given matrix B from the right
	// and stores the result in the columns of res
	// isSPD tells us whether the result should be SPD or not
	// retuns false if the dimensions didn't match (in such case res is unaltered), otherwise true
	bool MultiplyMatrixRight(const SparseSolver & B, SparseSolver & res,const bool isSPD) const;

	// Transposes the matrix 
	// and stores the result in the columns of res
	bool Transpose(SparseSolver & res) const;

	// will remove the rows startInd-endInd (including ends, zero-based) from the matrix
	// will discard any factors	
	bool DeleteRowRange(const int startInd, const int endInd);

	// will remove the columns startInd-endInd (including ends, zero-based) from the matrix
	// will discard any factors	
	bool DeleteColumnRange(const int startInd, const int endInd);

	// copies columns startInd through endInd (including, zero-based) from the matrix to res matrix
	// previous data in res, including factors, is erased.
	bool CopyColumnRange(const int startInd, const int endInd, SparseSolver & res) const;

	// copies rows startInd through endInd (including, zero-based) from the matrix to res matrix
	// previous data in res, including factors, is erased.
	bool CopyRowRange(const int startInd, const int endInd, SparseSolver & res) const;

	///Multiplies the given range of rows (zero-based, including startInd and endInd) with the given scalar.
	bool MultiplyRows(const int startInd, const int endInd, const taucsType val);

	///Multiplies all entries of the matrix with the given scalar. Will discard any factors.
	void MultiplyMatrixScalar(const taucsType val);

	/** Computes IdentityMatrix minus this matrix in-place. Will discard any factors.

		If an entry on the diagonal is missing, it will be created.
	*/
	void IdentityMinusMatrix();

	/** Computes IdentityMatrix plus this matrix in-place. Will discard any factors.

		If an entry on the diagonal is missing, it will be created.
	*/
	void IdentityPlusMatrix();

        // Alec:
        // Note on splice methods:
        // The following splice methods act like matlabs splicing and
        // reordering via index lists and the operator()
        // Both SpliceColumns and SpliceRows are ignorant of SPD, meaning that
        // they will not give "correct" splices for SPD matrices: splices of
        // This only matters if you are not setting all
        // SPD matrices are just splices of the lower triangle
        //
        // This only matters if you are not setting the upper half of the
        // matrix. To get correct splices you need to set both upper and lower
        // triangles before splicing.

        // Alec:
        // Copies certain cols in any order to a new matrix res
        // Like in Matlab B = A(:,column_indices)
        // Previous data in res is erased.
        // Returns false on error. As in column_indices are invalid
        // Result is always NOT SPD
        //
        // WARNING: If your matrix is only the lower half of an SPD matrix
        // then the result will only be a copy of the lower half of the 
        // spliced columns. See note above (Note on splice methods)
        bool SpliceColumns(
            const std::vector<size_t> column_indices, 
            SparseSolver & res) const;

        // Alec:
        // Copies certain rows in any order to a new matrix res
        // Like in Matlab B = A(:,row_indices)
        // Previous data in res is erased.
        // Returns false on error. As in row_indices are invalid
        // Result is always NOT SPD
        //
        // WARNING: If your matrix is only the lower half of an SPD matrix
        // then the result will only be a copy of the lower half of the 
        // spliced rows. See note above (Note on splice methods)
        bool SpliceRows(
            const std::vector<size_t> row_indices, 
            SparseSolver & res) const;

        // Wrapper than calls 
        //   SpliceColumns(column_indices, temp) 
        //   then 
        //   temp.SpliceRows(row_indices, res);
        bool SpliceColumnsThenRows(
            const std::vector<size_t> column_indices, 
            const std::vector<size_t> row_indices, 
            SparseSolver & res) const;

        // Vertical concatenation, works like matlab:
        // C = [A ; B];
        // If this matrix is A in the above the result is the concatenation of
        // this matrix on top of the supplied other matrix into the result
        // number of columns in A must equal number of columns in B
        //
        bool VerticalConcatenation(
            SparseSolver & other,
            SparseSolver & res);
        // Convert the matrix to IJV aka Coordinate format,
        // Outputs:
        //   vals             non-zero values, in no particular order
        //   row_indices      row indices corresponding to vals
        //   column_indices   column indices corresponding to vals
        //
        // All indices are 0-based
        void ConvertToIJV(
          std::vector<int> & row_indices,
          std::vector<int> & column_indices,
          std::vector<taucsType> & vals);

        // Convert lower triangle (including diagonal) of the matrix to IJV aka
        // Coordinate format,
        // Outputs:
        //   vals             non-zero values, in no particular order
        //   row_indices      row indices corresponding to vals
        //   column_indices   column indices corresponding to vals
        //
        // All indices are 0-based
        void ConvertLowerTriangleToIJV(
          std::vector<int> & row_indices,
          std::vector<int> & column_indices,
          std::vector<taucsType> & vals);

        // Convert the matrix to Compressed sparse column (CSC or CCS) format,
        // also known as Harwell Boeing format. As described:
        // http://netlib.org/linalg/html_templates/node92.html
        // or
        // http://en.wikipedia.org/wiki/Sparse_matrix
        //   #Compressed_sparse_column_.28CSC_or_CCS.29
        // Outputs:
        //   nnz              number of non zero entries
        //   num_rows         number of rows
        //   num_cols         number of cols
        //   vals             non-zero values, row indices running fastest,
        //                    size(vals) = nnz 
        //   row_indices      row indices corresponding to vals, 
        //                    size(row_indices) = nnz
        //   column_pointers  index in vals of first entry in each column, 
        //                    size(column_pointers) = num_cols+1
        //
        // All indices and pointers are 0-based
        void ConvertToHarwellBoeing(
          int & nnz,
          int & num_rows,
          int & num_cols,
          std::vector<taucsType> & vals,
          std::vector<int> & row_indices,
          std::vector<int> & column_pointers);

#ifdef PRINTMODE
        // Print the matrix in IJV form
        void PrintIJV();
#endif

protected:
#ifndef NO_TAUCS
	// Will create ATA matrix and AT matrix (in ccs format).
	void CreateATA();
	// Will create A matrix (in ccs format)
	void CreateA();

	void ClearFactorATA();
	void ClearFactorA();
	void ClearMatricesATA(); // clears m_A, m_ATA, m_AT
	void ClearMatricesA(); // clears the m_A matrix
#endif
};


#ifdef PRINTMODE
#include <cstdio>
inline void SparseSolver::PrintIJV(){
  std::map<int,taucsType>::iterator it;
  for (int c = 0; c < m_numCols; ++c) 
  {
    it = m_colsA[c].begin();
    while (it != m_colsA[c].end())
    {
      printf("%d %d %g\n",it->first,c, it->second);
      ++it;
    }
  }
}
#endif


inline void SparseSolver::ClearObject() {
#ifndef NO_TAUCS
	ClearFactorATA();
	ClearFactorA();
	ClearMatricesATA();
	ClearMatricesA();

	m_ATb.clear();
	m_etree = NULL;
#endif
	m_numCols = 0;
	m_numRows = 0;
	m_colsA.clear();
}

// adds a new entry to the sparse matrix A
// i and j are 0-based
// if the A^T A matrix had been factored, this will destroy the factor
inline void SparseSolver::AddIJV(const int i, const int j, const taucsType v) {
#ifndef NO_TAUCS
	ClearFactorATA();
	ClearFactorA();
	ClearMatricesATA();
	ClearMatricesA();
#endif

	m_colsA[j][i] = v;

	// update number of rows
	if ((i+1) > m_numRows)
		m_numRows = i+1;
}

// updates the entry in A by adding the value v to the current value
// i and j are 0-based
// if the A^T A matrix had been factored, this will destroy the factor
inline void SparseSolver::AddToIJV(const int i, const int j, const taucsType v) {
#ifndef NO_TAUCS
	ClearFactorATA();
	ClearFactorA();
	ClearMatricesATA();
	ClearMatricesA();
#endif

	// if no such entry existed
	if (m_colsA[j].find(i) == m_colsA[j].end())
		m_colsA[j][i] = 0;

	m_colsA[j][i] += v;

	// update number of rows
	if ((i+1) > m_numRows)
		m_numRows = i+1;
}

#ifndef NO_TAUCS
inline taucs_ccs_matrix * SparseSolver::GetATA_Copy() {
	if (! m_ATA) {
		CreateATA();
	}

	return MatrixCopy(m_ATA);
}
#endif

// Implementation
#define PRINTMODE

#ifndef NO_TAUCS
extern "C" {
    void chol_update(void *vF, int index, double val,void**vetree);
}
#endif

void SparseSolver::Reset(const int numRows, const int numCols)
{
	ClearObject();
	m_numRows = numRows;
	m_numCols = numCols;
	m_colsA.resize(m_numCols);
}


#ifndef NO_TAUCS
// Will create A^T A and its factorization, as well as store A^T
// returns true on success, false otherwise
bool SparseSolver::FactorATA() {
	if (m_SPD)
		return FactorA();

	ClearFactorATA();
	if (m_ATA == NULL)
		CreateATA();

	// factorization
	int	rc = taucs_linsolve(m_ATA, &m_factorATA,0,NULL,NULL,SIVANfactor,SIVANopt_arg);
	if (rc != TAUCS_SUCCESS) 
		return false;

	return true;
}


// Will create ATA matrix and AT matrix.
// returns true on success, false otherwise
void SparseSolver::CreateATA() {
	if (m_A == NULL)
		CreateA();

	ClearMatricesATA();

	if (m_SPD) {
		// don't need to multiply because m_A is invertible and symmetric
		return;
	}
	else {
		m_AT = MatrixTranspose(m_A);
		// Multiply to get ATA:
		m_ATA = Mul2NonSymmMatSymmResult(m_AT, m_A);
	}
}

// Will create A matrix (in ccs format)
void SparseSolver::CreateA() {
	ClearMatricesA();

	if (m_SPD) 
		m_A = CreateTaucsMatrixFromColumns(m_colsA, m_numRows, TAUCS_DOUBLE|TAUCS_SYMMETRIC|TAUCS_LOWER);
	else
		m_A = CreateTaucsMatrixFromColumns(m_colsA, m_numRows, TAUCS_DOUBLE);
}


// Will create *symbolic* factorization of ATA and return it; returns NULL on failure
void * SparseSolver::SymbolicFactorATA() {
	if (m_SPD) {
		if (m_A == NULL)
			CreateA();

		return taucs_ccs_factor_llt_symbolic(m_A);
	}
	else {
		if (m_ATA == NULL)
			CreateATA();

		return taucs_ccs_factor_llt_symbolic(m_ATA);
	}
}
// Will create the actual numerical factorization of the present ATA matrix
// with the same non-zero pattern as created by SymbolicFactorATA()
bool SparseSolver::FactorATA_UseSymbolic(void * symbFactor,
										 taucs_ccs_matrix ** PAP,
										 int * perm,
										 int * invperm) 
{	
	if (m_SPD) {
		if (m_A == NULL)
			CreateA();

		*PAP = taucs_ccs_permute_symmetrically(m_A, perm, invperm);
	}
	else {
		if (m_ATA == NULL)
			CreateATA();

		*PAP = taucs_ccs_permute_symmetrically(m_ATA, perm, invperm);
	}

	// compute numerical factorization
	if (taucs_ccs_factor_llt_numeric(*PAP, symbFactor) != 0)	{
		return false;
	}
	else {
		return true;	
	}
}

// Will create A and its factorization. If A is spd Cholesky factor, otherwise LU
// returns true on success, false otherwise
// Assumes A is square and invertible
bool SparseSolver::FactorA() {
	ClearFactorA();

	if (m_A == NULL)
		CreateA();

	// factorization
	int	rc;
	if (m_SPD)
		rc = taucs_linsolve(m_A, &m_factorA,0,NULL,NULL,SIVANfactor,SIVANopt_arg);
	else
		rc = taucs_linsolve(m_A, &m_factorA,0,NULL,NULL,SIVANfactorLU,SIVANopt_arg);

	if (rc != TAUCS_SUCCESS) 
		return false;

	return true;
}

// Solves the normal equations for Ax = b, namely A^T A x = A^T b
// it is assumed that enough space is allocated in x
// Uses the factorization provided in factor
// numRhs should be the the number of right-hand sides stored in b (and respectively,
// there should be enough space allocated in the solution vector x)
// returns true on success, false otherwise
bool SparseSolver::SolveATA_UseFactor(void * factor, 
									  taucs_ccs_matrix * PAP,
									  int * perm,
									  int * invperm,
									  const taucsType * b, taucsType * x, const int numRhs) {
	if (factor == NULL) {
		return false;
	}

	if (! m_SPD) {
		// prepare right-hand side
		if ((int)m_ATb.size() != m_numCols*numRhs)
			m_ATb.resize(m_numCols*numRhs);

		// multiply right-hand side
		for (int i = 0; i < numRhs; ++i) {
			MulMatrixVector(m_AT, b + i*m_numRows, (taucsType *)&(m_ATb.front()) + i*m_numCols);
		}
	}

	//	 solve the system:
	std::vector<taucsType>  PB(m_numCols*numRhs), PX(m_numCols*numRhs);

	// permute rhs
	for (int c = 0; c < numRhs; ++c) {
		taucsType * curPB = &PB[0] + c*m_numCols;
		const taucsType * curB = (m_SPD)? (b + c*m_numCols) : (&(m_ATb.front()) + c*m_numCols);

		for (int i=0; i<m_numCols; ++i)
			curPB[i] = curB[perm[i]];
	}

	// solve by back-substitution
	int rc;
	for (int c = 0; c < numRhs; ++c) {
		rc = taucs_supernodal_solve_llt(factor, &PX[0] + c*m_numCols, &PB[0] + c*m_numCols); 
		if (rc != TAUCS_SUCCESS) {
			return false;
		}
	}

	// re-permute x
	for (int c = 0; c < numRhs; ++c) {
		taucsType * curX = x + c*m_numCols;
		taucsType * curPX = &PX[0] + c*m_numCols;

		for (int i=0; i<m_numCols; ++i)
			curX[i] = curPX[invperm[i]];
	}

	return true;
}


// Solves the normal equations for Ax = b, namely A^T A x = A^T b
// it is assumed that enough space is allocated in x
// If a factorization of A^T A exists, it will be used, otherwise it will be created
// numRhs should be the the number of right-hand sides stored in b (and respectively,
// there should be enough space allocated in the solution vector x)
// returns true on success, false otherwise
bool SparseSolver::SolveATA(const taucsType * b, taucsType * x, const int numRhs) {
	if (m_SPD) {
		return SolveA(b, x, numRhs);
	}
	else {
		if (m_factorATA == NULL) {
			bool rc = FactorATA();
			if (!rc)
				return false;
		}
	}

	// prepare right-hand side
	if ((int)m_ATb.size() != m_numCols*numRhs)
		m_ATb.resize(m_numCols*numRhs);

	// multiply right-hand side
	for (int i = 0; i < numRhs; ++i) {
		MulMatrixVector(m_AT, b + i*m_numRows, (taucsType *)&(m_ATb.front()) + i*m_numCols);
	}

	int rc;
	// solve the system
	rc = taucs_linsolve(m_ATA,
						&m_factorATA,
						numRhs,
						x,
						(taucsType *)&(m_ATb.front()),
						SIVANsolve,
						SIVANopt_arg);
	return (rc == TAUCS_SUCCESS); 
}

// Same as SolveATA only that this solves the actual system Ax = b
// It is assumed that A is rectangular and invertible
bool SparseSolver::SolveA(const taucsType * b, taucsType * x, const int numRhs) {
	if (m_factorA == NULL) 
		if (!FactorA())
			return false;

	// solve the system
	int	rc = taucs_linsolve(m_A,
							&m_factorA,
							numRhs,
							x,
							(void *)b,
							SIVANsolve,
							SIVANopt_arg);

	return (rc == TAUCS_SUCCESS);
}

// The version of SolveATA with 3 right-hand sides
// returns true on success, false otherwise
bool SparseSolver::SolveATA3(const taucsType * bx, const taucsType * by, const taucsType * bz,
							       taucsType * x,        taucsType * y,        taucsType * z) 
{
	if (m_factorATA == NULL) {
		bool rc = FactorATA();
		if (!rc)
			return false;
	}

	// prepare right-hand side
	if ((int)m_ATb.size() != m_numCols*3)
		m_ATb.resize(m_numCols*3);

	// multiply right-hand side
	MulMatrixVector(m_AT, bx, (taucsType *)&(m_ATb.front()));
	MulMatrixVector(m_AT, by, (taucsType *)&(m_ATb.front()) + m_numCols);
	MulMatrixVector(m_AT, bz, (taucsType *)&(m_ATb.front()) + 2*m_numCols);
	

	// solve the system
	int	rc = taucs_linsolve(m_ATA,
							&m_factorATA,
							1,
							x,
							(taucsType *)&(m_ATb.front()),
							SIVANsolve,
							SIVANopt_arg);

	if (rc != TAUCS_SUCCESS)
		return false;

	rc = taucs_linsolve(m_ATA,
							&m_factorATA,
							1,
							y,
							(taucsType *)&(m_ATb.front()) + m_numCols,
							SIVANsolve,
							SIVANopt_arg);

	if (rc != TAUCS_SUCCESS)
		return false;

	rc = taucs_linsolve(m_ATA,
							&m_factorATA,
							1,
							z,
							(taucsType *)&(m_ATb.front()) +  2*m_numCols,
							SIVANsolve,
							SIVANopt_arg);

	return (rc == TAUCS_SUCCESS);
}

// The version of SolveATA with 2 right-hand sides
// returns true on success, false otherwise
bool SparseSolver::SolveATA2(const taucsType * bx, const taucsType * by,
								   taucsType * x,        taucsType * y)
{
	if (m_factorATA == NULL) {
		bool rc = FactorATA();
		if (!rc)
			return false;
	}

	// prepare right-hand side
	if ((int)m_ATb.size() != m_numCols*2)
		m_ATb.resize(m_numCols*2);

	// multiply right-hand side
	MulMatrixVector(m_AT, bx, (taucsType *)&(m_ATb.front()));
	MulMatrixVector(m_AT, by, (taucsType *)&(m_ATb.front()) + m_numCols);
	

	// solve the system
	int	rc = taucs_linsolve(m_ATA,
							&m_factorATA,
							1,
							x,
							(taucsType *)&(m_ATb.front()),
							SIVANsolve,
							SIVANopt_arg);

	if (rc != TAUCS_SUCCESS)
		return false;

	rc = taucs_linsolve(m_ATA,
							&m_factorATA,
							1,
							y,
							(taucsType *)&(m_ATb.front()) + m_numCols,
							SIVANsolve,
							SIVANopt_arg);

	return (rc == TAUCS_SUCCESS);
}

void SparseSolver::ClearFactorATA() {
	// release the factor
	taucs_linsolve(NULL,&m_factorATA,0, NULL,NULL,SIVANfactor,SIVANopt_arg);
	m_factorATA = NULL;
	m_etree = NULL;
}

void SparseSolver::ClearFactorA() {
	// release the factor
	if (m_SPD)
		taucs_linsolve(NULL,&m_factorA,0, NULL,NULL,SIVANfactor,SIVANopt_arg);
	else
		taucs_linsolve(NULL,&m_factorA,0, NULL,NULL,SIVANfactorLU,SIVANopt_arg);
	m_factorA = NULL;

}

void SparseSolver::ClearMatricesATA() {
	// free the matrices
	taucs_free(m_AT);
	m_AT = NULL;
	taucs_free(m_ATA);
	m_ATA = NULL;
}

void SparseSolver::ClearMatricesA() {
	// free the matrices
	taucs_free(m_A);
	m_A = NULL;
}
#endif

void SparseSolver::AddMatrix(SparseSolver & B)
{
  AddMatrix(0,0,B);
}

void SparseSolver::AddMatrix(const int i, const int j, SparseSolver & B)
{
  std::map<int,taucsType>::iterator it;
  for (int c = 0; c < B.m_numCols; ++c) 
  {
    it = B.m_colsA[c].begin();
    while (it != B.m_colsA[c].end())
    {
      // because AddToIJV is not O(1), this algorithm is O(Bnnz * Amaxrows),
      // ideally it would be O(Bnnz)
      AddToIJV(i+it->first,j+c,it->second);
      ++it;
    }
  }
}

bool SparseSolver::IsSymmetric()
{
  if(m_SPD)
  {
    return true;
  }

  if(GetNumColumns() != GetNumRows())
  {
    return false;
  }

  // Cheapskate way to test if symmetric
  SparseSolver AT;
  if(!Transpose(AT))
  {
    return false;
  }
  // loop over columns
  std::map<int,taucsType>::iterator it;
  std::map<int,taucsType>::iterator ATit;
  for (int c = 0; c < m_numCols; ++c) 
  {
    it = m_colsA[c].begin();
    ATit = AT.m_colsA[c].begin();
    while (it != m_colsA[c].end() && ATit !=AT.m_colsA[c].end())
    {
      if( it->first != ATit->first || it->second != ATit->second)
      {
        return false;
      }
      ++it;
      ++ATit;
    }
  }
  return true;

}

#ifndef NO_TAUCS
// allows to add an anchored vertex without destroying the factor, if there was one
// i is the anchor's number (i.e. the index of the mesh vertex that is anchored is i)
// w is the weight of the anchor in the original Ax=b system. HAS TO BE POSITIVE!!
// if there was no factor for A^T A, this will simply add an anchor row to the Ax=b system
// if there was a factor, the factor will be updated, approximately at the cost of
// a single solve.
void SparseSolver::AddAnchor(const int i, const taucsType w) {
	if (m_SPD) {
		// update the columns
		if (m_colsA[i].find(i) == m_colsA[i].end())
			m_colsA[i][i] = w*w;
		else
			m_colsA[i][i] += w*w;

		if (m_factorA != NULL) {
			// the matrix is SPD so we update the factor of A
			chol_update(m_factorA, i, w, &m_etree);
			// update the matrix
			m_A->taucs_values[ m_A->colptr[i] ] += w*w;
		}
	}
	else {
		m_colsA[i][m_numRows] = w;
		m_numRows++;

		if (m_factorATA != NULL) {
			// update the ATA factor
			chol_update(m_factorATA, i, w, &m_etree);
			// update ATA
			m_ATA->taucs_values[ m_ATA->colptr[i] ] += w*w;
			// update A and AT (for multiplying rhs)
			CreateA();
			taucs_free(m_AT);
            m_AT = MatrixTranspose(m_A);
		}
	}	 
}
#endif

void SparseSolver::MultiplyMatrixVector(const taucsType * v, taucsType * result, const int numCols) const {
	// make result all zero
	memset(result, 0, m_numRows * numCols * sizeof(taucsType));

	std::map<int,taucsType>::const_iterator iter;
	int  offset_v;
	int  offset_r;

	if (m_SPD) { // the matrix is symmetric
		for (int c = 0; c < numCols; ++c) {
			offset_v = m_numCols*c;
			for (int col = 0; col < m_numCols; ++col) {
				// going over column col of the matrix, multiplying
				// it by v[col] and setting the appropriate values
				// of vector result; also mirroring the other triangle
				for (iter = m_colsA[col].begin(); iter->first < col; ++iter) {
					result[iter->first + offset_v]	+= v[col + offset_v]*(iter->second);
					// mirroring
					result[col + offset_v]			+= v[iter->first + offset_v]*(iter->second);
				}
				if (iter->first == col)
					result[iter->first + offset_v]	+= v[col + offset_v]*(iter->second);
			}
		}
	}
	else {
		for (int c = 0; c < numCols; ++c) {
			offset_v = m_numCols*c;
			offset_r = m_numRows*c;
			for (int col = 0; col < m_numCols; ++col) {
				// going over column col of the matrix, multiplying
				// it by v[col] and setting the appropriate values
				// of vector result
				for (iter = m_colsA[col].begin(); iter != m_colsA[col].end(); ++iter) {
					result[iter->first + offset_r] += v[col + offset_v]*(iter->second);
				}
			}
		}
	}
}

// Multiplies the matrix by diagonal matrix D from the left, stores the result
// in the columns of res
void SparseSolver::MultiplyDiagMatrixMatrix(const taucsType * D, SparseSolver & res) const {
	res = SparseSolver(m_numRows, m_numCols, false);

	std::map<int,taucsType>::const_iterator iterA;
	std::map<int,taucsType>::iterator		iterRes;

	res.m_colsA = m_colsA;
	if (m_SPD) {
		for (int col = 0; col < m_numCols; ++col) {
			iterA = m_colsA[col].begin();
			iterRes = res.m_colsA[col].begin();
			for ( ; iterA != m_colsA[col].end() && iterA->first <= col; ++iterA, ++iterRes) {
				iterRes->second = (iterA->second) * D[iterA->first];
				// mirroring
				res.m_colsA[iterA->first][col] = (iterA->second) * D[col];
			}
		}
	}
	else {
		for (int col = 0; col < m_numCols; ++col) {
			iterA	= m_colsA[col].begin();
			iterRes = res.m_colsA[col].begin();
			for ( ; iterA != m_colsA[col].end(); ++iterA, ++iterRes) {
				iterRes->second = (iterA->second) * D[iterA->first];
			}
		}
	}
}

// Multiplies the matrix by the given matrix B from the right
// and stores the result in the columns of res
// isSPD tells us whether the result should be SPD or not
// retuns false if the dimensions didn't match (in such case res is unaltered), otherwise true
bool SparseSolver::MultiplyMatrixRight(const SparseSolver & B, SparseSolver & res,const bool isSPD) const {
	// for now ignore the option that one of the involved matrices is symmetric and doesn't have
	// stored under-diagonal values in m_colsA
	// we assume both matrices have all the values in m_colsA

	// (m x n) * (n x k)
	int m = m_numRows;
	int n = m_numCols;
	int k = B.m_numCols;

	if (B.m_numRows != n)
		return false;

	std::map<int,taucsType>::const_iterator iterBi, iterA;
	std::map<int,taucsType>::iterator iterRes;

	int			rowInd;
	taucsType	BiVal;

	res = SparseSolver(m, k, isSPD);
	for (int i = 0; i < k; ++i) { // go over the columns of res; res(i) = A*B(i)
		for (iterBi = B.m_colsA[i].begin(); iterBi != B.m_colsA[i].end(); ++iterBi) { // go over column B(i)
			// iterBi->second multiplies column A(iterBi->first)
			BiVal = iterBi->second;
			rowInd = iterBi->first;

			for (iterA = m_colsA[rowInd].begin(); iterA != m_colsA[rowInd].end(); ++iterA) { // go over column of A
				iterRes = res.m_colsA[i].find(iterA->first);
				if (iterRes == res.m_colsA[i].end()) // first occurence of this row number
					res.m_colsA[i][iterA->first] = (iterA->second)*BiVal;
				else
					iterRes->second += (iterA->second)*BiVal;
			}
		}
	}

	return true;
}

// Transposes the matrix 
// and stores the result in the columns of res
bool SparseSolver::Transpose(SparseSolver & res) const {
	res = SparseSolver(m_numCols, m_numRows, m_SPD);
	// if the matrix was symmetric to begin with, just copy the columns
	if (m_SPD) {
		res.m_colsA = m_colsA;
		return true;
	}
	// else need to transpose
	std::map<int,taucsType>::const_iterator iterA;
	for (int i=0; i < m_numCols; ++i) {
		for (iterA = m_colsA[i].begin(); iterA != m_colsA[i].end(); ++iterA) {
			res.m_colsA[iterA->first][i] = iterA->second;
		}
	}
	return true;
}

// will remove the rows startInd-endInd (including ends, zero-based) from the matrix
// will discard any factors	
bool SparseSolver::DeleteRowRange(const int startInd, const int endInd) {
	if (startInd < 0 || endInd >= m_numRows || startInd > endInd)
		return false;

#ifndef NO_TAUCS
	ClearFactorATA();
	ClearFactorA();
	ClearMatricesATA();
	ClearMatricesA();
#endif

	std::map<int,taucsType>::iterator it;
	std::map<int,taucsType>::iterator tempIt;
	int diff = endInd - startInd + 1;

	for (int c = 0; c < m_numCols; ++c) {
		it = m_colsA[c].begin();
		while (it != m_colsA[c].end()) {
			if (it->first >= startInd && it->first <= endInd) { // row that should be erased
				m_colsA[c].erase(it++);
                        } else {
				if (it->first > endInd) {// row number must be updated
					m_colsA[c][it->first - diff] = it->second;
					m_colsA[c].erase(it++);
				}
				else
                    ++it;
			}
		}
	}

	m_numRows -= diff;
	m_SPD = false;
	
	return true;
}

bool SparseSolver::MultiplyRows(const int startInd, const int endInd, const taucsType val)
{
    if (startInd < 0 || endInd >= m_numRows || startInd > endInd) return false;

#ifndef NO_TAUCS
    ClearFactorATA();
    ClearFactorA();
    ClearMatricesATA();
    ClearMatricesA();
#endif

    std::map<int,taucsType>::iterator it;

    for (int c = 0; c < m_numCols; ++c)
    {
        it = m_colsA[c].begin();
        while (it != m_colsA[c].end())
        {
            if (it->first >= startInd && it->first <= endInd) // row that should be multiplied
			{
                it->second *= val;
			}
            it++;
        }
    }

    return true;
}


// will remove the columns startInd-endInd (including ends, zero-based) from the matrix
// will discard any factors	
bool SparseSolver::DeleteColumnRange(const int startInd, const int endInd) {
	if (startInd < 0 || endInd >= m_numCols || startInd > endInd)
		return false;

#ifndef NO_TAUCS
	ClearFactorATA();
	ClearFactorA();
	ClearMatricesATA();
	ClearMatricesA();
#endif

	m_colsA.erase(m_colsA.begin() + startInd, m_colsA.begin() + (endInd+1));

	m_numCols -= (endInd - startInd + 1);

	return true;
}

// copies columns startInd through endInd (including, zero-based) from the matrix to res matrix
// previous data in res, including factors, is erased.
bool SparseSolver::CopyColumnRange(const int startInd, const int endInd, SparseSolver & res) const {
	if (startInd < 0 || endInd >= m_numCols || startInd > endInd)
		return false;

	int k = endInd - startInd + 1;

	res = SparseSolver(m_numRows, k, false);
	for (int c = startInd; c <= endInd; ++c)
		res.m_colsA[c - startInd] = m_colsA[c];

	return true;
}

// copies rows startInd through endInd (including, zero-based) from the matrix to res matrix
// previous data in res, including factors, is erased.
bool SparseSolver::CopyRowRange(const int startInd, const int endInd, SparseSolver & res) const {
	if (startInd < 0 || endInd >= m_numRows || startInd > endInd)
		return false;

	res.m_colsA = m_colsA;
	for (int c = 0; c < m_numCols; ++c) {
		std::map<int,taucsType>::iterator it = res.m_colsA[c].begin();
		std::map<int,taucsType>::iterator tempIt;
		for ( ; it != res.m_colsA[c].end(); ++it) {
			if (it->first < startInd || it->first > endInd) {
				tempIt = it;
				it++;
				res.m_colsA[c].erase(tempIt);
				continue;
			}
		}
	}

	return true;
}


void SparseSolver::MultiplyMatrixScalar(const taucsType val)
{
#ifndef NO_TAUCS
    ClearFactorATA();
    ClearFactorA();
    ClearMatricesATA();
    ClearMatricesA();
#endif

	//Over all entries
    for (int c=0;c<m_numCols;c++)
    {
		for(std::map< int, taucsType >::iterator it=m_colsA[c].begin();it!=m_colsA[c].end();it++)
        {
			it->second *= val; //Multiply
        }
    }
}

void SparseSolver::IdentityMinusMatrix()
{
#ifndef NO_TAUCS
    ClearFactorATA();
    ClearFactorA();
    ClearMatricesATA();
    ClearMatricesA();
#endif

	//Over the diagonal
    for (int c=0;c<m_numCols;c++)
    {
		std::map< int, taucsType >::iterator it = m_colsA[c].find(c);
		if (it != m_colsA[c].end())
		{
			it->second = 1.0 - it->second;
		}
		else
		{
			AddIJV(c, c, 1.0);
		}
    }
}


void SparseSolver::IdentityPlusMatrix()
{
#ifndef NO_TAUCS
    ClearFactorATA();
    ClearFactorA();
    ClearMatricesATA();
    ClearMatricesA();
#endif

	//Over the diagonal
    for (int c=0;c<m_numCols;c++)
    {
		std::map< int, taucsType >::iterator it = m_colsA[c].find(c);
		if (it != m_colsA[c].end())
		{
			it->second += 1.0;
		}
		else
		{
			AddIJV(c, c, 1.0);
		}
    }
}

//
// Memory: 1*spliced = 1*O(column_indices.size * numRows)
// Runtime: O(size of spliced) = O(column_indices.size * numRows)
bool SparseSolver::SpliceColumns(
    const std::vector<size_t> column_indices,
    SparseSolver & res) const
{
  // Check that all indices in column_indices are valid
  for(
      std::vector<size_t>::const_iterator old_index = column_indices.begin();
      old_index != column_indices.end();
      old_index++)
  {
    if((*old_index) >= (size_t)m_numCols)
    {
      return false;
    }
  }

  // create result
  res = SparseSolver(m_numRows,column_indices.size(),false);

  // splice (and reorder) columns
  size_t new_index = 0;
  for(
      std::vector<size_t>::const_iterator old_index = column_indices.begin();
      old_index != column_indices.end();
      old_index++,new_index++)
  {
    res.m_colsA[new_index] = m_colsA[(*old_index)];
  }

  return true;
}

// Cheating way of handling SpliceRows. 
// SpliceRows = Transpose -> SpliceColumns -> Transpose
// 
// Sparse (not dense), but not perfect implementation:
// Memory: 1*original + 2*spliced
// Runtime: O(
//   transpose of original + 
//   splicecolumns of original + 
//   transpose of spliced)
//   = O(
//   size of original + 
//   size of spliced
//   size of spliced)
//   = O( size of original )
//
// A good implementation should be:
// Memory: 1*spliced = 1*O(row_indices.size * numCols)
// Runtime: O(size of spliced) = O(row_indices.size * numCols)
bool SparseSolver::SpliceRows(
  const std::vector<size_t> row_indices,
  SparseSolver & res) const 
{
  // Check that all indices in row_indices are valid
  for(
      std::vector<size_t>::const_iterator old_index = row_indices.begin();
      old_index != row_indices.end();
      old_index++)
  {
    if((*old_index) >= (size_t)m_numRows)
    {
      return false;
    }
  }

  // transpose (now columns are rows, rows are columns)
  SparseSolver transpose;
  if(!Transpose(transpose))
  {
    return false;
  }

  // Splice columns of transpose
  SparseSolver spliced_transpose;
  if(!transpose.SpliceColumns(row_indices, spliced_transpose))
  {
    return false;
  }

  // Transpose again (rows back to rows, columns back to columns)
  if(!spliced_transpose.Transpose(res))
  {
    return false;
  }

  return true;
}

bool SparseSolver::SpliceColumnsThenRows(
  const std::vector<size_t> column_indices, 
  const std::vector<size_t> row_indices, 
  SparseSolver & res) const 
{
  SparseSolver temp;
  if(!SpliceColumns(column_indices,temp))
  {
    return false;
  }
  if(!temp.SpliceRows(row_indices,res))
  {
    return false;
  }
  return true;
}
    

bool SparseSolver::VerticalConcatenation(
    SparseSolver & other,
    SparseSolver & res)
{
  // if other if skinnier than we'll padd with zeros
  if(other.GetNumColumns() > m_numCols)
  {
    return false;
  }
  res = SparseSolver(m_numRows+other.GetNumRows(), m_numCols,false);
  // copy all columns from this matrix to result
  res.m_colsA = m_colsA;
  // copy all entries in other to result with rows offset by "m_numRows"
  std::map<int,taucsType>::iterator it;
  for (int c = 0; c < other.m_numCols; ++c)
  {
    it = other.m_colsA[c].begin();
    while (it != other.m_colsA[c].end())
    {
      res.AddIJV(it->first+m_numRows,c,it->second);
      ++it;
    }
  }

  return true;
}

void SparseSolver::ConvertToIJV(
  std::vector<int> & row_indices,
  std::vector<int> & column_indices,
  std::vector<taucsType> & vals)
{
  // Get number of non-zeros
  size_t nnz = 0;
  // loop over columns
  for (int c = 0; c < m_numCols; ++c) 
  {
    // add up the number of non-zeros in each column
    nnz += m_colsA[c].size();
  }
  // allocate space in outputs
  vals.resize(nnz);
  row_indices.resize(nnz);
  column_indices.resize(nnz);

  size_t index = 0;
  
  std::map<int,taucsType>::iterator it;
  // loop over columns
  for (int c = 0; c < m_numCols; ++c) 
  {
    // loop over rows
    it = m_colsA[c].begin();
    while (it != m_colsA[c].end())
    {
      row_indices[index] = it->first;
      column_indices[index] = c;
      vals[index] = it->second;
      // increment index
      index++;
      
      ++it;
    }
  }
}


void SparseSolver::ConvertLowerTriangleToIJV(
  std::vector<int> & row_indices,
  std::vector<int> & column_indices,
  std::vector<taucsType> & vals)
{
  // Get number of ALL non-zeros
  size_t nnz = 0;
  // loop over columns
  for (int c = 0; c < m_numCols; ++c) 
  {
    // add up the number of non-zeros in each column
    nnz += m_colsA[c].size();
  }
  // allocate (too much) space in outputs
  vals.resize(nnz);
  row_indices.resize(nnz);
  column_indices.resize(nnz);

  size_t index = 0;
  
  std::map<int,taucsType>::iterator it;
  // loop over columns
  for (int c = 0; c < m_numCols; ++c) 
  {
    // loop over rows
    it = m_colsA[c].begin();
    while (it != m_colsA[c].end())
    {
      // if row is less than or equal to column then we're in the lower
      // triangle
      if(it->first>= c)
      {
        row_indices[index] = it->first;
        column_indices[index] = c;
        vals[index] = it->second;
        // increment index
        index++;
      }
      
      ++it;
    }
  }
  // resize to actually number of non-zeros in lower triangle
  vals.resize(index);
  row_indices.resize(index);
  column_indices.resize(index);

}

void SparseSolver::ConvertToHarwellBoeing(
  int & nnz,
  int & num_rows,
  int & num_cols,
  std::vector<taucsType> & vals,
  std::vector<int> & row_indices,
  std::vector<int> & column_pointers)
{
  num_rows = m_numRows;
  num_cols = m_numCols;
  nnz = 0;
  // loop over columns
  for (int c = 0; c < m_numCols; ++c) 
  {
    // add up the number of non-zeros in each column
    nnz += m_colsA[c].size();
  }

  vals.resize(nnz);
  row_indices.resize(nnz);
  column_pointers.resize(num_cols+1);

  // loop over columns
  int column_pointer = 0;
  int val_index = 0;
  std::map<int,taucsType>::iterator it;
  int c;
  for (c = 0; c < m_numCols; ++c) 
  {
    column_pointers[c] = column_pointer;
    column_pointer += m_colsA[c].size();
    // over rows in this column
    it = m_colsA[c].begin();
    while (it != m_colsA[c].end())
    {
      vals[val_index] = it->second;
      row_indices[val_index] = it->first;
      val_index++;
      ++it;
    }
  }
  // by convention column_pointers[num_cols] = nnz
  column_pointers[c] = column_pointer;
}
