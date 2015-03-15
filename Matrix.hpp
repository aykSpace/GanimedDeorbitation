// Created by Yap Chun Wei
// Version 1.0
// Version 1.1
	// - Removed Identity matrix implementation because not well implemented.
// Version 1.2
	// - Add RowVector and ColVector classes.
	// - Fix errors involving reverse iterators.
	// - Changed function empty().
// Version 1.3
	// - Added memberspaces row and col to StoragePolicy and Expression Templates.
	// - Replace rowsize, colsize, rowbegin, rowend, colbegin, and colend with
	//   row::size, col::size, row(i).begin, row(i).end, col(i).begin and col(i).end.
	// - Added row::insert, row::erase, row::push_back, row::pop_back, 
	//   col::insert, col::erase, col::push_back, col::pop_back.
	// - Added pop_back for RowVector and ColVector.
// Version 1.4
	// - Change memberspaces row and col to Row and Col to prevent conflicts causing row(i)
	//	 not to work when RowVector or ColVector is declared.
// Version 1.5 (16 April 2006)
	// - Fix code so that it is no longer dependent on STLport and will work on Visual Studio 2005, gcc 3.3.2.
	// - Because of the fix, the constancy of the row and col iterators may not really work.
	// - Added row::push_front, row::pop_front, col::push_front, col::pop_front.
	// - Added push_front and pop_front for RowVector and ColVector.

// Variable Suffixes
	// _ - Class members
	// T - Typename (typename)
	// P - Policy (template<typename> class)
	// Type - New type by typedef

#ifndef MATRIX_HPP
#define MATRIX_HPP

namespace YMatrix // Forward declaration of namespace
{
// Forward declaration of class Matrix.
// This matrix is the only matrix that user should construct in their code.
template<typename ValueT,template<typename> class MatrixTypeP,template<typename,typename> class MathP>
class Matrix;
template<typename ValueT,template<typename> class MatrixTypeP,template<typename,typename> class MathP>
class RowVector;
template<typename ValueT,template<typename> class MatrixTypeP,template<typename,typename> class MathP>
class ColVector;

template<typename ValueT> class DenseMatrix;
template<typename ExprT> class MET;
template<typename ValueT, typename MatrixT>  class MathMatrix;
template<typename ValueT, typename MatrixT>  class MathMETMatrix;
template<typename ValueT, typename MatrixT>  class NonMathMatrix;

template<typename MatrixT> struct Transpose;
template<typename MatrixT> struct Diagonal;
template<typename MatrixT> struct Covariance;
template<typename MatrixT> struct Power;
template<typename MatrixT> struct Mean;
template<typename MatrixT> struct Median;
template<typename MatrixT> struct Sum;
template<typename MatrixT> struct CumulativeSum;
template<typename MatrixT> struct Inverse;
template<typename MatrixT> struct Det;
}

#include <algorithm>
#include <iostream>
#include <iterator>
#include <functional>
#include <cmath>
#include <numeric>
#include <vector>
using namespace std;

namespace YMatrix
{
// Forward declaration of class Matrix.
// This matrix is the only matrix that user should construct in their code.
template<typename ValueT,template<typename> class MatrixTypeP,template<typename,typename> class MathP>
class Matrix;
template<typename ValueT,template<typename> class MatrixTypeP,template<typename,typename> class MathP>
class RowVector;
template<typename ValueT,template<typename> class MatrixTypeP,template<typename,typename> class MathP>
class ColVector;

template<typename ValueT>
class DenseMatrix
// PURPOSE: Store elements in a dense matrix format.
// PURPOSE: Provide functions to access and modify elements.
// REQUIRE: Not to be constructed directly by user.
// PROMISE: Elements will be stored as a dense matrix and can be accessed and modified.
// NOTES: Three critical variables are exposed to the user: rows_, cols_ and matrix_.
// NOTES: This is to allow user to eliminate overheads in performance critical routines
// NOTES: by using the variables directly instead of accessing them through accessor functions.
// NOTES: rows_ and cols_ contains the number of rows and Columns of the matrix respectively.
// NOTES: matrix_ contains the entire matrix data.
// NOTES: The function ResizeIterators() is also purposely exposed to allow users who modify
// NOTES: rows_,cols_ and matrix_ directly to resize and update all internal iterators.
{	
protected: // Typedefs for internal use.
	typedef DenseMatrix<ValueT> Self;
	typedef ValueT*				VectorType;
	typedef ValueT*				VectorTypeIterator;
	typedef ValueT				MatrixStorageType;
	typedef VectorType			MatrixType;

public: // Memberspaces.
	class Row
	{
	private: // Typedefs for internal use.
		typedef Row	Self;
		typedef DenseMatrix<ValueT>	EnclosingClass;

	public:	// STL-like typedefs.
		typedef size_t								size_type;
		typedef VectorTypeIterator					iterator;
		typedef const iterator						const_iterator;
		typedef VectorTypeIterator					reverse_iterator;
		typedef const reverse_iterator				const_reverse_iterator;

	// This is a special case where there is no structure for row iterators because of the nature
	// of storage for the matrix. By using ValueT* as the type for row iterators, performance
	// of the matrix is higher than using a structure for row iterators. Thus row iterators
	// is used whenever possible instead of Column iterators to maintain performance.

	public: // Constructors, destructor, copy constructors and assignment.
		// Constructor requiring address of enclosing class i.e. address of DenseMatrix<ValueT>.
		Row(EnclosingClass& ec) : enclosingClass_(ec) {};
		// Destructor.
		~Row() {};

	public: // Functions.
		// Change current row (Do not change for all matrix types).
		Self&		operator()(size_type rowNo)			{currentIndex_ = rowNo; return *this;}
		const Self&	operator()(size_type rowNo) const	{currentIndex_ = rowNo; return *this;}

		// Return the total number of rows in matrix.
		size_type size() const {return enclosingClass_.rows_;}

		// Insert row at specified row
		template<typename ForwardIterator>
		void insert(size_type rowNo, const ForwardIterator& first)
		{
			enclosingClass_.resize(max(enclosingClass_.rows_, rowNo)+1, enclosingClass_.cols_);
			for (size_type i=enclosingClass_.rows_-1; i>rowNo; --i)
			{
				copy ((*this)(i-1).begin(), (*this)(i-1).end(), (*this)(i).begin());
			}
			copy (first, first+enclosingClass_.cols_, (*this)(rowNo).begin());
		}

		// Erase row at specified row
		void erase(size_type rowNo)
		{
			if (rowNo < enclosingClass_.rows_)
			{
				for (size_type i=rowNo+1; i<enclosingClass_.rows_; ++i)
				{
					copy ((*this)(i).begin(), (*this)(i).end(), (*this)(i-1).begin());
				}
				enclosingClass_.resize(enclosingClass_.rows_-1, enclosingClass_.cols_);
			}
		}

		// Add row at end of matrix
		template<typename ForwardIterator>
		void push_back(const ForwardIterator& first)
		{
			this->insert(enclosingClass_.rows_, first);
		}

		// Add row at beginning of matrix
		template<typename ForwardIterator>
		void push_front(const ForwardIterator& first)
		{
			this->insert(0, first);
		}

		// Remove row at end of matrix
		void pop_back()
		{
			this->erase(enclosingClass_.rows_-1);
		}

		// Remove row at beginning of matrix
		void pop_front()
		{
			this->erase(0);
		}

		// Iterators access (Do not change for all matrix types).
		// Must be called by row(i) and not row.
		iterator 		begin() 		{return enclosingClass_.rowIteratorStart_[currentIndex_];} 
		const_iterator 	begin() const	{return enclosingClass_.rowIteratorStart_[currentIndex_];}
		iterator 		end() 			{return enclosingClass_.rowIteratorFinish_[currentIndex_];}   
		const_iterator 	end() const 	{return enclosingClass_.rowIteratorFinish_[currentIndex_];}

		reverse_iterator		rbegin()		{return reverse_iterator(enclosingClass_.rowIteratorFinish_[currentIndex_]);}
		const_reverse_iterator	rbegin() const	{return const_reverse_iterator(enclosingClass_.rowIteratorFinish_[currentIndex_]);}
		reverse_iterator		rend()			{return reverse_iterator(enclosingClass_.rowIteratorStart_[currentIndex_]);}
		const_reverse_iterator	rend() const	{return const_reverse_iterator(enclosingClass_.rowIteratorStart_[currentIndex_]);}

	private: // Member variables (Do not change for all matrix types).
		EnclosingClass& enclosingClass_;
		mutable size_type currentIndex_;

	private: // Undefined constructors and copy constructors.
		Row();
		Row(const Self&);
	} row;
	friend class Row;

	class Col
	{
		// Column iterator.
		struct ColIterator
		{
		private: // Typedefs for internal use.
			typedef ColIterator			Self;
			typedef ColIterator			Nonconst_self;
			typedef const ColIterator	Const_self;

		public: // STL-like typedefs.
			typedef random_access_iterator_tag	iterator_category;
			typedef ValueT						value_type;
			typedef ValueT&						reference;
			typedef ValueT*						pointer;
			typedef size_t						size_type;
			typedef ptrdiff_t					difference_type;

		public:	// Constructors, destructor, copy constructors and assignment.
			// Default constructor
			ColIterator() {};
			// Constructor requiring vector type iterator and total number of Columns
			ColIterator(VectorTypeIterator, size_type sizeCol) : current_(x), cols_(sizeCol) {};
			// Destructor
			~ColIterator() {};
			// Copy constructor
			ColIterator(const Nonconst_self& x) : current_(x.current_), cols_(x.cols_) {};
			// Self = Vector type iterator
			Self& operator=(VectorTypeIterator x) {current_ = x; return *this;}

		public: // Functions (Do not change for all types of iterators).
			difference_type operator-(const Self& x) const {return this->Subtract(x);}
			Self& operator++() {this->Increment(); return *this;}
			Self operator++(int)  {Self tmp(*this); ++*this; return tmp;}
			Self& operator--() {this->Decrement(); return *this;}
			Self operator--(int) {Self tmp(*this); --*this; return tmp;}
			Self& operator+=(difference_type n) {this->Advance(n); return *this;}
			Self operator+(difference_type n) const {Self tmp(*this); return tmp += n;}
			Self& operator-=(difference_type n) {return *this += -n;}
			Self operator-(difference_type n) const {Self tmp(*this); return tmp -= n;}

		public:	// Functions
			void SetCols(size_type sizeCol) {cols_ = sizeCol;}
			// Dereference iterator.
			reference operator*() const {return *(this->current_);}
			reference operator->() const {return (this->current_);}
			reference operator[](difference_type n) const {return *(*this + n);}
			// Difference between two iterators.
			difference_type Subtract(const Self& x) const {return (current_ - x.current_)/cols_;} 
			// Iterator increment.
			void Increment() {current_+=cols_;}
			// Iterator decrement.
			void Decrement() {current_-=cols_;} 
			// Iterator advance.
			void Advance(difference_type n) {current_ += n*cols_;} 

		public:	// Friend functions.
			friend Self operator+(ptrdiff_t n, const Self& x) {return x + n;}
			friend bool operator==(const Self& x, const Self& y) {return (x.cols_==y.cols_ && x.current_==y.current_);}
			friend bool operator<(const Self& x, const Self& y) {return (x.cols_==y.cols_ && x.current_ < y.current_);}
			friend bool operator!=(const Self& x, const Self& y) {return !(x==y);}
			friend bool operator>(const Self& x, const Self& y) {return y < x;}
			friend bool operator>=(const Self& x, const Self& y) {return !(x < y);}
			friend bool operator<=(const Self& x, const Self& y) {return !(y < x);}
			
		public:	// Member variables.
			size_type cols_;
			VectorTypeIterator current_;
		};
	
	private: // Typedefs for internal use.
		typedef Col	Self;
		typedef DenseMatrix<ValueT>	EnclosingClass;

	public:	// STL-like typedefs.
		typedef size_t				size_type;
		typedef ColIterator			iterator;
		typedef const ColIterator	const_iterator;
		typedef ColIterator			reverse_iterator;
		typedef const ColIterator	const_reverse_iterator;

	public: // Constructors, destructor, copy constructors and assignment.
		// Constructor requiring address of enclosing class i.e. address of DenseMatrix<ValueT>.
		Col(EnclosingClass& ec) : enclosingClass_(ec) {};
		// Destructor.
		~Col() {};

	public: // Functions.
		// Change current row (Do not change for all matrix types).
		Self&		operator()(size_type colNo)			{currentIndex_ = colNo; return *this;}
		const Self&	operator()(size_type colNo) const	{currentIndex_ = colNo; return *this;}

		// Return the total number of columns in matrix.
		size_type size() const {return enclosingClass_.cols_;}

		// Insert column at specified column
		template<typename ForwardIterator>
		void insert(size_type colNo, ForwardIterator first)
		{
			enclosingClass_.resize(enclosingClass_.rows_, max(enclosingClass_.cols_, colNo)+1);
			for (size_type i=enclosingClass_.cols_-1; i>colNo; --i)
			{
				copy ((*this)(i-1).begin(), (*this)(i-1).end(), (*this)(i).begin());
			}
			copy (first, first+enclosingClass_.rows_, (*this)(colNo).begin());
		}

		// Erase column at specified column
		void erase(size_type colNo)
		{
			for (size_type i=colNo+1; i<enclosingClass_.cols_; ++i)
			{
				copy ((*this)(i).begin(), (*this)(i).end(), (*this)(i-1).begin());
			}
			enclosingClass_.resize(enclosingClass_.rows_, enclosingClass_.cols_-1);
		}

		// Add column at end of matrix
		template<typename ForwardIterator>
		void push_back(const ForwardIterator first)
		{
			this->insert(enclosingClass_.cols_, first);
		}

		// Add column at beginning of matrix
		template<typename ForwardIterator>
		void push_front(const ForwardIterator first)
		{
			this->insert(0, first);
		}

		// Remove column at end of matrix
		void pop_back()
		{
			this->erase(enclosingClass_.cols_-1);
		}

		// Remove column at beginning of matrix
		void pop_front()
		{
			this->erase(0);
		}

		// Iterators access (Do not change for all matrix types).
		// Must be called by col(i) and not col.
		iterator 		begin() 		{return enclosingClass_.colIteratorStart_[currentIndex_];} 
		const_iterator 	begin() const	{return const_iterator(enclosingClass_.colIteratorStart_[currentIndex_]);}
		iterator 		end() 			{return enclosingClass_.colIteratorFinish_[currentIndex_];}   
		const_iterator 	end() const 	{return const_iterator(enclosingClass_.colIteratorFinish_[currentIndex_]);}

		reverse_iterator		rbegin()		{return reverse_iterator(enclosingClass_.colIteratorFinish_[currentIndex_]);}
		const_reverse_iterator	rbegin() const	{return const_reverse_iterator(enclosingClass_.colIteratorFinish_[currentIndex_]);}
		reverse_iterator		rend()			{return reverse_iterator(enclosingClass_.colIteratorStart_[currentIndex_]);}
		const_reverse_iterator	rend() const	{return const_reverse_iterator(enclosingClass_.colIteratorStart_[currentIndex_]);}

	private: // Member variables (Do not change for all matrix types).
		EnclosingClass& enclosingClass_;
		mutable size_type currentIndex_;

	private: // Undefined constructors and copy constructors.
		Col();
		Col(const Self&);
	} col;
	friend class Col;
	
public:	// STL-like typedefs.
	typedef ValueT     			value_type;
	typedef ValueT&      		reference;
	typedef const ValueT& 		const_reference;
	typedef ValueT*       		pointer;
	typedef const ValueT* 		const_pointer;
	typedef ptrdiff_t			difference_type;
	typedef size_t				size_type;

public: // Constructors, destructor, copy constructors and assignment.
	// Default constructor.
	DenseMatrix() : rows_(0), cols_(0), matrix_(new MatrixStorageType[rows_*cols_]),
		row(*this), col(*this),
		rowIteratorStart_(0), rowIteratorFinish_(0), colIteratorStart_(0), colIteratorFinish_(0)
	{
	};
	// Constructor requiring number of rows, Columns and initial value for elements.
	DenseMatrix(size_type sizerow, size_type sizeCol, value_type x=value_type()) :
		rows_(sizerow), cols_(sizeCol), matrix_(new MatrixStorageType[rows_*cols_]),
		row(*this), col(*this),
		rowIteratorStart_(0), rowIteratorFinish_(0), colIteratorStart_(0), colIteratorFinish_(0)
	{
		ResizeIterators();
		for (size_type i=0; i<rows_; ++i) fill(row(i).begin(), row(i).end(), x);
	}
	// Constructor requiring number of rows, Columns and a 1D array containing initial values for elements.
	DenseMatrix(size_type sizerow, size_type sizeCol, value_type* am) :
		rows_(sizerow), cols_(sizeCol), matrix_(new MatrixStorageType[rows_*cols_]),
		row(*this), col(*this),
		rowIteratorStart_(0), rowIteratorFinish_(0), colIteratorStart_(0), colIteratorFinish_(0)
	{
		ResizeIterators();
		*this = am;
	}
	// Constructor requiring number of rows, Columns and a 2D array containing initial values for elements.
	DenseMatrix(size_type sizerow, size_type sizeCol, value_type** am) :
		rows_(sizerow), cols_(sizeCol), matrix_(new MatrixStorageType[rows_*cols_]),
		row(*this), col(*this),
		rowIteratorStart_(0), rowIteratorFinish_(0), colIteratorStart_(0), colIteratorFinish_(0)
	{
		ResizeIterators();
		*this = am;
	}
	// Destructor.
	virtual ~DenseMatrix() 
	{
		if (matrix_) delete[] matrix_;
		if (rowIteratorStart_) delete[] rowIteratorStart_;
		if (rowIteratorFinish_) delete[] rowIteratorFinish_;
		if (colIteratorStart_) delete[] colIteratorStart_;
		if (colIteratorFinish_) delete[] colIteratorFinish_;
	}
	// Copy constructor.
	DenseMatrix(const Self& dm) : 
		row(*this), col(*this),
		rowIteratorStart_(0), rowIteratorFinish_(0), colIteratorStart_(0), colIteratorFinish_(0)
	{
		size_type totalElements = dm.row.size() * dm.col.size();
		matrix_ = new MatrixStorageType[totalElements];
		for (size_type i=0; i<totalElements; ++i) matrix_[i] = dm.matrix_[i];
		rows_ = dm.rows_;
		cols_ = dm.cols_;
		ResizeIterators();
	}
	// Self = Same type of matrix.
	Self& operator=(Self dm) {swap(dm); return *this;}
	// Self = Any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	Self& operator=(const Matrix<ValueTT,MatrixTypePT,MathPT>& m) 
	{
		size_type sizerow = m.row.size();
		size_type sizeCol = m.col.size();
		Self tempMatrix(sizerow, sizeCol);
		for (size_type i=0; i<sizerow; ++i)
		{
			copy (m.row(i).begin(), m.row(i).end(), tempMatrix.row(i).begin());
		}
		swap(tempMatrix);
		return *this;
	}
	Self& operator=(const value_type* am)
	// PURPOSE: Copy a 1D array into the matrix.
	// REQUIRE: The number of array elements must be equal or larger than the number of
	// REQUIRE: matrix elements.
	// REQUIRE: The sequence of array elements must be concatenation of rows to be put
	// REQUIRE: into the matrix.
	// PROMISE: The matrix will be filled with elements from the array.
	{
		size_type totalElements = row.size() * col.size();
		for (size_type i=0; i<totalElements; ++i)
		{
			matrix_[i] = am[i];
		}
		return *this;
	}
	Self& operator=(value_type** am)
	// PURPOSE: Copy a 2D array into the matrix.
	// REQUIRE: The size of the array must be the same as the matrix.
	// PROMISE: The matrix will be filled with elements from the array.
	{
		size_type sizerow = row.size();
		size_type sizeCol = col.size();
		for (size_type i=0; i<sizerow; ++i)
		{
			//copy (am[i], am[i]+sizeCol, row(i).begin()); 
			std::copy(am[i], am[i]+sizeCol, stdext::checked_array_iterator<value_type*>(row(i).begin(), sizeCol));
			//std::copy(_copy.distances, _copy.distances + _copy.sizeOfmassives, stdext::checked_array_iterator<double*> (distances,_copy.sizeOfmassives));
		}
		return *this;
	}

public: // Functions
	void resize(size_type sizerow, size_type sizeCol, value_type x=value_type())
	// PURPOSE: Resize matrix, filling new elements with supplied or default values and
	// PURPOSE: removing extra elements.
	// REQUIRE: None.
	// PROMISE: Matrix will be resized to the desired size.
	{
		if (sizerow != rows_ || sizeCol != cols_)
		{
			Self tempMatrix(sizerow, sizeCol, x);
			for (size_type i=min(sizerow, rows_); i--;)
			{
				for (size_type j=min(sizeCol, cols_); j--; )
				{
					tempMatrix(i,j) = (*this)(i,j);
				}
			}
			swap(tempMatrix);
		}
	}
	void swap(Self& dm)
	// PURPOSE: Swap the contents of another matrix with the present matrix.
	// REQUIRE: The matrix to be swapped must be of the same type.
	// PROMISE: The entire contents of both matrices will be swapped.
	{
		std::swap(matrix_, dm.matrix_);
		std::swap(rows_, dm.rows_);
		std::swap(cols_, dm.cols_);
		std::swap(rowIteratorStart_, dm.rowIteratorStart_);
		std::swap(rowIteratorFinish_, dm.rowIteratorFinish_);
		std::swap(colIteratorStart_, dm.colIteratorStart_);
		std::swap(colIteratorFinish_, dm.colIteratorFinish_);
	}

public: // Friend functions.
	friend ostream& operator<<(ostream& os, const Self& dm)
	// PURPOSE: Print the elements of the matrix to output stream.
	// REQUIRE: None.
	// PROMISE: All the elements of the matrix will be output to the output stream.
	{
		for (size_t i=0; i<dm.rows_; ++i)
		{
			for (size_t j=0; j<dm.cols_; ++j)
			{
				os << dm(i,j) << "\t";
			}
			os << "\n";
		}
		os << "\n";
		return os;
	}

public: // Comparison operators.
	bool operator==(const Self& m) const 
	{
		if (rows_!=m.rows_ && cols_!=m.cols_) return false;
		else
		{
			size_type totalElements = rows_ * cols_;
			for (size_type i=0; i<totalElements; ++i)
			{
				// Element-wise comparison.
				if (matrix_[i] != m.matrix_[i]) return false;
			}
		}
		return true;
	}
	bool operator!=(const Self& m) const {return !(*this == m);} 

public: // Functions.
	// Update matrix details.
	// Not used by this matrix because normally nothing needs to be updated.
	void Update() {}; 
	
public: // Member variables.
	// Exposed to user to allow user to use it directly for performance critical routines.
	size_type rows_;	// Total number of rows in matrix.
	size_type cols_;	// Total number of Columns in matrix.
	MatrixType matrix_;	// Matrix elements storage.
	
public:
	void ResizeIterators()
	// PURPOSE: Resize all internal iterators.
	// REQUIRE: None.
	// PROMISE: All internal iterators will be resized according to the matrix size.
	{
		if (rowIteratorStart_) delete[] rowIteratorStart_;
		if (rowIteratorFinish_) delete[] rowIteratorFinish_;
		if (colIteratorStart_) delete[] colIteratorStart_;
		if (colIteratorFinish_) delete[] colIteratorFinish_;

		rowIteratorStart_ = new typename Row::iterator[rows_];
		rowIteratorFinish_ = new typename Row::iterator[rows_];
		colIteratorStart_ = new typename Col::iterator[cols_];
		colIteratorFinish_ = new typename Col::iterator[cols_];

		UpdateIterators();
	}

protected:
	void UpdateIterators()
	// PURPOSE: Update all internal iterators.
	// REQUIRE: None.
	// PROMISE: All internal iterators will be updated.
	{
		for (size_type i=0; i<rows_; ++i)
		{
			rowIteratorStart_[i] = matrix_ + (i*cols_);
			rowIteratorFinish_[i] = matrix_ + ((i+1)*cols_);
		}
		size_type total = rows_*cols_;
		for (size_type j=0; j<cols_; ++j)
		{
			colIteratorStart_[j].SetCols(cols_);
			colIteratorStart_[j] = matrix_ + j;

			colIteratorFinish_[j].SetCols(cols_);
			colIteratorFinish_[j] = matrix_ + j + total;
		}
	}

public:	// Member access (Do not change for all matrix types).
	const_reference		operator()(size_type posrow, 
								   size_type posCol) const	{return rowIteratorStart_[posrow][difference_type(posCol)];}
	reference			operator()(size_type posrow, 
								   size_type posCol)		{return rowIteratorStart_[posrow][difference_type(posCol)];}
	const_reference		at(size_type posrow, 
						   size_type posCol) const			{RangeCheck(posrow, posCol); return (*this)(posrow,posCol);}
	reference			at(size_type posrow, 
						   size_type posCol)				{RangeCheck(posrow, posCol); return (*this)(posrow,posCol);}
protected: // Supporting functions for member access (Do not change for all matrix types).
	void RangeCheck(size_type posrow, size_type posCol) {if (posrow>=row.size() || posCol>=col.size()) throw("Out of range");}

public: // Functions (Do not change for all matrix types).
	// Erase matrix.
	void clear() {resize(0, 0);} 
	// Check if matrix is empty.
	bool empty() const {return (row.size()==0 || col.size()==0);} 

protected: // Essential member variables (Do not change for all matrix types).
	typename Row::iterator* rowIteratorStart_;	// Starting point of row iterators.
	typename Row::iterator* rowIteratorFinish_;	// Ending point of row iterators.
	typename Col::iterator* colIteratorStart_;	// Starting point of column iterators.
	typename Col::iterator* colIteratorFinish_;	// Ending point of column iterators.
};

template<typename ValueT, typename MatrixT> 
class MathMatrix
// PURPOSE: Provide mathematical functions for matrices.
// REQUIRE: ValueT must be a numerical type.
// REQUIRE: Not to be constructed directly by user.
// PROMISE: Mathematical manipulation can be done on the matrices.
{
private: // Typedefs for internal use.
	typedef MathMatrix<ValueT,MatrixT>	Self;
	typedef ValueT						value_type;

protected: // Constructors, destructor, copy constructors and assignment.
	// No default constructor.
	// Constructor requiring address of derived (child) matrix i.e. address of Matrix<ValueT,MatrixTypeP,MathP>.
	MathMatrix(MatrixT& matrix) : childMatrix_(matrix) {};
	// Destructor.
	virtual ~MathMatrix() {};

public: // Mathematical operators
	/*********************************** Addition *************************************/
	// Self += Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT& operator+=(const Matrix<ValueT,MatrixTypePT,MathPT>& mm) // Add two matrices
	{
		size_t sizerow = childMatrix_.row.size();
		if (sizerow!=mm.row.size() && childMatrix_.col.size()!=mm.col.size()) throw("Matrix size mismatch!");
		for (size_t i=0; i<sizerow; ++i)
		{
			 transform(childMatrix_.row(i).begin(), childMatrix_.row(i).end(), mm.row(i).begin(), childMatrix_.row(i).begin(), plus<value_type>());
		}
		childMatrix_.Update();
		return childMatrix_;
	}
	// Self + Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT operator+(const Matrix<ValueT,MatrixTypePT,MathPT>& v) const {MatrixT w(childMatrix_); w += v; return w;} // Add two matrices

	/*********************************** Subtraction *************************************/
	// Self -= Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT& operator-=(const Matrix<ValueT,MatrixTypePT,MathPT>& mm) // Subtract two matrices
	{
		size_t sizerow = childMatrix_.row.size();
		if (sizerow!=mm.row.size() && childMatrix_.col.size()!=mm.col.size()) throw("Matrix size mismatch!");
		for (size_t i=0; i<sizerow; ++i)
		{
			 transform(childMatrix_.row(i).begin(), childMatrix_.row(i).end(), mm.row(i).begin(), childMatrix_.row(i).begin(), minus<value_type>());
		}
		childMatrix_.Update();
		return childMatrix_;
	}
	// Self - Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT operator-(const Matrix<ValueT,MatrixTypePT,MathPT>& v) const {MatrixT w(childMatrix_); w -= v; return w;} // Subtract two matrices

	/*********************************** Multiplication *************************************/
	// Self *= Any matrix (Dot product).
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT& operator*=(const Matrix<ValueT,MatrixTypePT,MathPT>& x)
	{
		childMatrix_ = childMatrix_ * x;
		return childMatrix_;
	}
	// Self * Any matrix (Dot product).
	MatrixT operator*(const MatrixT& x) const
	{
		// Code is put in here instead of in operator*= to speed up matrix multiplication
		// by avoid unnecessary creation of temporaries.
		// This technique is not useful for addition and substraction because they didn't
		// require creation of a tempMatrix in their code.
		size_t sizerow1 = childMatrix_.row.size();
		size_t sizeCol1 = childMatrix_.col.size();
		size_t sizerow2 = x.row.size();
		size_t sizeCol2 = x.col.size();

		if (sizeCol1!=sizerow2) throw("Matrix size mismatch!");
		MatrixT tempMatrix(sizerow1, sizeCol2);
		for (size_t i=0; i<sizerow1; ++i)
		{
			for (size_t j=0; j<sizeCol2; ++j)
			{
				for (size_t k=0; k<sizeCol1; ++k)
				{
					tempMatrix(i, j) += childMatrix_(i,k) * x(k,j);
				}
			}
		}
		return tempMatrix;
	}
	// Self *= Any value (Scalar product).
	MatrixT& operator*=(value_type v)
	{
		size_t sizerow = childMatrix_.row.size();
		for (size_t i=0; i<sizerow; ++i)
		{
			 transform(childMatrix_.row(i).begin(), childMatrix_.row(i).end(), childMatrix_.row(i).begin(), std::bind2nd(multiplies<value_type>(), v));
		}
		childMatrix_.Update();
		return childMatrix_;
	}
	// Self * Any value (Scalar product).
	MatrixT operator*(value_type v) const {MatrixT w(childMatrix_); w *= v; return w;}

private: // Member variables.
	MatrixT& childMatrix_;

private: // Undefined constructors and copy constructors.
	MathMatrix();
	MathMatrix(const Self&);
};

template<typename ExprT>
class MET
// PURPOSE: Matrix Expression Template (MET).
// PURPOSE: Provide expression templates to speed up mathematical operator calculations like
// PURPOSE: MatrixA = MatrixB + MatrixC - MatrixD + MatrixE.
// REQUIRE: Not to be constructed directly by user.
// PROMISE: Addition and subtraction may be speed up.
{
private: // Typedefs for internal use.
	typedef MET<ExprT>			Self;
	typedef typename ExprT::value_type	value_type;

public: // Memberspaces.
	class row
	{
	private: // Typedefs for internal use.
		typedef MET<ExprT> EnclosingClass;
		friend class EnclosingClass;

	private:	// Constructors, destructor, copy constructors and assignment.
		// Constructor requiring address of enclosing class i.e. address of MET<ExprT>.
		row(EnclosingClass& p) : enclosingClass_(p) {};

	public: // Functions.
		size_t size() const
		// PURPOSE: Return total number of rows of the resultant matrix as a result
		// PURPOSE: of operations performed by the expression.
		// REQUIRE: All matrices in expression must be of the same size.
		// PROMISE: The total number of rows will be returned.
		{return enclosingClass_.expr_.row.size();}

	private: // Member variables.
		EnclosingClass& enclosingClass_;
	} row;

	class col
	{
	private: // Typedefs for internal use.
		typedef MET<ExprT> EnclosingClass;
		friend class EnclosingClass;

	private:	// Constructors, destructor, copy constructors and assignment.
		// Constructor requiring address of enclosing class i.e. address of MET<ExprT>.
		col(EnclosingClass& p) : enclosingClass_(p) {};

	public: // Functions.
		size_t size() const 
		// PURPOSE: Return total number of Columns of the resultant matrix as a result
		// PURPOSE: of operations performed by the expression.
		// REQUIRE: All matrices in expression must be of the same size.
		// PROMISE: The total number of Columns will be returned.
		{return enclosingClass_.expr_.col.size();}

	private: // Member variables.
		EnclosingClass& enclosingClass_;
	} col;
	
public:	// Constructors, destructor, copy constructors and assignment.
	// No default constructor.
	// Constructor requiring an expression.
	MET(ExprT expr) : expr_(expr), row(*this), col(*this) {};
	// Copy constructor.
	MET(const Self& s) : expr_(s.expr_), row(*this), col(*this) {};
	// Destructor.
	~MET() {};

public: // Functions.
	// Get the resultant value of element(i,j) as a result of operations performed
	// by the expression.
	value_type operator()(size_t i, size_t j) const {return (expr_(i, j));}

private: // Member variables.
	ExprT expr_;	// Expression.
};

template<typename ValueT, typename MatrixT> 
class MathMETMatrix
// PURPOSE: Provide mathematical functions for matrices.
// PURPOSE: Provide expression templates to speed up mathematical operator calculations.
// REQUIRE: ValueT must be a numerical type.
// REQUIRE: Not to be constructed directly by user.
// REQUIRE: User must ensure matrices must be of the same size for addition and subtraction.
// PROMISE: Mathematical manipulation can be done on the matrices.
// PROMISE: Mathematical operator calculations may speed up or slow down depending on
// PROMISE: the formula. Test with both MathMETMatrix and MathMatrix to find out.
{
private: // Typedefs for internal use.
	typedef MathMETMatrix<ValueT,MatrixT>	Self;
	typedef ValueT							value_type;

protected: // Supporting classes for expression template.
	template<typename OperatorT, typename LeftT, typename RightT>
	class METElementBinaryOp
	// PURPOSE: Provide holder for binary operations involving individual elements of matrix
	// REQUIRE: None.
	// PROMISE: None.
	{
	public: // Typedefs for internal use.
		typedef METElementBinaryOp<OperatorT,LeftT,RightT> Self;
		typedef typename OperatorT::value_type value_type;

	public: // Memberspaces.
		class row
		{
		private: // Typedefs for internal use.
			typedef METElementBinaryOp<OperatorT,LeftT,RightT> EnclosingClass;
			friend class EnclosingClass;

		private:	// Constructors, destructor, copy constructors and assignment.
			// Constructor requiring address of enclosing class i.e. address of METElementBinaryOp<OperatorT,LeftT,RightT>.
			row(EnclosingClass& p) : enclosingClass_(p) {};

		public: // Functions.
			size_t size() const 
			// PURPOSE: Return total number of Columns of the resultant matrix as a result
			// PURPOSE: of operations performed by the expression.
			// REQUIRE: All matrices in expression must be of the same size.
			// PROMISE: The total number of Columns will be returned.
			{return enclosingClass_.left_->row.size();}

		private: // Member variables.
			EnclosingClass& enclosingClass_;
		} row;

		class col
		{
		private: // Typedefs for internal use.
			typedef METElementBinaryOp<OperatorT,LeftT,RightT> EnclosingClass;
			friend class EnclosingClass;

		private:	// Constructors, destructor, copy constructors and assignment.
			// Constructor requiring address of enclosing class i.e. address of METElementBinaryOp<OperatorT,LeftT,RightT>.
			col(EnclosingClass& p) : enclosingClass_(p) {};

		public: // Functions.
			size_t size() const 
			// PURPOSE: Return total number of Columns of the resultant matrix as a result
			// PURPOSE: of operations performed by the expression.
			// REQUIRE: All matrices in expression must be of the same size.
			// PROMISE: The total number of Columns will be returned.
			{return enclosingClass_.left_->col.size();}

		private: // Member variables.
			EnclosingClass& enclosingClass_;
		} col;

	public:	// Constructors, destructor, copy constructors and assignment.
		// Constructor requiring address of left and right expression.
		METElementBinaryOp(const LeftT* left, const RightT* right) : row(*this), col(*this) {left_ = left; right_ = right;}
		// Copy constructor.
		METElementBinaryOp(const Self& s) : left_(s.left_), right_(s.right_), row(*this), col(*this) {};
		// Destructor.
		~METElementBinaryOp() {};

	public: // Functions.
		// Get the resultant value of element(i,j) as a result of operations performed
		// by the operator.
		value_type operator()(size_t i, size_t j) const {return OperatorT::Evaluate((*left_)(i,j), (*right_)(i,j));}

	private: // Member variables.
		const LeftT* left_;
		const RightT* right_;
	};

	struct METAdd
	// PURPOSE: Addition operator of expression template.
	// PURPOSE: Provide static function to add two elements of matrix/matrices together.
	// REQUIRE: None.
	// PROMISE: The result of the addition will be returned.
	{
		typedef ValueT value_type;
		static value_type Evaluate(value_type x, value_type y) {return x+y;}	
	};

	struct METSubtract
	// PURPOSE: Subtraction operator of expression template.
	// PURPOSE: Provide static function to subtract two elements of matrix/matrices together.
	// REQUIRE: None.
	// PROMISE: The result of the subtraction will be returned.
	{
		typedef ValueT value_type;
		static value_type Evaluate(value_type x, value_type y) {return x-y;}		
	};

	struct METMultiply
	// PURPOSE: Multiplication operator of expression template.
	// PURPOSE: Provide static function to multiply two elements of matrix/matrices together.
	// REQUIRE: None.
	// PROMISE: The result of the multiplication will be returned.
	{
		typedef ValueT value_type;
		static value_type Evaluate(value_type x, value_type y) {return x*y;}
	};

protected: // Constructors, destructor, copy constructors and assignment.
	// No default constructor.
	// Constructor requiring address of derived (child) matrix i.e. address of Matrix<ValueT,MatrixTypeP,MathP>.
	MathMETMatrix(MatrixT& matrix) : childMatrix_(matrix) {};
	// Destructor
	virtual ~MathMETMatrix() {};

public: // Mathematical operators.
	/*********************************** Addition *************************************/
	// Self += Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT& operator+=(const Matrix<ValueT,MatrixTypePT,MathPT>& x) 
	{
		childMatrix_ = childMatrix_ + x;
		return childMatrix_;
	}
	// Self += Any expression.
	template<typename ExprT>
	MatrixT& operator+=(const MET<ExprT>& x) 
	{
		childMatrix_ = childMatrix_ + x;
		return childMatrix_;
	}
	// Self + Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MET<METElementBinaryOp<METAdd,MatrixT,Matrix<ValueT,MatrixTypePT,MathPT> > > 
	operator+(const Matrix<ValueT,MatrixTypePT,MathPT>& y)
	{
		typedef METElementBinaryOp<METAdd,MatrixT,Matrix<ValueT,MatrixTypePT,MathPT> > ExprT;
		return MET<ExprT>(ExprT(&childMatrix_,&y));
	}
	// Any expression + Any matrix.
	template<typename A,
		     template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	friend MET<METElementBinaryOp<METAdd,MET<A>,Matrix<ValueT,MatrixTypePT,MathPT> > > 
	operator+(const MET<A>& x, const Matrix<ValueT,MatrixTypePT,MathPT>& y)
	{
		typedef METElementBinaryOp<METAdd,MET<A>,Matrix<ValueT,MatrixTypePT,MathPT> > ExprT;
		return MET<ExprT>(ExprT(&x,&y));
	}
	// Any matrix + Any expression.
	template<typename A,
		     template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	friend MET<METElementBinaryOp<METAdd,Matrix<ValueT,MatrixTypePT,MathPT>,MET<A> > > 
	operator+(const Matrix<ValueT,MatrixTypePT,MathPT>& x, const MET<A>& y)
	{
		typedef METElementBinaryOp<METAdd,Matrix<ValueT,MatrixTypePT,MathPT>,MET<A> > ExprT;
		return MET<ExprT>(ExprT(&x,&y));
	}
	// Any expression + Any expression.
	template<typename A, typename B>
	friend MET<METElementBinaryOp<METAdd,MET<A>,MET<B> > > 
	operator+(const MET<A>& x, const MET<B>& y)
	{
		typedef METElementBinaryOp<METAdd,MET<A>,MET<B> > ExprT;
		return MET<ExprT>(ExprT(&x,&y));
	}

	/*********************************** Subtraction *************************************/
	// Self -= Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT& operator-=(const Matrix<ValueT,MatrixTypePT,MathPT>& x) 
	{
		childMatrix_ = childMatrix_ - x;
		return childMatrix_;
	}
	// Self -= Any expression.
	template<typename ExprT>
	MatrixT& operator-=(const MET<ExprT>& x) 
	{
		childMatrix_ = childMatrix_ - x;
		return childMatrix_;
	}
	// Self - Any matrix.
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MET<METElementBinaryOp<METSubtract,MatrixT,Matrix<ValueT,MatrixTypePT,MathPT> > > 
	operator-(const Matrix<ValueT,MatrixTypePT,MathPT>& y)
	{
		typedef METElementBinaryOp<METSubtract,MatrixT,Matrix<ValueT,MatrixTypePT,MathPT> > ExprT;
		return MET<ExprT>(ExprT(&childMatrix_,&y));
	}
	// Any expression - Any matrix.
	template<typename A,
		     template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	friend MET<METElementBinaryOp<METSubtract,MET<A>,Matrix<ValueT,MatrixTypePT,MathPT> > > 
	operator-(const MET<A>& x, const Matrix<ValueT,MatrixTypePT,MathPT>& y)
	{
		typedef METElementBinaryOp<METSubtract,MET<A>,Matrix<ValueT,MatrixTypePT,MathPT> > ExprT;
		return MET<ExprT>(ExprT(&x,&y));
	}
	// Any matrix - Any expression.
	template<typename A,
		     template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	friend MET<METElementBinaryOp<METSubtract,Matrix<ValueT,MatrixTypePT,MathPT>,MET<A> > > 
	operator-(const Matrix<ValueT,MatrixTypePT,MathPT>& x, const MET<A>& y)
	{
		typedef METElementBinaryOp<METSubtract,Matrix<ValueT,MatrixTypePT,MathPT>,MET<A> > ExprT;
		return MET<ExprT>(ExprT(&x,&y));
	}
	// Any expression - Any expression.
	template<typename A, typename B>
	friend MET<METElementBinaryOp<METSubtract,MET<A>,MET<B> > > 
	operator-(const MET<A>& x, const MET<B>& y)
	{
		typedef METElementBinaryOp<METSubtract,MET<A>,MET<B> > ExprT;
		return MET<ExprT>(ExprT(&x,&y));
	}

	/*********************************** Multiplication *************************************/
	// Self *= Any matrix (Dot product).
	template<template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	MatrixT& operator*=(const Matrix<ValueT,MatrixTypePT,MathPT>& x)
	{
		childMatrix_ = childMatrix_ * x;
		return childMatrix_;
	}
	// Self * Any matrix (Dot product).
	MatrixT operator*(const MatrixT& x) const
	{
		// Code is put in here instead of in operator*= to speed up matrix multiplication
		// by avoid unnecessary creation of temporaries.
		size_t sizerow1 = childMatrix_.row.size();
		size_t sizeCol1 = childMatrix_.col.size();
		size_t sizerow2 = x.row.size();
		size_t sizeCol2 = x.col.size();

		if (sizeCol1!=sizerow2) throw("Matrix size mismatch!");
		MatrixT tempMatrix(sizerow1, sizeCol2);
		for (size_t i=0; i<sizerow1; ++i)
		{
			for (size_t j=0; j<sizeCol2; ++j)
			{
				for (size_t k=0; k<sizeCol1; ++k)
				{
					tempMatrix(i, j) += childMatrix_(i,k) * x(k,j);
				}
			}
		}
		return tempMatrix;
	}
	// Self *= Any value (Scalar product).
	MatrixT& operator*=(value_type v)
	{
		size_t sizerow = childMatrix_.row.size();
		size_t sizeCol = childMatrix_.col.size();
		for (size_t i=0; i<sizerow; ++i)
		{
			for (size_t j=0; j<sizeCol; ++j)
			{
				childMatrix_(i,j) *= v;
			}
		}
		childMatrix_.Update();
		return childMatrix_;
	}
	// Self * Any value (Scalar product).
	MatrixT operator*(value_type v) const {MatrixT w(childMatrix_); w *= v; return w;}

private: // Member variables.
	MatrixT& childMatrix_;

private: // Undefined constructors and copy constructors.
	MathMETMatrix();
	MathMETMatrix(const Self&);
};

template<typename ValueT, typename MatrixT> 
class NonMathMatrix
// PURPOSE: Indicate matrix does not have mathematical functions.
// REQUIRE: Not to be constructed directly by user.
// PROMISE: Mathematical functions cannot be performed on matrix.
{
private: // Typedefs for internal use.
	typedef NonMathMatrix<ValueT,MatrixT>	Self;
	typedef ValueT							value_type;

protected: // Constructors, destructor, copy constructors and assignment.
	// No default constructor.
	// Constructor requiring address of derived (child) matrix i.e. address of Matrix<ValueT,MatrixTypeP,MathP>.
	NonMathMatrix(MatrixT& matrix) : childMatrix_(matrix) {};
	// Destructor.
	virtual ~NonMathMatrix() {};

private: // Member variables.
	MatrixT& childMatrix_;

private: // Undefined constructors and copy constructors.
	NonMathMatrix();
	NonMathMatrix(const Self&);
};

template<typename ValueT,
		 template<typename> class MatrixTypeP,
		 template<typename,typename> class MathP>
class Matrix : public MatrixTypeP<ValueT>, public MathP<ValueT,Matrix<ValueT, MatrixTypeP, MathP> >
// PURPOSE: Provide a unified interface for all types of matrices.
// REQUIRE: User should be mindful of the matrix sizes if mathematical calculations is required.
// PROMISE: Each type of matrices can be constructed similarly and yet have different capabilities.
{
private: // Typedefs for internal use.
	typedef Matrix<ValueT,MatrixTypeP,MathP>	Self;
	typedef MathP<ValueT,Self>					MathBase;
	typedef MatrixTypeP<ValueT>					StorageBase;

public: // STL-like typedefs.
	typedef typename StorageBase::value_type     	value_type;
	typedef typename StorageBase::reference      	reference;
	typedef typename StorageBase::const_reference 	const_reference;
	typedef typename StorageBase::pointer     		pointer;
	typedef typename StorageBase::const_pointer		const_pointer;
	typedef typename StorageBase::difference_type	difference_type;
	typedef typename StorageBase::size_type			size_type;

public: // Constructors, destructor, copy constructors and assignment.
	// Default constructor.
	Matrix() : MathBase(*this) {};
	// Constructor requiring number of rows, Columns and initial value for elements.
	Matrix(size_type sizerow, size_type sizeCol, value_type x=value_type())
		: MathBase(*this), StorageBase(sizerow, sizeCol, x) {};
	// Constructor requiring number of rows, Columns and a 1D array containing initial values for elements.
	Matrix(size_type sizerow, size_type sizeCol, value_type* am)
		: MathBase(*this), StorageBase(sizerow, sizeCol, am) {};
	// Constructor requiring number of rows, Columns and a 2D array containing initial values for elements.
	Matrix(size_type sizerow, size_type sizeCol, value_type** am)
		: MathBase(*this), StorageBase(sizerow, sizeCol, am) {};
	// Destructor.
	virtual ~Matrix() {};
	// Copy constructor.
	Matrix(const Self& m) : MathBase(*this), StorageBase(m) {};
	// Copy constructor taking any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	Matrix(const Matrix<ValueTT,MatrixTypePT,MathPT>& m) : MathBase(*this) {StorageBase::operator=(m);}
	// Self = Same type of matrix.
	Self& operator=(Self m) {swap(m); return *this;}
	// Self = Any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	Self& operator=(const Matrix<ValueTT,MatrixTypePT,MathPT>& m) {StorageBase::operator=(m); return *this;}
	Self& operator=(const value_type* am) 
	// PURPOSE: Copy a 1D array into the matrix.
	// REQUIRE: The number of array elements must be equal or larger than the number of
	// REQUIRE: matrix elements.
	// REQUIRE: The sequence of array elements must be concatenation of rows to be put
	// REQUIRE: into the matrix.
	// PROMISE: The matrix will be filled with elements from the array.
	{StorageBase::operator=(am); return *this;}
	Self& operator=(value_type** am) 
	// PURPOSE: Copy a 2D array into the matrix.
	// REQUIRE: The size of the array must be the same as the matrix.
	// PROMISE: The matrix will be filled with elements from the array.
	{StorageBase::operator=(am); return *this;}
	// Self = Any expression.
	template<typename ExprT>
	Self& operator=(const MET<ExprT>& expr)
	{
		size_t sizerow = expr.row.size();
		size_t sizeCol = expr.col.size();
		
		// Profiling indicates creation of this temporary matrix is a time-consuming step
		// because of the ResizeIterators() function during the creation
		// However this temporary matrix is essential to ensure safety during assignment
		Self tempMatrix(sizerow,sizeCol);
		
		for (size_t i=0; i<sizerow; ++i)
		{
			for (size_t j=0; j<sizeCol; ++j)
			{
				tempMatrix(i,j) = expr(i,j);		
			}
		}
		swap(tempMatrix);
		Update();
		return *this;
	}

public: // Functions.
	// Swap two rows in matrix.
	void swaprows(size_type i1, size_type i2)
	{
		typename Self::Row::iterator iti1 = row(i1).begin();
		typename Self::Row::iterator iti2 = row(i2).begin();
		for (size_type j=0; j<cols_; ++j, ++iti1, ++iti2) std::swap(*iti1, *iti2);
		Update();
	}
	// Swap two Columns in matrix.
	void swapcols(size_type j1, size_type j2)
	{
		typename Self::Col::iterator itj1 = col(j1).begin();
		typename Self::Col::iterator itj2 = col(j2).begin();
		for (size_type i=0; i<rows_; ++i, ++itj1, ++itj2) std::swap(*itj1, *itj2);
		Update();
	}
};

template<typename ValueT,
		 template<typename> class MatrixTypeP,
		 template<typename,typename> class MathP>
class RowVector : public Matrix<ValueT,MatrixTypeP,MathP>
// PURPOSE: Provide a specialised form of a matrix. 
// REQUIRE: User should be mindful of the vector sizes if mathematical calculations is required.
// PROMISE: Each type of vector can be constructed similarly and yet have different capabilities.
{
private: // Typedefs for internal use.
	typedef RowVector<ValueT,MatrixTypeP,MathP>	Self;
	typedef Matrix<ValueT,MatrixTypeP,MathP>	Base;

public: // STL-like typedefs.
	typedef typename Base::value_type     		value_type;
	typedef typename Base::reference      		reference;
	typedef typename Base::const_reference 		const_reference;
	typedef typename Base::pointer     			pointer;
	typedef typename Base::const_pointer		const_pointer;
	typedef typename Base::difference_type		difference_type;
	typedef typename Base::size_type			size_type;

public: // Vector-only STL-like typedefs.
	typedef typename Base::Row::iterator				iterator;
	typedef typename Base::Row::const_iterator			const_iterator;
	typedef typename Base::Row::reverse_iterator		reverse_iterator;
	typedef typename Base::Row::const_reverse_iterator	const_reverse_iterator;

public: // Constructors, destructor, copy constructors and assignment.
	// Default constructor.
	RowVector() : Base() {};
	// Destructor.
	virtual ~RowVector() {};
	// Copy constructor.
	RowVector(const Self& m) : Base(m) {};
	// Copy constructor taking any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	RowVector(const Matrix<ValueTT,MatrixTypePT,MathPT>& m) : Base(m) {};
	// Self = Same type of RowVector.
	Self& operator=(Self m) {swap(m); return *this;}
	// Self = Any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	Self& operator=(const Matrix<ValueTT,MatrixTypePT,MathPT>& m) {Base::operator=(m); return *this;}
	Self& operator=(const value_type* am) 
	// PURPOSE: Copy a 1D array into the RowVector.
	// REQUIRE: The number of array elements must be equal or larger than the number of
	// REQUIRE: RowVector elements.
	// PROMISE: The RowVector will be filled with elements from the array.
	{Base::operator=(am); return *this;}
	// Self = Any expression.
	template<typename ExprT>
	Self& operator=(MET<ExprT> expr)
	{
		size_t sizeCol = expr.col.size();
		
		// Profiling indicates creation of this temporary matrix is a time-consuming step
		// because of the ResizeIterators() function during the creation
		// However this temporary matrix is essential to ensure safety during assignment
		Self tempMatrix(1,sizeCol);
		
		for (size_t j=0; j<sizeCol; ++j)
		{
			tempMatrix(0,j) = expr(0,j);		
		}
		swap(tempMatrix);
		Update();
		return *this;
	}

public: // Vector-only functions.
	// Constructor requiring size of vector and initial value for elements.
	RowVector(size_type n, value_type x=value_type())
		: Base(1, n, x) {};
	// Constructor requiring size of vector and a 1D array containing initial values for elements.
	RowVector(size_type n, value_type* am)
		: Base(1, n, am) {};
	// Copy constructor taking any RowVector of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	RowVector(const RowVector<ValueTT,MatrixTypePT,MathPT>& m) {Base::operator=(m);}
	// Self = Any RowVector of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	Self& operator=(const RowVector<ValueTT,MatrixTypePT,MathPT>& m) {Base::operator=(m); return *this;}

	const_reference operator[](size_type pos) const {return rowIteratorStart_[0][difference_type(pos)];}
	reference		operator[](size_type pos)		{return rowIteratorStart_[0][difference_type(pos)];}
	const_reference	at(size_type pos) const			{RangeCheck(0, pos); return (*this)[pos];}
	reference		at(size_type pos)				{RangeCheck(0, pos); return (*this)[pos];}
	void push_back(const_reference x) {value_type v[1]; v[0] = x; this->col.push_back(v);}
	void push_front(const_reference x) {value_type v[1]; v[0] = x; this->col.push_front(v);}
	void pop_back() {this->col.pop_back();}
	void pop_front() {this->col.pop_front();}
	void resize(size_type n, value_type x=value_type()) {Base::resize(1, n, x);}
	size_type size() const {return Base::col.size();}


public: // Vector-only iterators access.
	iterator& 			begin() 		{return rowIteratorStart_[0];}  
	const_iterator& 	begin() const	{return rowIteratorStart_[0];}
	iterator& 			end() 			{return rowIteratorFinish_[0];}    
	const_iterator& 	end() const 	{return rowIteratorFinish_[0];}

	reverse_iterator&		rbegin()		{return reverse_iterator(rowIteratorFinish_[0]);}
	const_reverse_iterator&	rbegin() const	{return const_reverse_iterator(rowIteratorFinish_[0]);}
	reverse_iterator&		rend()			{return reverse_iterator(rowIteratorStart_[0]);}
	const_reverse_iterator&	rend() const	{return const_reverse_iterator(rowIteratorStart_[0]);}
};

template<typename ValueT,
		 template<typename> class MatrixTypeP,
		 template<typename,typename> class MathP>
class ColVector : public Matrix<ValueT,MatrixTypeP,MathP>
// PURPOSE: Provide a specialised form of a matrix. 
// REQUIRE: User should be mindful of the vector sizes if mathematical calculations is required.
// PROMISE: Each type of vector can be constructed similarly and yet have different capabilities.
{
private: // Typedefs for internal use.
	typedef ColVector<ValueT,MatrixTypeP,MathP>	Self;
	typedef Matrix<ValueT,MatrixTypeP,MathP>	Base;

public: // STL-like typedefs.
	typedef typename Base::value_type     		value_type;
	typedef typename Base::reference      		reference;
	typedef typename Base::const_reference 		const_reference;
	typedef typename Base::pointer     			pointer;
	typedef typename Base::const_pointer		const_pointer;
	typedef typename Base::difference_type		difference_type;
	typedef typename Base::size_type			size_type;
	
public: // Vector-only STL-like typedefs.
	typedef typename Base::Col::iterator				iterator;
	typedef typename Base::Col::const_iterator			const_iterator;
	typedef typename Base::Col::reverse_iterator		reverse_iterator;
	typedef typename Base::Col::const_reverse_iterator	const_reverse_iterator;


public: // Constructors, destructor, copy constructors and assignment.
	// Default constructor.
	ColVector() : Base() {};
	// Destructor.
	virtual ~ColVector() {};
	// Copy constructor.
	ColVector(const Self& m) : Base(m) {};
	// Copy constructor taking any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	ColVector(const Matrix<ValueTT,MatrixTypePT,MathPT>& m) : Base(m) {};
	// Self = Same type of ColVector.
	Self& operator=(Self m) {swap(m); return *this;}
	// Self = Any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	Self& operator=(const Matrix<ValueTT,MatrixTypePT,MathPT>& m) {Base::operator=(m); return *this;}
	Self& operator=(const value_type* am) 
	// PURPOSE: Copy a 1D array into the ColVector.
	// REQUIRE: The number of array elements must be equal or larger than the number of
	// REQUIRE: ColVector elements.
	// PROMISE: The ColVector will be filled with elements from the array.
	{Base::operator=(am); return *this;}
	// Self = Any expression.
	template<typename ExprT>
	Self& operator=(MET<ExprT> expr)
	{
		size_t sizerow = expr.row.size();
		
		// Profiling indicates creation of this temporary matrix is a time-consuming step
		// because of the ResizeIterators() function during the creation
		// However this temporary matrix is essential to ensure safety during assignment
		Self tempMatrix(sizerow,sizeCol);
		
		for (size_t i=0; i<sizerow; ++i)
		{
			tempMatrix(i,0) = expr(i,0);		
		}
		swap(tempMatrix);
		Update();
		return *this;
	}

public: // Vector-only functions.
	// Constructor requiring size of vector and initial value for elements.
	ColVector(size_type n, value_type x=value_type())
		: Base(n, 1, x) {};
	// Constructor requiring size of vector and a 1D array containing initial values for elements.
	ColVector(size_type n, value_type* am)
		: Base(n, 1, am) {};
	// Constructor requiring size of vector and a 2D array containing initial values for elements.
	ColVector(size_type n, value_type** am)
		: Base(n, 1, am) {};
	// Copy constructor taking any matrix of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	ColVector(const ColVector<ValueTT,MatrixTypePT,MathPT>& m) {Base::operator=(m);}
	// Self = Any ColVector of interconvertible ValueT.
	template<typename ValueTT,template<typename> class MatrixTypePT,template<typename,typename> class MathPT>
	Self& operator=(const ColVector<ValueTT,MatrixTypePT,MathPT>& m) {Base::operator=(m); return *this;}

	const_reference operator[](size_type pos) const {return rowIteratorStart_[pos][0];}
	reference		operator[](size_type pos)		{return rowIteratorStart_[pos][0];}
	const_reference	at(size_type pos) const			{RangeCheck(pos, 0); return (*this)[pos];}
	reference		at(size_type pos)				{RangeCheck(pos, 0); return (*this)[pos];}
	void push_back(const_reference x) {value_type v[1]; v[0] = x; this->row.push_back(v);}
	void push_front(const_reference x) {value_type v[1]; v[0] = x; this->row.push_front(v);}
	void pop_back() {this->row.pop_back();}
	void pop_front() {this->row.pop_front();}
	void resize(size_type n, value_type x=value_type()) {Base::resize(n, 1, x);}
	size_type size() const {return row.size();}

public: // Vector-only iterators access.
	iterator& 			begin() 		{return colIteratorStart_[0];} 
	const_iterator& 	begin() const	{return const_iterator(colIteratorStart_[0]);}
	iterator& 			end() 			{return colIteratorFinish_[0];}   
	const_iterator& 	end() const 	{return const_iterator(colIteratorFinish_[0]);}

	reverse_iterator&		rbegin()		{return reverse_iterator(colIteratorFinish_[0]);}
	const_reverse_iterator&	rbegin() const	{return const_reverse_iterator(colIteratorFinish_[0]);}
	reverse_iterator&		rend()			{return reverse_iterator(colIteratorStart_[0]);}
	const_reverse_iterator&	rend() const	{return const_reverse_iterator(colIteratorStart_[0]);}
};

template<typename MatrixT>
struct Transpose
// PURPOSE: Get the transpose of a matrix, i.e. swap rows with Columns and vice versa.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: The transpose of the matrix will be generated.
{
	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT transposedMatrix;
		(*this)(matrix, transposedMatrix);
		return transposedMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& transposedMatrix)
	{
		typename MatrixT::size_type rows = matrix.row.size();
		transposedMatrix.resize(matrix.col.size(), rows);
		for (typename MatrixT::size_type i=0; i<rows; ++i)
		{
			copy (matrix.row(i).begin(), matrix.row(i).end(), transposedMatrix.col(i).begin());
		}
	}
};

template<typename MatrixT>
struct Diagonal
// PURPOSE: (1) Get the main diagonal of a matrix or 
// PURPOSE: (2) Put a vector into the main diagonal of a matrix.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: Returns a (1) single column matrix or
// PROMISE: (2) a square matrix.
{
	typedef typename MatrixT::value_type value_type;

	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT tempMatrix;
		(*this)(matrix, tempMatrix);
		return tempMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& tempMatrix)
	{
		const size_t maxRows = matrix.row.size();
		const size_t maxCols = matrix.col.size();

		if (maxRows==1 && maxCols!=1)
		{
			// Matrix is a row vector
			tempMatrix.resize(maxCols, maxCols, value_type(0));
			for (size_t i=0; i<maxCols; ++i) tempMatrix(i,i) = matrix(0,i);
		}
		else if (maxRows!=1 && maxCols==1)
		{
			// Matrix is a column vector
			tempMatrix.resize(maxRows, maxRows, value_type(0));
			for (size_t i=0; i<maxRows; ++i) tempMatrix(i,i) = matrix(i,0);
		}
		else
		{
			// Matrix is a matrix
			size_t rows = min(maxRows, maxCols);
			tempMatrix.resize(rows, 1, value_type(0));
			for (size_t i=0; i<rows; ++i) tempMatrix(i,0) = matrix(i,i);
		}
	}
};

template<typename MatrixT>
struct Covariance
// PURPOSE: Get the covariance of a matrix.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: The covariance of the matrix will be generated.
{
	typedef typename MatrixT::value_type value_type;
	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT covMatrix;
		(*this)(matrix, covMatrix);
		return covMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& covMatrix)
	{
		const value_type ZERO(0);
		const value_type ONE(1);
		
		size_t maxRows = matrix.row.size();
		size_t maxCols = matrix.col.size();
		
		if (maxRows == 1) covMatrix = MatrixT(1,1,ZERO);
		else
		{
			MatrixT means;
			YMatrix::Mean<MatrixT>()(matrix, means);
			MatrixT temp = matrix;
			for (size_t j=0; j<maxCols; ++j)
			{
				for (size_t i=0; i<maxRows; ++i) temp(i,j) -= means(0,j);
			}
			covMatrix = Transpose<MatrixT>()(temp) * temp * (ONE/value_type(maxRows-1));
		}
	}
};

template<typename MatrixT>
struct Power
// PURPOSE: Returns a matrix with all its elements raised to a defined power.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: A matrix with all the original elements raised to a defined power will be returned.
{
	typedef typename MatrixT::value_type value_type;

	MatrixT operator()(const MatrixT& matrix, const value_type& power)
	{
		MatrixT powerMatrix;
		(*this)(matrix, powerMatrix, power);
		return powerMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& powerMatrix, const value_type& power)
	{
		const size_t maxRows = matrix.row.size();
		const size_t maxCols = matrix.col.size();

		powerMatrix.resize(maxRows, maxCols);
		for (size_t i=0; i<maxRows; ++i)
		{
			for (size_t j=0; j<maxCols; ++j) powerMatrix(i,j) = pow(matrix(i,j), power);	
		}
	}
};

template<typename MatrixT>
struct Mean
// PURPOSE: Returns a row vector with the mean of each column.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: A 'row vector' matrix with the mean of each column will be returned.
{
	typedef typename MatrixT::value_type value_type;

	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT meanMatrix;
		(*this)(matrix, meanMatrix);
		return meanMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& meanMatrix)
	{
		const size_t maxRows = matrix.row.size();
		Sum<MatrixT>()(matrix, meanMatrix);
		meanMatrix *= value_type(1) / value_type(maxRows);
	}
};

template<typename MatrixT>
struct Median
// PURPOSE: Returns a row vector with the median of each column.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: A 'row vector' matrix with the median of each column will be returned.
{
	typedef typename MatrixT::value_type value_type;

	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT medianMatrix;
		(*this)(matrix, medianMatrix);
		return medianMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& medianMatrix)
	{
		const size_t maxCols = matrix.col.size();
		medianMatrix.resize(1, maxCols);
		for (size_t i=0; i<maxCols; ++i)
		{
			vector<value_type> v(matrix.col(i).begin(), matrix.col(i).end());
			nth_element(v.begin(), v.begin() + v.size()/2, v.end());
			medianMatrix(0,i) = (*(v.begin() + v.size()/2));
		}
	}
};

template<typename MatrixT>
struct Sum
// PURPOSE: Returns a row vector with the sum over each column.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: A 'row vector' matrix with the sum over each column will be returned.
{
	typedef typename MatrixT::value_type value_type;

	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT sumMatrix;
		(*this)(matrix, sumMatrix);
		return sumMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& sumMatrix)
	{
		const size_t maxRows = matrix.row.size();
		const size_t maxCols = matrix.col.size();

		sumMatrix.resize(1, maxCols, value_type(0));
		for (size_t j=0; j<maxCols; ++j)
		{
			for (size_t i=0; i<maxRows; ++i) sumMatrix(0,j) += matrix(i,j);
		}
	}
};

template<typename MatrixT>
struct CumulativeSum
// PURPOSE: Returns the cumulative sum of elements.
// REQUIRE: A matrix conforming to the template class Matrix<ValueT,MatrixTypeP,MathP>.
// PROMISE: A matrix of the cumulative sum of elements will be returned.
{
	typedef typename MatrixT::value_type value_type;

	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT cumsumMatrix;
		(*this)(matrix, cumsumMatrix);
		return cumsumMatrix;
	}
	void operator()(const MatrixT& matrix, MatrixT& cumsumMatrix)
	{
		const size_t maxRows = matrix.row.size();
		const size_t maxCols = matrix.col.size();

		cumsumMatrix.resize(maxRows, maxCols, value_type(0));
		for (size_t j=0; j<maxCols; ++j)
		{
			cumsumMatrix(0,j) = matrix(0,j);
			for (size_t i=1; i<maxRows; ++i) cumsumMatrix(i,j) = cumsumMatrix(i-1,j) + matrix(i,j);
		}
	}
};

template<typename MatrixT>
struct Inverse
{
	typedef typename MatrixT::value_type value_type;

	MatrixT operator()(const MatrixT& matrix)
	{
		MatrixT InvertMatrix;
		(*this)(matrix, InvertMatrix);
		return InvertMatrix;
	}

	void operator()(const MatrixT& matrix, MatrixT& inverseMatrix)
	{
		const int maxRows = matrix.row.size();
		const int maxCols = matrix.col.size();
		inverseMatrix = matrix;
		if (maxRows == maxCols && maxRows >= 2)
		{
			value_type temp;
			
			value_type **E = new value_type *[maxRows];
			for (int i = 0; i < maxRows; ++i)
				E[i] = new value_type[maxRows];

			for (int i = 0; i < maxRows; i++)
			{
				for (int j = 0; j < maxRows; j++)
				{
					E[i][j] = 0.0;
					if (i == j)
						E[i][j] = 1.0;
				}
			}

			for (int k = 0; k < maxRows; k++)
			{
				temp = inverseMatrix(k,k);
				for (int j = 0; j < maxRows; j++)
				{
					inverseMatrix(k,j) /= temp;
					E[k][j] /= temp;
				}
				for (int i = k + 1; i < maxRows; i++)
				{
					temp = inverseMatrix(i,k);
					for (int j = 0; j < maxRows; j++)
					{
						inverseMatrix(i,j) -= inverseMatrix(k,j) * temp;
						E[i][j] -= E[k][j] * temp;
					}
				}
			}

			for (int k = maxRows - 1; k > 0; k--)
			{
				for (int i = k - 1; i >= 0; i--)
				{
					temp = inverseMatrix(i,k);
					for (int j = 0; j < maxRows; j++)
					{
						inverseMatrix(i,j) -= inverseMatrix(k,j) * temp;
						E[i][j] -= E[k][j] * temp;
					}
				}
			}

			for (int i = 0; i < maxRows; i++)
			{
				for (int j = 0; j < maxRows; j++)
					inverseMatrix(i,j) = E[i][j];
			}				
			for (int i = 0; i < maxRows; i++)
				delete [] E[i];
			delete [] E;
		}
		else
			inverseMatrix.empty();
	}
};

template<typename MatrixT>
struct Det
{
	
	double operator()(const MatrixT& matrix)
	{
		double det;
		(*this)(matrix,det);
		return det;
	}

	void operator()(const MatrixT& matrix, double& det)
	{
		const int maxRows = matrix.row.size();
		const int maxCols = matrix.col.size();
		MatrixT mtr = matrix;
		if (maxRows == maxCols)
		{
			for (int step = 0; step < maxRows - 1; ++step)
			{
				for (int row = step + 1; row < maxRows; ++row)
				{
					double coeff = -mtr(row,step)/mtr(step,step);
					for (int col = step; col < maxCols; col++)
						mtr(row,col) += mtr(step,col)*coeff;
				}
			}
			det = 1;
			for (int i = 0; i < maxRows; ++i)
			{
				det *= mtr(i,i);
			}
		}
		else
			throw logic_error("Number of cols != number of rows");
	}
	
};


} // Namespace YMatrix end.

#endif

