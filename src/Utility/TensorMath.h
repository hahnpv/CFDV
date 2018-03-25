#pragma once 

	///////////////////////////////////////////////////////////////////////////
	/// TENSOR MATHEMATICS													///
	///////////////////////////////////////////////////////////////////////////


	// Addition of two tensors - so long as indices align, and dimensions match,
	// it does not matter, we just add the two data structures.
template<class T, int ndim>
Tensor<T, ndim> operator+( Tensor<T, ndim> & left, Tensor<T, ndim> & right)
{
//	Tensor<T, ndim> result = left;
	Tensor<T, ndim> result(left.imax(), left.jmax(), left.kmax(), left.lmax());

	for (int i = 0; i < result.ref->size; i++)
	{
		result.ref->data[i] = left.ref->data[i] + right.ref->data[i];
	}

	return result;
}

// Subtraction of two tensors - so long as indices align, and dimensions match,
// it does not matter, we just add the two data structures.
template<class T, int ndim>
Tensor<T, ndim> operator-(Tensor<T, ndim> & left, Tensor<T, ndim> & right)
{
	//	Tensor<T, ndim> result = left;
	Tensor<T, ndim> result(left.imax(), left.jmax(), left.kmax(), left.lmax());

	for (int i = 0; i < result.ref->size; i++)
	{
		result.ref->data[i] = left.ref->data[i] - right.ref->data[i];
	}

	return result;
}

	// Multiplication of a scalar
template<class T, int ndim>
Tensor<T, ndim> operator*(Tensor<T, ndim> & left, T c)		// T * c
{
	for (int i = 0; i < left.ref->size; i++)
	{
		left.ref->data[i] *= c;
	}	

	return left;
}

template<class T, int ndim>
Tensor<T, ndim> operator*(T c, Tensor<T, ndim> & left)		// necessary for c*T
{
	return left * c;
}


	// -= operator
template<class T, int ndim>
Tensor<T, ndim> & operator-=(Tensor<T, ndim> & lhs, Tensor<T, ndim> & rhs)
{
	lhs.ref = rhs.ref;
	for (int i = 0; i < lhs.ref->size; i++)
	{
		lhs.ref->data[i] *= -1.;
	}
	return lhs;
}


//   d(r, s) * (flow.F(i, s) + flow.G(i, s))
// = d(r, s) * FG(i, s)
// = x(r, i)
// check continuum notes.
// need the concept of indices at this point ... 
// methodology, I think, would be to loop over indices to find the matching indice. 
// then set up multiplication accordingly.