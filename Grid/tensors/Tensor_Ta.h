/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/tensors/Tensor_Ta.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef GRID_MATH_TA_H
#define GRID_MATH_TA_H


NAMESPACE_BEGIN(Grid);

/////////////////////////////////////////////// 
// Ta function for scalar, vector, matrix
/////////////////////////////////////////////// 
/*
  accelerator_inline ComplexF Ta( const ComplexF &arg){    return arg;}
  accelerator_inline ComplexD Ta( const ComplexD &arg){    return arg;}
  accelerator_inline RealF Ta( const RealF &arg){    return arg;}
  accelerator_inline RealD Ta( const RealD &arg){    return arg;}
*/

template<class vtype> accelerator_inline iScalar<vtype> Ta(const iScalar<vtype>&r)
{
  iScalar<vtype> ret;
  ret._internal = Ta(r._internal);
  return ret;
}
template<class vtype,int N> accelerator_inline iVector<vtype,N> Ta(const iVector<vtype,N>&r)
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = Ta(r._internal[i]);
  }
  return ret;
}
template<class vtype,int N> accelerator_inline iMatrix<vtype,N> Ta(const iMatrix<vtype,N> &arg)
{
  iMatrix<vtype,N> ret;

  double factor = (1.0/(double)N);
  ret= (arg - adj(arg))*0.5;
  ret=ret - (trace(ret)*factor);
  return ret;
}

template<class vtype> accelerator_inline iScalar<vtype> SpTa(const iScalar<vtype>&r)
{
  iScalar<vtype> ret;
  ret._internal = SpTa(r._internal);
  return ret;
}
template<class vtype,int N> accelerator_inline iVector<vtype,N> SpTa(const iVector<vtype,N>&r)
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = SpTa(r._internal[i]);
  }
  return ret;
}
template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr>
accelerator_inline iMatrix<vtype,N> SpTa(const iMatrix<vtype,N> &arg)
{
  // Generalises Ta to Sp2n
  // Applies the following projections
  // P_{antihermitian} P_{antihermitian-Sp-algebra} P_{traceless}
  // where the ordering matters
  // P_{traceless} subtracts the trace
  // P_{antihermitian-Sp-algebra} provides the block structure of the algebra based on U = exp(T) i.e. anti-hermitian generators
  // P_{antihermitian} does in-adj(in) / 2
  iMatrix<vtype,N> ret(arg);
  double factor = (1.0/(double)N);
  vtype nrm;
  nrm = 0.5;
    
  ret = arg - (trace(arg)*factor);
    
  for(int c1=0;c1<N/2;c1++)
  {
      for(int c2=0;c2<N/2;c2++)
      {
          ret._internal[c1][c2] = nrm*(conjugate(ret._internal[c1+N/2][c2+N/2]) + ret._internal[c1][c2]); // new[up-left] = old[up-left]+old*[down-right]
          ret._internal[c1][c2+N/2] = nrm*(ret._internal[c1][c2+N/2] - conjugate(ret._internal[c1+N/2][c2])); // new[up-right] = old[up-right]-old*[down-left]
      }
      for(int c2=N/2;c2<N;c2++)
      {
          ret._internal[c1+N/2][c2-N/2] = -conjugate(ret._internal[c1][c2]);  //  reconstructs lower blocks
          ret._internal[c1+N/2][c2] = conjugate(ret._internal[c1][c2-N/2]);   //  from upper blocks
      }
  }
    
  ret = (ret - adj(ret))*0.5;

  return ret;
}

/////////////////////////////////////////////// 
// ProjectOnGroup function for scalar, vector, matrix 
// Projects on orthogonal, unitary group
/////////////////////////////////////////////// 

template<class vtype> accelerator_inline iScalar<vtype> ProjectOnGroup(const iScalar<vtype>&r)
{
  iScalar<vtype> ret;
  ret._internal = ProjectOnGroup(r._internal);
  return ret;
}
template<class vtype,int N> accelerator_inline iVector<vtype,N> ProjectOnGroup(const iVector<vtype,N>&r)
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = ProjectOnGroup(r._internal[i]);
  }
  return ret;
}
template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr> 
accelerator_inline iMatrix<vtype,N> ProjectOnGroup(const iMatrix<vtype,N> &arg)
{
  // need a check for the group type?
  iMatrix<vtype,N> ret(arg);
  vtype nrm;
  vtype inner;
  for(int c1=0;c1<N;c1++){

    // Normalises row c1
    zeroit(inner);	
    for(int c2=0;c2<N;c2++)
      inner += innerProduct(ret._internal[c1][c2],ret._internal[c1][c2]);

    nrm = sqrt(inner);
    nrm = 1.0/nrm;
    for(int c2=0;c2<N;c2++)
      ret._internal[c1][c2]*= nrm;
      
    // Remove c1 from rows c1+1...N-1
    for (int b=c1+1; b<N; ++b){
      decltype(ret._internal[b][b]*ret._internal[b][b]) pr;
      zeroit(pr);
      for(int c=0; c<N; ++c)
	pr += conjugate(ret._internal[c1][c])*ret._internal[b][c];
	  
      for(int c=0; c<N; ++c){
	ret._internal[b][c] -= pr * ret._internal[c1][c];
      }
    }
  }

  // Normalise last row
  {
    int c1 = N-1;
    zeroit(inner);	
    for(int c2=0;c2<N;c2++)
      inner += innerProduct(ret._internal[c1][c2],ret._internal[c1][c2]);

    nrm = sqrt(inner);
    nrm = 1.0/nrm;
    for(int c2=0;c2<N;c2++)
      ret._internal[c1][c2]*= nrm;
  }
  // assuming the determinant is ok
  return ret;
}

// re-do for sp2n

// Ta cannot be defined here for Sp2n because I need the generators from the Sp class
// It is defined in gauge impl types

template<class vtype> accelerator_inline iScalar<vtype> ProjectOnSpGroup(const iScalar<vtype>&r)
{
  iScalar<vtype> ret;
  ret._internal = ProjectOnSpGroup(r._internal);
  return ret;
}
template<class vtype,int N> accelerator_inline iVector<vtype,N> ProjectOnSpGroup(const iVector<vtype,N>&r)
{
  iVector<vtype,N> ret;
  for(int i=0;i<N;i++){
    ret._internal[i] = ProjectOnSpGroup(r._internal[i]);
  }
  return ret;
}


// int N is 2n in Sp(2n)
template<class vtype,int N, typename std::enable_if< GridTypeMapper<vtype>::TensorLevel == 0 >::type * =nullptr>
accelerator_inline iMatrix<vtype,N> ProjectOnSpGroup(const iMatrix<vtype,N> &arg)
{
  // need a check for the group type?
  iMatrix<vtype,N> ret(arg);
  vtype nrm;
  vtype inner;
  
  for(int c1=0;c1<N/2;c1++)
  {
      
    for (int b=0; b<c1; b++)                  // remove the b-rows from U_c1
    {
      decltype(ret._internal[b][b]*ret._internal[b][b]) pr;
      decltype(ret._internal[b][b]*ret._internal[b][b]) prn;
      zeroit(pr);
      zeroit(prn);
          
      for(int c=0; c<N; c++)
      {
        pr += conjugate(ret._internal[c1][c])*ret._internal[b][c];        // <U_c1 | U_b >
        prn += conjugate(ret._internal[c1][c])*ret._internal[b+N/2][c];   // <U_c1 | U_{b+N} >
      }
       

      for(int c=0; c<N; c++)
      {
        ret._internal[c1][c] -= (conjugate(pr) * ret._internal[b][c] + conjugate(prn) * ret._internal[b+N/2][c] );    //  U_c1 -= (  <U_c1 | U_b > U_b + <U_c1 | U_{b+N} > U_{b+N}  )
      }
    }
    
    zeroit(inner);
    for(int c2=0;c2<N;c2++)
    {
      inner += innerProduct(ret._internal[c1][c2],ret._internal[c1][c2]);
    }
      
    nrm = sqrt(inner);
    nrm = 1.0/nrm;
    for(int c2=0;c2<N;c2++)
    {
      ret._internal[c1][c2]*= nrm;
    }
      
    for(int c2=0;c2<N/2;c2++)
    {
      ret._internal[c1+N/2][c2+N/2] = conjugate(ret._internal[c1][c2]);          // down right in the new matrix = (up-left)* of the old matrix
    }
      
    for(int c2=N/2;c2<N;c2++)
    {
      ret._internal[c1+N/2][c2-N/2] = -conjugate(ret._internal[c1][c2]);;     // down left in the new matrix = -(up-right)* of the old
    }
  }
  return ret;
}

NAMESPACE_END(Grid);

#endif
