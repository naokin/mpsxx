/**
 * \file  MPSblas.h
 * 
 * Here functions are implemented that perform basic multi-linear algebra computations matrix product states and matrix product operators.
 * The naming of the blas library has been taken over for many functions, i.e. matrix-vector gemv is MPO/MPS contraction.
 * Standard operators *, + , - are overloaded, --> MPO * MPO, MPO * MPS.
 *
 * @mainpage 
 * \section CONTENTS
 *
 *  - MPS/MPO blas equivalent for convenient MPS/MPO algebra.
 *  
 *  - MPO generator for Hubbard and Molecular Hamiltonians
 *
 *  - Simple DMRG code implemented in terms of MPS/MPO language
 *
 * \section DEPENDENCY
 *
 *  - Boost Library (http://www.boost.org/)
 *
 *  - BTAS Library (https://github.com/naokin/btas.git)
 *
 *  - GNU GCC 4.7.0 or Later / Intel C/C++ Compiler 13.0 or Later
 *
 *  - CBLAS & LAPACK Libraries / Intel MKL Library
 */
#ifndef _BTAS_MPSBLAS_H
#define _BTAS_MPSBLAS_H 1

#include <iostream>
#include <iomanip>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

using std::cout;
using std::endl;
using std::ostream;

#include <btas/blas_cxx_interface.h>

#include <btas/TVector.h>

#include <btas/DENSE/DArray.h>
#include <btas/QSPARSE/QSDArray.h>
#include <btas/QSPARSE/QSDcontract.h>
#include "btas/QSPARSE/QSDArray.h"

using namespace btas;

namespace mpsxx {

   enum MPS_DIRECTION {

      Left,//!left to right
      Right //!right to left

   };

   enum MPS_TYPE { WAVEFUNCTION, LEFTCANONICAL, RIGHTCANONICAL };

   //!typedefine MPX as a std::vector< QSDArray<N> >  for functions the same for both MPS and MPO
   template<size_t N,class Q>
      using MPX = std::vector< QSDArray<N,Q> >;

   //!typedefine MPS as an std::vector< QSDArray<3> > 
   template<class Q>
      using MPS = std::vector< QSDArray<3,Q> >;

   //!typedefine MPO as an std::vector< QSDArray<4> > 
   template<class Q>
      using MPO = std::vector< QSDArray<4,Q> >;

   /**
    * given the length, total quantumnumber and physical quantumnumbers, and some cutoff block dimension
    * @param L length of the chain
    * @param qt total quantumnumber
    * @param qp array of physical quantumnumbers on the sites
    * @param qr std::vector of Qshapes length L, will contain the right bond quantumnumbers on exit, input will be destroyed
    * @param dr std::vector of Dshapes length L, will contain the right bond dimensions corresponding to the numbers on exit, input will be destroyed
    * @param D cutoff block dimension (max dimension of each quantumsector)
    */
   template<class Q>
      void calc_qdim(int L,const Q &qt,const Qshapes<Q> &qp,std::vector< Qshapes<Q> > &qr,std::vector<Dshapes> &dr,int D){

         //shape of the physical index
         Dshapes dp(qp.size(),1);

         qr[0] = qp;
         dr[0] = dp;

         for(int i = 1;i < L - 1;++i){

            qr[i] = qr[i-1] * qp;

            for(unsigned int j = 0;j < dr[i - 1].size();++j)
               for(unsigned int k = 0;k < dp.size();++k)
                  dr[i].push_back(dr[i-1][j]*dp[k]);

            //sort the list of quantumnumbers
            Q srtq;
            int srtd;

            for(int j = 0;j < qr[i].size();++j){

               for(int k = j + 1;k < qr[i].size();++k){

                  if(qr[i][k] < qr[i][j]){

                     srtq = qr[i][j];
                     srtd = dr[i][j];

                     qr[i][j] = qr[i][k];
                     dr[i][j] = dr[i][k];

                     qr[i][k] = srtq;
                     dr[i][k] = srtd;

                  }

               }

            }

            //remove quantumnumbers that occur multiple times
            int j = 0;

            while(j < qr[i].size()){

               int k = j + 1;

               while(k < qr[i].size()){

                  //this removes redundant
                  if( qr[i][k] == qr[i][j] ){

                     //erase the redundant quantumnumber
                     qr[i].erase(qr[i].begin() + k);

                     //add the dimension to the right block
                     dr[i][j] += dr[i][k];

                     //erase the redundant dimension
                     dr[i].erase(dr[i].begin() + k);

                  }
                  else
                     ++k;

                  //if dimension is too large, set to D
                  if(dr[i][j] > D)
                     dr[i][j] = D;

               }

               ++j;

            }

         }

         qr[L-1] = Qshapes<Q>(1,qt);
         dr[L-1] = Dshapes(1,1);

         Qshapes<Q> tmpq;
         Dshapes tmpd;

         for(int i = L - 2;i >= 0;--i){

            tmpq.clear();

            for(int j = 0;j < qr[i+1].size();++j)
               for(int k = qp.size() - 1;k >= 0;--k)
                  tmpq.push_back(qr[i + 1][j] * (-qp[k]));

            tmpd.clear();

            for(int j = 0;j < dr[i+1].size();++j)
               for(int k = dp.size() - 1;k >= 0;--k)
                  tmpd.push_back(dr[i+1][j]*dp[k]);

            //sort the list of temporary quantumnumbers
            Q srtq;
            int srtd;

            for(int j = 0;j < tmpq.size();++j){

               for(int k = j + 1;k < tmpq.size();++k){

                  if(tmpq[k] < tmpq[j]){

                     srtq = tmpq[j];
                     srtd = tmpd[j];

                     tmpq[j] = tmpq[k];
                     tmpd[j] = tmpd[k];

                     tmpq[k] = srtq;
                     tmpd[k] = srtd;

                  }

               }

            }

            int j = 0;

            while(j < tmpq.size()){

               int k = j + 1;

               while(k < tmpq.size()){

                  //this removes redundant
                  if( tmpq[k] == tmpq[j] ){

                     //erase the redundant quantumnumber
                     tmpq.erase(tmpq.begin() + k);

                     //add the dimension to the right block
                     tmpd[j] += tmpd[k];

                     //erase the redundant dimension
                     tmpd.erase(tmpd.begin() + k);

                  }
                  else
                     ++k;

                  //if dimension is too large, set to D
                  if(tmpd[j] > D)
                     tmpd[j] = D;

               }

               ++j;

            }

            //remove irrelevant quantum blocks from below: i.e. which are not present in both tmpq and qr[i]
            for(int j = 0;j < qr[i].size();++j){

               int flag = 0;

               //if its present: set flag to 1
               for(int k = 0;k < tmpq.size();++k)
                  if(qr[i][j] == tmpq[k])
                     flag = 1;

               //if not: erase element
               if(flag == 0){

                  qr[i].erase(qr[i].begin() + j);
                  dr[i].erase(dr[i].begin() + j);

                  --j;

               }

            }

            //now replace the dimensions
            for(unsigned int k = 0;k < qr[i].size();++k){

               //is there a quantumnumber in tmpq equal to qr[i][k]?
               for(unsigned int l = 0;l < tmpq.size();++l){

                  //if there is, take the smallest dimension
                  if(qr[i][k] == tmpq[l]){

                     if(dr[i][k] > tmpd[l])
                        dr[i][k] = tmpd[l];

                  }

               }

            }

         }

      }

   /**
    * create an MPS chain of length L initialized randomly on total Q number qt, with physical quantumnumber qp
    * @param L length of the chain
    * @param qt total quantumnumber
    * @param qp Qshapes object containing the physical quantumnumbers
    * @param D maximal dimension of the quantum blocks
    * @param f_random_generator predefined function which calls a random number generator
    * @return the MPS chain randomly filled and with correct quantumnumbers and dimensions
    */
   template<class Q>
      MPS<Q> create(int L,const Q &qt,const Qshapes<Q> &qp,int D,const function<double(void)>& f_random_generator){ 

         //shape of the physical index
         Dshapes dp(qp.size(),1);

         std::vector< Qshapes<Q> > qr(L);
         std::vector<Dshapes> dr(L);

         calc_qdim(L,qt,qp,qr,dr,D);

         //now allocate the tensors!
         TVector<Qshapes<Q>,3> qshape;
         TVector<Dshapes,3> dshape;

         //first 0
         Qshapes<Q> ql(1,Q::zero());
         Dshapes dl(ql.size(),1);

         qshape = make_array(ql,qp,-qr[0]);
         dshape = make_array(dl,dp,dr[0]);

         //construct an MPS
         MPS<Q> A(L);

         A[0].resize(Q::zero(),qshape,dshape);
         A[0].generate(f_random_generator);

         //then the  middle ones
         for(int i = 1;i < L;++i){

            ql = qr[i - 1];
            dl = dr[i - 1];

            qshape = make_array(ql,qp,-qr[i]);
            dshape = make_array(dl,dp,dr[i]);

            A[i].resize(Q::zero(),qshape,dshape);
            A[i].generate(f_random_generator);

         }

         return A;

      }

   /**
    * create an MPS chain of length L initialized on a contant number on total Q number qt, with physical quantumnumber qp
    * @param L length of the chain
    * @param qt total quantumnumber
    * @param qp Qshapes object containing the physical quantumnumbers
    * @param D maximal dimension of the quantum blocks
    * @param value the number the mps is initialized onto, standard 0
    * @return the MPS chain randomly filled and with correct quantumnumbers and dimensions
    */
   template<class Q>
      MPS<Q> create(int L,const Q &qt,const Qshapes<Q> &qp,int D,double value = 0.0){ 

         //shape of the physical index
         Dshapes dp(qp.size(),1);

         std::vector< Qshapes<Q> > qr(L);
         std::vector<Dshapes> dr(L);

         calc_qdim(L,qt,qp,qr,dr,D);

         //now allocate the tensors!
         TVector<Qshapes<Q>,3> qshape;
         TVector<Dshapes,3> dshape;

         //first 0
         Qshapes<Q> ql(1,Q::zero());
         Dshapes dl(ql.size(),1);

         qshape = make_array(ql,qp,-qr[0]);
         dshape = make_array(dl,dp,dr[0]);

         //construct an MPS
         MPS<Q> A(L);

         A[0].resize(Q::zero(),qshape,dshape);
         A[0] = value;

         //then the  middle ones
         for(int i = 1;i < L;++i){

            ql = qr[i - 1];
            dl = dr[i - 1];

            qshape = make_array(ql,qp,-qr[i]);
            dshape = make_array(dl,dp,dr[i]);

            A[i].resize(Q::zero(),qshape,dshape);
            A[i] = value;

         }

         return A;

      }

   /**
    * @param L length of the chain
    * @param qp physical quantumnumbers
    * @param occ std::vector of length L ints containing the local physical quantumnumber on every site in the product state
    * @return create an product state chain of length L with physical indices qp and
    */
   template<class Q>
      MPS<Q> product_state(int L,const Qshapes<Q> &qp,const std::vector<int> &occ){ 

         //shape of the physical index
         Dshapes dp(qp.size(),1);

         //now allocate the tensors!
         TVector<Qshapes<Q>,3> qshape;
         TVector<Dshapes,3> dshape;

         MPS<Q> A(L);

         Qshapes<Q> qz;
         qz.push_back(Q::zero());

         Qshapes<Q> qr;
         qr.push_back(qp[occ[0]]);

         Dshapes dz;
         dz.push_back(1);

         qshape = make_array(qz,qp,-qr);
         dshape = make_array(dz,dp,dz);

         A[0].resize(Q::zero(),qshape,dshape,1.0);

         Q tmpq;

         for(int i = 1;i < L;++i){

            tmpq = qr[0] * qp[occ[i]];

            qr.clear();
            qr.push_back(tmpq);

            qshape = make_array(-A[i - 1].qshape(2),qp,-qr);
            dshape = make_array(dz,dp,dz);

            A[i].resize(Q::zero(),qshape,dshape,1.0);

         }

         return A;

      }

   /**
    * scale the MPX with a constant factor
    * @param alpha scalingfactor
    * @param mpx the MPX to be scaled
    */
   template<size_t N,class Q>
      void scal(double alpha,MPX<N,Q> &mpx){

         int L = mpx.size();

         int sign;

         if(alpha > 0)
            sign = 1;
         else
            sign = -1;

         alpha = pow(fabs(alpha),1.0/(double)L);

         QSDscal(sign * alpha,mpx[0]);

         for(int i = 1;i < mpx.size();++i)
            QSDscal(alpha,mpx[i]);

      }

   /**
    * MPS/O equivalent of the axpy blas function: Y <- alpha X + Y
    * taking the direct sum of the individual tensors in the chain
    * @param alpha double scaling factor
    * @param X input MPX
    * @param Y output MPX: alpha * X will be added to the input Y and put in output Y
    */
   template<size_t N,class Q>
      void axpy(double alpha,const MPX<N,Q> &X,MPX<N,Q> &Y){

         //first check if we can sum these two:
         if(X.size() != Y.size())
            BTAS_THROW(false, "Error: input MP objects do not have the same length!");

         int L = X.size();

         QSDArray<N> tmp1;
         QSDArray<N> tmp2;

         IVector<N-1> left;

         for(int i = 0;i < N-1;++i)
            left[i] = i;

         //first left: scale the B term
         QSDscal(1.0/alpha,Y[0]);

         QSDdsum(X[0],Y[0],left,tmp1);

         //merge the column quantumnumbers together
         TVector<Qshapes<Q>,1> qmerge;
         TVector<Dshapes,1> dmerge;

         qmerge[0] = tmp1.qshape(N-1);
         dmerge[0] = tmp1.dshape(N-1);

         QSTmergeInfo<1> info(qmerge,dmerge);

         //then merge
         Y[0].clear();
         QSTmerge(tmp1,info,Y[0]);

         //rescale again
         QSDscal(alpha,Y[0]);

         IVector<N-2> middle;

         for(int i = 1;i < N-1;++i)
            middle[i - 1] = i;

         //row and column addition in the middle of the chain
         for(int i = 1;i < L - 1;++i){

            tmp1.clear();
            QSDdsum(X[i],Y[i],middle,tmp1);

            //merge the row quantumnumbers together
            qmerge[0] = tmp1.qshape(0);
            dmerge[0] = tmp1.dshape(0);

            info.reset(qmerge,dmerge);

            //then merge
            tmp2.clear();
            QSTmerge(info,tmp1,tmp2);

            //column quantumnumbers
            qmerge[0] = tmp2.qshape(N-1);
            dmerge[0] = tmp2.dshape(N-1);

            info.reset(qmerge,dmerge);

            //then merge
            Y[i].clear();
            QSTmerge(tmp2,info,Y[i]);

         }

         IVector<N-1> right;

         for(int i = 0;i < N-1;++i)
            right[i] = i + 1;

         //finally the right
         tmp1.clear();
         QSDdsum(X[L-1],Y[L-1],right,tmp1);

         //merge the row quantumnumbers together
         qmerge[0] = tmp1.qshape(0);
         dmerge[0] = tmp1.dshape(0);

         info.reset(qmerge,dmerge);

         //then merge
         Y[L - 1].clear();
         QSTmerge(info,tmp1,Y[L-1]);

      }

   /**
    * construct new MPX AB that is the sum of A + B: this is done by making a larger MPX object with larger bond dimension,
    * taking the direct sum of the individual tensors in the chain
    * @param A input MPX
    * @param B input MPX
    * @return the MPX result
    */
   template<size_t N,class Q>
      MPX<N,Q> operator+(const MPX<N,Q> &A,const MPX<N,Q> &B){

         MPX<N,Q> AB(B);
         axpy(1.0,A,AB);

         return AB;

      }

   /**
    * construct new MPX AB that is the difference of A and B: A - B this is done by making an MPX object with larger bond dimension,
    * taking the direct sum of the individual tensors in the chain
    * @param A input MPX
    * @param B input MPX
    * @return the MPX result
    */
   template<size_t N,class Q>
     MPX<N,Q> operator-(const MPX<N,Q> &A,const MPX<N,Q> &B){

        MPX<N,Q> AB(A);
        axpy(-1.0,B,AB);

        return AB;

     }

   /**
    * Compress an MP object by performing an SVD
    * @param mpx is the input MPX, will be lost/overwritten by the compressed MPX
    * @param dir direction of the canonicalization, from left to right if Left, right to left if Right
    * @param D if > 0   this specifies the number of states to be kept
    *          if == 0  all the states are kept
    *          if < 0 all singular values > 10^-D are kept
    * @return the total discarded weight
    */
   template<size_t N,class Q>
      double compress(MPX<N,Q> &mpx,const MPS_DIRECTION &dir,int D){

         int L = mpx.size();//length of the chain

         double dweight = 0.0;

         if(dir == Left) {

            SDArray<1> S;//singular values
            QSDArray<2> V;//V^T
            QSDArray<N> U;//U --> unitary left normalized matrix

            for(int i = 0;i < L - 1;++i){

               //redistribute the norm over the chain: for stability reasons
               double nrm = sqrt(QSDdotc(mpx[i],mpx[i]));

               QSDscal(1.0/nrm,mpx[i]);

               scal(nrm,mpx);

               //then svd
               dweight += QSDgesvd(RightArrow,mpx[i],S,U,V,D);

               //copy unitary to mpx
               QSDcopy(U,mpx[i]);

               //paste S and V together
               SDdidm(S,V);

               //and multiply with mpx on the next site
               U = mpx[i + 1];

               //when compressing dimensions will change, so reset:
               mpx[i + 1].clear();

               QSDcontract(1.0,V,shape(1),U,shape(0),0.0,mpx[i + 1]);

            }

            //redistribute the norm over the chain
            double nrm = sqrt(QSDdotc(mpx[L-1],mpx[L-1]));

            QSDscal(1.0/nrm,mpx[L-1]);

            scal(nrm,mpx);

         }
         else{//right

            SDArray<1> S;//singular values
            QSDArray<N> V;//V^T --> unitary right normalized matrix
            QSDArray<2> U;//U

            for(int i = L - 1;i > 0;--i){

               //redistribute the norm over the chain: for stability reasons
               double nrm = sqrt(QSDdotc(mpx[i],mpx[i]));

               QSDscal(1.0/nrm,mpx[i]);

               scal(nrm,mpx);

               //then SVD: 
               dweight += QSDgesvd(RightArrow,mpx[i],S,U,V,D);

               //copy unitary to mpx
               QSDcopy(V,mpx[i]);

               //paste U and S together
               SDdimd(U,S);

               //and multiply with mpx on the next site
               V = mpx[i - 1];

               //when compressing dimensions will change, so reset:
               mpx[i - 1].clear();

               QSDcontract(1.0,V,shape(N-1),U,shape(0),0.0,mpx[i - 1]);

            }

            double nrm = sqrt(QSDdotc(mpx[0],mpx[0]));

            QSDscal(1.0/nrm,mpx[0]);

            scal(nrm,mpx);

         }

         return dweight;

      }

   /**
    * clean up the MPX, i.e. make sure the right quantumblocks are connected, remove unnecessary quantumnumbers and blocks
    * @param mpx input MPX, will be changed 'cleaned' on exit
    */
   template<size_t N,class Q>
      void clean(MPX<N,Q> &mpx){

         Dshapes dr;

         //from left to right
         for(int i = 0;i < mpx.size() - 1;++i){

            dr = mpx[i].dshape()[N - 1];

            std::vector<Q> qrem;

            for(int j = 0;j < dr.size();++j)
               if(dr[j] == 0)
                  qrem.push_back(mpx[i].qshape()[N - 1][j]);//what is the quantumnumber with 0 dimension?

            if(qrem.size() != 0){

               //remove the zero blocks from site i
               for(int j = 0;j < qrem.size();++j){

                  //find the index corresponding to quantumnumber qrem[j]
                  Qshapes<Q> qr = mpx[i].qshape()[N - 1];

                  for(int k = 0;k < qr.size();++k)
                     if(qr[k] == qrem[j])
                        mpx[i].erase(N - 1,k);

               }

               for(int j = 0;j < qrem.size();++j){

                  //remove the corresponding blocks on the 0 leg of the next site
                  Qshapes<Q> ql = mpx[i + 1].qshape()[0];

                  for(int k = 0;k < ql.size();++k)
                     if(ql[k] == -qrem[j])
                        mpx[i + 1].erase(0,k);

               }

            }

         }

         //and back from right to left
         for(int i = mpx.size() - 1;i > 0;--i){

            dr = mpx[i].dshape()[0];//actually dl now

            std::vector<Q> qrem;

            for(int j = 0;j < dr.size();++j)
               if(dr[j] == 0)
                  qrem.push_back(mpx[i].qshape()[0][j]);//what is the quantumnumber with 0 dimension?

            if(qrem.size() != 0){

               //remove the zero blocks from site i
               for(int j = 0;j < qrem.size();++j){

                  //find the index corresponding to quantumnumber qrem[j]
                  Qshapes<Q> qr = mpx[i].qshape()[0];

                  for(int k = 0;k < qr.size();++k)
                     if(qr[k] == qrem[j])
                        mpx[i].erase(0,k);

               }

               for(int j = 0;j < qrem.size();++j){

                  //remove the corresponding blocks on the (nlegs-1) leg of the previous site
                  Qshapes<Q> ql = mpx[i - 1].qshape()[N-1];

                  for(int k = 0;k < ql.size();++k)
                     if(ql[k] == -qrem[j])
                        mpx[i - 1].erase(N-1,k);

               }

            }

         }

      }

   /**
    * @param dir go from left to right (Left) or right to left (Right) for contraction
    * @param A input MPS
    * @param O input MPO
    * @param B input MPS
    * @return the number containing < A | O | B >
    */
   template<class Q>
      double inprod(const MPS_DIRECTION &dir,const MPS<Q> &A,const MPO<Q> &O,const MPS<Q> &B){

         //first check if we can sum these two:
         if(A.size() != B.size() || A.size() != O.size())
            BTAS_THROW(false, "Error: input objects do not have the same length!");

         int L = A.size();

         if(dir == Left){

            enum {j,k,l,m,n,o};

            //from left to right
            QSDArray<5> loc;

            QSDindexed_contract(1.0,O[0],shape(m,n,k,o),A[0],shape(j,k,l),0.0,loc,shape(j,m,n,l,o));

            //merge 2 rows together
            TVector<Qshapes<Q>,2> qmerge;
            TVector<Dshapes,2> dmerge;

            for(int i = 0;i < 2;++i){

               qmerge[i] = loc.qshape(i);
               dmerge[i] = loc.dshape(i);

            }

            QSTmergeInfo<2> info(qmerge,dmerge);

            QSDArray<4> tmp;
            QSTmerge(info,loc,tmp);

            //this will contain the right going part
            QSDArray<3> EO;

            QSDindexed_contract(1.0,B[0].conjugate(),shape(j,k,l),tmp,shape(j,k,m,n),0.0,EO,shape(m,n,l));

            QSDArray<4> I1;
            QSDArray<4> I2;

            for(int i = 1;i < L;++i){

               I1.clear();

               QSDindexed_contract(1.0,EO,shape(j,k,l),A[i],shape(j,m,n),0.0,I1,shape(k,l,n,m));

               I2.clear();

               QSDindexed_contract(1.0,I1,shape(k,l,n,m),O[i],shape(k,j,m,o),0.0,I2,shape(l,j,n,o));

               EO.clear();

               QSDindexed_contract(1.0,I2,shape(l,j,n,o),B[i].conjugate(),shape(l,j,k),0.0,EO,shape(n,o,k));

               //bad style: if no blocks remain, return zero
               if(EO.begin() == EO.end())
                  return 0.0;

            }

            return (*(EO.find(shape(0,0,0))->second))(0,0,0);

         }
         else{

            enum {j,k,l,m,n,o};

            //from right to left
            QSDArray<5> loc;

            QSDindexed_contract(1.0,O[L - 1],shape(j,k,l,m),A[L - 1],shape(o,l,n),0.0,loc,shape(o,j,k,n,m));

            //merge 2 columns together
            TVector<Qshapes<Q>,2> qmerge;
            TVector<Dshapes,2> dmerge;

            for(int i = 0;i < 2;++i){

               qmerge[i] = loc.qshape(3 + i);
               dmerge[i] = loc.dshape(3 + i);

            }

            QSTmergeInfo<2> info(qmerge,dmerge);

            QSDArray<4> tmp;
            QSTmerge(loc,info,tmp);

            //this will contain the left going part
            QSDArray<3> EO;
            QSDindexed_contract(1.0,tmp,shape(j,k,l,m),B[L-1].conjugate(),shape(n,l,m),0.0,EO,shape(j,k,n));

            QSDArray<4> I1;
            QSDArray<4> I2;

            for(int i = L - 2;i >= 0;--i){

               I1.clear();

               QSDindexed_contract(1.0,A[i],shape(j,k,l),EO,shape(l,m,n),0.0,I1,shape(j,k,m,n));

               I2.clear();

               QSDindexed_contract(1.0,O[i],shape(l,o,k,m),I1,shape(j,k,m,n),0.0,I2,shape(j,l,o,n));

               EO.clear();

               QSDindexed_contract(1.0,B[i].conjugate(),shape(k,o,n),I2,shape(j,l,o,n),0.0,EO,shape(j,l,k));

               //bad style: if no blocks remain, return zero
               if(EO.begin() == EO.end())
                  return 0.0;

            }

            return (*(EO.find(shape(0,0,0))->second))(0,0,0);

         }

      }

   /**
    * MPO/S equivalent of the blas gemv function: Y <- alpha * A X + beta Y
    * @param alpha scaling factor of the input MPO
    * @param A input MPO
    * @param X input MPS
    * @param beta scaling factor of the output MPS
    * @param Y output MPS, its content will change on exit.
    */
   template<class Q>
      void gemv(double alpha,const MPO<Q> &A,const MPS<Q> &X,double beta,MPS<Q> &Y){

         //first check if length is the same
         if(A.size() != X.size())
            BTAS_THROW(false, "Error: input objects do not have the same length!");

         if(fabs(beta) < 1.0e-15){

            int L = A.size();

            Y.resize(L);

            enum {j,k,l,m,n,o};

            QSDArray<5> tmp;
            QSDArray<4> mrows;

            for(int i = 0;i < L;++i){

               //clear the tmp object first
               tmp.clear();

               QSDindexed_contract(1.0,A[i],shape(j,k,l,m),X[i],shape(n,l,o),0.0,tmp,shape(n,j,k,o,m));

               //merge 2 rows together
               TVector<Qshapes<Q>,2> qmerge;
               TVector<Dshapes,2> dmerge;

               for(int r = 0;r < 2;++r){

                  qmerge[r] = tmp.qshape(r);
                  dmerge[r] = tmp.dshape(r);

               }

               QSTmergeInfo<2> info(qmerge,dmerge);

               //clear the mrows object first
               mrows.clear();

               //then merge
               QSTmerge(info,tmp,mrows);

               //merge 2 columns together
               for(int r = 2;r < 4;++r){

                  qmerge[r - 2] = mrows.qshape(r);
                  dmerge[r - 2] = mrows.dshape(r);

               }

               info.reset(qmerge,dmerge);

               QSTmerge(mrows,info,Y[i]);

            }

            if( fabs(alpha - 1.0) > 1.0e-15)
               scal(alpha,Y);

         }
         else{//beta != 0.0:

            int L = A.size();

            //first check if we can sum these two:
            if(L != Y.size())
               BTAS_THROW(false, "Error: input objects do not have the same length!");

            scal(beta/alpha,Y);

            enum {j,k,l,m,n,o};

            QSDArray<5> tmp;
            QSDArray<4> mrows;

            QSDindexed_contract(1.0,A[0],shape(j,k,l,m),X[0],shape(n,l,o),0.0,tmp,shape(n,j,k,o,m));

            //merge 2 rows together
            TVector<Qshapes<Q>,2> qmerge1;
            TVector<Dshapes,2> dmerge1;

            for(int r = 0;r < 2;++r){

               qmerge1[r] = tmp.qshape(r);
               dmerge1[r] = tmp.dshape(r);

            }

            QSTmergeInfo<2> info1(qmerge1,dmerge1);

            //clear the mrows object first
            mrows.clear();

            //then merge
            QSTmerge(info1,tmp,mrows);

            //merge 2 columns together
            for(int r = 2;r < 4;++r){

               qmerge1[r - 2] = mrows.qshape(r);
               dmerge1[r - 2] = mrows.dshape(r);

            }

            info1.reset(qmerge1,dmerge1);

            QSDArray<3> Ax;
            QSTmerge(mrows,info1,Ax);

            IVector<2> left;

            for(int i = 0;i < 2;++i)
               left[i] = i;

            QSDArray<3> tmp1;
            QSDArray<3> tmp2;
            QSDdsum(Ax,Y[0],left,tmp1);

            //merge the column quantumnumbers together
            TVector<Qshapes<Q>,1> qmerge2;
            TVector<Dshapes,1> dmerge2;

            qmerge2[0] = tmp1.qshape(2);
            dmerge2[0] = tmp1.dshape(2);

            QSTmergeInfo<1> info2(qmerge2,dmerge2);

            //then merge
            Y[0].clear();
            QSTmerge(tmp1,info2,Y[0]);

            IVector<1> middle;
            middle[0] = 1;

            for(int i = 1;i < L - 1;++i){

               //clear the tmp object first
               tmp.clear();

               QSDindexed_contract(1.0,A[i],shape(j,k,l,m),X[i],shape(n,l,o),0.0,tmp,shape(n,j,k,o,m));

               //merge 2 rows together
               for(int r = 0;r < 2;++r){

                  qmerge1[r] = tmp.qshape(r);
                  dmerge1[r] = tmp.dshape(r);

               }

               info1.reset(qmerge1,dmerge1);

               //clear the mrows object first
               mrows.clear();

               //then merge
               QSTmerge(info1,tmp,mrows);

               //merge 2 columns together
               for(int r = 2;r < 4;++r){

                  qmerge1[r - 2] = mrows.qshape(r);
                  dmerge1[r - 2] = mrows.dshape(r);

               }

               info1.reset(qmerge1,dmerge1);

               //this makes the AX
               QSTmerge(mrows,info1,Ax);

               tmp1.clear();
               QSDdsum(Ax,Y[i],middle,tmp1);

               //merge the row quantumnumbers together
               qmerge2[0] = tmp1.qshape(0);
               dmerge2[0] = tmp1.dshape(0);

               info2.reset(qmerge2,dmerge2);

               //then merge
               tmp2.clear();
               QSTmerge(info2,tmp1,tmp2);

               //column quantumnumbers
               qmerge2[0] = tmp2.qshape(2);
               dmerge2[0] = tmp2.dshape(2);

               info2.reset(qmerge2,dmerge2);

               //then merge
               Y[i].clear();
               QSTmerge(tmp2,info2,Y[i]);

            }

            //last site
            //clear the tmp object first
            tmp.clear();

            QSDindexed_contract(1.0,A[L - 1],shape(j,k,l,m),X[L - 1],shape(n,l,o),0.0,tmp,shape(n,j,k,o,m));

            //merge 2 rows together
            for(int r = 0;r < 2;++r){

               qmerge1[r] = tmp.qshape(r);
               dmerge1[r] = tmp.dshape(r);

            }

            info1.reset(qmerge1,dmerge1);

            //clear the mrows object first
            mrows.clear();

            //then merge
            QSTmerge(info1,tmp,mrows);

            //merge 2 columns together
            for(int r = 2;r < 4;++r){

               qmerge1[r - 2] = mrows.qshape(r);
               dmerge1[r - 2] = mrows.dshape(r);

            }

            info1.reset(qmerge1,dmerge1);

            //this makes the AX
            QSTmerge(mrows,info1,Ax);

            IVector<2> right;

            for(int i = 0;i < 2;++i)
               right[i] = i + 1;

            //finally the right
            tmp1.clear();
            QSDdsum(Ax,Y[L - 1],right,tmp1);

            //merge the row quantumnumbers together
            qmerge2[0] = tmp1.qshape(0);
            dmerge2[0] = tmp1.dshape(0);

            info2.reset(qmerge2,dmerge2);

            //then merge
            QSTmerge(info2,tmp1,Y[L-1]);

            if( fabs(alpha - 1.0) > 1.0e-15)
               scal(alpha,Y);

         }

      }

   /**
    * MPO/S equivalent of a matrix vector multiplication. Let an MPO act on an MPS and return the new MPS
    * @param O input MPO
    * @param A input MPS
    * @return the new MPS object created by the multiplication
    */
   template<class Q>
      MPS<Q> operator*(const MPO<Q> &O,const MPS<Q> &A){

         MPS<Q> OA(O.size());

         gemv(1.0,O,A,0.0,OA);

         return OA;

      }

   /**
    * print the total bond dimensions
    */
   template<size_t N,class Q>
      void print_dim(const MPX<N,Q> &mpx){

         for(int i = 0;i < mpx.size();++i){

            int dim = 0;

            for(int j = 0;j < mpx[i].dshape(N-1).size();++j)
               dim += mpx[i].dshape(N-1)[j];

            cout << i << "\t" << dim << endl;

         }

      }

   /**
    * Let MPO A act on MPS X and compress at the same time to finite dimension D
    * @param A input MPO
    * @param X input MPS
    * @param Y output MPS, its content will change on exit.
    * @param D dimension of the compression
    * @return the discarded weigth of the compression
    */
   template<class Q>
      double gemv_compress(const MPO<Q> &A,const MPS<Q> &X,MPS<Q> &Y,int D){

         //first check if length is the same
         if(A.size() != X.size())
            BTAS_THROW(false, "Error: input objects do not have the same length!");

         int L = A.size();

         Y.resize(L);

         enum {j,k,l,m,n,o};

         double dweight = 0.0;

         QSDArray<5> tmp5;
         QSDArray<4> mrows;
         QSDArray<4> tmp4;

         SDArray<1> S;//singular values
         QSDArray<2> V;//V^T
         QSDArray<3> U;//U --> unitary left normalized matrix

         //first the most left tensor
         QSDindexed_contract(1.0,A[0],shape(j,k,l,m),X[0],shape(n,l,o),0.0,tmp5,shape(n,j,k,m,o));

         //merge 2 rows together
         TVector<Qshapes<Q>,2> qmerge;
         TVector<Dshapes,2> dmerge;

         for(int r = 0;r < 2;++r){

            qmerge[r] = tmp5.qshape(r);
            dmerge[r] = tmp5.dshape(r);

         }

         QSTmergeInfo<2> info(qmerge,dmerge);

         //clear the mrows object first
         mrows.clear();

         //then merge
         QSTmerge(info,tmp5,mrows);

         for(int i = 0;i < L - 1;++i){

            //merge 2 columns together
            for(int r = 2;r < 4;++r){

               qmerge[r - 2] = mrows.qshape(r);
               dmerge[r - 2] = mrows.dshape(r);

            }

            info.reset(qmerge,dmerge);

            QSTmerge(mrows,info,Y[i]);

            //svd
            S.clear();
            U.clear();
            V.clear();

            dweight += QSDgesvd(RightArrow,Y[i],S,U,V,0);

            //copy unitary to mpx
            QSDcopy(U,Y[i]);

            //paste S and V together
            SDdidm(S,V);

            //now expand V back
            U.clear();

            QSTexpand(V,info,U);

            //paste the next terms, A[i] and X[i] to U
            tmp4.clear();

            QSDindexed_contract(1.0,U,shape(j,k,l),X[i + 1],shape(l,m,n),0.0,tmp4,shape(j,k,m,n));

            mrows.clear();

            QSDindexed_contract(1.0,tmp4,shape(j,k,m,n),A[i + 1],shape(k,l,m,o),0.0,mrows,shape(j,l,o,n));

         }

         //merge 2 columns together
         for(int r = 2;r < 4;++r){

            qmerge[r - 2] = mrows.qshape(r);
            dmerge[r - 2] = mrows.dshape(r);

         }

         info.reset(qmerge,dmerge);

         QSTmerge(mrows,info,Y[L - 1]);

         print_dim(Y);

         dweight += compress(Y,Right,D);

         return dweight;

      }

   /**
    * MPO equivalent of a matrix matrix multiplication gemm in blas. MPO action on MPO gives new MPO: alpha A-B + C|MPS>
    * or new MPO is C <- alpha * A*B + beta * C
    * @param alpha scaling factor
    * @param A input MPO
    * @param B input MPO
    * @param beta scaling factor
    * @param C output MPO, input of this matrix will be added to A*B
    */
   template<class Q>
      void gemm(double alpha,const MPO<Q> &A,const MPO<Q> &B,double beta,MPO<Q> &C){

         //first check if we can sum these two:
         if(A.size() != B.size())
            BTAS_THROW(false, "Error: input objects do not have the same length!");

         int L = A.size();

         if( fabs(beta) < 1.0e-15 ){//only A * B

            enum {j,k,l,m,n,o,p};

            QSDArray<6> tmp;
            QSDArray<5> mrows;

            for(int i = 0;i < L;++i){

               //clear the tmp object first
               tmp.clear();

               QSDindexed_contract(1.0,A[i],shape(n,o,k,p),B[i],shape(j,k,l,m),0.0,tmp,shape(n,j,o,l,p,m));

               //merge 2 rows together
               TVector<Qshapes<Q>,2> qmerge;
               TVector<Dshapes,2> dmerge;

               for(int r = 0;r < 2;++r){

                  qmerge[r] = tmp.qshape(r);
                  dmerge[r] = tmp.dshape(r);

               }

               QSTmergeInfo<2> info(qmerge,dmerge);

               //clear the mrows object first
               mrows.clear();

               //then merge
               QSTmerge(info,tmp,mrows);

               //merge 2 columns together
               for(int r = 3;r < 5;++r){

                  qmerge[r - 3] = mrows.qshape(r);
                  dmerge[r - 3] = mrows.dshape(r);

               }

               info.reset(qmerge,dmerge);

               C[i].clear();
               QSTmerge(mrows,info,C[i]);

            }

            if(fabs(alpha - 1.0) > 1.0e-15)
               scal(alpha,C);


         }
         else{//beta != 0.0

            //first check if we can sum these two:
            if(L != C.size())
               BTAS_THROW(false, "Error: input objects do not have the same length!");

            scal(beta/alpha,C);

            //first left

            enum {j,k,l,m,n,o,p};

            QSDArray<6> tmp;
            QSDArray<5> mrows;

            //clear the tmp object first
            tmp.clear();

            QSDindexed_contract(1.0,A[0],shape(n,o,k,p),B[0],shape(j,k,l,m),0.0,tmp,shape(n,j,o,l,p,m));

            //merge 2 rows together
            TVector<Qshapes<Q>,2> qmerge1;
            TVector<Dshapes,2> dmerge1;

            for(int r = 0;r < 2;++r){

               qmerge1[r] = tmp.qshape(r);
               dmerge1[r] = tmp.dshape(r);

            }

            QSTmergeInfo<2> info1(qmerge1,dmerge1);

            //clear the mrows object first
            mrows.clear();

            //then merge
            QSTmerge(info1,tmp,mrows);

            //merge 2 columns together
            for(int r = 3;r < 5;++r){

               qmerge1[r - 3] = mrows.qshape(r);
               dmerge1[r - 3] = mrows.dshape(r);

            }

            info1.reset(qmerge1,dmerge1);

            //make AB
            QSDArray<4> AB;
            QSTmerge(mrows,info1,AB);

            IVector<3> left;

            for(int i = 0;i < 3;++i)
               left[i] = i;

            QSDArray<4> tmp1;
            QSDArray<4> tmp2;
            QSDdsum(AB,C[0],left,tmp1);

            //merge the column quantumnumbers together
            TVector<Qshapes<Q>,1> qmerge2;
            TVector<Dshapes,1> dmerge2;

            qmerge2[0] = tmp1.qshape(3);
            dmerge2[0] = tmp1.dshape(3);

            QSTmergeInfo<1> info2(qmerge2,dmerge2);

            //then merge
            C[0].clear();
            QSTmerge(tmp1,info2,C[0]);

            IVector<2> middle;

            for(int i = 0;i < 2;++i)
               middle[i] = i + 1;

            for(int i = 1;i < L - 1;++i){

               //clear the tmp object first
               tmp.clear();

               QSDindexed_contract(1.0,A[i],shape(n,o,k,p),B[i],shape(j,k,l,m),0.0,tmp,shape(n,j,o,l,p,m));

               //merge 2 rows together
               for(int r = 0;r < 2;++r){

                  qmerge1[r] = tmp.qshape(r);
                  dmerge1[r] = tmp.dshape(r);

               }

               info1.reset(qmerge1,dmerge1);

               //clear the mrows object first
               mrows.clear();

               //then merge
               QSTmerge(info1,tmp,mrows);

               //merge 2 columns together
               for(int r = 3;r < 5;++r){

                  qmerge1[r - 3] = mrows.qshape(r);
                  dmerge1[r - 3] = mrows.dshape(r);

               }

               info1.reset(qmerge1,dmerge1);

               AB.clear();
               QSTmerge(mrows,info1,AB);

               tmp1.clear();
               QSDdsum(AB,C[i],middle,tmp1);

               //merge the row quantumnumbers together
               qmerge2[0] = tmp1.qshape(0);
               dmerge2[0] = tmp1.dshape(0);

               info2.reset(qmerge2,dmerge2);

               //then merge
               tmp2.clear();
               QSTmerge(info2,tmp1,tmp2);

               //column quantumnumbers
               qmerge2[0] = tmp2.qshape(3);
               dmerge2[0] = tmp2.dshape(3);

               info2.reset(qmerge2,dmerge2);

               //then merge
               C[i].clear();
               QSTmerge(tmp2,info2,C[i]);

            }

            //clear the tmp object first
            tmp.clear();

            QSDindexed_contract(1.0,A[L - 1],shape(n,o,k,p),B[L - 1],shape(j,k,l,m),0.0,tmp,shape(n,j,o,l,p,m));

            //merge 2 rows together
            for(int r = 0;r < 2;++r){

               qmerge1[r] = tmp.qshape(r);
               dmerge1[r] = tmp.dshape(r);

            }

            info1.reset(qmerge1,dmerge1);

            //clear the mrows object first
            mrows.clear();

            //then merge
            QSTmerge(info1,tmp,mrows);

            //merge 2 columns together
            for(int r = 3;r < 5;++r){

               qmerge1[r - 3] = mrows.qshape(r);
               dmerge1[r - 3] = mrows.dshape(r);

            }

            info1.reset(qmerge1,dmerge1);

            AB.clear();
            QSTmerge(mrows,info1,AB);

            IVector<3> right;

            for(int i = 0;i < 3;++i)
               right[i] = i + 1;

            //finally the right
            tmp1.clear();
            QSDdsum(AB,C[L - 1],right,tmp1);

            //merge the row quantumnumbers together
            qmerge2[0] = tmp1.qshape(0);
            dmerge2[0] = tmp1.dshape(0);

            info2.reset(qmerge2,dmerge2);

            //then merge
            QSTmerge(info2,tmp1,C[L-1]);

            if( fabs(alpha - 1.0) > 1.0e-15)
               scal(alpha,C);

         }

      }

   /**
    * MPO equivalent of a matrix matrix multiplication. MPO action on MPO gives new MPO: O1-O2|MPS>
    * @param O1 input MPO
    * @param O2 input MPO
    * @return the new MPO object created by the multiplication
    */
   template<class Q>
      MPO<Q> operator*(const MPO<Q> &O1,const MPO<Q> &O2){

         MPO<Q> O12(O1.size());

         gemm(1.0,O1,O2,0.0,O12);

         return O12;

      }


   /**
    * the contraction of two MPS's
    * @return the overlap of two MPS objects
    * @param dir if dir = mpsxx::Left contract from left to right, if not contract from right to left
    * @param X input MPS
    * @param Y input MPS
    */
   template<class Q>
      double dot(const MPS_DIRECTION &dir,const MPS<Q> &X,const MPS<Q> &Y){

         int L = X.size();

         if(L != Y.size())
            cout << "Error: input MPS objects do not have the same length!" << endl;

         if(X[L-1].qshape(2) != Y[L-1].qshape(2))
            cout << "Error: input MPS objects do not have the same total quantumnumbers!" << endl;

         QSDArray<2> E;

         //going from left to right
         if(dir == Left){

            QSDcontract(1.0,X[0],shape(0,1),Y[0].conjugate(),shape(0,1),0.0,E);

            //this will contain an intermediate
            QSDArray<3> I;

            for(unsigned int i = 1;i < L;++i){

               //construct intermediate, i.e. past X to E
               QSDcontract(1.0,E,shape(0),X[i],shape(0),0.0,I);

               //clear structure of E
               E.clear();

               //construct E for site i by contracting I with Y
               QSDcontract(1.0,I,shape(0,1),Y[i].conjugate(),shape(0,1),0.0,E);

               I.clear();

               //bad style: if no blocks remain, return zero
               if(E.begin() == E.end())
                  return 0.0;

            }

         }
         else{ //going from right to left

            enum {j,k,l,m,n,o};

            QSDindexed_contract(1.0,X[L-1],shape(j,k,l),Y[L-1].conjugate(),shape(m,k,l),0.0,E,shape(j,m));

            //this will contain an intermediate
            QSDArray<3> I;

            for(int i = L - 2;i >= 0;--i){

               //construct intermediate, i.e. paste X to E
               QSDindexed_contract(1.0,X[i],shape(j,k,l),E,shape(l,m),0.0,I,shape(j,k,m));

               //clear structure of E
               E.clear();

               //construct E for site i by contracting I with Y
               QSDindexed_contract(1.0,Y[i].conjugate(),shape(j,k,l),I,shape(m,k,l),0.0,E,shape(m,j));

               I.clear();

               //bad style: if no blocks remain, return zero
               if(E.begin() == E.end())
                  return 0.0;

            }

         }

         return (*(E.find(shape(0,0))->second))(0,0);

      }

   /**
    * the contraction of two MPS's: easier notation, always from left to right
    * @return the overlap of two MPS objects
    * @param X input MPS
    * @param Y input MPS
    */
   template<class Q>
      double operator*(const MPS<Q> &X,const MPS<Q> &Y){

         return dot(mpsxx::Left,X,Y);

      }

   /**
    * @return the norm of the state
    */
   template<class Q>
      double nrm2(const MPS<Q> &mps){

         return sqrt(dot(Left,mps,mps));

      }

   /**
    * @return the distance between 2 mps's ||X - Y||_2
    */
   template<class Q>
      double dist(const MPS<Q> &X,const MPS<Q> &Y){

         return dot(Left,X,X) + dot(Left,Y,Y) - 2.0 * dot(Left,X,Y);

      }

   /**
    * normalize the MPS
    */
   template<class Q>
      void normalize(MPS<Q> &mps){

         double nrm = nrm2(mps);

         scal(1.0/nrm,mps);

      }

   /**
    * @return the MPO that is the result of the expontential of the input operator MPO O: output wil be e^O = 1 + O + 1/2 O^2 + ...
    * @param O input MPO
    * @param cutoff vector of size the number of terms in the expansion that will be kept, and containing the dimension for svd for every order
    */
   template<class Q>
      MPO<Q> exp(const MPO<Q> &O,const std::vector<int> &cutoff){

         std::vector< MPO<Q> > term(cutoff.size());

         //form the list of contributing terms in the expansion
         term[0] = O;

         compress(term[0],Left,0);
         compress(term[0],Right,cutoff[0]);

         for(int i = 1;i < cutoff.size();++i){

            term[i] = O*term[i - 1];
            compress(term[i],Left,0);
            compress(term[i],Right,cutoff[i]);
            scal(1.0/(i + 1.0),term[i]);

         }

         //now sum all the terms together:
         for(int i = 1;i < cutoff.size();++i){

            axpy(1.0,term[i],term[0]);
            compress(term[0],Left,0);
            compress(term[0],Right,0);

         }

         return term[0];

      }

   /**
    * @return the MPS that is the result of the expontential of the operator MPO O acting on input MPS A: output will be A + OA + 1/2 (O^2)A + ...
    * @param O input MPO
    * @param A input MPS
    * @param cutoff vector of size the number of terms in the expansion that will be kept, and containing the dimension for svd for every order
    */
   template<class Q>
      MPS<Q> exp(const MPO<Q> &O,const MPS<Q> &A,const std::vector<int> &cutoff){

         std::vector< MPS<Q> > term(cutoff.size());

         //form the list of contributing terms in the expansion
         term[0] = O*A;

         compress(term[0],Left,0);
         compress(term[0],Right,cutoff[0]);

         for(int i = 1;i < cutoff.size();++i){

            term[i] = O*term[i - 1];
            compress(term[i],Left,0);
            compress(term[i],Right,cutoff[i]);
            scal(1.0/(i + 1.0),term[i]);

         }

         //now sum all the terms together:
         for(int i = 1;i < cutoff.size();++i){

            axpy(1.0,term[i],term[0]);
            compress(term[0],Left,0);
            compress(term[0],Right,0);

         }

         axpy(1.0,A,term[0]);

         return term[0];

      }

   /**
    * @return the MPS that is the result of the exponential of the operator MPO O acting on input MPS A: output will be A + OA + 1/2 (O^2)A + ... 
    * expanded until order 'order' and compressed to dimension D
    * @param O input MPO
    * @param A input MPS
    * @param order order of the expansion
    * @param D compression dimensions
    */
   template<class Q>
      MPS<Q> exp(const MPO<Q> &O,const MPS<Q> &A,int order,int D){

         MPS<Q> eOA(A);
         MPS<Q> tmp;

         gemv(1.0/(double)order,O,A,1.0,eOA);

         compress(eOA,Left,0);
         compress(eOA,Right,D);

         for(int n = order - 1;n > 0;--n){

            tmp = std::move(eOA);
            eOA = A;

            gemv(1.0/(double)n,O,tmp,1.0,eOA);

            compress(eOA,Left,0);
            compress(eOA,Right,D);

         }

         return eOA;

      }


   /**
    * @param mpx will be written to file
    * @param filename name of the file
    * save the MPX object to a file in binary format.
    */
   template<size_t N,class Q>
      void save_mpx(const MPX<N,Q> &mpx,const char *filename){

         for(int i = 0;i < mpx.size();++i){

            char name[50];

            sprintf(name,"%s/%d.mpx",filename,i);

            std::ofstream fout(name);
            boost::archive::binary_oarchive oar(fout);

            oar << mpx[i];

         }

      }

   /**
    * @param mpx will be constructed from file
    * @param filename name of the file
    * load the MPX object from a file in binary format.
    */
   template<size_t N,class Q>
      void load_mpx(MPX<N,Q> &mpx,const char *filename){

         for(int i = 0;i < mpx.size();++i){

            char name[50];

            sprintf(name,"%s/%d.mpx",filename,i);

            std::ifstream fin(name);
            boost::archive::binary_iarchive iar(fin);
            iar >> mpx[i];

         }

      }

}

#endif 
