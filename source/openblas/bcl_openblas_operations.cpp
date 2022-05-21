// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "openblas/bcl_openblas_operations.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_default_flag_types.h"
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_operations.h"
#include "linal/bcl_linal_operations_interface.h"
#include "linal/bcl_linal_vector_interface.h"
#include "util/bcl_util_assert.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically
#include "OpenBlas/cblas.h"

namespace bcl
{
  namespace openblas
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new OperationsInterface
    template< typename t_DataType>
    Operations< t_DataType> *Operations< t_DataType>::Clone() const
    {
      return new Operations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Operations< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    // This namespace contains openblas functions for float and double seperately
    namespace
    {
      //! @brief dot product of two vectors of doubles
      double CblasDotProduct
      (
        const linal::VectorInterface< double> &VECTOR_A,
        const linal::VectorInterface< double> &VECTOR_B
      )
      {
        return cblas_ddot( VECTOR_A.GetSize(), VECTOR_A.Begin(), 1, VECTOR_B.Begin(), 1);
      }

      //! @brief dot product of two vectors of float
      float CblasDotProduct
      (
        const linal::VectorInterface< float> &VECTOR_A,
        const linal::VectorInterface< float> &VECTOR_B
      )
      {
        return cblas_sdot( VECTOR_A.GetSize(), VECTOR_A.Begin(), 1, VECTOR_B.Begin(), 1);
      }

      //! @brief outer product of two vectors
      //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
      //! @param VECTOR_U left hand side vector with m elements
      //! @param VECTOR_V right hand side vector with n elements
      //! @return Matrix with m*n elements (outer product of vector u and v)
      void CblasOuterProduct
      (
        const linal::VectorConstInterface< double> &VECTOR_U,
        const linal::VectorConstInterface< double> &VECTOR_V,
        linal::MatrixInterface< double> &MATRIX
      )
      {
        cblas_dger
        (
          CblasRowMajor, MATRIX.GetNumberRows(), MATRIX.GetNumberCols(), 1.0,
          VECTOR_U.Begin(), 1, VECTOR_V.Begin(), 1, MATRIX.Begin(), MATRIX.GetNumberCols()
        );
      }

      //! @brief outer product of two vectors
      //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
      //! @param VECTOR_U left hand side vector with m elements
      //! @param VECTOR_V right hand side vector with n elements
      //! @return Matrix with m*n elements (outer product of vector u and v)
      void CblasOuterProduct
      (
        const linal::VectorConstInterface< float> &VECTOR_U,
        const linal::VectorConstInterface< float> &VECTOR_V,
        linal::MatrixInterface< float> &MATRIX
      )
      {
        cblas_sger
        (
          CblasRowMajor, MATRIX.GetNumberRows(), MATRIX.GetNumberCols(), 1.0,
          VECTOR_U.Begin(), 1, VECTOR_V.Begin(), 1, MATRIX.Begin(), MATRIX.GetNumberCols()
        );
      }

      //! @brief norm of a vector of double
      double CblasNrm2( const linal::VectorInterface< double> &VECTOR)
      {
        return cblas_dnrm2( VECTOR.GetSize(), VECTOR.Begin(), 1);
      }

      //! @brief norm of a vector of float
      float CblasNrm2( const linal::VectorInterface< float> &VECTOR)
      {
        return cblas_snrm2( VECTOR.GetSize(), VECTOR.Begin(), 1);
      }

      //! @brief matrix-matrix multiplication
      //! @return stores the results of MATRIX_AxMATRIX_B in MATRIX_C
      void CblasGemm
      (
        linal::MatrixInterface< double> &MATRIX_C,
        const linal::MatrixConstInterface< double> &MATRIX_A,
        const linal::MatrixConstInterface< double> &MATRIX_B
      )
      {
        cblas_dgemm
        (
          CblasRowMajor, CblasNoTrans, CblasNoTrans, MATRIX_A.GetNumberRows(),
          MATRIX_B.GetNumberCols(), MATRIX_A.GetNumberCols(), 1.0, MATRIX_A.Begin(),
          MATRIX_A.GetNumberCols(), MATRIX_B.Begin(), MATRIX_B.GetNumberCols(),
          1.0, MATRIX_C.Begin(), MATRIX_C.GetNumberCols()
        );
      }

      //! @brief matrix-matrix multiplication
      //! @return stores the results of MATRIX_AxMATRIX_B in MATRIX_C
      void CblasGemm
      (
        linal::MatrixInterface< float> &MATRIX_C,
        const linal::MatrixConstInterface< float> &MATRIX_A,
        const linal::MatrixConstInterface< float> &MATRIX_B
      )
      {
        cblas_sgemm
        (
          CblasRowMajor, CblasNoTrans, CblasNoTrans, MATRIX_A.GetNumberRows(),
          MATRIX_B.GetNumberCols(), MATRIX_A.GetNumberCols(), 1.0, MATRIX_A.Begin(),
          MATRIX_A.GetNumberCols(), MATRIX_B.Begin(), MATRIX_B.GetNumberCols(),
          1.0, MATRIX_C.Begin(), MATRIX_C.GetNumberCols()
        );
      }

      //! @brief matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void CblasGemv
      (
        linal::VectorInterface< double> &STORAGE,
        const linal::MatrixConstInterface< double> &MATRIX,
        const linal::VectorConstInterface< double> &FEATURE
      )
      {
        cblas_dgemv
        (
          CblasRowMajor, CblasNoTrans, MATRIX.GetNumberRows(),
          MATRIX.GetNumberCols(), 1, MATRIX.Begin(),
          MATRIX.GetNumberCols(), FEATURE.Begin(), 1, 1,
          STORAGE.Begin(), 1
        );
      }

      //! @brief matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void CblasGemv
      (
        linal::VectorInterface< float> &STORAGE,
        const linal::MatrixConstInterface< float> &MATRIX,
        const linal::VectorConstInterface< float> &FEATURE
      )
      {
        cblas_sgemv
        (
          CblasRowMajor, CblasNoTrans, MATRIX.GetNumberRows(),
          MATRIX.GetNumberCols(), 1, MATRIX.Begin(),
          MATRIX.GetNumberCols(), FEATURE.Begin(), 1, 1,
          STORAGE.Begin(), 1
        );
      }
      //! @brief matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void CblasGevm
      (
        linal::VectorInterface< double> &STORAGE,
        const linal::MatrixConstInterface< double> &MATRIX,
        const linal::VectorConstInterface< double> &FEATURE
      )
      {
        cblas_dgemv
        (
          CblasRowMajor, CblasTrans, MATRIX.GetNumberRows(),
          MATRIX.GetNumberCols(), 1, MATRIX.Begin(),
          MATRIX.GetNumberCols(), FEATURE.Begin(), 1, 1,
          STORAGE.Begin(), 1
        );
      }

      //! @brief matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void CblasGevm
      (
        linal::VectorInterface< float> &STORAGE,
        const linal::MatrixConstInterface< float> &MATRIX,
        const linal::VectorConstInterface< float> &FEATURE
      )
      {
        cblas_sgemv
        (
          CblasRowMajor, CblasTrans, MATRIX.GetNumberRows(),
          MATRIX.GetNumberCols(), 1, MATRIX.Begin(),
          MATRIX.GetNumberCols(), FEATURE.Begin(), 1, 1,
          STORAGE.Begin(), 1
        );
        //BCL_Debug(STORAGE);
      }
    } // ends namespace

    //! @brief dot product of two vectors
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return dot product of VECTOR_A * VECTOR_B
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::DotProduct
    (
      const linal::VectorInterface< t_DataType> &VECTOR_A,
      const linal::VectorInterface< t_DataType> &VECTOR_B
    ) const
    {
      BCL_Assert
      (
        VECTOR_A.GetSize() == VECTOR_B.GetSize(), "Can only compute dot products between vectors of same size"
      );
      return CblasDotProduct( VECTOR_A, VECTOR_B);
    }

    //! @brief outer product of two vectors
    //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
    //! @param VECTOR_U left hand side vector with m elements
    //! @param VECTOR_V right hand side vector with n elements
    //! @return Matrix with m*n elements (outer product of vector u and v)
    template< typename t_DataType>
    linal::Matrix< t_DataType> Operations< t_DataType>::OuterProduct
    (
      const linal::VectorConstInterface< t_DataType> &VECTOR_U,
      const linal::VectorConstInterface< t_DataType> &VECTOR_V
    )
    {
      linal::Matrix< t_DataType> outer_product( VECTOR_U.GetSize(), VECTOR_V.GetSize());
      CblasOuterProduct( VECTOR_U, VECTOR_V, outer_product);
      return outer_product;
    }

    //! @brief norm of a vector
    //! @param VECTOR vector
    //! @return the Euclidean norm of the VECTOR
    template< typename t_DataType>
    t_DataType Operations< t_DataType>::Norm( const linal::VectorInterface< t_DataType> &VECTOR) const
    {
      return CblasNrm2( VECTOR);
    }

    //! @brief matrix * vector
    //! @param MATRIX matrix
    //! @param VECTOR vector
    //! @return product of MATRIX * VECTOR
    template< typename t_DataType>
    linal::Vector< t_DataType> Operations< t_DataType>::Multiply
    (
      const linal::MatrixConstInterface< t_DataType> &MATRIX,
      const linal::VectorInterface< t_DataType> &VECTOR
    ) const
    {
      BCL_Assert
      (
        VECTOR.GetSize() == MATRIX.GetNumberCols(), "Can only multiply a vector and a matrix of same size"
      );
      linal::Vector< t_DataType> vector_b( MATRIX.GetNumberRows());
      VectorPlusEqualsMatrixTimesVector( vector_b, MATRIX, VECTOR);
      return vector_b;
    }

    //! @brief matrix-matrix multiplication
    //! @param MATRIX_A first matrix
    //! @param MATRIX_B second matrix
    //! @return product of MATRIX_A * MATRIX_B
    template< typename t_DataType>
    linal::Matrix< t_DataType> Operations< t_DataType>::Multiply
    (
      const linal::MatrixConstInterface< t_DataType> &MATRIX_A,
      const linal::MatrixConstInterface< t_DataType> &MATRIX_B
    ) const
    {
      BCL_Assert
      (
        MATRIX_A.GetNumberCols() == MATRIX_B.GetNumberRows(), "Multiply matrix and matrix of wrong size"
      );
      linal::Matrix< t_DataType> matrix_c( MATRIX_A.GetNumberRows(), MATRIX_B.GetNumberCols());

      CblasGemm( matrix_c, MATRIX_A, MATRIX_B);
      return matrix_c;
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void Operations< t_DataType>::VectorEqualsVectorTimesMatrix
    (
      linal::VectorInterface< t_DataType> &STORAGE,
      const linal::VectorConstInterface< t_DataType> &FEATURE,
      const linal::MatrixConstInterface< t_DataType> &MATRIX
    ) const
    {
      const size_t num_rows( MATRIX.GetNumberRows());
      const size_t num_cols( MATRIX.GetNumberCols());
      BCL_Assert
      (
        num_cols == STORAGE.GetSize() && num_rows == FEATURE.GetSize(),
        "non-matching dimensions! " + util::Format()( num_rows) + "X" + util::Format()( num_cols)
        + " * " + util::Format()( FEATURE.GetSize()) + " * " + util::Format()( STORAGE.GetSize())
      );
      // initialize the storage
      STORAGE = 0;
      CblasGevm( STORAGE, MATRIX, FEATURE);
    }

    //! @brief matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void Operations< t_DataType>::VectorPlusEqualsMatrixTimesVector
    (
      linal::VectorInterface< t_DataType> &STORAGE,
      const linal::MatrixConstInterface< t_DataType> &MATRIX,
      const linal::VectorConstInterface< t_DataType> &FEATURE
    ) const
    {
      const size_t num_rows( MATRIX.GetNumberRows());
      const size_t num_cols( MATRIX.GetNumberCols());
      BCL_Assert
      (
        num_cols == FEATURE.GetSize() && num_rows == STORAGE.GetSize(),
        "non-matching dimensions! " + util::Format()( num_rows) + "X" + util::Format()( num_cols)
        + " * " + util::Format()( FEATURE.GetSize()) + " * " + util::Format()( STORAGE.GetSize())
      );
      CblasGemv( STORAGE, MATRIX, FEATURE);
    }

    //! @brief instance of this operations class in Operations enumerator
    template< typename t_DataType>
    const typename linal::Operations< t_DataType>::EnumType Operations< t_DataType>::e_OpenBlas
    (
      linal::GetOperationsNonConst< t_DataType>().AddEnum
      (
        "OpenBlas", util::ShPtr< linal::OperationsInterface< t_DataType> >( new Operations< t_DataType>())
      )
    );

    namespace
    {
      //! @brief Initialize the logger from the command line flag
      void UpdateCurrentOperationFromCommandLineFlag()
      {
        if( GetFlagOpenblas()->GetFlag())
        {
          linal::GetOperationsNonConst< float>().SetDefaultOperationsType( Operations< float>::e_OpenBlas);
          linal::GetOperationsNonConst< double>().SetDefaultOperationsType( Operations< double>::e_OpenBlas);
          const int thread_num
          (
            util::ConvertStringToNumericalValue< int>( GetFlagOpenblas()->GetFirstParameter()->GetValue())
          );
          openblas_set_num_threads( thread_num);
          // switch was successful
          BCL_MessageStd( "OpenBlas operations can be used now");
        }
        else
        {
          linal::GetOperationsNonConst< float>().SetDefaultOperationsType( linal::Operations< float>::Operator( "CPU"));
          linal::GetOperationsNonConst< double>().SetDefaultOperationsType( linal::Operations< double>::Operator( "CPU"));
        }
      }
    }

    //! command line flag to be used to set Logger over the command line
    util::ShPtr< command::FlagInterface> &GetFlagOpenblas()
    {
      static util::ShPtr< command::FlagInterface> s_openblas_flag;
      if( !s_openblas_flag.IsDefined())
      {
        s_openblas_flag = util::ShPtr< command::FlagInterface>
        (
          new command::FlagStatic
          (
            "openblas",
            "use the openblas operations for matrix operations",
            util::ShPtrVector< command::ParameterInterface>::Create
            (
              util::ShPtr< command::ParameterInterface>
              (
                new command::Parameter
                (
                  "thread number",
                  "number of cores to be used by openblas operation",
                  command::ParameterCheckRanged< int>( 1, 1024),
                  "1"
                )
              )
            ),
            &UpdateCurrentOperationFromCommandLineFlag
          )
        );
        command::GetAppDefaultFlags().AddDefaultFlag( s_openblas_flag, command::e_OpenBlas);
      }
      return s_openblas_flag;
    }

    template< typename t_DataType>
    bool Operations< t_DataType>::s_IsFlagInstantiated = GetFlagOpenblas().IsDefined();

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Operations< float>;
    template class BCL_API Operations< double>;

  } // namespace openblas
} // namespace bcl
