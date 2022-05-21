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
#include "linal/bcl_linal.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_distance_geometry.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_principal_component_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
    //! @param SYMMETRIC_MATRIX containing pairwise distances for target
    //! @return function value of the given argument
    storage::Vector< VectorND< double, 3> > DistanceGeometry::operator()
    (
      const storage::SymmetricMatrix< double> &SYMMETRIC_MATRIX
    ) const
    {
      // Set result dimension
      const size_t result_dimension( 3);

      // Get input matrix size
      const size_t matrix_size( SYMMETRIC_MATRIX.GetSize());

      // Create matrix to be filled as the metric matrix
      Matrix< double> metric_matrix( matrix_size, matrix_size, 0.0);

      // TODO: Use math::Sqr()

      // Make a metric matrix using M = [m_ij] of dimension n, where m_ij = 1/2( p_0i^2 + p_0j^2 - p_ij^2)
      for( size_t i( 0); i < matrix_size; ++i)
      {
        for( size_t j( 0); j < matrix_size; ++j)
        {
          double p_0i = SYMMETRIC_MATRIX( i, 0);
          double p_0j = SYMMETRIC_MATRIX( 0, j);
          double p_ij = SYMMETRIC_MATRIX( i, j);
          metric_matrix( i, j) = ( pow( p_0i, 2) + pow( p_0j, 2) - pow( p_ij, 2)) / 2.0;
        }
      }

      // Create matrix to be filled along the diagonal with eigen values and empty eigen vectors matrix
      Matrix< double> eigen_values_matrix( metric_matrix);
      Matrix< double> eigen_vectors_matrix( matrix_size, matrix_size, 0.0);

      // Calculate eigen values
      SymmetricEigenSolver< double> eigensolver( metric_matrix, true);

      // Place diagonal containing values in a vector
      Vector< double> sorted_eigen_values( eigensolver.GetSortedEigenvalues());

      // Get sorted eigen_vectors matrix
      Matrix< double> sorted_eigen_vectors_matrix( eigensolver.GetSortedEigenvectors());

      // Create vector for checked and zeroed eigenvalues beyond the targeted dimension
      Matrix< double> clean_eigenvalues( 1, result_dimension, 0.0);

      // Take top three eigenvectors and check that they are all positive or zero, if negative set to zero
      for( size_t index( 0); index < result_dimension; ++index)
      {
        if( sorted_eigen_values( index) > 0)
        {
          // Take square root of each eigenvalue
          clean_eigenvalues( 0, index) = sqrt( sorted_eigen_values( index));
        }
      };

      // Multiply eigenvector by corresponding eigenvalue
      Matrix< double> result_matrix( matrix_size, result_dimension, 0.0);
      for( size_t col( 0); col < result_dimension; ++col)
      {
        for( size_t row( 0); row < matrix_size; ++row)
        {
          result_matrix( row, col) = clean_eigenvalues( 0, col) * sorted_eigen_vectors_matrix( row, col);
        }
      }

      // Create vector of size N to take all coordinates
      storage::Vector< VectorND< double, 3> > result_vector( matrix_size);

      // Iterate through n matrix placing first three positions into vector of VectorND of size result_dimension (3)
      VectorND< double, 3> temp_vect_nd;
      for( size_t i( 0); i < matrix_size; ++i)
      {
        for( size_t index( 0); index < result_dimension; ++index)
        {
          temp_vect_nd( index) = result_matrix( i, index);
        }
        result_vector( i) = temp_vect_nd;
      }

// TODO: Take first three using constructor from pointer

      // Return vector with coordinates
      return result_vector;

      // TODO: Add check later on for size mismatch ??
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DistanceGeometry::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DistanceGeometry::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // write members
      return OSTREAM;
    }

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix2x2.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Matrix2x2< float>;
    template class BCL_API Matrix2x2< double>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix3x3.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Matrix3x3< float>;
    template class BCL_API Matrix3x3< double>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix_const_interface.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API MatrixConstInterface< double>;
    template class BCL_API MatrixConstInterface< float>;
    template class BCL_API MatrixConstInterface< int>;
    template class BCL_API MatrixConstInterface< unsigned int>;
    template class BCL_API MatrixConstInterface< unsigned long>;
    template class BCL_API MatrixConstInterface< unsigned long long>;
    template class BCL_API MatrixConstInterface< char>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Matrix< double>;
    template class BCL_API Matrix< float>;
    template class BCL_API Matrix< int>;
    template class BCL_API Matrix< unsigned int>;
    template class BCL_API Matrix< unsigned long>;
    template class BCL_API Matrix< unsigned long long>;
    template class BCL_API Matrix< char>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix_interface.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API MatrixInterface< double>;
    template class BCL_API MatrixInterface< float>;
    template class BCL_API MatrixInterface< int>;
    template class BCL_API MatrixInterface< unsigned int>;
    template class BCL_API MatrixInterface< unsigned long>;
    template class BCL_API MatrixInterface< unsigned long long>;
    template class BCL_API MatrixInterface< char>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix_inversion_cholesky.hpp"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    // register the class with the enumerated class instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface>
    MatrixInversionCholesky< t_DataType>::s_Instance
    (
      util::Enumerated< MatrixInversionInterface< t_DataType> >::AddInstance
      (
        new MatrixInversionCholesky< t_DataType>()
      )
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API MatrixInversionCholesky< double>;
    template class BCL_API MatrixInversionCholesky< float>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix_inversion_gauss_jordan.hpp"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

      // register the class with the enumerated class instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface>
    MatrixInversionGaussJordan< t_DataType>::s_PivotInstance
    (
      util::Enumerated< MatrixInversionInterface< t_DataType> >::AddInstance
      (
        new MatrixInversionGaussJordan< t_DataType>( true)
      )
    );

    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface>
    MatrixInversionGaussJordan< t_DataType>::s_NoPivotInstance
    (
      util::Enumerated< MatrixInversionInterface< t_DataType> >::AddInstance
      (
        new MatrixInversionGaussJordan< t_DataType>( false)
      )
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API MatrixInversionGaussJordan< double>;
    template class BCL_API MatrixInversionGaussJordan< float>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix_inversion_interface.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API MatrixInversionInterface< double>;
    template class BCL_API MatrixInversionInterface< float>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_matrix_inversion_moore_penrose.hpp"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

      // register the class with the enumerated class instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface>
    MatrixInversionMoorePenrose< t_DataType>::s_Instance
    (
      util::Enumerated< MatrixInversionInterface< t_DataType> >::AddInstance
      (
        new MatrixInversionMoorePenrose< t_DataType>()
      )
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API MatrixInversionMoorePenrose< double>;
    template class BCL_API MatrixInversionMoorePenrose< float>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_operations.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_interface.h"
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Operations< t_DataType>::Operations() :
      util::Enumerate< util::ShPtr< OperationsInterface< t_DataType> >, Operations< t_DataType> >( false)
    {
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

    //! @brief get the default operations class; linal unless -opencl was given with a valid platform
    template< typename t_DataType>
    const OperationsInterface< t_DataType> &Operations< t_DataType>::GetDefaultOperations() const
    {
      return **m_DefaultType;
    }

    //! @brief set the default operations type from given flag
    template< typename t_DataType>
    const typename Operations< t_DataType>::EnumType &Operations< t_DataType>::SetDefaultOperationsType( const EnumType &DEFAULT)
    {
      m_DefaultType = DEFAULT;
      return m_DefaultType;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Operations< float>;
    template class BCL_API Operations< double>;

  } // namespace linal

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< linal::OperationsInterface<  float> >, linal::Operations<  float> >;
    template class BCL_API Enumerate< ShPtr< linal::OperationsInterface< double> >, linal::Operations< double> >;

  } // namespace util
} // namespace bcl
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
#include "linal/bcl_linal_operations_cpu.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_const_interface.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new OperationsCPU
    template< typename t_DataType>
    OperationsCPU< t_DataType>  *OperationsCPU< t_DataType>::Clone() const
    {
      return new OperationsCPU( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &OperationsCPU< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief dot product of two vectors
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return dot product of VECTOR_A * VECTOR_B
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::DotProduct
    (
      const VectorConstInterface< t_DataType> &VECTOR_A,
      const VectorConstInterface< t_DataType> &VECTOR_B
    ) const
    {
      // check for identical size
      if( VECTOR_A.GetSize() != VECTOR_B.GetSize())
      {
        return util::GetUndefined< t_DataType>();
      }

      // inner product
      return std::inner_product( VECTOR_A.Begin(), VECTOR_A.End(), VECTOR_B.Begin(), t_DataType( 0));
    }

    //! @brief outer product of two vectors
    //! u x v = A while u has m elements, v has n elements, A is a m*n matrix (m rows, n cols)
    //! @param VECTOR_U left hand side vector with m elements
    //! @param VECTOR_V right hand side vector with n elements
    //! @return Matrix with m*n elements (outer product of vector u and v)
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::OuterProduct
    (
      const VectorConstInterface< t_DataType> &VECTOR_U,
      const VectorConstInterface< t_DataType> &VECTOR_V
    )
    {
      Matrix< t_DataType> outer_product( VECTOR_U.GetSize(), VECTOR_V.GetSize());

      // ptr to first element of first row in matrix
      t_DataType *ptr_matrix( outer_product.Begin());

      // each result row is the product of the current VECTOR_U element multiplied with VECTOR_V
      for( const t_DataType *ptr_u( VECTOR_U.Begin()), *ptr_u_end( VECTOR_U.End()); ptr_u != ptr_u_end; ++ptr_u)
      {
        std::transform
        (
          VECTOR_V.Begin(), VECTOR_V.End(),
          ptr_matrix,
          std::bind2nd( std::multiplies< t_DataType>(), *ptr_u)
        );
        // ptr_matrix points to first element in next row
        ptr_matrix += VECTOR_V.GetSize();
      }

      // end
      return outer_product;
    }

    //! @brief norm of a vector
    //! @param VECTOR vector
    //! @return the euclidean norm of the VECTOR
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Norm( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Norm();
    }

    //! @brief matrix-matrix multiplication
    //! @param MATRIX_A matrix to be multiplied
    //! @param MATRIX_B matrix to be multiplied
    //! @return resulting linal::Matrix< t_DataType>
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::Multiply
    (
      const MatrixConstInterface< t_DataType> &MATRIX_A,
      const MatrixConstInterface< t_DataType> &MATRIX_B
    ) const
    {
      return MATRIX_A * MATRIX_B;
    }

    //! @brief matrix-vector multiplication
    //! @param MATRIX matrix to be multiplied
    //! @param VECTOR vector to be multiplied
    //! @return resulting linal::Vector< t_DataType>
    template< typename t_DataType>
    Vector< t_DataType> OperationsCPU< t_DataType>::Multiply
    (
      const MatrixConstInterface< t_DataType> &MATRIX,
      const VectorConstInterface< t_DataType> &VECTOR
    ) const
    {
      return MATRIX * VECTOR;
    }

    //! @brief calculate pair-wise distances between lists of vectors
    //! @param LIST_VECTORS list of vectors
    //! @return triangular matrix of values of the distances
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::DistanceMatrix
    (
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS
    ) const
    {
      const size_t num_elements( LIST_VECTORS.GetSize());
      Matrix< t_DataType> results( num_elements, num_elements, util::GetUndefined< t_DataType>());
      // iterate over rows
      for( size_t row_num( 0); row_num < num_elements; ++row_num)
      {
        const VectorConstInterface< t_DataType> &vector_a( *LIST_VECTORS( row_num));
        // iterate over cols
        for( size_t col_num( row_num); col_num < num_elements; ++col_num)
        {
          // assign to matrix position
          results( row_num, col_num) = Distance( vector_a, *LIST_VECTORS( col_num));
        }
      }
      return results;
    }

    //! @brief calculate pair-wise distances between lists of vectors
    //! @param LIST_VECTORS_A list of vectors
    //! @param LIST_VECTORS_B list of vectors
    //! @return matrix of values of the distances
    template< typename t_DataType>
    Matrix< t_DataType> OperationsCPU< t_DataType>::DistanceMatrix
    (
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_A,
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS_B
    ) const
    {
      const size_t num_rows( LIST_VECTORS_A.GetSize());
      const size_t num_cols( LIST_VECTORS_B.GetSize());
      Matrix< t_DataType> results( num_rows, num_cols, util::GetUndefined< t_DataType>());
      // iterate over rows
      for( size_t row_num( 0); row_num < num_rows; ++row_num)
      {
        const VectorConstInterface< t_DataType> &vector_a( *LIST_VECTORS_A( row_num));
        // iterate over cols
        for( size_t col_num( 0); col_num < num_cols; ++col_num)
        {
          // assign to matrix position
          results( row_num, col_num) = Distance( vector_a, *LIST_VECTORS_B( col_num));
        }
      }
      return results;
    }

    //! @brief returns sum of all elements in vector
    //! @param VECTOR the vector from which to get the sum
    //! @return the sum
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Sum( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Sum();
    }

    //! @brief reduction kernel
    //! @param VECTOR the vector to reduce
    //! @return the reduced result
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Min( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Min();
    }

    //! @brief gets max
    //! @param VECTOR the vector in which to find the max
    //! @return the max value
    template< typename t_DataType>
    t_DataType OperationsCPU< t_DataType>::Max( const VectorConstInterface< t_DataType> &VECTOR) const
    {
      return VECTOR.Max();
    }

    //! @brief gets the min and max for each column in a matrix
    //! @param MATRIX the matrix input
    //! @return a storage vector of math ranges with min and max for each column in the matrix
    template< typename t_DataType>
    storage::Vector< math::Range< t_DataType> > OperationsCPU< t_DataType>::MinMax( const MatrixConstInterface< t_DataType> &MATRIX) const
    {
      storage::Vector< math::Range< t_DataType> > ranges;
      ranges.AllocateMemory( MATRIX.GetNumberCols());
      const size_t number_data_pts( MATRIX.GetNumberRows());

      // find ranges of each column
      t_DataType tmp( 0);
      for( size_t col( 0), number_features( MATRIX.GetNumberCols()); col < number_features; ++col)
      {
        // set initial min and max
        t_DataType min( std::numeric_limits< t_DataType>::infinity());
        t_DataType max( -std::numeric_limits< t_DataType>::infinity());

        // iterate through all values in that feature col
        for( size_t row( 0); row < number_data_pts; ++row)
        {
          // set as min or max if qualifies
          tmp = MATRIX( row, col);
          if( tmp < min)
          {
            min = tmp;
          }
          else if( tmp > max)
          {
            max = tmp;
          }
        }
        // add to range vector
        ranges.PushBack( math::Range< t_DataType>( min, max));
      }
      return ranges;
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void OperationsCPU< t_DataType>::VectorEqualsVectorTimesMatrix
    (
      VectorInterface< t_DataType> &STORAGE,
      const VectorConstInterface< t_DataType> &FEATURE,
      const MatrixConstInterface< t_DataType> &MATRIX
    ) const
    {
      linal::VectorEqualsVectorTimesMatrix( STORAGE, FEATURE, MATRIX);
    }

    //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
    //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
    template< typename t_DataType>
    void OperationsCPU< t_DataType>::VectorPlusEqualsMatrixTimesVector
    (
      VectorInterface< t_DataType> &STORAGE,
      const MatrixConstInterface< t_DataType> &MATRIX,
      const VectorConstInterface< t_DataType> &FEATURE
    ) const
    {
      linal::VectorPlusEqualsMatrixTimesVector( STORAGE, MATRIX, FEATURE);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &OperationsCPU< t_DataType>::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_DataType>
    std::ostream &OperationsCPU< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

    //! @brief instance of this operations class in Operations enumerator
    template< typename t_DataType>
    const typename Operations< t_DataType>::EnumType OperationsCPU< t_DataType>::e_CPU
    (
      Operations< t_DataType>::GetEnums().SetDefaultOperationsType
      (
        Operations< t_DataType>::GetEnums().AddEnum
        (
          "CPU", util::ShPtr< OperationsInterface< t_DataType> >( new OperationsCPU< t_DataType>())
        )
      )
    );

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API OperationsCPU< float>;
    template class BCL_API OperationsCPU< double>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_operations_interface.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_si_ptr_vector.h"

using bcl::command::Parameter;
using bcl::command::FlagStatic;

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief creates a matrix from a SiPtrVector of VectorInterfaces
    //! @param LIST_VECTORS the SiPtrVector of VectorInterfaces to be transferred into a matrix
    //! @return the matrix created from the SiPtrList of VectorInterfaces
    template< typename t_DataType>
    Matrix< t_DataType> OperationsInterface< t_DataType>::CopyListToMatrix
    (
      const util::SiPtrVector< const VectorConstInterface< t_DataType> > &LIST_VECTORS
    )
    {
      // check if list is empty
      if( LIST_VECTORS.IsEmpty())
      {
        return Matrix< t_DataType>();
      }

      // getting matrix dimensions
      const size_t num_rows( LIST_VECTORS.GetSize());
      const size_t num_cols( LIST_VECTORS.FirstElement()->GetSize());

      // allocate matrix
      Matrix< t_DataType> matrix( num_rows, num_cols, t_DataType( 0));
      t_DataType *row_itr( matrix.Begin());

      // iterate over list and rows
      for
      (
          typename util::SiPtrVector< const VectorConstInterface< t_DataType> >::const_iterator
          list_itr( LIST_VECTORS.Begin()),
          list_itr_end( LIST_VECTORS.End());
          list_itr != list_itr_end;
          ++list_itr, row_itr += num_cols
      )
      {
        // get current vector and verify size consistency
        const VectorConstInterface< t_DataType> &current_vector( **list_itr);
        if( current_vector.GetSize() != num_cols)
        {
          BCL_MessageCrt( "Supplied list of vectors have different lengths!");
          return Matrix< t_DataType>();
        }

        // copy vector into row
        std::copy( current_vector.Begin(), current_vector.End(), row_itr);
      }

      // end
      return matrix;
    }

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API OperationsInterface< float>;
    template class BCL_API OperationsInterface< double>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_vector_2d.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Vector2D::s_Instance
    (
      GetObjectInstances().AddInstance( new Vector2D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Vector2D::Vector2D() :
      m_X( 0.0),
      m_Y( 0.0)
    {
    }

    //! @brief constructor from a single value
    //! @param VALUE common value for all the elements
    Vector2D::Vector2D( const double &VALUE) :
      m_X( VALUE),
      m_Y( VALUE)
    {
    }

    //! @brief constructor from two values
    //! @param X, Y, Z two elements
    Vector2D::Vector2D( const double &X, const double &Y) :
      m_X( X),
      m_Y( Y)
    {
    }

    //! @brief constructor from a pointer to two elements
    //! @param PTR_DATA a pointer to two elements
    Vector2D::Vector2D( const double *PTR_DATA) :
      m_X( PTR_DATA[ 0]),
      m_Y( PTR_DATA[ 1])
    {
    }

    //! @brief constructor from VectorInterface
    //! @param VECTOR vector
    Vector2D::Vector2D( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "given vector has wrong length");
      std::copy( VECTOR.Begin(), VECTOR.End(), &m_X);
    }

    //! copy constructor
    Vector2D *Vector2D::Clone() const
    {
      return new Vector2D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Vector2D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns true if all the coordinates are defined
    //! @return boolean true if all coordinates are defined - false otherwise
    bool Vector2D::IsDefined() const
    {
      return util::IsDefined( m_X) && util::IsDefined( m_Y);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief set all elements to elements from another vector
    //! @param VECTOR vector, which must contain 2 elements
    //! @return reference to this, after copying vector
    Vector2D &Vector2D::operator =( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "Tried to copy non-2d vector onto vector2d");
      m_X = VECTOR( 0);
      m_Y = VECTOR( 1);
      return *this;
    }

    //! @brief set all elements to a vector 2d
    //! @param VECTOR the vector to copy into this vector
    //! @return reference to this, after copying vector
    Vector2D &Vector2D::operator =( const Vector2D &VECTOR)
    {
      m_X = VECTOR.m_X;
      m_Y = VECTOR.m_Y;
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! change all values
    Vector2D &Vector2D::Set( const double &X, const double &Y)
    {
      m_X = X;
      m_Y = Y;
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Vector2D::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_X, ISTREAM);
      io::Serialize::Read( m_Y, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Vector2D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Y, OSTREAM, 0);

      // return the stream
      return OSTREAM;
    }

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_vector_2d_operations.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    //! @brief calculates the absolute difference between two linal::Vector2Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return absolute difference between two linal::Vector2D (points)
    Vector2D AbsoluteDifference( const Vector2D &VECTOR_A, const Vector2D &VECTOR_B)
    {
      Vector2D returned_vector( VECTOR_B - VECTOR_A);
      math::Absolute( returned_vector);
      return returned_vector;
    }

    //! @brief calculates the unit vector starting from one linal::Vector2D to another
    //! @param ORIGIN vector of origin
    //! @param TARGET target vector
    //! @return the unit vector between ORIGIN and TARGET
    Vector2D UnitVector( const Vector2D &ORIGIN, const Vector2D &TARGET)
    {
      Vector2D unit_vector( TARGET);
      unit_vector -= ORIGIN;
      unit_vector.Normalize();
      return unit_vector;
    }

  } // namespace linal
} // namespace bcl
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
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "linal/bcl_linal_vector_3d.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Vector3D::s_Instance
    (
      GetObjectInstances().AddInstance( new Vector3D())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Vector3D::Vector3D() :
      m_X( 0.0),
      m_Y( 0.0),
      m_Z( 0.0)
    {
    }

    //! @brief constructor from a single value
    //! @param VALUE common value for all the elements
    Vector3D::Vector3D( const double &VALUE) :
      m_X( VALUE),
      m_Y( VALUE),
      m_Z( VALUE)
    {
    }

    //! @brief constructor from three values
    //! @param X, Y, Z three elements
    Vector3D::Vector3D( const double &X, const double &Y, const double &Z) :
      m_X( X),
      m_Y( Y),
      m_Z( Z)
    {
    }

    //! @brief constructor from a pointer to three elements
    //! @param PTR_DATA a pointer to three elements
    Vector3D::Vector3D( const double *PTR_DATA) :
      m_X( PTR_DATA[ 0]),
      m_Y( PTR_DATA[ 1]),
      m_Z( PTR_DATA[ 2])
    {
    }

    //! @brief constructor from VectorInterface
    //! @param VECTOR vector
    Vector3D::Vector3D( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "given vector has wrong length");
      std::copy( VECTOR.Begin(), VECTOR.End(), &m_X);
    }

    //! copy constructor
    Vector3D *Vector3D::Clone() const
    {
      return new Vector3D( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Vector3D::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns true if all the coordinates are defined
    //! @return boolean true if all coordinates are defined - false otherwise
    bool Vector3D::IsDefined() const
    {
      return util::IsDefined( m_X) && util::IsDefined( m_Y) && util::IsDefined( m_Z);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief set all elements to elements from another vector
    //! @param VECTOR vector, which must contain 3 elements
    //! @return reference to this, after copying vector
    Vector3D &Vector3D::operator =( const VectorConstInterface< double> &VECTOR)
    {
      BCL_Assert( VECTOR.GetSize() == s_Size, "Tried to copy non-3d vector onto vector3d");
      m_X = VECTOR( 0);
      m_Y = VECTOR( 1);
      m_Z = VECTOR( 2);
      return *this;
    }

    //! @brief set all elements to a vector 3d
    //! @param VECTOR the vector to copy into this vector
    //! @return reference to this, after copying vector
    Vector3D &Vector3D::operator =( const Vector3D &VECTOR)
    {
      m_X = VECTOR.m_X;
      m_Y = VECTOR.m_Y;
      m_Z = VECTOR.m_Z;
      return *this;
    }

    //! @brief add a vector to this one
    //! @param VECTOR the vector to add
    //! @return reference to this
    Vector3D &Vector3D::operator +=( const Vector3D &VECTOR)
    {
      m_X += VECTOR.m_X;
      m_Y += VECTOR.m_Y;
      m_Z += VECTOR.m_Z;
      return *this;
    }

    //! @brief add a scalar to this vector
    //! @param SCALAR the scalar to add
    //! @return reference to this
    Vector3D &Vector3D::operator +=( const double &SCALAR)
    {
      m_X += SCALAR;
      m_Y += SCALAR;
      m_Z += SCALAR;
      return *this;
    }

    //! @brief subtract a vector from this one
    //! @param VECTOR the vector to subtract
    //! @return reference to this
    Vector3D &Vector3D::operator -=( const Vector3D &VECTOR)
    {
      m_X -= VECTOR.m_X;
      m_Y -= VECTOR.m_Y;
      m_Z -= VECTOR.m_Z;
      return *this;
    }

    //! @brief subtract a scalar from this vector
    //! @param SCALAR the value to subtract
    //! @return reference to this
    Vector3D &Vector3D::operator -=( const double &SCALAR)
    {
      m_X -= SCALAR;
      m_Y -= SCALAR;
      m_Z -= SCALAR;
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! change all values
    Vector3D &Vector3D::Set( const double &X, const double &Y, const double &Z)
    {
      m_X = X;
      m_Y = Y;
      m_Z = Z;
      return *this;
    }

    //! translate vector with translation
    Vector3D &Vector3D::Translate( const Vector3D &TRANSLATE)
    {
      m_X += TRANSLATE.m_X;
      m_Y += TRANSLATE.m_Y;
      m_Z += TRANSLATE.m_Z;
      return *this;
    }

    //! translate vector with translation
    Vector3D &Vector3D::Translate( const double X, const double Y, const double Z)
    {
      m_X += X;
      m_Y += Y;
      m_Z += Z;
      return *this;
    }

    //! rotate vector with rotation matrix
    Vector3D &Vector3D::Rotate( const math::RotationMatrix3D &ROTATE)
    {
      const double *row_x( ROTATE.m_RotationMatrix3D[ 0]);
      const double *row_y( ROTATE.m_RotationMatrix3D[ 1]);
      const double *row_z( ROTATE.m_RotationMatrix3D[ 2]);

      // this is a simple, unrolled version of operator *= ROTATE.m_RotationMatrix3D
      Set
      (
        m_X * row_x[ 0] + m_Y * row_y[ 0] + m_Z * row_z[ 0],
        m_X * row_x[ 1] + m_Y * row_y[ 1] + m_Z * row_z[ 1],
        m_X * row_x[ 2] + m_Y * row_y[ 2] + m_Z * row_z[ 2]
      );
      return *this;
    }

    //! transform vector with transformation matrix
    Vector3D &Vector3D::Transform( const math::TransformationMatrix3D &TRANSFORM)
    {
      const double *row_x( TRANSFORM.m_TransformationMatrix3D[ 0]);
      const double *row_y( TRANSFORM.m_TransformationMatrix3D[ 1]);
      const double *row_z( TRANSFORM.m_TransformationMatrix3D[ 2]);
      const double *row_t( TRANSFORM.m_TransformationMatrix3D[ 3]); // translation row

      Set
      (
        m_X * row_x[ 0] + m_Y * row_y[ 0] + m_Z * row_z[ 0] + row_t[ 0],
        m_X * row_x[ 1] + m_Y * row_y[ 1] + m_Z * row_z[ 1] + row_t[ 1],
        m_X * row_x[ 2] + m_Y * row_y[ 2] + m_Z * row_z[ 2] + row_t[ 2]
      );
      return *this;
    }

    //! return random translation vector equally distributed in a sphere of RADIUS
    Vector3D &Vector3D::SetRandomTranslation( const double &RADIUS)
    {
      //random angle between 0 and 2pi
      const double phi( 2 * math::g_Pi * random::GetGlobalRandom().Double());

      //random cos angle between -1 and 1
      const double theta( std::acos( random::GetGlobalRandom().Random( double( -1.0), double( 1.0))));

      //random distance
      const double distance( math::Pow( random::GetGlobalRandom().Random( RADIUS * RADIUS * RADIUS), 1.0 / 3.0));

      //product of sin( theta) * distance
      const double sin_theta_times_distance( std::sin( theta) * distance);

      //set the three values
      return Set(
                  std::cos( phi) * sin_theta_times_distance,
                  std::sin( phi) * sin_theta_times_distance,
                  std::cos( theta) * distance
                );
    }

    //! return random translation vector equally distributed in a ellipse of RADII
    Vector3D &Vector3D::SetRandomTranslation( const Vector3D &RADII)
    {
      // rejection sampling method. More efficient than parametric method due to absence of trig function calls
      // On average the loop terminates after 2 pi / 3 = 2.28 iterations (volume of cube / volume of ellipsoid)
      Vector3D radii_sq( math::Sqr( RADII.X()), math::Sqr( RADII.Y()), math::Sqr( RADII.Z()));
      if( RADII.X() && RADII.Y() && RADII.Z())
      {
        do
        {
          m_X = random::GetGlobalRandom().Double( math::Range< double>( -RADII.X(), RADII.X()));
          m_Y = random::GetGlobalRandom().Double( math::Range< double>( -RADII.Y(), RADII.Y()));
          m_Z = random::GetGlobalRandom().Double( math::Range< double>( -RADII.Z(), RADII.Z()));
        } while
          (
            math::Sqr( m_X) / radii_sq.X()
            + math::Sqr( m_Y) / radii_sq.Y()
            + math::Sqr( m_Z) / radii_sq.Z()
            > 1.0
          );
      }
      else
      {
        radii_sq.X() = std::max( radii_sq.X(), 1e-38);
        radii_sq.Y() = std::max( radii_sq.Y(), 1e-38);
        radii_sq.Z() = std::max( radii_sq.Z(), 1e-38);
        do
        {
          m_X = RADII.X() ? random::GetGlobalRandom().Double( math::Range< double>( -RADII.X(), RADII.X())) : 0.0;
          m_Y = RADII.Y() ? random::GetGlobalRandom().Double( math::Range< double>( -RADII.Y(), RADII.Y())) : 0.0;
          m_Z = RADII.Z() ? random::GetGlobalRandom().Double( math::Range< double>( -RADII.Z(), RADII.Z())) : 0.0;
        } while
          (
            math::Sqr( m_X) / radii_sq.X()
            + math::Sqr( m_Y) / radii_sq.Y()
            + math::Sqr( m_Z) / radii_sq.Z()
            > 1.0
          );
      }

      //set the three values
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! print in a user readable format
    //! @return the string containing the information
    std::string Vector3D::ToString() const
    {
      return std::string( util::Format()( m_X) + " " + util::Format()( m_Y) + " " + util::Format()( m_Z));
    }

    //! @brief normalize the vector (such that inner product == 1)
    Vector3D &Vector3D::Normalize()
    {
      const double sqnorm( SquareNorm());
      if( sqnorm != 1.0)
      {
        const double norm( math::Sqrt( sqnorm));
        m_X /= norm;
        m_Y /= norm;
        m_Z /= norm;
      }
      return *this;
    }

    //! @brief normalize the vector (such that inner product == 1)
    //! Overrides function of the same name in VectorInterface for performance reasons
    Vector3D Vector3D::Normalized() const
    {
      Vector3D c( *this);
      return c.Normalize();
    }

    //! @brief norm = length of vector
    //! @return length of vector
    double Vector3D::Norm() const
    {
      return math::Sqrt( SquareNorm());
    }

    //! @brief square norm = square length of vector
    //! @return square length of vector
    double Vector3D::SquareNorm() const
    {
      return m_X * m_X + m_Y * m_Y + m_Z * m_Z;
    }

    //! @brief sum up all elements
    //! @return sum of all elements
    double Vector3D::Sum() const
    {
      return m_X + m_Y + m_Z;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Vector3D::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_X, ISTREAM);
      io::Serialize::Read( m_Y, ISTREAM);
      io::Serialize::Read( m_Z, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Vector3D::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_X, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_Y, OSTREAM, 0) << '\t';
      io::Serialize::Write( m_Z, OSTREAM, 0);

      // return the stream
      return OSTREAM;
    }

  } // namespace linal
} // namespace bcl
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
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_message.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    //! @brief calculates the absolute difference between two linal::Vector3Ds
    //! @param VECTOR_A first vector
    //! @param VECTOR_B second vector
    //! @return absolute difference between two linal::Vector3D (points)
    Vector3D AbsoluteDifference( const Vector3D &VECTOR_A, const Vector3D &VECTOR_B)
    {
      Vector3D returned_vector( VECTOR_B - VECTOR_A);
      returned_vector.X() = math::Absolute( returned_vector.X());
      returned_vector.Y() = math::Absolute( returned_vector.Y());
      returned_vector.Z() = math::Absolute( returned_vector.Z());
      return returned_vector;
    }

    //! @brief calculates the unit vector starting from one linal::Vector3D to another
    //! @param ORIGIN vector of origin
    //! @param TARGET target vector
    //! @return the unit vector between ORIGIN and TARGET
    Vector3D UnitVector( const Vector3D &ORIGIN, const Vector3D &TARGET)
    {
      Vector3D unit_vector( TARGET);
      unit_vector -= ORIGIN;
      unit_vector.Normalize();
      return unit_vector;
    }

    //! @brief dihedral angle between four points (A->B -x-> C->D)
    //! @brief see http://en.wikipedia.org/wiki/Dihedral_angle for reference
    //! @param VECTOR_A first vector (point)
    //! @param VECTOR_B second vector (point)
    //! @param VECTOR_C third vector (point)
    //! @param VECTOR_D fourth vector (point)
    //! @return dihedral angle between four points
    double Dihedral
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const Vector3D &VECTOR_D
    )
    {
      // calculate the two cross products (b1xb2 and b2xb3)
      const Vector3D cross_b1_b2( CrossProduct( VECTOR_A, VECTOR_B, VECTOR_B, VECTOR_C));
      const Vector3D cross_b2_b3( CrossProduct( VECTOR_B, VECTOR_C, VECTOR_C, VECTOR_D));

      // calculate the vectors b1 and b2
      const Vector3D b1( VECTOR_B - VECTOR_A);

      // get the distance b -c
      const double distance( Distance( VECTOR_C, VECTOR_B));

      // calculate dihedral angle
      const double dihedral( std::atan2( distance * b1 * cross_b2_b3, cross_b1_b2 * cross_b2_b3));

      return util::IsDefined( dihedral) ? dihedral : 0.0;
    }

    //! @brief calculates coordinates using dihedral angle information (point X in C->B->A->X)
    //! @param VECTOR_A first point
    //! @param VECTOR_B second point
    //! @param VECTOR_C third point
    //! @param DISTANCE_XA distance between X and VECTOR_A
    //! @param ANGLE_XAB angle between X, VECTOR_A, and VECTOR_B
    //! @param DIHEDRAL_XABC dihedral angle between X, VECTOR_A, VECTOR_B, and VECTOR_C
    //! @return point X in C->B->A->X
    Vector3D CoordinatesDihedral
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const double DISTANCE_XA,
      const double ANGLE_XAB,
      const double DIHEDRAL_XABC
    )
    {
      const Vector3D a( UnitVector( VECTOR_B, VECTOR_A));
      const Vector3D b( UnitVector( VECTOR_C, VECTOR_B));
      const Vector3D c( CrossProduct( a, b).Normalize());
      const Vector3D d( CrossProduct( a, c).Normalize());
      const Vector3D x
      (
        DISTANCE_XA *
        (
          a * cos( math::g_Pi - ANGLE_XAB) - c * sin( math::g_Pi - ANGLE_XAB) * sin( DIHEDRAL_XABC) +
          d * sin( math::g_Pi - ANGLE_XAB) * cos( DIHEDRAL_XABC)
        )
      );

      return VECTOR_A + x;
    }

    //! @brief calculates coordinates using angle information (point X in B,C->A->X)
    //! @param VECTOR_A first point
    //! @param VECTOR_B second point
    //! @param VECTOR_C third point
    //! @param DISTANCE_XA distance between X and VECTOR_A
    //! @param ANGLE_XAB angle between X, VECTOR_A, and VECTOR_B
    //! @param ANGLE_XAC angle between X, VECTOR_A, and VECTOR_C
    //! @param SIDE true if on a side
    //! @param VECTOR_SIDE vector of the side
    //! @return point X in B,C->A->X
    Vector3D CoordinatesAngle
    (
      const Vector3D &VECTOR_A,
      const Vector3D &VECTOR_B,
      const Vector3D &VECTOR_C,
      const double DISTANCE_XA,
      const double ANGLE_XAB,
      const double ANGLE_XAC,
      const bool SIDE,
      const Vector3D &VECTOR_SIDE
    )
    {
      // in plane components and direction
      const Vector3D a( UnitVector( VECTOR_B, VECTOR_A));
      Vector3D b( UnitVector( VECTOR_C, VECTOR_A));

      // determine whether a and b are colinear. In which case, use a different vector C
      while( math::Absolute( ProjAngleCosinus( a, b)) >= 0.99)
      {
        BCL_MessageDbg( "Lines were colinear, CoordinatesAngle returning non-unique vector");
        b.Rotate( math::RotationMatrix3D().SetRand( math::g_Pi));
      }

      const Vector3D c( cos( math::g_Pi - ANGLE_XAB) * a);
      const Vector3D d( cos( math::g_Pi - ANGLE_XAC) * b);

      const double bac_angl( ProjAngle( c, d));
      const double a_dist( c.Norm());
      const double b_dist( d.Norm());
      const double c_dist
      (
        math::Sqrt
        (
          std::max
          (
            0.0,
            ( math::Sqr( a_dist) + math::Sqr( b_dist) - 2 * a_dist * b_dist * cos( bac_angl))
            / math::Sqr( sin( bac_angl)) -
            math::Sqr( a_dist)
          )
        )
      );

      //in a,b plane
      const Vector3D e( CrossProduct( c, d).Normalize());
      const Vector3D f( CrossProduct( e, c).Normalize() * c_dist);
      const Vector3D g( c + f);

      //orthogonal to a,b plane
      Vector3D h( e * math::Sqrt( std::max( 0.0, 1.0 - g.SquareNorm())));

      //side
      if( util::IsDefined( VECTOR_SIDE( 0)))
      {
        const Vector3D i( VECTOR_SIDE - VECTOR_A);
        const double dihe_1( Dihedral( i, Vector3D(), a, b));
        const double dihe_2( Dihedral( e, Vector3D(), a, b));
        if( dihe_2 && ( ( dihe_1 / dihe_2 > 0.0) != SIDE))
        {
          h *= ( -1.0);
        }
      }

      //sum and length
      const Vector3D x( Vector3D( g + h).Normalize() * DISTANCE_XA);

      return VECTOR_A + x;
    }

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_vector_const_interface.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API VectorConstInterface< double>;
    template class BCL_API VectorConstInterface< float>;
    template class BCL_API VectorConstInterface< int>;
    template class BCL_API VectorConstInterface< unsigned int>;
    template class BCL_API VectorConstInterface< unsigned long>;
    template class BCL_API VectorConstInterface< unsigned long long>;
    template class BCL_API VectorConstInterface< char>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_vector.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Vector< double>;
    template class BCL_API Vector< float>;
    template class BCL_API Vector< int>;
    // these are necessary to ensure that we get both the cl_uint type, which is used in
    // some opencl classes, as well as size_t.  It is possible that cl_uint and size_t are the same
    // type on some platforms.
    template class BCL_API Vector< unsigned int>;
    template class BCL_API Vector< unsigned long>;
    template class BCL_API Vector< unsigned long long>;
    template class BCL_API Vector< char>;

  } // namespace linal
} // namespace bcl
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
#include "linal/bcl_linal_vector_interface.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API VectorInterface< double>;
    template class BCL_API VectorInterface< float>;
    template class BCL_API VectorInterface< int>;
    template class BCL_API VectorInterface< unsigned int>;
    template class BCL_API VectorInterface< unsigned long>;
    template class BCL_API VectorInterface< unsigned long long>;
    template class BCL_API VectorInterface< char>;

  } // namespace linal
} // namespace bcl
