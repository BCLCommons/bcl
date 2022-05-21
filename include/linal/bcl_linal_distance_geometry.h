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

#ifndef BCL_LINAL_DISTANCE_GEOMETRY_H_
#define BCL_LINAL_DISTANCE_GEOMETRY_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_vector.h"
#include "bcl_linal_vector_nd.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_symmetric_matrix.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DistanceGeometry
    //! @brief Uses function that performs distance geometry calculation and returns a vector of points (3 dimensional)
    //! @details This method is based on the following: Distance Geometry of Molecules, A Sara, B Elena, G Simone, I
    //!          Samuli, P Chara, Sabak Anna, Weenink Jan Willem, Instructor: Simon Kokkendorff,
    //!          14th ECMI modelling week Lund 2000,
    //!
    //! @see @link example_linal_distance_geometry.cpp @endlink
    //! @author teixeipl
    //! @date Jun 8, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DistanceGeometry :
      public math::FunctionInterfaceSerializable< storage::SymmetricMatrix< double>, storage::Vector< VectorND< double, 3> > >
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DistanceGeometry()
      {
      }

      //! @brief Clone function
      //! @return pointer to new DistanceGeometry
      DistanceGeometry *Clone() const
      {
        return new DistanceGeometry( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        // Get BCL Standardized Class name
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
      //! @param SYMMETRIC_MATRIX containing pairwise distances for target
      //! @return function value of the given argument
      virtual storage::Vector< VectorND< double, 3> > operator()
      (
        const storage::SymmetricMatrix< double> &SYMMETRIC_MATRIX
      ) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class DistanceGeometry

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_DISTANCE_GEOMETRY_H_
