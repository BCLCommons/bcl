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

#ifndef BCL_CLUSTER_DISTANCES_EUCLIDEAN_H_
#define BCL_CLUSTER_DISTANCES_EUCLIDEAN_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DistancesEuclidean
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_cluster_distances_euclidean.cpp @endlink
    //! @author alexanns
    //! @date Jun 2, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_PrecisionType>
    class DistancesEuclidean :
      public math::FunctionInterface
      <
        storage::VectorND< 2, util::SiPtr< const linal::Vector< t_PrecisionType> > >, t_PrecisionType
      >
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
      DistancesEuclidean()
      {
      }

      //! @brief Clone function
      //! @return pointer to new DistancesFeatureVectors
      DistancesEuclidean *Clone() const
      {
        return new DistancesEuclidean( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() which takes two t_DataTypes and returns the distance between them
      //! @param OBJECTS is a VectorND which has the two objects whose distance is needed
      //! @return returns a t_PrecisionType which is the distance between the two objects
      t_PrecisionType operator()
      (
        const storage::VectorND< 2, util::SiPtr< const linal::Vector< t_PrecisionType> > > &OBJECTS
      ) const
      {
        const linal::VectorConstInterface< t_PrecisionType> &features_a( *OBJECTS.First());
        const linal::VectorConstInterface< t_PrecisionType> &features_b( *OBJECTS.Second());

        return linal::Distance( features_a, features_b);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class DistancesEuclidean

    // instantiate s_Instance
    template< typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> DistancesEuclidean< t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new DistancesEuclidean< t_PrecisionType>())
    );

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_DISTANCES_EUCLIDEAN_H_
