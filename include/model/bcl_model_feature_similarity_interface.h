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

#ifndef BCL_MODEL_FEATURE_SIMILARITY_INTERFACE_H_
#define BCL_MODEL_FEATURE_SIMILARITY_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_const_interface.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureSimilarityMeasuresInterface
    //! @brief interface for similarity measures
    //!
    //! @tparam t_DataType can be float, double, int
    //!
    //! @author loweew
    //! @date Oct 11, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class FeatureSimilarityMeasuresInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone the FeatureSimilarityMeasuresInterface
      //! @return pointer to new FeatureSimilarityMeasuresInterface
      virtual FeatureSimilarityMeasuresInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief operator makes similarity measurement
      //! @param INPUT the input matrix to compare against the input matrix from the constructor
      //! @return the matrix of similarity coefficients
      virtual linal::Matrix< t_DataType> operator()( const linal::MatrixConstInterface< t_DataType> &INPUT) const = 0;

    ///////////////
    // operators //
    ///////////////

    }; // template class FeatureSimilarityMeasuresInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_FEATURE_SIMILARITY_INTERFACE_H_
