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

#ifndef BCL_SCORE_AA_ASSIGNMENT_MEAN_SIMILARITY_MATRIX_H_
#define BCL_SCORE_AA_ASSIGNMENT_MEAN_SIMILARITY_MATRIX_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "math/bcl_math_function_interface_serializable.h"

namespace bcl
{
  namespace score
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAAssignmentMeanSimilarityMatrix
    //! @brief TODO: Similarity matrix to be used by the Pearson-based correlated mutation calculation
    //! @details TODO: This similarity matrix is the mean matrix produced in (Lena 2011) after a training period to determine the optimal matrix
    //!
    //! @see @link example_score_AAAssignmentMeanSimilarityMatrix.cpp @endlink
    //! @author teixeipl
    //! @date Jul 30, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAAssignmentMeanSimilarityMatrix :
      public math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const biol::AABase> >, double>
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::ObjectInterface *s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AAAssignmentMeanSimilarityMatrix();

      //! @brief Clone function
      //! @return pointer to new AAAssignmentMeanSimilarityMatrix
      AAAssignmentMeanSimilarityMatrix *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class AAAssignmentMeanSimilarityMatrix

  } // namespace score
  
} // namespace bcl

#endif // BCL_SCORE_AA_ASSIGNMENT_MEAN_SIMILARITY_MATRIX_H_
