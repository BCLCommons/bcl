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

#ifndef BCL_CONTACT_AA_CORRELATION_FROM_FILE_H_
#define BCL_CONTACT_AA_CORRELATION_FROM_FILE_H_

// include the namespace header
#include "bcl_contact.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_contact_correlation_matrix.h"
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "descriptor/bcl_descriptor_base_pair.h"
#include "descriptor/bcl_descriptor_type.h"
#include "score/bcl_score.h"
#include "score/bcl_score_aa_assignments.h"
#include "score/bcl_score_assignment_with_gap.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AACorrelationFromFile
    //! @brief AACorrelationFromFiles the output of descriptors into one vector
    //!
    //! @see @link example_contact_aa_correlation_from_file.cpp @endlink
    //! @author teixeipl, mendenjl
    //! @date Feb 11, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AACorrelationFromFile :
      public descriptor::BasePair< biol::AABase, float>
    {
    private:

    //////////
    // data //
    //////////
      //! AAAssignment score
      score::AAAssignment m_Score;

      //! Vector provides position indexing from the target sequence to the correlation matrix based on the alignment
      storage::Map< util::SiPtr< const biol::AABase>, size_t, biol::AALessThanSeqID> m_SeqToMatrixColMap;

      //! Correlation matrix storing data from file
      CorrelationMatrix m_CorrelationMatrix;

      //! Pointer to assignment gap
      score::AssignmentWithGap< biol::AABase> m_Assignment;

      //! Correlation matrix file suffix if explicitly specified
      std::string m_FileSuffix;

      //! Extension for correlation matrix BCL file that contains data
      std::string m_FileExtension;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor that takes a filename string
      AACorrelationFromFile( const std::string &SUFFIX = "");

      //! @brief virtual copy constructor
      AACorrelationFromFile *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT_A,
        const iterate::Generic< const biol::AABase> &ELEMENT_B,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      descriptor::CachePreference GetNormalCachePreference() const;

      //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
      //! since the results are often in the cache
      void LoadFiles();

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

      //! @brief adds the given pair score with the given weight to the scoring function
      //! @param SCORE the pair score to add
      //! @param SCORE_WEIGHT the weight for the pair score
      //! @param SCORE_FCT the scoring function (will be modified)
      static void AddScore
      (
        const score::AAAssignment &SCORE,
        double SCORE_WEIGHT,
        function::BinarySum< const biol::AABase, const biol::AABase, double> &SCORE_FCT
      );

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    }; // class AACorrelationFromFile

  } // namespace contact
} // namespace bcl

#endif // BCL_CONTACT_AA_CORRELATION_FROM_FILE_H_
