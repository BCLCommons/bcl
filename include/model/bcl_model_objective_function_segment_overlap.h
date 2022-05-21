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

#ifndef BCL_MODEL_OBJECTIVE_FUNCTION_SEGMENT_OVERLAP_H_
#define BCL_MODEL_OBJECTIVE_FUNCTION_SEGMENT_OVERLAP_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_objective_function_categorical_max.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ObjectiveFunctionSegmentOverlap
    //! @brief Segment overlap is a common method for benchmarking secondary structure prediction methods
    //! The version implemented here is the updated definition from Zemla et al. - PROTEINS: Structure, Function, and
    //! Genetics, 34, 1999, pp. 220-223.
    //!
    //! @see @link example_model_objective_function_segment_overlap.cpp @endlink
    //! @author mendenjl
    //! @date Mar 18, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ObjectiveFunctionSegmentOverlap :
      public ObjectiveFunctionCategoricalMax
    {

    private:

    //////////
    // data //
    //////////

      //! bool to output the segment overlaps for each protein
      bool m_PerIdOutput;

      //! bool to output sub-classification segment overlaps, e.g. for secondary structure, SOV-H, SOV-E, SOV-C
      bool m_OutputSubclassOverlaps;

      //! bool; whether to weight by per-element (true) or per-sequence (false)
      bool m_WeightPerElement;

      //! Sequence Identification label (Normally ProteinId)
      util::ObjectDataLabel m_SequenceIdLabel;

      //! Element Identification label (Normally AASeqID)
      util::ObjectDataLabel m_ElementIdLabel;

      //! Positions of missing id elements (first missing id in any contiguous block)
      //! Missing ids always terminate the current segment
      storage::Vector< storage::Vector< size_t> > m_MissingIds;

      //! Boundary ids of all sequences. 1 number for each sequence (+1 for the end of the dataset)
      storage::Vector< size_t> m_BoundaryIds;

      //! Counts of each class, Sequences X # result cols
      linal::Matrix< size_t> m_ClassCounts;

      //! names of each sequence
      storage::Vector< std::string> m_SequenceNames;

      //! Sizes of each sequence
      storage::Vector< size_t> m_SequenceSizes;

      //! Typedef for the container used for the segments of a given sequence
      //! Triplet is ordered as SubClass, Start, End
      typedef storage::Vector< storage::Triplet< size_t, size_t, size_t> > SegmentContainer;

      //! Experimental segments; Sequences X Classes X Segments X (Class, Start, End)
      storage::Vector< storage::Vector< SegmentContainer> > m_ActualSegments;

      //! Allowed permutations of result classes. If an alternate permutation is used, it will be noted in the output
      storage::Vector< storage::Vector< size_t> > m_AllowedPermutations;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ObjectiveFunctionSegmentOverlap();

      //! @brief Clone function
      //! @return pointer to new Evaluator
      ObjectiveFunctionSegmentOverlap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_Name( "SegmentOverlap");
        return s_Name;
      }

      //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
      //! @param DATA monitoring dataset results, non-scaled
      //! @param IDS ids; can be used by the objective function
      void SetData
      (
        const FeatureDataSet< float> &DATA,
        const FeatureDataSet< char> &IDS = FeatureDataSet< char>()
      );

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
      //! @param EXPERIMENTAL feature dataset with experimental values
      //! @param PREDICTED feature dataset with predicted values
      //! @return objective function value based on the given data of experimental and predicted values
      float operator()
      (
        const FeatureDataSetInterface< float> &EXPERIMENTAL,
        const FeatureDataSetInterface< float> &PREDICTED
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief get the segments for a given subclass and sequence
      //! @param SUBCLASS subclass index; if just computing secondary structure using on method, this is always 0; but
      //!        when using multiple methods, it may range from 0-number classifications
      //! @param SEQUENCE sequence index
      //! @param CHOSEN_SUBCLASS_MATRIX matrix of chosen sub-classes
      //! @return vector of segments for the given sequence for the given subclass
      SegmentContainer GetSegments
      (
        const size_t &SUBCLASS,
        const size_t &SEQUENCE,
        const linal::MatrixConstInterface< size_t> &CHOSEN_CLASSES,
        const size_t &PERMUTATION = size_t( 0)
      ) const;

      //! @brief compute segment overlap from predicted to experimental
      //! @param SUBCLASS subclass index; if just computing secondary structure using on method, this is always 0; but
      //!        when using multiple methods, it may range from 0-number classifications
      //! @param SEQUENCE_ID sequence index
      //! @param RESTRICT_CLASS optional value; if set, only consider sub
      //! @param RESTRICT_PERMUTATION optional value; if set, only consider the indicated permutation
      //! @return the segment overlap and the permutation index used
      storage::Pair< float, size_t> GetSegmentOverlap
      (
        const size_t &SUBCLASS,
        const size_t &SEQUENCE,
        const size_t &RESTRICT_CLASS = util::GetUndefined< size_t>(),
        const size_t &RESTRICT_PERMUTATION = util::GetUndefined< size_t>()
      ) const;

      //! @brief create the segment string, useful for checking results
      //! @param SUBCLASS subclass index; if just computing secondary structure using on method, this is always 0; but
      //!        when using multiple methods, it may range from 0-number classifications
      //! @param SEQUENCE_ID sequence index
      //! @return string with format: XXXXX ID,  PSEC,  OSEC
      std::string GetSegmentString( const size_t &SUBCLASS, const size_t &SEQUENCE) const;

    };

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_OBJECTIVE_FUNCTION_SEGMENT_OVERLAP_H_
