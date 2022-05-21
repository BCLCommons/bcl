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

#ifndef BCL_BIOL_ALIGN_BY_AA_DATA_H_
#define BCL_BIOL_ALIGN_BY_AA_DATA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_base.h"
#include "align/bcl_align_pairwise_aligner_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignByAAData
    //! @brief implements align::EngineInterface for aligning amino acids in a model to a full sequence/model by using AAData
    //! @details for two given alignments, assuming one of them is a reference alignment, meaning it includes all the
    //! amino acids in the other alignment, it creates an alignment of these two alignments by matching the assignments
    //! by addresses of AAData. This class is useful when aligning a model to the native structure.
    //!
    //! @see @link example_biol_align_by_aa_data.cpp @endlink
    //! @author karakam
    //! @date Jan 13, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AlignByAAData :
      public align::PairwiseAlignerInterface< AABase>
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
      AlignByAAData();

      //! @brief Clone function
      //! @return pointer to new AlignByAAData
      AlignByAAData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief SetScoringFunction sets the score::Assignment scoring function
      //! @param SCORE is the score::Assignment to be used for scoring an assignment
      void SetScoringFunction( const score::AssignmentWithGap< AABase> &SCORE)
      {
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
      //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
      //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
      storage::Pair< align::AlignmentNode< AABase>, double> AlignPair
      (
        util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
        util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B
      ) const;

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
      //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
      //! @param COMPARISON method for comparing the aa data to see if it is equal or not
      //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
      storage::Pair< align::AlignmentNode< AABase>, double> AlignPair
      (
        util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
        util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B,
        const util::BinaryFunctionInterface< AABase, AABase, bool> &COMPARISON
      ) const;

      //! @brief Align is the function which does the actually aligns data together
      //! @param NODE is the alignment node, which must alread contain the two parent alignment interfaces
      //! @param COMPARISON method for comparing the aa data to see if it is equal or not
      //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
      double AlignPairWithNode
      (
        align::AlignmentInterface< AABase> &NODE,
        const util::BinaryFunctionInterface< AABase, AABase, bool> &COMPARISON
      ) const;

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
      //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
      //! @param TEMPLATE_ALIGNMENT is the template alignment, ALIGNMENT_A and ALIGNMENT_B will be aligned according to this
      //!        alignment - the first sequence must correspond to ALIGNMENT_A and the second to ALIGNMENT_B
      //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
      storage::Pair< align::AlignmentNode< AABase>, double> AlignPair
      (
        util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
        util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B,
        const util::ShPtr< align::AlignmentInterface< AABase> > &TEMPLATE_ALIGNMENT
      ) const;

      //! @brief Align is the function which does the actually aligns data together
      //! @param NODE is the alignment node, which must alread contain the two parent alignment interfaces
      //! @param COMPARISON method for comparing the aa data to see if it is equal or not
      //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
      double AlignPairWithNode
      (
        align::AlignmentInterface< AABase> &NODE,
        const util::ShPtr< align::AlignmentInterface< AABase> > &TEMPLATE_ALIGNMENT
      ) const;

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

    }; // class AlignByAAData

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_ALIGN_BY_AA_DATA_H_
