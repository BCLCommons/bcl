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
#include "biol/bcl_biol_align_by_pdb_id.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AlignByPdbID::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignByPdbID())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AlignByPdbID::AlignByPdbID()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AlignByPdbID
    AlignByPdbID *AlignByPdbID::Clone() const
    {
      return new AlignByPdbID( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AlignByPdbID::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
    //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByPdbID::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B
    ) const
    {
      // make sure both have depth of at least 1
      if( ALIGNMENT_A->GetDepth() == 0 || ALIGNMENT_B->GetDepth() == 0)
      {
        BCL_MessageCrt( "The provided alignments should have depths of at least 1!");
        return storage::Pair< align::AlignmentNode< AABase>, double>( align::AlignmentNode< AABase>(), util::GetUndefinedDouble());
      }

      // initialize return type
      align::AlignmentNode< AABase> alignment( ALIGNMENT_A, ALIGNMENT_B);
      storage::Pair< align::AlignmentNode< AABase>, double> alignment_and_score( alignment, util::GetUndefinedDouble());

      // initialize iterators
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a( ALIGNMENT_A->GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a_end( ALIGNMENT_A->GetAssignments().End());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b( ALIGNMENT_B->GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b_end( ALIGNMENT_B->GetAssignments().End());

      // initialize empty pointer to be used for gaps
      const util::SiPtr< const AABase> empty_ptr;

      // while there are still elements left in ALIGNMENT_B
      while( itr_a != itr_a_end && itr_b != itr_b_end)
      {
        util::SiPtr< const AABase> amino_acid_a( ( *itr_a)->GetMembers().FirstElement());
        util::SiPtr< const AABase> amino_acid_b( ( *itr_b)->GetMembers().FirstElement());

        // match the chains
        if( amino_acid_a->GetChainID() == amino_acid_b->GetChainID())
        {
          // match pdb id
          if( amino_acid_a->GetPdbID() == amino_acid_b->GetPdbID())
          {
            // match pdb icode
            if( amino_acid_a->GetPdbICode() == amino_acid_b->GetPdbICode())
            {
              // move both iterators
              ++itr_a;
              ++itr_b;
            }
            // pdb i code smaller
            else if( amino_acid_a->GetPdbICode() < amino_acid_b->GetPdbICode())
            {
              ++itr_a;
              amino_acid_b = empty_ptr;
            }
            // pdb i code larger
            else
            {
              ++itr_b;
              amino_acid_a = empty_ptr;
            }
          }
          // pdb id smaller
          else if( amino_acid_a->GetPdbID() < amino_acid_b->GetPdbID())
          {
            ++itr_a;
            amino_acid_b = empty_ptr;
          }
          // pdb id larger
          else
          {
            ++itr_b;
            amino_acid_a = empty_ptr;
          }
        }
        // chain id smaller
        else if( amino_acid_a->GetChainID() < amino_acid_b->GetChainID())
        {
          ++itr_a;
          amino_acid_b = empty_ptr;
        }
        // chain id larger
        else
        {
          ++itr_b;
          amino_acid_a = empty_ptr;
        }

        // insert assignment of itr_a and itr_b
        alignment_and_score.First().Append
        (
          util::ShPtr< align::Assignment< AABase> >( new align::Assignment< AABase>( amino_acid_a, amino_acid_b))
        );
      } // end of while loop

      // iterate while itr_a reaches end
      while( itr_a != itr_a_end)
      {
        // insert assignment of itr_a with gap
        alignment_and_score.First().Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( ( *itr_a)->GetMembers().FirstElement(), empty_ptr)
          )
        );

        // move the next residue
        ++itr_a;
      }

      // iterate while itr_b reaches end
      while( itr_b != itr_b_end)
      {
        // insert assignment of itr_b with gap
        alignment_and_score.First().Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( empty_ptr, ( *itr_b)->GetMembers().FirstElement())
          )
        );

        // move the next residue
        ++itr_b;
      }

      // end
      return alignment_and_score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AlignByPdbID::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AlignByPdbID::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
