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
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_align_by_aa_data.h"

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
    const util::SiPtr< const util::ObjectInterface> AlignByAAData::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignByAAData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AlignByAAData::AlignByAAData()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AlignByAAData
    AlignByAAData *AlignByAAData::Clone() const
    {
      return new AlignByAAData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AlignByAAData::GetClassIdentifier() const
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
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByAAData::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B
    ) const
    {
      // end
      return AlignPair( ALIGNMENT_A, ALIGNMENT_B, AACompareDataPtr());
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
    //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
    //! @param COMPARISON method for comparing the aa data to see if it is equal or not
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByAAData::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B,
      const util::BinaryFunctionInterface< AABase, AABase, bool> &COMPARISON
    ) const
    {
      // initialize return type
      align::AlignmentNode< AABase> node( ALIGNMENT_A, ALIGNMENT_B);
      const double score( AlignPairWithNode( node, COMPARISON));
      return storage::Pair< align::AlignmentNode< AABase>, double>( node, score);
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param NODE is the alignment node, which must alread contain the two parent alignment interfaces
    //! @param COMPARISON method for comparing the aa data to see if it is equal or not
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    double AlignByAAData::AlignPairWithNode
    (
      align::AlignmentInterface< AABase> &NODE,
      const util::BinaryFunctionInterface< AABase, AABase, bool> &COMPARISON
    ) const
    {
      iterate::Generic< const align::AlignmentInterface< AABase> >
        itr_child_alignments( NODE.GetChildAlignmentsIterator());
      BCL_Assert
      (
        itr_child_alignments.GetSize() == size_t( 2),
        "Expected two child alignments already in the node"
      );
      const align::AlignmentInterface< AABase> &alignment_first( *itr_child_alignments);
      const align::AlignmentInterface< AABase> &alignment_second( *++itr_child_alignments);

      const align::AlignmentInterface< AABase> &alignment_a
      (
        alignment_first.GetSize() >= alignment_second.GetSize()
        ? alignment_first
        : alignment_second
      );
      const align::AlignmentInterface< AABase> &alignment_b
      (
        alignment_first.GetSize() >= alignment_second.GetSize()
        ? alignment_second
        : alignment_first
      );

      // make sure both have depth of at least 1
      if( alignment_a.GetDepth() == 0 || alignment_b.GetDepth() == 0)
      {
        BCL_MessageCrt( "The provided alignments should have depths of at least 1!");
        return util::GetUndefinedDouble();
      }
      double score( util::GetUndefinedDouble());

      // initialize iterators
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a( alignment_a.GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a_end( alignment_a.GetAssignments().End());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b( alignment_b.GetAssignments().Begin());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b_end( alignment_b.GetAssignments().End());

      // initialize empty pointer to be used for gaps
      const util::SiPtr< const AABase> empty_ptr;

      // while there are still elements left in ALIGNMENT_B
      while( itr_b != itr_b_end)
      {
        // check that a is not at the end yet
        if( itr_a == itr_a_end)
        {
          BCL_MessageCrt( "The provided alignment a should have all aadata of b!");
          NODE.ResetAssignments();
          score = util::GetUndefined< double>();
          return score;
        }

        const align::Assignment< AABase> &assignment_a( **itr_a), &assignment_b( **itr_b);
        BCL_MessageVrb( "Trying to align itr_a=" + assignment_a.ToString() + " and itr_b=" + assignment_b.ToString());

        // iterate until the AAData pointer by both pointers are the same
        if( !COMPARISON( *( assignment_a.GetMembers().FirstElement()), *( assignment_b.GetMembers().FirstElement())))
        {
          BCL_MessageDbg( "Assignment comparison returned false, not aligning");

          // insert assignment of itr_a with gap
          NODE.Append
          (
            util::ShPtr< align::Assignment< AABase> >
            (
              new align::Assignment< AABase>( assignment_a.GetMembers().FirstElement(), empty_ptr)
            )
          );

          // move to next amino acid
          ++itr_a;
        }
        else // aa data pointer agrees
        {
          BCL_MessageDbg( "Assignment comparison returned true, aligning");

          // insert assignment of itr_a and itr_b
          NODE.Append
          (
            util::ShPtr< align::Assignment< AABase> >
            (
              new align::Assignment< AABase>
              (
                assignment_a.GetMembers().FirstElement(), assignment_b.GetMembers().FirstElement()
              )
            )
          );

          // move both iterators
          ++itr_a;
          ++itr_b;
        }
      } // end of while loop

      // iterate while itr_a reaches end
      while( itr_a != itr_a_end)
      {
        // insert assignment of itr_a with gap
        NODE.Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( ( *itr_a)->GetMembers().FirstElement(), empty_ptr)
          )
        );

        // move the next residue
        ++itr_a;
      }

      // end
      return score;
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A is the Alignment to be aligned with ALIGNMENT_B
    //! @param ALIGNMENT_B is the Alignment to be aligned with ALIGNMENT_A
    //! @param TEMPLATE_ALIGNMENT is the template alignment, ALIGNMENT_A and ALIGNMENT_B will be aligned according to this
    //!        alignment - the first sequence must correspond to ALIGNMENT_A and the second to ALIGNMENT_B
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    storage::Pair< align::AlignmentNode< AABase>, double> AlignByAAData::AlignPair
    (
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_A,
      util::ShPtr< align::AlignmentInterface< AABase> > &ALIGNMENT_B,
      const util::ShPtr< align::AlignmentInterface< AABase> > &TEMPLATE_ALIGNMENT
    ) const
    {
      // initialize return type
      align::AlignmentNode< AABase> alignment( ALIGNMENT_A, ALIGNMENT_B);
      const double score( AlignPairWithNode( alignment, TEMPLATE_ALIGNMENT));
      return storage::Pair< align::AlignmentNode< AABase>, double>( alignment, score);
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param NODE is the alignment node, which must alread contain the two parent alignment interfaces
    //! @param COMPARISON method for comparing the aa data to see if it is equal or not
    //! @return returns a pair of the Alignment of ALIGNMENT_A and ALIGNMENT_B and a double which is the score
    double AlignByAAData::AlignPairWithNode
    (
      align::AlignmentInterface< AABase> &NODE,
      const util::ShPtr< align::AlignmentInterface< AABase> > &TEMPLATE_ALIGNMENT
    ) const
    {
      iterate::Generic< const align::AlignmentInterface< AABase> >
        itr_child_alignments( NODE.GetChildAlignmentsIterator());
      BCL_Assert
      (
        itr_child_alignments.GetSize() == size_t( 2),
        "Expected two child alignments already in the node"
      );
      const align::AlignmentInterface< AABase> &alignment_a( *itr_child_alignments);
      const align::AlignmentInterface< AABase> &alignment_b( *++itr_child_alignments);

      double score( util::GetUndefinedDouble());

      // initialize iterators
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a( alignment_a.GetAssignments().Begin());
      const util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_a_end( alignment_a.GetAssignments().End());
      util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b( alignment_b.GetAssignments().Begin());
      const util::ShPtrList< align::Assignment< AABase> >::const_iterator itr_b_end( alignment_b.GetAssignments().End());

      // initialize empty pointer to be used for gaps
      const util::SiPtr< const AABase> empty_ptr;

      // iterate through the template alignment
      for
      (
        util::ShPtrList< align::Assignment< AABase> >::const_iterator
          template_itr( TEMPLATE_ALIGNMENT->GetAssignments().Begin()),
          template_itr_end( TEMPLATE_ALIGNMENT->GetAssignments().End());
        template_itr != template_itr_end; ++template_itr
      )
      {
        // get the AA's
        util::SiPtr< const AABase> sp_aa_a( empty_ptr);
        util::SiPtr< const AABase> sp_aa_b( empty_ptr);
        const util::SiPtr< const AABase> &sp_aa_template_a( ( *template_itr)->GetMembers().FirstElement());
        const util::SiPtr< const AABase> &sp_aa_template_b( ( *template_itr)->GetMembers().LastElement());

        // if template is >= alignment A
        if
        (
          itr_a != itr_a_end &&
          sp_aa_template_a.IsDefined() &&
          !AALessThanSeqID()( sp_aa_template_a, ( *itr_a)->GetMembers().FirstElement())
        )
        {
          // iterate along the alignment until it is no longer smaller
          while( itr_a != itr_a_end && AALessThanSeqID()( ( *itr_a)->GetMembers().FirstElement(), sp_aa_template_a))
          {
            ++itr_a;
          }

          // if the data match
          if( itr_a != itr_a_end && AACompareData()( sp_aa_template_a, ( *itr_a)->GetMembers().FirstElement()))
          {
            // use this assignment
            sp_aa_a = ( *itr_a)->GetMembers().FirstElement();

            // increment
            ++itr_a;
          }
        }

        // if template is >= alignment B
        if
        (
          itr_b != itr_b_end &&
          sp_aa_template_b.IsDefined() &&
          !AALessThanSeqID()( sp_aa_template_b, ( *itr_b)->GetMembers().FirstElement())
        )
        {
          // iterate along the alignment until it is no longer smaller
          while( itr_b != itr_b_end && AALessThanSeqID()( ( *itr_b)->GetMembers().FirstElement(), sp_aa_template_b))
          {
            ++itr_b;
          }

          // if the data match
          if( itr_b != itr_b_end && AACompareData()( sp_aa_template_b, ( *itr_b)->GetMembers().FirstElement()))
          {
            // use this assignment
            sp_aa_b = ( *itr_b)->GetMembers().FirstElement();

            // increment
            ++itr_b;
          }
        }

        // pushback the assignment
        NODE.Append
        (
          util::ShPtr< align::Assignment< AABase> >
          (
            new align::Assignment< AABase>( sp_aa_a, sp_aa_b)
          )
        );
      }
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AlignByAAData::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AlignByAAData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
