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
#include "score/bcl_score_accessibility_hydrophobic_moment_magnitude.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_statistics.h"
#include "util/bcl_util_colors.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AccessibilityHydrophobicMomentMagnitude::s_Instance
    (
      GetObjectInstances().AddInstance( new AccessibilityHydrophobicMomentMagnitude())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AccessibilityHydrophobicMomentMagnitude::AccessibilityHydrophobicMomentMagnitude() :
      AccessibilityHydrophobicMoment()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param ENVIRONMENT the type of environment the accessibility was measured in that should be scored
    //! @param WINDOW_SIZES for sstype, number of restraints included in each window the moment will be calculated for
    AccessibilityHydrophobicMomentMagnitude::AccessibilityHydrophobicMomentMagnitude
    (
      const restraint::AccessibilityAA::EnvironmentType &ENVIRONMENT,
      const storage::Map< biol::SSType, size_t> WINDOW_SIZES
    ) :
      AccessibilityHydrophobicMoment( ENVIRONMENT, WINDOW_SIZES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AccessibilityHydrophobicMomentMagnitude
    AccessibilityHydrophobicMomentMagnitude *AccessibilityHydrophobicMomentMagnitude::Clone() const
    {
      return new AccessibilityHydrophobicMomentMagnitude( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AccessibilityHydrophobicMomentMagnitude::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an accessibility profile that has been assigned and scores it
    //! @param ASSIGNMENT the profile assignment that is going to be scored
    //! @return double that is the score of the accessibility profile assignment
    double AccessibilityHydrophobicMomentMagnitude::operator()
    (
      const restraint::AccessibilityProfileAssignment &ASSIGNMENT
    ) const
    {
      // exposures to hold the magnitudes of hydrophobic moments from structure
      // accessilities to hold the magnitudes of the hydrophobic moments from experiment
      storage::List< double> exposures, accessibilities;

      // iterate over the sses and their associated restraints
      for
      (
        storage::Map
        <
          util::SiPtr< const assemble::SSE>, storage::List< restraint::AccessibilityAAAssignment>
        >::const_iterator
          itr( ASSIGNMENT.GetSSEAssignments().Begin()), itr_end( ASSIGNMENT.GetSSEAssignments().End());
        itr != itr_end;
        ++itr
      )
      {
        // make sure the sse siptr is defined and get a reference
        BCL_Assert( itr->first.IsDefined(), "SiPtr is not defined");
        const assemble::SSE &current_sse( *itr->first);
        BCL_MessageDbg( "scoring sse " + util::Format()( current_sse.GetIdentification()));

        // the list of assignments associated with the current sse
        const storage::List< restraint::AccessibilityAAAssignment> &current_sse_assignments( itr->second);

        // determine window size
        size_t window_size( util::GetUndefinedSize_t());
        storage::Map< biol::SSType, size_t>::const_iterator size_itr( m_WindowSizes.Find( current_sse.GetType()));
        if( size_itr != m_WindowSizes.End())
        {
          window_size = size_itr->second;
        }
        // true if the window size is not defined or is larger than the number of assignments
        if( !util::IsDefined( window_size) || window_size > current_sse_assignments.GetSize())
        {
          // set window size to the number of assignments for this sse
          window_size = current_sse_assignments.GetSize();
        }

        // get the list of moments for overlapping windows of the assignments
        const storage::List< AccessibilityHydrophobicMoment::Window> windows
        (
          CalculateHydrophobicMomentWindows( current_sse_assignments, window_size, m_AccessibilityType)
        );

        // no windows created
        if( windows.IsEmpty())
        {
          // go to next sse
          continue;
        }

        const size_t number_of_windows( windows.GetSize());
        BCL_MessageDbg( "number_of_windows " + util::Format()( number_of_windows));

        // iterate through the list of overlapping windows of assignments
        for
        (
          storage::List< AccessibilityHydrophobicMoment::Window>::const_iterator
            window_itr( windows.Begin()), window_itr_end( windows.End());
          window_itr != window_itr_end;
          ++window_itr
        )
        {
          const double exposure_magnitude     ( window_itr->GetCalculatedMoment().Norm());
          const double accessibility_magnitude( window_itr->GetExperimentMoment().Norm());
          exposures.PushBack( exposure_magnitude);
          accessibilities.PushBack( accessibility_magnitude);
        }
      }

      BCL_MessageDbg( "exposures are\n" + util::Format()( exposures));
      BCL_MessageDbg( "accessibilities are\n" + util::Format()( accessibilities));

      // calculate score as negative of spearman correlation
      double score
      (
        -math::Statistics::CorrelationSpearman< storage::List< double>::const_iterator>
        (
          exposures.Begin(), exposures.End(), accessibilities.Begin(), accessibilities.End()
        )
      );
      BCL_MessageDbg( "CorrelationSpearman\n" + util::Format()( score));

      // if correlation is not defined
      if( !util::IsDefined( score))
      {
        // set score to zero
        score = 0;
      }

      BCL_MessageDbg( "returning score \n" + util::Format()( score));

      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AccessibilityHydrophobicMomentMagnitude::Read( std::istream &ISTREAM)
    {
      // read members
      AccessibilityHydrophobicMoment::Read( ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AccessibilityHydrophobicMomentMagnitude::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      AccessibilityHydrophobicMoment::Write( OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief write detailed scheme and values to OSTREAM
    //! @param ARGUMENT Argument to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
     std::ostream &AccessibilityHydrophobicMomentMagnitude::WriteDetailedSchemeAndValues
    (
      const restraint::AccessibilityProfileAssignment &ASSIGNMENT,
      std::ostream &OSTREAM
    ) const
     {
       OSTREAM << "AccessibilityHydrophobicMomentMagnitude::WriteDetailedSchemeAndValues\n";
       OSTREAM << "nr sses is " << ASSIGNMENT.GetSSEAssignments().GetSize() << '\n';

       // exposures to hold the magnitudes of hydrophobic moments from structure
       // accessilities to hold the magnitudes of the hydrophobic moments from experiment
       storage::List< double> exposures, accessibilities;
       for
       (
         storage::Map
         <
           util::SiPtr< const assemble::SSE>,
           storage::List< restraint::AccessibilityAAAssignment>,
           assemble::SSELessThanNoOverlap
         >::const_iterator
           sse_itr( ASSIGNMENT.GetSSEAssignments().Begin()), sse_itr_end( ASSIGNMENT.GetSSEAssignments().End());
         sse_itr != sse_itr_end;
         ++sse_itr
       )
       {
         // make sure the sse siptr is defined and get a reference
         BCL_Assert( sse_itr->first.IsDefined(), "SiPtr is not defined");
         const assemble::SSE &current_sse( *sse_itr->first);
         BCL_MessageDbg( "scoring sse " + util::Format()( current_sse.GetIdentification()));

         // the list of assignments associated with the current sse
         const storage::List< restraint::AccessibilityAAAssignment> &current_sse_assignments( sse_itr->second);

         // determine window size
         size_t window_size( util::GetUndefinedSize_t());
         storage::Map< biol::SSType, size_t>::const_iterator size_itr( m_WindowSizes.Find( current_sse.GetType()));
         if( size_itr != m_WindowSizes.End())
         {
           window_size = size_itr->second;
         }
         // true if the window size is not defined or is larger than the number of assignments
         if( !util::IsDefined( window_size) || window_size > current_sse_assignments.GetSize())
         {
           // set window size to the number of assignments for this sse
           window_size = current_sse_assignments.GetSize();
         }
         storage::List< AccessibilityHydrophobicMoment::Window> window_list
         (
           AccessibilityHydrophobicMoment::CalculateHydrophobicMomentWindows
           (
             current_sse_assignments, window_size, m_AccessibilityType
           )
         );

         std::string sse_id( sse_itr->first->GetIdentification());
         {
           util::StringReplacement space_replacer( util::StringReplacement::e_Any, " ", "");
           space_replacer.ReplaceEachIn( sse_id);
         }
         {
           util::StringReplacement space_replacer( util::StringReplacement::e_Any, "<==>", "");
           space_replacer.ReplaceEachIn( sse_id);
         }

         if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
         {
           io::OFStream write;
           io::File::MustOpenOFStream
           (
             write, "show_hydrophobic_moment_" + sse_id + "_" + m_AccessibilityType.GetString() + ".py"
           );
           AccessibilityHydrophobicMoment::ShowHydrophobicMomentWindows
           (
             window_list,
             write,
             *sse_itr->first,
             std::string( "sse_id"),
             util::GetColors().e_Cyan
           );
           io::File::CloseClearFStream( write);
         }
         OSTREAM << "for sse " << current_sse.GetIdentification() << " nr assignments is "
                 << ASSIGNMENT.GetSSEAssignments().GetSize() << " and nr windows is " << window_list.GetSize() << '\n';

         // iterate over the window list and print out the windows
         // iterate through the list of overlapping windows of assignments
         for
         (
           storage::List< AccessibilityHydrophobicMoment::Window>::const_iterator
             window_itr( window_list.Begin()), window_itr_end( window_list.End());
           window_itr != window_itr_end;
           ++window_itr
         )
         {
           OSTREAM << window_itr->GetIdentification() << '\n';

           const double exposure_magnitude     ( window_itr->GetCalculatedMoment().Norm());
           const double accessibility_magnitude( window_itr->GetExperimentMoment().Norm());
           exposures.PushBack( exposure_magnitude);
           accessibilities.PushBack( accessibility_magnitude);

           OSTREAM << "exposure_magnitude " << exposure_magnitude
                   << "\naccessibility_magnitude " << accessibility_magnitude << '\n';
         }

       }

       const double score
       (
         -math::Statistics::CorrelationSpearman
         (
           exposures.Begin(), exposures.End(), accessibilities.Begin(), accessibilities.End()
         )
       );

       OSTREAM << "score: " << score << '\n';

       return OSTREAM;
     }

  } // namespace score
} // namespace bcl
