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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "opti/bcl_opti_tracker_with_history.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_evolution_population.hpp"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_tracker_with_history.cpp
  //! @brief this example tests the implementation of TrackerWithHistory
  //!
  //! @author geanesar
  //! @date Nov 1, 2014
  //! @remarks status complete 
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiTrackerWithHistory :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiTrackerWithHistory
    ExampleOptiTrackerWithHistory *Clone() const
    {
      return new ExampleOptiTrackerWithHistory( *this);
    }

  //////////
  // data //
  //////////

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

      // Construction
      opti::TrackerWithHistory< int, int> tracker;
      opti::TrackerWithHistory< int, int> tracker_limit( size_t( 2));
      opti::TrackerWithHistory< int, int> tracker_infinite( size_t( 0));
      opti::TrackerWithHistory< int, int> tracker_full( opti::e_LargerIsBetter, 1, 2);

      // Copying
      util::ShPtr< opti::TrackerWithHistory< int, int> > new_tracker( tracker_full.Clone());
      BCL_ExampleCheck( new_tracker.IsDefined(), true);

      // Input/output
      BCL_ExampleCheck( tracker.GetClassIdentifier(), "bcl::opti::TrackerWithHistory<int,int>");

      util::ShPtr< storage::Pair< int, int> > init( new storage::Pair< int, int>( 1, 1));
      util::ShPtr< storage::Pair< int, int> > to_track( new storage::Pair< int, int>( 0, 0));
      util::ShPtr< storage::Pair< int, int> > final( new storage::Pair< int, int>( 1, 1));

      tracker_full.Track( init);
      
      // These should not increase the number of histories
      tracker_infinite.Track( init);
      tracker_infinite.Track( init);

      BCL_ExampleCheck( tracker_full.GetCurrent()->First(), init->First());

      tracker_full.SetPhase( opti::e_Iteration);
      tracker_infinite.SetPhase( opti::e_Iteration);

      tracker_full.Track( to_track, opti::e_Improved);
      tracker_infinite.Track( to_track, opti::e_Improved);
      tracker_full.Track( to_track, opti::e_Improved);
      tracker_infinite.Track( to_track, opti::e_Improved);
      tracker_full.Track( final, opti::e_Improved);
      tracker_infinite.Track( final, opti::e_Improved);

      // Initial should still be the same
      BCL_ExampleCheck( tracker_full.GetInitial()->First(), init->First());
      
      // Check history sizes are correct
      BCL_ExampleCheck( tracker_full.GetHistory().GetSize(), 2);
      BCL_ExampleCheck( tracker_infinite.GetHistory().GetSize(), 4);

      // Check last entries
      BCL_ExampleCheck( ( *tracker_full.GetHistory().Begin())->First(), final->First());
      BCL_ExampleCheck( ( *tracker_infinite.GetHistory().Begin())->First(), final->First());

      // Changing max history size
      tracker_infinite.SetMaxHistorySize( 2);
      BCL_ExampleCheck( tracker_infinite.GetHistory().GetSize(), 2);
      BCL_ExampleCheck( ( *tracker_infinite.GetHistory().Begin())->First(), final->First());

      // Reset/clear
      tracker_full.Reset();

      BCL_ExampleCheck( tracker_full.GetHistory().GetSize(), 0);
      BCL_ExampleCheck( tracker_full.GetInitial().IsDefined(), false);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; // class ExampleOptiTrackerWithHistory

   const ExampleClass::EnumType ExampleOptiTrackerWithHistory::s_Instance
   (
     GetExamples().AddEnum( ExampleOptiTrackerWithHistory())
   );

} // namespace bcl
