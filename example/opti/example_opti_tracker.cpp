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
#include "opti/bcl_opti_tracker.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_tracker.cpp
  //! @brief this example tests the implementation of opti::Tracker
  //!
  //! @author fischea
  //! @date Dec 14, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiTracker :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiTracker
    ExampleOptiTracker *Clone() const
    {
      return new ExampleOptiTracker( *this);
    }

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
      // create improvement types
      const opti::ImprovementType type_larger( opti::e_LargerIsBetter);
      const opti::ImprovementType type_smaller( opti::e_SmallerIsBetter);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor
      util::ShPtr< opti::Tracker< double, double> > sp_tracker( new opti::Tracker< double, double>( type_larger));

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 0 &&
        sp_tracker->GetBestIteration() == 0 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 0,
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 0. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 0. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0."
      );

      // test cloning function
      util::ShPtr< opti::Tracker< double, double> > sp_tracker_clone( sp_tracker->Clone());

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == sp_tracker_clone->GetIteration() &&
        sp_tracker->GetBestIteration() == sp_tracker_clone->GetBestIteration() &&
        sp_tracker->GetIterationsSinceLastImprovement() == sp_tracker_clone->GetIterationsSinceLastImprovement(),
        "Members of clone don't have the same values like the members of the original tracker: " +
        util::Format()( *sp_tracker) + "\n" + util::Format()( *sp_tracker_clone)
      );

    ////////////////////////////////
    // data access and operations //
    ////////////////////////////////

      // test get class identifier
      BCL_ExampleCheck
      (
        sp_tracker->GetClassIdentifier(),
        ( GetStaticClassName< opti::Tracker< double, double> >())
      );

      // argument result pairs for testing
      util::ShPtr< storage::Pair< double, double> > sp_pair_1( new storage::Pair< double, double>( 1.0, 2.0));
      util::ShPtr< storage::Pair< double, double> > sp_pair_2( new storage::Pair< double, double>( 2.0, 3.0));
      util::ShPtr< storage::Pair< double, double> > sp_pair_3( new storage::Pair< double, double>( 3.0, 1.0));
      util::ShPtr< storage::Pair< double, double> > sp_pair_4( new storage::Pair< double, double>( 4.0, 4.0));

      // track first pair
      sp_tracker->Track( sp_pair_1);
      util::ShPtr< storage::Pair< double, double> > sp_pair_current( sp_tracker->GetCurrent());
      util::ShPtr< storage::Pair< double, double> > sp_pair_best( sp_tracker->GetBest());

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 1 &&
        sp_tracker->GetBestIteration() == 1 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 0 &&
        sp_pair_current->First() == sp_pair_1->First() &&
        sp_pair_current->Second() == sp_pair_1->Second() &&
        sp_pair_best->First() == sp_pair_1->First() &&
        sp_pair_best->Second() == sp_pair_1->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 1. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 1. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_1) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_1)
      );

      // track second pair
      sp_tracker->Track( sp_pair_2);
      sp_pair_current = sp_tracker->GetCurrent();
      sp_pair_best = sp_tracker->GetBest();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 2 &&
        sp_tracker->GetBestIteration() == 2 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 0 &&
        sp_pair_current->First() == sp_pair_2->First() &&
        sp_pair_current->Second() == sp_pair_2->Second() &&
        sp_pair_best->First() == sp_pair_2->First() &&
        sp_pair_best->Second() == sp_pair_2->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 2. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 2. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_2) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_2)
      );

      // track third pair
      sp_tracker->Track( sp_pair_3);
      sp_pair_current = sp_tracker->GetCurrent();
      sp_pair_best = sp_tracker->GetBest();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 3 &&
        sp_tracker->GetBestIteration() == 2 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 1 &&
        sp_pair_current->First() == sp_pair_3->First() &&
        sp_pair_current->Second() == sp_pair_3->Second() &&
        sp_pair_best->First() == sp_pair_2->First() &&
        sp_pair_best->Second() == sp_pair_2->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 3. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 2. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 1. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_3) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_2)
      );

      // track fourth pair
      sp_tracker->Track( sp_pair_4);
      sp_pair_current = sp_tracker->GetCurrent();
      sp_pair_best = sp_tracker->GetBest();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 4 &&
        sp_tracker->GetBestIteration() == 4 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 0 &&
        sp_pair_current->First() == sp_pair_4->First() &&
        sp_pair_current->Second() == sp_pair_4->Second() &&
        sp_pair_best->First() == sp_pair_4->First() &&
        sp_pair_best->Second() == sp_pair_4->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 4. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 4. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_4) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_4)
      );

      // test tracker with improvement type "smaller is better"
      sp_tracker = util::ShPtr< opti::Tracker< double, double> >( new opti::Tracker< double, double>( type_smaller));

      // track first pair
      sp_tracker->Track( sp_pair_1);
      sp_pair_current = sp_tracker->GetCurrent();
      sp_pair_best = sp_tracker->GetBest();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 1 &&
        sp_tracker->GetBestIteration() == 1 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 0 &&
        sp_pair_current->First() == sp_pair_1->First() &&
        sp_pair_current->Second() == sp_pair_1->Second() &&
        sp_pair_best->First() == sp_pair_1->First() &&
        sp_pair_best->Second() == sp_pair_1->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 1. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 1. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_1) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_1)
      );

      // track second pair
      sp_tracker->Track( sp_pair_2);
      sp_pair_current = sp_tracker->GetCurrent();
      sp_pair_best = sp_tracker->GetBest();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 2 &&
        sp_tracker->GetBestIteration() == 1 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 1 &&
        sp_pair_current->First() == sp_pair_2->First() &&
        sp_pair_current->Second() == sp_pair_2->Second() &&
        sp_pair_best->First() == sp_pair_1->First() &&
        sp_pair_best->Second() == sp_pair_1->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 2. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 1. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 1. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_2) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_1)
      );

      // track third pair
      sp_tracker->Track( sp_pair_3);
      sp_pair_current = sp_tracker->GetCurrent();
      sp_pair_best = sp_tracker->GetBest();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 3 &&
        sp_tracker->GetBestIteration() == 3 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 0 &&
        sp_pair_current->First() == sp_pair_3->First() &&
        sp_pair_current->Second() == sp_pair_3->Second() &&
        sp_pair_best->First() == sp_pair_3->First() &&
        sp_pair_best->Second() == sp_pair_3->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 3. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 3. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_3) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_3)
      );

      // track fourth pair
      sp_tracker->Track( sp_pair_4);
      sp_pair_current = sp_tracker->GetCurrent();
      sp_pair_best = sp_tracker->GetBest();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 4 &&
        sp_tracker->GetBestIteration() == 3 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 1 &&
        sp_pair_current->First() == sp_pair_4->First() &&
        sp_pair_current->Second() == sp_pair_4->Second() &&
        sp_pair_best->First() == sp_pair_3->First() &&
        sp_pair_best->Second() == sp_pair_3->Second(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 4. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 3. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0. " +
        "Current is " + util::Format()( sp_pair_current) + " and should be " + util::Format()( sp_pair_4) +
        "Best is " + util::Format()( sp_pair_best) + " and should be " + util::Format()( sp_pair_3)
      );

      // test reset
      sp_tracker->Reset();

      BCL_Example_Check
      (
        sp_tracker->GetIteration() == 0 &&
        sp_tracker->GetBestIteration() == 0 &&
        sp_tracker->GetIterationsSinceLastImprovement() == 0 &&
        !sp_tracker->GetCurrent().IsDefined() &&
        !sp_tracker->GetBest().IsDefined(),
        "Iteration count is " + util::Format()( sp_tracker->GetIteration()) + " and should be 0. " +
        "Best iteration is " + util::Format()( sp_tracker->GetBestIteration()) + " and should be 0. " +
        "Iteration since last improvement is " + util::Format()( sp_tracker->GetIterationsSinceLastImprovement()) +
        " and should be 0. " +
        "Current is " + util::Format()( sp_tracker->GetCurrent()) + " and should be undefined." +
        "Best is " + util::Format()( sp_tracker->GetBest()) + " and should be undefined."
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write tracker
      WriteBCLObject( *sp_tracker);

      // read tracker
      opti::Tracker< double, double> tracker_read( opti::e_SmallerIsBetter);
      ReadBCLObject( tracker_read);

      // compare read in object to written object
      BCL_Example_Check
      (
        sp_tracker->GetIteration() == tracker_read.GetIteration() &&
        sp_tracker->GetBestIteration() == tracker_read.GetBestIteration() &&
        sp_tracker->GetIterationsSinceLastImprovement() == tracker_read.GetIterationsSinceLastImprovement(),
        "Members of read in object don't have the same values like the members of the original tracker: " +
        util::Format()( *sp_tracker) + "\n" + util::Format()( tracker_read)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

    }; // class ExampleOptiTracker

    const ExampleClass::EnumType ExampleOptiTracker::s_Instance
    (
      GetExamples().AddEnum( ExampleOptiTracker())
    );

} // namespace bcl
