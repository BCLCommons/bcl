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
#include "opti/bcl_opti_criterion_elapsed_time.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_convergence_argument.cpp
  //! @brief tests the implementation of the CriterionElapsedTime
  //!
  //! @author butkiem1, fischea
  //! @date Jan 4, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiCriterionElapsedTime :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiCriterionElapsedTime
    ExampleOptiCriterionElapsedTime *Clone() const
    {
      return new ExampleOptiCriterionElapsedTime( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opti::CriterionElapsedTime< double, double> criterion_default;

      // constructor with parameters
      opti::CriterionElapsedTime< double, double> criterion( util::Time( 1, 0));

      // test clone
      util::ShPtr< opti::CriterionElapsedTime< double, double> > sp_criterion( criterion.Clone());

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        criterion.GetClassIdentifier(),
        ( GetStaticClassName< opti::CriterionElapsedTime< double, double> >())
      );

    ////////////////
    // operations //
    ////////////////

      // create tracker for testing
      opti::Tracker< double, double> tracker( opti::e_SmallerIsBetter);

      // criterion should not be met since time is not up yet
      BCL_ExampleCheck( criterion.CriteriaMet( tracker), false);

      // delay for a little longer than a second because some machines delay updating the system clock for a few cycles
      // with MinGW-32 bit has approximately millisecond accuracy, while Linux has microsecond accuracy
      util::Time::Delay( util::Time( 1, 2000));

      // time should be up, criterion should be met
      BCL_ExampleCheck( criterion.CriteriaMet( tracker), true);

    //////////////////////
    // input and output //
    //////////////////////

      // write criterion
      WriteBCLObject( *sp_criterion);

      // read criterion
      opti::CriterionElapsedTime< double, double> criterion_read;
      ReadBCLObject( criterion_read);

      // compare read in object to written object
      // criterion should not be met since time is not up yet
      BCL_ExampleCheck( criterion_read.CriteriaMet( tracker), false);

      // delay for a little longer than a second because some machines delay updating the system clock for a few cycles
      // with MinGW-32 bit has approximately millisecond accuracy, while Linux has microsecond accuracy
      util::Time::Delay( util::Time( 1, 2000));

      // time should be up, criterion should be met
      BCL_ExampleCheck( criterion_read.CriteriaMet( tracker), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleOptiCriterionElapsedTime

  const ExampleClass::EnumType ExampleOptiCriterionElapsedTime::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiCriterionElapsedTime())
  );

} // namespace bcl
