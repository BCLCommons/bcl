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
#include "restraint/bcl_restraint_epr_decay.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_epr_decay.cpp
  //! @brief tests the implementation of EPRDecay, which stores the results of EPR decay measurements.
  //!
  //! @author fischea
  //! @date Nov 11, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleRestraintEPRDecay :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleEPRDecay
    ExampleRestraintEPRDecay *Clone() const
    {
      return new ExampleRestraintEPRDecay( *this);
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
      restraint::EPRDecay default_decay;

      // construct from spin-labeling sites
      restraint::EPRDecay decay( 'A', 10, 'B', 212);
      restraint::EPRDecay::SLPair sl_pair( decay.GetSpinLabelingSites());
      BCL_ExampleCheck
      (
        sl_pair.First().First() == 'A' && sl_pair.First().Second() == 10 && sl_pair.Second().First() == 'B' &&
        sl_pair.Second().Second() == 212,
        true
      );

    /////////////////
    // data access //
    /////////////////

      // test getter for class name identifier
      BCL_ExampleCheck( default_decay.GetClassIdentifier(), ( GetStaticClassName< restraint::EPRDecay>()));

      // add some measurements
      decay.AddMeasurement( 10.0, 11.0);
      const storage::Vector< storage::Pair< double, double> > &measurements( decay.GetMeasurements());
      BCL_ExampleCheck( measurements.GetSize(), 1);
      BCL_ExampleCheck( measurements( 0).First(), 10.0);
      BCL_ExampleCheck( measurements( 0).Second(), 11.0);
      decay.AddMeasurement( 2010.0, 11.0);
      decay.AddMeasurement( 99.0, -12.0);
      BCL_ExampleCheck( measurements.GetSize(), 3);

    ///////////////
    // operators //
    ///////////////

      return 0;
    }

  }; // class ExampleEPRDecay

  //! single instance of this class
  const ExampleClass::EnumType ExampleRestraintEPRDecay::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintEPRDecay())
  );

} // namespace bcl
