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
#include "mc/bcl_mc_metropolis.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_temperature_default.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_metropolis.cpp
  //!
  //! @author karakam, fischea
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMcMetropolis :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief clone function
    //! @return pointer to a new ExampleMcOptimizationLoopHash
    ExampleMcMetropolis *Clone() const
    {
      return new ExampleMcMetropolis( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the temperature control and other initialization variable for Metropolis
      const util::ShPtr< mc::TemperatureInterface> sp_temperature( new mc::TemperatureDefault( 100.0));
      const bool keep_history( true);
      const double min_change( 0.001);

      // construct from parameters
      mc::Metropolis< double> param_metropolis( sp_temperature, keep_history, min_change);

      // construct from serialization
      const std::string parameters( param_metropolis.GetString( true));
      mc::Metropolis< double> serialized_metropolis;
      serialized_metropolis.AssertRead( parameters);

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck
      (
        serialized_metropolis.GetClassIdentifier(), ( GetStaticClassName< mc::Metropolis< double> >())
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

  }; //end ExampleMcMetropolis

  //! single instance of this class
  const ExampleClass::EnumType ExampleMcMetropolis::s_Instance
  (
    GetExamples().AddEnum( ExampleMcMetropolis())
  );

} // namespace bcl
