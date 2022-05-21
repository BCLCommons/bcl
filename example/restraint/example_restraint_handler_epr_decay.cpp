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
#include "restraint/bcl_restraint_handler_epr_decay.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_handler_epr_decay.cpp
  //! @brief tests the implementation of HandlerEPRDecay, which reads in files containing EPR decay measurements.
  //!
  //! @author fischea
  //! @date Nov 11, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleHandlerEPRDecay :
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
    //! @return pointer to a new ExampleHandlerEPRDecay
    ExampleHandlerEPRDecay *Clone() const
    {
      return new ExampleHandlerEPRDecay( *this);
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

    //////////////////////
    // data preparation //
    //////////////////////

      // create an input stream to the example file holding the EPR decay data
      const std::string epr_decay_filename( AddExampleInputPathToFilename( e_Restraint, "2lzm_epr_decay"));
      io::IFStream epr_decay_file;
      BCL_ExampleMustOpenInputFile( epr_decay_file, epr_decay_filename);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      restraint::HandlerEPRDecay default_decay_hander;

    /////////////////
    // data access //
    /////////////////

      // test getter for class name identifier
      BCL_ExampleCheck( default_decay_hander.GetClassIdentifier(), ( GetStaticClassName< restraint::HandlerEPRDecay>()));

    ////////////////
    // operations //
    ////////////////

      // read in the example EPR decay file
      storage::Vector< restraint::EPRDecay> restraints( default_decay_hander.ReadRestraints( epr_decay_file));
      io::File::CloseClearFStream( epr_decay_file);

      // check if all restraints were read in correctly
      const size_t number_spin_labeling_sites( 2);
      BCL_ExampleCheck( restraints.GetSize(), number_spin_labeling_sites);
      const size_t number_measurements( restraints( 0).GetMeasurements().GetSize() + restraints( 1).GetMeasurements().GetSize());
      const size_t correct_number_measurements( 322);
      BCL_ExampleCheck( number_measurements, correct_number_measurements);

      return 0;
    }

  }; // class ExampleHandlerEPRDecay

  //! single instance of this class
  const ExampleClass::EnumType ExampleHandlerEPRDecay::s_Instance
  (
    GetExamples().AddEnum( ExampleHandlerEPRDecay())
  );

} // namespace bcl
