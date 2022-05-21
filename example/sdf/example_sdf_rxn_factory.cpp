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
#include "sdf/bcl_sdf_rxn_factory.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "chemistry/bcl_chemistry_reaction_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_rxn_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_rxn_factory.cpp
  //! @details Demonstrates how the RXNFactory class is used to read and write SmallMolecule
  //! and SmallMoleculeEnsembles from files
  //!
  //! @author geanesar, combss, mendenjl
  //! @date Jun 02, 2014
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfRXNFactory :
    public ExampleInterface
  {
  public:

    ExampleSdfRXNFactory *Clone() const
    {
      return new ExampleSdfRXNFactory( *this);
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

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Chemistry, "multi_reaction_1.rxn"));

      sdf::RXNHandler handler_default;
      sdf::RXNHandler handler( read);

      io::File::CloseClearFStream( read);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construction
      chemistry::ReactionComplete new_reaction
      (
        sdf::RXNFactory::MakeReactionComplete( handler)
      );

    //////////////////
    // output check //
    //////////////////

      BCL_ExampleCheck( new_reaction.GetNumberReactants(), 2);
      BCL_ExampleCheck( new_reaction.GetNumberProducts(), 2);

//      // test the stream operators
//      BCL_MessageStd( "Outputting factory_a: " + util::Format()( factory_a));

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleSdfRXNFactory

  const ExampleClass::EnumType ExampleSdfRXNFactory::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfRXNFactory())
  );
} // namespace bcl

