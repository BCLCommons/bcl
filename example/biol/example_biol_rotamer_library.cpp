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
#include "biol/bcl_biol_rotamer_library.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_rotamer_library.cpp
  //! @brief this example tests the implementation of the rotamer library for amino acid side chains
  //!
  //! @author fischea
  //! @date Nov 13, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleBiolRotamerLibrary :
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
    //! @return pointer to a new ExampleBiolRotamerLibrary
    ExampleBiolRotamerLibrary *Clone() const
    {
      return new ExampleBiolRotamerLibrary( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////
    // data preparation //
    //////////////////////

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from an existing rotamer library
      const std::string library_file_path( AddExampleInputPathToFilename( e_Biology, "rotamer_library"));
      util::ShPtr< biol::RotamerLibrary> sp_library( biol::RotamerLibrary::CreateRotamerLibrary( library_file_path));

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( sp_library->GetClassIdentifier(), ( GetStaticClassName< biol::RotamerLibrary>()));

      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  }; // class ExampleBiolRotamerLibrary

  //! single instance of this class
  const ExampleClass::EnumType ExampleBiolRotamerLibrary::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolRotamerLibrary())
  );

} // namespace bcl
