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
#include "mc/bcl_mc_mutates.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_mutates.h"

// external includes - sorted alphabetically

namespace bcl
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_mutates.cpp
  //! @brief This example tests the implementation of the class mc::Mutates, which initializes all mutates for
  //! structure prediction and adds them to the enumerate.
  //!
  //! @author fischea
  //! @date Aug 17, 2016
  //! @remarks status complete
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleMcMutates :
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
    //! @return pointer to a new ExampleMcMutates
    ExampleMcMutates *Clone() const
    {
      return new ExampleMcMutates( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! get the single instance of this class
      mc::Mutates &mutates( mc::Mutates::GetInstance());

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( mutates.GetClassIdentifier(), ( GetStaticClassName< mc::Mutates>()));

    ////////////////
    // operations //
    ////////////////

      // initialize mutates and check if they have been added to the enumerator
      mutates.Initialize();
      BCL_ExampleCheck( fold::GetMutates().GetEnumCount() != 0, true);
      BCL_MessageStd
      (
        "Final mutate enums: " + util::Format()( fold::GetMutates().GetEnumCount())
      );

      return 0;
    }

  }; // class ExampleMcMutates

  //! single instance of this class
  const ExampleClass::EnumType ExampleMcMutates::s_Instance
  (
     GetExamples().AddEnum( ExampleMcMutates())
  );

} // namespace bcl
