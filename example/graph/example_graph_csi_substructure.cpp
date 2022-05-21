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
#include "graph/bcl_graph_csi_substructure.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_csi_substructure.cpp
  //!
  //! @author
  //! @date Mar 21, 2012
  //! @remarks status empty
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphCSISubstructure :
    public ExampleInterface
  {
  public:

    ExampleGraphCSISubstructure *Clone() const
    {
      return new ExampleGraphCSISubstructure( *this);
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

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphLargestCommonSubstructure

  const ExampleClass::EnumType ExampleGraphCSISubstructure::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphCSISubstructure())
  );

} // namespace bcl
