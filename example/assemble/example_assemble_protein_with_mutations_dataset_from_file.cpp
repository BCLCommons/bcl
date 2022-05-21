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
#include "assemble/bcl_assemble_protein_with_mutations_dataset_from_file.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_protein_with_mutations_dataset_from_file.cpp
  //!
  //! @author mendenjl
  //! @date Jan 21, 2019
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleProteinWithMutationsDatasetFromFile :
    public ExampleInterface
  {
  public:

    ExampleAssembleProteinWithMutationsDatasetFromFile *Clone() const
    {
      return new ExampleAssembleProteinWithMutationsDatasetFromFile( *this);
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
    //////////////////
    // construction //
    //////////////////

      // testing this class will be performed once some of the descriptors have been converted to the new design

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleProteinWithMutationsDatasetFromFile

  const ExampleClass::EnumType ExampleAssembleProteinWithMutationsDatasetFromFile::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleProteinWithMutationsDatasetFromFile())
  );

} // namespace bcl
