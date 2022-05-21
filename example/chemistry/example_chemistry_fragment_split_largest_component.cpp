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
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_split_largest_component.cpp
  //!
  //! @author kothiwsk, brownbp1
  //! @date Oct 23, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentSplitLargestComponent :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentSplitLargestComponent *Clone() const
    { return new ExampleChemistryFragmentSplitLargestComponent( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      chemistry::FragmentSplitLargestComponent largest;

    /////////////////
    // data access //
    /////////////////

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "indinavir_sulfate.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      chemistry::FragmentEnsemble largest_component_ensemble( largest( ensemble.GetMolecules().FirstElement()));

      const std::string largest_component_file( AddExampleOutputPathToFilename( largest_component_ensemble, "indinavir_sulfate.largest_component.sdf"));

      io::OFStream output;

      BCL_ExampleMustOpenOutputFile( output, largest_component_file);
      largest_component_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename_component( largest_component_file + ".correct");

      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( largest_component_file, correct_filename_component),
          true,
          "largest componenet found"
        )
      )
      {
        remove( largest_component_file.c_str());
      }

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

    static const ExampleClass::EnumType ExampleChemistryFragmentSplitLargestComponent_Instance;

  }; //end ExampleChemistryFragmentSplitLargestComponent

  const ExampleClass::EnumType ExampleChemistryFragmentSplitLargestComponent::ExampleChemistryFragmentSplitLargestComponent_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentSplitLargestComponent())
  );

} // namespace bcl
