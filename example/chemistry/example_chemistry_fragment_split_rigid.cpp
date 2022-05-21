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
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_fragment_split_rigid.cpp
  //!
  //! @author mendenjl
  //! @date Jan 15, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryFragmentSplitRigid :
    public ExampleInterface
  {
  public:

    ExampleChemistryFragmentSplitRigid *Clone() const
    { return new ExampleChemistryFragmentSplitRigid( *this);}

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

      chemistry::FragmentSplitRigid rigid;

    /////////////////
    // data access //
    /////////////////

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      chemistry::FragmentEnsemble rigid_ensemble( rigid( ensemble.GetMolecules().FirstElement()));
      const std::string rigid_split( AddExampleOutputPathToFilename( ensemble, "split_rigid.sdf"));

      io::OFStream output;

      BCL_ExampleMustOpenOutputFile( output, rigid_split);
      rigid_ensemble.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename_rigid( rigid_split + ".correct");

      if
      (
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( rigid_split, correct_filename_rigid),
          true,
          "graph was split correctly"
        )
      )
      {
        remove( rigid_split.c_str());
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

    static const ExampleClass::EnumType ExampleChemistryFragmentSplitRigid_Instance;

  }; //end ExampleChemistryFragmentSplitRigid

  const ExampleClass::EnumType ExampleChemistryFragmentSplitRigid::ExampleChemistryFragmentSplitRigid_Instance
  (
    GetExamples().AddEnum( ExampleChemistryFragmentSplitRigid())
  );

} // namespace bcl
