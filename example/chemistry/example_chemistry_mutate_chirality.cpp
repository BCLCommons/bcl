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
#include "chemistry/bcl_chemistry_mutate_chirality.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_mutate_chirality.cpp
  //! @author kothiwsk
  //! @date Mar 31, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMutateChirality :
    public ExampleInterface
  {
  public:

    ExampleChemistryMutateChirality *Clone() const
    {
      return new ExampleChemistryMutateChirality( *this);
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

      // create input stream for reading a fragment ensemble
      io::IFStream input;

      // create input stream for reading molecule of interest which needs to generated
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "taxol.sdf"));
      // read in molecule
      chemistry::FragmentEnsemble molecules( input, sdf::e_Remove);
      // close stream
      io::File::CloseClearFStream( input);

      chemistry::FragmentComplete &first_molecule( molecules.GetMolecules().FirstElement());

      chemistry::MutateChirality first_mutate( first_molecule);
      chemistry::FragmentEnsemble new_conformations;
      new_conformations.PushBack( *first_mutate( first_molecule).GetArgument());

      for( size_t counts( 0); counts < size_t( 5); ++counts)
      {
        new_conformations.PushBack( *first_mutate( new_conformations.GetMolecules().LastElement()).GetArgument());
      }

      const std::string mutated_conformations( AddExampleOutputPathToFilename( new_conformations, "mutated_chirality.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, mutated_conformations);
      new_conformations.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename( mutated_conformations + ".correct");
      std::string correct_win_filename( mutated_conformations + ".win.correct");

      // check if the generated file is right
      if( io::File::FilesMatch( mutated_conformations, correct_filename))
      {
        BCL_ExampleCheck( io::File::FilesMatch( mutated_conformations, correct_filename), true);
        remove( mutated_conformations.c_str());
      }
      else if( BCL_ExampleCheck( io::File::FilesMatch( mutated_conformations, correct_win_filename), true))
      {
        remove( mutated_conformations.c_str());
      }

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
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryMutateChirality

  const ExampleClass::EnumType ExampleChemistryMutateChirality::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMutateChirality())
  );

} // namespace bcl
