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
#include "chemistry/bcl_chemistry_sample_conformations.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_sample_conformations.cpp
  //!
  //! @author kothiwsk
  //! @date Oct 14, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistrySampleConformations :
    public ExampleInterface
  {
  public:

    ExampleChemistrySampleConformations *Clone() const
    {
      return new ExampleChemistrySampleConformations( *this);
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
    int RunTestWithImplementation
    (
      const std::string &ARG,
      const chemistry::FragmentEnsemble &MOLECULES
    ) const
    {
      random::GetGlobalRandom().SetSeedFromCommandlineFlag();
      util::Implementation< chemistry::RotamerLibraryInterface> rotamer_lib;
      std::stringstream error_stream;
      if( !rotamer_lib.TryRead( util::ObjectDataLabel( ARG), error_stream))
      {
        return 0;
      }

      // fragment ensemble for storing molecules that are assembled
      chemistry::FragmentEnsemble conformations;

      chemistry::SampleConformations conformation_sampling( *rotamer_lib, "SymmetryRMSD", 0.25, 5, 10, false, 0.01);
      std::string read_label( conformation_sampling.GetString());
      read_label[ read_label.find( "cluster=0") + 8] = '1';

      conformation_sampling.AssertRead( util::ObjectDataLabel( read_label));
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
          itr( MOLECULES.Begin()), itr_end( MOLECULES.End());
        itr != itr_end;
        ++itr
      )
      {
        conformations.Append( conformation_sampling( *itr).First());
      }

      // write out assembled fragments to an example file.
      const std::string conformations_io( AddExampleOutputPathToFilename( conformations, "sample_conformations.sdf"));
      io::OFStream output;
      BCL_ExampleMustOpenOutputFile( output, conformations_io);
      conformations.WriteMDL( output);
      io::File::CloseClearFStream( output);

      std::string correct_filename( conformations_io + ".correct");
      std::string correct_filename_mac( conformations_io + ".correct.mac");
      std::string correct_filename_win( conformations_io + ".correct.win");
      if( io::File::FilesMatchWithinAbsoluteTolerance( conformations_io, correct_filename, 0.1))
      {
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( conformations_io, correct_filename, 0.1),
          true
        );
        remove( conformations_io.c_str());
      }
      else if( io::File::FilesMatchWithinAbsoluteTolerance( conformations_io, correct_filename_mac, 0.1))
      {
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( conformations_io, correct_filename_mac, 0.1),
          true
        );
        remove( conformations_io.c_str());
      }
      else
      {
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( conformations_io, correct_filename_win, 0.1),
          true
        );
      }
      return 0;
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create input stream for reading rotamer library
      io::IFStream input;

      // create input stream for reading molecule of interest
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));

      // read in the first molecule molecules
      chemistry::FragmentEnsemble molecules( input, sdf::e_Saturate);

      // close stream
      io::File::CloseClearFStream( input);

      RunTestWithImplementation( "cod", molecules);

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

  }; //end ExampleChemistrySampleConformations

  const ExampleClass::EnumType ExampleChemistrySampleConformations::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistrySampleConformations())
  );

} // namespace bcl
