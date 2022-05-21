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
#include "molecule/bcl_app_molecule_compare.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "chemistry/bcl_chemistry.h"
#include "io/bcl_io_file.h"
#include "sched/bcl_sched_schedulers.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_molecule_compare.cpp
  //!
  //! @author mendenjl
  //! @date May 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppMoleculeCompare :
    public ExampleInterface
  {
  public:

    ExampleAppMoleculeCompare *Clone() const
    {
      return new ExampleAppMoleculeCompare( *this);
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

      ApplicationExampleHelper molecule_compare_helper( app::MoleculeCompare::MoleculeCompare_Instance);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // check that flags are needed
      BCL_ExampleCheck( molecule_compare_helper.CheckCommandString( false), false);

      const std::string generated_matrix_200_200_filename
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "tanimoto.csd200.csd200.dat")
      );

      const std::string generated_matrix_1115_5_filename
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "tanimoto.csd1115.test5.dat")
      );

      // file containing the correct tanimoto matrix
      const std::string correct_matrix_200_200_filename
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "tanimoto.csd200.csd200.dat.correct")
      );
      const std::string correct_matrix_1115_5_filename
      (
        AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "tanimoto.csd1115.test5.dat.correct")
      );

      // file containing 1st 200 molecules from the csd
      const std::string fragment_file_csd200
      (
        AddExampleInputPathToFilename( e_Chemistry, "csd200first_with_atom_types.sdf")
      );
      // file containing 1st 1115 molecules from the csd
      const std::string fragment_file_csd1115
      (
        AddExampleInputPathToFilename( e_Chemistry, "csd_first_1115_simple.sdf")
      );

      // file containing five arbitrary structures from the csd
      const std::string fragment_file_test5
      (
        AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf")
      );

      // create a command line that will generate descriptors from a stored sdf
      molecule_compare_helper.ResetFlagsAndParameters();

      const double tanimoto_coefficient( 0.25);
      molecule_compare_helper.AddParameter( fragment_file_csd200);
      molecule_compare_helper.AddParameter( fragment_file_csd200);
      molecule_compare_helper.SetFlag( "output", generated_matrix_200_200_filename);
      molecule_compare_helper.SetFlag
      (
        "method",
        "LargestCommonSubstructureTanimoto(atom comparison=ElementType,bond comparison=BondOrderOrAromaticWithRingness,"
        "min = " + util::Format()( tanimoto_coefficient) + ")"
      );

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( molecule_compare_helper.CheckCommandString( true), true);

      // if threading is available, make use of it; this way we also test that threads are working properly
      if( sched::Scheduler( "PThread").IsDefined())
      {
        molecule_compare_helper.SetFlag( "scheduler", storage::Vector< std::string>::Create( "PThread", "8"));
      }

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_compare_helper.RunCommand(), 0))
      {
        // if the application ran successfully, check that files match
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( generated_matrix_200_200_filename, correct_matrix_200_200_filename),
            true,
            "csi matrix was generated correctly"
          )
        )
        {
          io::IFStream input;
          BCL_ExampleMustOpenInputFile( input, generated_matrix_200_200_filename);
          linal::Matrix< double> generated_matrix;
          io::Serialize::Read( generated_matrix, input);
          double min_nonzero_tanimoto_coefficient( 1.0), max_tanimoto_coefficient( 0.0);
          for
          (
            double *itr_matrix( generated_matrix.Begin()), *itr_matrix_end( generated_matrix.End());
            itr_matrix != itr_matrix_end;
            ++itr_matrix
          )
          {
            if( *itr_matrix < min_nonzero_tanimoto_coefficient && *itr_matrix != double( 0.0))
            {
              min_nonzero_tanimoto_coefficient = *itr_matrix;
            }
            else if( *itr_matrix > max_tanimoto_coefficient)
            {
              max_tanimoto_coefficient = *itr_matrix;
            }
          }
          BCL_ExampleIndirectCheck( min_nonzero_tanimoto_coefficient >= tanimoto_coefficient, true, "min_tanimoto flag");
          BCL_ExampleIndirectCheck
          (
            max_tanimoto_coefficient,
            double( 1.0),
            "tanimoto of identical molecules should be 1"
          );
          remove( generated_matrix_200_200_filename.c_str());
        }
      }
      molecule_compare_helper.SetFlag
      (
        "",
        storage::Vector< std::string>::Create( fragment_file_csd1115, fragment_file_test5)
      );
      molecule_compare_helper.SetFlag( "output", generated_matrix_1115_5_filename);

      molecule_compare_helper.SetFlag
      (
        "method",
        "LargestCommonSubstructureTanimoto(atom comparison=AtomType,bond comparison=BondOrder,"
        "min = " + util::Format()( tanimoto_coefficient) + ")"
      );

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( molecule_compare_helper.RunCommand(), 0))
      {
        // if the application ran successfully, check that files match
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( generated_matrix_1115_5_filename, correct_matrix_1115_5_filename),
            true,
            "csi matrix was generated correctly"
          )
        )
        {
          remove( generated_matrix_1115_5_filename.c_str());
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end Exampleexample_name

  const ExampleClass::EnumType ExampleAppMoleculeCompare::s_Instance
  (
    GetExamples().AddEnum( ExampleAppMoleculeCompare())
  );

} // namespace bcl
