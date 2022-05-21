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
#include "chemistry/bcl_chemistry_molecule_environment.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_environment_bender.h"
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_molecule_environment.cpp
  //!
  //! @author vuot2
  //! @date   09/01/2016
  //! @remarks status complete
  //! @remarks reviewed by
  //!
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryMoleculeEnvironment :
    public ExampleInterface
  {
  public:

    ExampleChemistryMoleculeEnvironment *Clone() const
    {
      return new ExampleChemistryMoleculeEnvironment( *this);
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

    //! @brief reads in filename of sdf and output the molecule environment
    chemistry::MoleculeEnvironment ReadFile( const std::string &NAME) const
    {
      // read sdf file
      io::IFStream input_sdf;
      const std::string filename
      (
        AddExampleInputPathToFilename( e_Chemistry, NAME)
      );
      io::File::TryOpenIFStream( input_sdf, filename);

      // load information into fragment complete
      chemistry::FragmentComplete mol
      (
        sdf::FragmentFactory::MakeFragment( input_sdf)
      );
      io::File::CloseClearFStream( input_sdf);

      chemistry::MoleculeEnvironment mol_env
      (
        chemistry::AtomEnvironmentBender::e_Element,
        mol
      );
      return mol_env;
    }

    int Run() const
    {
      // Create MoleculeEnvironment objects
      chemistry::MoleculeEnvironment taxol
      (
        ReadFile( "taxol.sdf")
      );

      chemistry::MoleculeEnvironment cyclohexane
      (
        ReadFile( "cyclohexane.sdf")
      );

      chemistry::MoleculeEnvironment cyclohexaneol
      (
        ReadFile( "cyclohexane1_2diol.sdf")
      );

      chemistry::MoleculeEnvironment benzene
      (
        ReadFile( "1-Phenylethanol.sdf")
      );
      // Checks general class function
      BCL_ExampleCheck( taxol.GetClassIdentifier(), "bcl::chemistry::MoleculeEnvironment");

      // Checks Tanimoto score of two identical molecular atom environments
      BCL_ExampleCheckWithinTolerance( taxol.TanimotoScore( taxol), 1.0, 0.000001);

      // Checks the Tanimoto index taxol and cyclohexane
      BCL_ExampleCheckWithinTolerance( taxol.TanimotoScore( cyclohexane), 0.030303, 0.000001);

      // Checks the Tanimoto index cyclohexane and cyclohexane-1,2diol
      BCL_ExampleCheckWithinTolerance
      (
        cyclohexane.TanimotoScore( cyclohexaneol), 0.4, 0.000001
      );

      // Checks the Buser score of taxol w itself
      BCL_ExampleCheckWithinTolerance( taxol.BuserScore( taxol), 1.0, 0.000001);

      // Check the Buser score of taxol and cyclohexane
      BCL_ExampleCheckWithinTolerance( taxol.BuserScore( cyclohexane), 0.138003, 0.00001);

      // Check the Buser score of cyclohexane and cyclohexanediol
      BCL_ExampleCheckWithinTolerance( cyclohexane.BuserScore( cyclohexaneol), 0.607625, 0.000001);

      // Check the string representation of molecule atom environments
      std::string cyclo_str [ 8] =
      {
        "[_C ][-C -C ]",
        "[_C ][-C -C ]",
        "[_C ][-C -C ]",
        "[_C ][-C -C ]",
        "[_C ][-C -C -O ]",
        "[_C ][-C -C -O ]",
        "[_O ][-C ]",
        "[_O ][-C ]"
      };
      storage::Vector< std::string> cyclo_v( 8, cyclo_str);
      auto mol_env( cyclohexaneol.GetMoleculeEnvironment());
      auto str_iter( cyclo_v.Begin());
      for
      (
        auto iter = mol_env.Begin(), end_iter = mol_env.End();
        iter != end_iter; ++str_iter, ++iter
      )
      {
        BCL_ExampleCheck( iter->UnHash(), *str_iter);
      }

      auto mol_e( benzene.GetMoleculeEnvironment());
      for
      (
          auto iter = mol_e.Begin(),
          end_iter = mol_e.End(); iter != end_iter; ++iter
      )
      {
        std::cout << iter->UnHash() << std::endl;
      }
      /*
      [ _C][ -C][ -C - O]
      [ _C][ -C ~C ~C][ -C - O ~C ~C]
      [ _C][ -C - C - O][ ~C ~C]
      [_C ][~C ~C ][-C ~C ~C ]
      [_C ][~C ~C ][-C ~C ~C ]
      [_C ][~C ~C ][~C ~C ]
      [_C ][~C ~C ][~C ~C ]
      [_C ][~C ~C ][~C ~C ]
      [ _O][ -C][ -C - C]
      */

      // Generate all the hashes of the database
      /*/ / from a vector of strings of filenames
      storage::Vector< std::string> filenames;
      //filenames.PushBack( AddExampleInputPathToFilename( e_Chemistry, "cyclohexane.sdf") );
      filenames.PushBack( AddExampleInputPathToFilename( e_Chemistry, "fixChemCart/ChemCart_mostof.sdf") );
      filenames.PushBack( AddExampleInputPathToFilename
        (
          e_Chemistry, "Vanderbilt_Discovery_Collection_--_Set_1_cleaned.sdf"
        ));
      filenames.PushBack( AddExampleInputPathToFilename
        (
          e_Chemistry, "emolecules_cleaned_3D_1-532_pass_pains_reos_pass_lipinski_unique_labeled_cleaned.sdf"
        ));
      // creates a fragment feed from filename vector
      FragmentFeed feed( filenames, sdf::e_Maintain);
      std::cout << "feed_size = " << feed->GetSize() << std::endl;
      //size_t index(0);

      // read in from the stream until we reach the end of the file or the last index
      for( ; feed.NotAtEnd(); ++feed)
      {
        //std::cout << "index = " << index << std::endl;
        //std::cout << "" << std::endl;
        chemistry::MoleculeEnvironment environment_elem
        ( chemistry::AtomEnvironmentBender::e_AtomRC, *feed);
        const chemistry::MoleculeEnvironment::t_MoleculeEnv &mol_env( environment_elem.GetMoleculeEnvironment());
        for
        (
          chemistry::MoleculeEnvironment::t_MoleculeEnv::const_iterator iter = mol_env.Begin(),
          end_iter = mol_env.End(); iter != end_iter; ++iter
        )
        {
          std::cout << iter->UnHash() << std::endl;
        }
        //++index;
      } */
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryMoleculeEnvironment

  const ExampleClass::EnumType ExampleChemistryMoleculeEnvironment::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryMoleculeEnvironment())
  );

} // namespace bcl
