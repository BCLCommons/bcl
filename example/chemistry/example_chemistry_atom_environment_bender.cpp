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
#include "chemistry/bcl_chemistry_atom_environment_bender.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_atom_environment_bender.cpp
  //!
  //! @author vuot2
  //! @date Aug 23, 2017
  //! @remarks status complete
  //! @remarks
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryAtomEnvironmentBender :
    public ExampleInterface
  {

  public:
    ExampleChemistryAtomEnvironmentBender *Clone() const
    {
      return new ExampleChemistryAtomEnvironmentBender( *this);
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
      // read sdf file
      io::IFStream input_sdf;
      const std::string filename( AddExampleInputPathToFilename( e_Chemistry, "azetin3ol.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, filename);

      // load information into small_mol_conformation
      chemistry::FragmentEnsemble molecules( input_sdf);

      // close the input stream
      io::File::CloseClearFStream( input_sdf);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // iterates through all the fragment complete objects in the fragment ensemble
      for( chemistry::FragmentEnsemble::iterator itr = molecules.Begin(), itr_end = molecules.End(); itr != itr_end; ++itr)
      {
        chemistry::FragmentComplete molecule( *itr);

        // Atom index is 2
        chemistry::AtomEnvironmentBender element_AE( 2, chemistry::AtomEnvironmentBender::e_Element, molecule);
        //chemistry::AtomEnvironmentBender elemRC_AE ( 2, chemistry::AtomEnvironmentBender::e_ElemRC, molecule);
        chemistry::AtomEnvironmentBender atom_AE( 2, chemistry::AtomEnvironmentBender::e_Atom, molecule);
        //chemistry::AtomEnvironmentBender atomRC_AE ( 2, chemistry::AtomEnvironmentBender::e_AtomRC, molecule);

        // Check the class parameter functions
        BCL_ExampleCheck( element_AE.GetClassIdentifier(), "bcl::chemistry::AtomEnvironmentBender");
        BCL_ExampleCheck( element_AE.GetAlias(), "Element");

        // Check the Unhash function
        BCL_ExampleCheck( element_AE.UnHash(), "[_C ][-C -C -O ]");
        //BCL_ExampleCheck( elemRC_AE.UnHash(), "[_C.c ][-O -C.c -C.c ][=N.c ]");
        BCL_ExampleCheck( atom_AE.UnHash(), "[_C_TeTeTeTe ][-C_TrTrTrPi -C_TeTeTeTe -O_Te2Te2TeTe ]");
        //BCL_ExampleCheck( atomRC_AE.UnHash(), "[_C_TeTeTeTe.c ][-O_Te2Te2TeTe -C_TrTrTrPi.c -C_TeTeTeTe.c ][=N_Tr2TrTrPi.c ]");

        // Check the constructor from the string representation of the molecule
        chemistry::AtomEnvironmentBender element_AE2( chemistry::AtomEnvironmentBender::e_Element, element_AE.UnHash());
        BCL_ExampleCheck( element_AE, element_AE2);

        //chemistry::AtomEnvironmentBender elemRC_AE2( chemistry::AtomEnvironmentBender::e_ElemRC, elemRC_AE.UnHash());
        //BCL_ExampleCheck( elemRC_AE, elemRC_AE2);

        chemistry::AtomEnvironmentBender atom_AE2( chemistry::AtomEnvironmentBender::e_Atom, atom_AE.UnHash());
        BCL_ExampleCheck( atom_AE, atom_AE2);

        //chemistry::AtomEnvironmentBender atomRC_AE2( chemistry::AtomEnvironmentBender::e_AtomRC, atomRC_AE.UnHash());
        //BCL_ExampleCheck( atomRC_AE, atomRC_AE2);
      }

      // Test the aromatic bond
      // read sdf file
      BCL_ExampleMustOpenInputFile( input_sdf, AddExampleInputPathToFilename( e_Chemistry, "benzene_prepared.sdf"));

      // load information into small_mol_conformation
      chemistry::FragmentEnsemble molecule2( input_sdf);

      // close the input stream
      io::File::CloseClearFStream( input_sdf);
      // iterates through all the fragment complete objects in the fragment ensemble
      for( chemistry::FragmentEnsemble::iterator itr = molecule2.Begin(), itr_end = molecule2.End(); itr != itr_end; ++itr)
      {
        chemistry::FragmentComplete molecule( *itr);

        chemistry::AtomEnvironmentBender element_AE( 2, chemistry::AtomEnvironmentBender::e_Element, molecule);
        std::cout << "benzene" << element_AE.UnHash() << std::endl;

        chemistry::AtomEnvironmentBender atom_AE( 2, chemistry::AtomEnvironmentBender::e_Atom, molecule);
        std::cout << "benzene" << atom_AE.UnHash() << std::endl;
      }
      return 0;
    } // Run
    /*
    //! @brief check if the StringToAE works correctly
    // Print out all the atom environments of a molecule
    static void CheckAtomEnv( const chemistry::AtomEnvironmentBender &MOLECULE)
    {
      chemistry::AtomEnvironmentBender::t_AtomEnvironment atom_env( MOLECULE.GetAtomEnvironment());

      //check that the number and the identities of the atoms in the atom environment are correct
      std::cout << "The atom of interest is: " << MOLECULE.GetAtomOfInterestIndex() << std::endl;
      int i = 0;
      for( chemistry::AtomEnvironmentBender::t_AtomEnvironment::const_iterator itr = atom_env.Begin(), itr_end = atom_env.End(); itr != itr_end; ++itr, ++i)
      {
        std::cout << i << ": ";
        for( storage::Map< size_t, size_t>::const_iterator map_itr = itr->Begin(), map_itr_end = itr->End(); map_itr != map_itr_end; ++map_itr)
        {
          std::cout << map_itr->first << "(" << map_itr->second << ") ";
        }
        std::cout <<std::endl;
      }
      // check the the string of unhashed atom environment
      std::cout << "the uncoded string elem is: " << MOLECULE.UnHash() << std::endl;
    }
    */

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryAtomEnvironmentBender

  const ExampleClass::EnumType ExampleChemistryAtomEnvironmentBender::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryAtomEnvironmentBender())
  );

} // namespace bcl
