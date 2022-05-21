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
#include "chemistry/bcl_chemistry_configuration_set.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_configuration_set.cpp
  //!
  //! @author kothiwsk
  //! @date Mar 09, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConfigurationSet :
    public ExampleInterface
  {
  public:

    ExampleChemistryConfigurationSet *Clone() const
    {
      return new ExampleChemistryConfigurationSet( *this);
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

      // load a ensembles directly from input streams
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_E.sdf"));
      chemistry::FragmentEnsemble ensemble_a( input);
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_Z.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "azulene.sdf"));
      chemistry::FragmentEnsemble azulene_ensemble( input);
      ensemble.Append( azulene_ensemble);
      chemistry::FragmentConfigurationShared azulene_config( *azulene_ensemble.Begin());
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "Testosterone_configurations.sdf"));
      chemistry::FragmentEnsemble testosterone_isomers( input);
      ensemble.Append( testosterone_isomers);
      chemistry::FragmentConfigurationShared first_testosterone_isomer( *testosterone_isomers.Begin());
      io::File::CloseClearFStream( input);

      storage::List< chemistry::FragmentConfigurationShared> configuration;

      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        configuration.PushBack( chemistry::FragmentConfigurationShared( *itr));
      }

      // constructor
      chemistry::ConfigurationSet configuration_set
      (
        iterate::Generic< const chemistry::ConfigurationInterface>( configuration.Begin(), configuration.End())
      );

      BCL_ExampleCheck( configuration_set.GetConstitutions().GetSize(), size_t( 3));
      BCL_ExampleCheck( configuration_set.GetConfigurationsMap().GetSize(), size_t( 3));
      BCL_ExampleCheck( configuration_set.Find( azulene_config)->IsDefined(), true);
      BCL_ExampleCheck( configuration_set.Find( first_testosterone_isomer)->IsDefined(), true);

      // get the constitution stored for testosterone
      const util::SiPtr< const chemistry::ConstitutionInterface> testosterone_constitution
      (
        ( *configuration_set.Find( first_testosterone_isomer))->GetConstitution()
      );
      BCL_ExampleCheck
      (
        configuration_set.GetConfigurationsMap().Find( testosterone_constitution)->second.GetSize(),
        size_t( 44)
      );
      BCL_ExampleCheck( configuration_set.GetConfigurations().GetSize(), size_t( 46));
      BCL_ExampleCheck( configuration_set.GetSize(), configuration_set.GetConfigurations().GetSize());

      configuration_set.Insert( chemistry::FragmentConfigurationShared( ensemble_a.GetMolecules().FirstElement()));

      BCL_ExampleCheck( configuration_set.GetConstitutions().GetSize(), size_t( 3));
      BCL_ExampleCheck( configuration_set.GetConfigurationsMap().GetSize(), size_t( 3));
      BCL_ExampleCheck( configuration_set.GetConfigurations().GetSize(), size_t( 47));
      BCL_ExampleCheck( configuration_set.GetSize(), configuration_set.GetConfigurations().GetSize());

      chemistry::AtomVector< chemistry::AtomComplete> atoms_complete
      (
        configuration.Begin()->GetAtomInfo(),
        configuration.Begin()->GetBondInfo()
      );
      storage::Vector< size_t> new_indices( configuration.Begin()->GetNumberAtoms());
      for( size_t i( 0), n( new_indices.GetSize()); i < n; ++i)
      {
        new_indices( i) = i;
      }

      bool found_first_molecule_after_reordering( true);
      for( size_t i( 0); i < 1000; ++i)
      {
        new_indices.Shuffle();
        atoms_complete.Reorder( new_indices);
        chemistry::FragmentConfigurationShared perturbed = chemistry::FragmentConfigurationShared( chemistry::FragmentComplete( atoms_complete, ""));
        chemistry::ConfigurationSet::const_iterator itr( configuration_set.Find( perturbed));
        if( itr != configuration_set.Begin())
        {
          found_first_molecule_after_reordering = false;
          break;
        }
      }
      BCL_ExampleIndirectCheck( found_first_molecule_after_reordering, true, "Find");

      bool find_respects_chirality( true);
      for( size_t i( 0); i < 1000; ++i)
      {
        storage::Vector< sdf::AtomInfo> atom_info( configuration.Begin()->GetAtomInfo());
        atom_info( 0).SetChirality( chemistry::e_SChirality);
        chemistry::AtomVector< chemistry::AtomComplete> atoms_complete
        (
          atom_info,
          configuration.Begin()->GetBondInfo()
        );
        new_indices.Shuffle();
        atoms_complete.Reorder( new_indices);
        chemistry::FragmentConfigurationShared perturbed = chemistry::FragmentConfigurationShared( chemistry::FragmentComplete( atoms_complete, ""));
        chemistry::ConfigurationSet::const_iterator itr( configuration_set.Find( perturbed));
        if( itr != configuration_set.End())
        {
          find_respects_chirality = false;
          break;
        }
      }
      BCL_ExampleIndirectCheck( find_respects_chirality, true, "Find");

      bool find_respects_atom_types( true);
      for( size_t i( 0); i < 1000; ++i)
      {
        storage::Vector< sdf::AtomInfo> atom_info( configuration.Begin()->GetAtomInfo());
        atom_info( 0).SetAtomType( chemistry::GetAtomTypes().Al_DiDi);
        chemistry::AtomVector< chemistry::AtomComplete> atoms_complete
        (
          atom_info,
          configuration.Begin()->GetBondInfo()
        );
        new_indices.Shuffle();
        atoms_complete.Reorder( new_indices);
        chemistry::FragmentConfigurationShared perturbed = chemistry::FragmentConfigurationShared( chemistry::FragmentComplete( atoms_complete, ""));
        chemistry::ConfigurationSet::const_iterator itr( configuration_set.Find( perturbed));
        if( itr != configuration_set.End())
        {
          find_respects_atom_types = false;
          break;
        }
      }
      BCL_ExampleIndirectCheck( find_respects_atom_types, true, "Find");

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

  }; //end ExampleChemistryConfigurationSet

  const ExampleClass::EnumType ExampleChemistryConfigurationSet::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConfigurationSet())
  );

} // namespace bcl
