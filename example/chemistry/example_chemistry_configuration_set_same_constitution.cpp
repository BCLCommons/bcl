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
#include "chemistry/bcl_chemistry_configuration_set_same_constitution.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_configuration_set_same_constitution.cpp
  //!
  //! @author kothiwsk, mendenjl
  //! @date Mar 07, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConfigurationSetSameConstitution :
    public ExampleInterface
  {
  public:

    ExampleChemistryConfigurationSetSameConstitution *Clone() const
    {
      return new ExampleChemistryConfigurationSetSameConstitution( *this);
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

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "Testosterone_configurations.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      storage::List< chemistry::FragmentConfigurationShared> configuration;

      // constructor
      chemistry::ConfigurationSetSameConstitution configuration_set;

      // create a dummy isomorphism object; this is only necessary because we are using this class directly,
      // instead of using it through configuration set, as is normally intended
      // Ordinarily, the isomorphism vector is used to cache isomorphisms between CSIs at different layers for
      // performance and data integrity reasons
      storage::Vector< size_t> dummy_isomorphism;
      chemistry::ConstitutionSet constitutions;

      // create a constitution for the configurations to share
      util::ShPtr< chemistry::FragmentConstitutionShared> constitution
      (
        new chemistry::FragmentConstitutionShared( *ensemble.Begin())
      );

      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        configuration.PushBack( chemistry::FragmentConfigurationShared( *itr));

        // call constitutions with a dummy insert object to get the isomorphism of this molecule onto the core molecule
        constitutions.Insert( configuration.LastElement().GetConstitution(), dummy_isomorphism);
        configuration_set.Insert( configuration.LastElement(), constitution, dummy_isomorphism);
      }

      BCL_ExampleCheck( configuration_set.GetConfigurations().GetSize(), size_t( 44));

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

  }; //end ExampleChemistryConfigurationSetSameConstitution

  const ExampleClass::EnumType ExampleChemistryConfigurationSetSameConstitution::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConfigurationSetSameConstitution())
  );

} // namespace bcl
