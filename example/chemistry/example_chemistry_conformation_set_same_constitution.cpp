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
#include "chemistry/bcl_chemistry_conformation_set_same_constitution.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedrals.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_set_same_constitution.cpp
  //!
  //! @author kothiwsk
  //! @date Mar 14, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationSetSameConstitution :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationSetSameConstitution *Clone() const
    {
      return new ExampleChemistryConformationSetSameConstitution( *this);
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
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "pentadiene_conformers.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      // declare a list to store configurations
      storage::List< chemistry::FragmentConfigurationShared> configuration;

      // constructor for conformation comparison
      chemistry::ConformationComparisonByDihedrals conformation_comparison;

      util::ShPtr< chemistry::FragmentConstitutionShared> constitution
      (
        new chemistry::FragmentConstitutionShared( *ensemble.Begin())
      );

      // declare a ConformationSetSameConstitution object
      chemistry::ConformationSetSameConstitution conformation_set( conformation_comparison, 1.0);
      size_t count( 0);
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
          itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        BCL_MessageStd( "molecule number " + util::Format()( count));

        // create a dummy isomorphism object; this is only necessary because we are using this class directly,
        // instead of using it through conformation set, as is normally intended
        // Ordinarily, the isomorphism vector is used to cache isomorphisms between CSIs at different layers for
        // performance and data integrity reasons
        storage::Vector< size_t> dummy_isomorphism;
        conformation_set.Insert( *itr, constitution, dummy_isomorphism);
      }

      BCL_ExampleCheck( conformation_set.GetConformations().GetSize(), size_t( 3));

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

  }; //end ExampleChemistryConformationSetSameConstitution

  const ExampleClass::EnumType ExampleChemistryConformationSetSameConstitution::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationSetSameConstitution())
  );

} // namespace bcl
