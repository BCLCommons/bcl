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
#include "chemistry/bcl_chemistry_conformation_set_same_configuration.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedrals.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_set_same_configuration.cpp
  //!
  //! @author kothiwsk
  //! @date Mar 12, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationSetSameConfiguration :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationSetSameConfiguration *Clone() const
    {
      return new ExampleChemistryConformationSetSameConfiguration( *this);
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

      // declare a list to store configuration
      storage::List< chemistry::FragmentConfigurationShared> configuration;

      // constructor for conformation comparison
      chemistry::ConformationComparisonByDihedrals conformation_comparison;

      // initialize ConformationSetSameConfiguration
      chemistry::ConformationSetSameConfiguration conformation_set( conformation_comparison, 1.0);

      // insert conformers into ConformationSetSameConfiguration
      size_t count( 0);
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
          itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        conformation_set.Insert( *itr);
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

  }; //end ExampleChemistryConformationSetSameConfiguration

  const ExampleClass::EnumType ExampleChemistryConformationSetSameConfiguration::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationSetSameConfiguration())
  );

} // namespace bcl
