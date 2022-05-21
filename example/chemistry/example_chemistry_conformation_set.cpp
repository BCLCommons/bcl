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
#include "chemistry/bcl_chemistry_conformation_set.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedrals.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_conformation_set.cpp
  //!
  //! @author kothiwsk
  //! @date Mar 13, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConformationSet :
    public ExampleInterface
  {
  public:

    ExampleChemistryConformationSet *Clone() const
    {
      return new ExampleChemistryConformationSet( *this);
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

      // load an ensemble directly from an input stream
      io::IFStream input;
      BCL_ExampleMustOpenInputFile
      (
        input,
        AddExampleInputPathToFilename( e_Chemistry, "conformation_sampling_fragments.sdf")
      );

      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      chemistry::ConformationComparisonByDihedrals conformation_comparison;
      // constructor
      chemistry::ConformationSet conformation_set( conformation_comparison, 1.0);
      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator
          itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        conformation_set.Insert( *itr);
      }

      BCL_ExampleCheck( conformation_set.GetConfigurations().GetSize(), size_t( 102));
      BCL_ExampleCheck( conformation_set.GetConformations().GetSize(), size_t( 123));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryConformationSet

  const ExampleClass::EnumType ExampleChemistryConformationSet::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConformationSet())
  );

} // namespace bcl
