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
#include "chemistry/bcl_chemistry_constitution_set.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_constitution_set.cpp
  //!
  //! @author kothiwsk, mendenjl
  //! @date Feb 25, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryConstitutionSet :
    public ExampleInterface
  {
  public:

    ExampleChemistryConstitutionSet *Clone() const
    {
      return new ExampleChemistryConstitutionSet( *this);
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
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_E.sdf"));
      chemistry::FragmentEnsemble ensemble( input);
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "1_3_pentadiene_Z.sdf"));
      ensemble.Append( chemistry::FragmentEnsemble( input));
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "azulene.sdf"));
      ensemble.Append( chemistry::FragmentEnsemble( input));
      io::File::CloseClearFStream( input);

      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "Testosterone_configurations.sdf"));
      ensemble.Append( chemistry::FragmentEnsemble( input));
      io::File::CloseClearFStream( input);

      storage::List< chemistry::FragmentConstitutionShared> constitution;

      for
      (
        storage::List< chemistry::FragmentComplete>::const_iterator itr( ensemble.Begin()), itr_end( ensemble.End());
        itr != itr_end;
        ++itr
      )
      {
        constitution.PushBack( chemistry::FragmentConstitutionShared( *itr));
      }

      // constructor
      chemistry::ConstitutionSet constitution_set
      (
        iterate::Generic< const chemistry::ConstitutionInterface>( constitution.Begin(), constitution.End())
      );

      // check the number of unique constitutions in the given library
      BCL_ExampleCheck( constitution_set.GetConstitutions().GetSize(), size_t( 3));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryConstitutionSet

  const ExampleClass::EnumType ExampleChemistryConstitutionSet::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryConstitutionSet())
  );

} // namespace bcl
