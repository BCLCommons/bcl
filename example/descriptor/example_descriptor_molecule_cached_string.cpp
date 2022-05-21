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
#include "descriptor/bcl_descriptor_molecule_cached_string.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_molecule_cached_string.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 14, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMoleculeCachedString :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMoleculeCachedString *Clone() const
    {
      return new ExampleDescriptorMoleculeCachedString( *this);
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
      // create input stream for reading a smallmolecule ensemble
      io::IFStream input;
      BCL_ExampleMustOpenInputFile( input, AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf"));
      // read in ensemble
      chemistry::FragmentEnsemble ensemble( input);
      // close stream
      io::File::CloseClearFStream( input);

      // make a property that looks up the mGluR5 Category
      descriptor::MoleculeCachedString mglur5_category;
      BCL_ExampleCheck( mglur5_category.GetAlias(), "Cached");
      // create the label for the property
      util::ObjectDataLabel mglur5_category_label( "Cached(mGluR5 Category,size=12)");

      BCL_ExampleAssert( mglur5_category.TryRead( util::ObjectDataLabel( "mGluR5 Category"), util::GetLogger()), true);
      // check the resulting label
      BCL_ExampleCheck
      (
        (
          util::Implementation< descriptor::Base< chemistry::AtomConformationalInterface, char> >( mglur5_category).GetLabel()
        ),
        mglur5_category_label
      );

      // make a second property that looks for a property that will be added containing common delimiters
      descriptor::MoleculeCachedString test_delimiters_ignored;
      BCL_ExampleAssert
      (
        test_delimiters_ignored.TryRead( util::ObjectDataLabel( "Delimiters"), util::GetLogger()),
        true
      );
      BCL_ExampleAssert
      (
        test_delimiters_ignored.TryRead( util::ObjectDataLabel( "Cached(Delimiters,size=100)"), util::GetLogger()),
        true
      );

      // get the first molecule out of the ensemble
      chemistry::ConformationInterface &first_molecule( *ensemble.Begin());

      // check that the mglur5 category retriever works
      BCL_ExampleCheck( std::string( mglur5_category.SumOverObject( first_molecule).Begin(), 11), "Potentiator");

      // add a property called Delimiters to the first molecule
      std::string delimiters( "\"' \t()[]\r\n,+-;:<>{}!@#$%^&*()=`\\|");
      first_molecule.StoreProperty( "Delimiters", delimiters);

      // test that delimiters are ignored
      BCL_ExampleCheck
      (
        util::TrimString( std::string( test_delimiters_ignored.SumOverObject( first_molecule).Begin(), delimiters.size())),
        delimiters
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMoleculeCachedString

  const ExampleClass::EnumType ExampleDescriptorMoleculeCachedString::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMoleculeCachedString())
  );

} // namespace bcl
