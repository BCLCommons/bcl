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
#include "chemistry/bcl_chemistry_collector_valence.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_chemistry_collector_valence.cpp
  //!
  //! @author mendenjl, sliwosgr
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleChemistryCollectorValence :
    public ExampleInterface
  {
  public:

    ExampleChemistryCollectorValence *Clone() const
    {
      return new ExampleChemistryCollectorValence( *this);
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

      // default constructor
      chemistry::CollectorValence cvalence_default;

      // clone
      const util::ShPtr< chemistry::CollectorValence> cvalence_clone( cvalence_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // verify that clone returned a chemistry::CollectorValence
      BCL_ExampleIndirectCheck( cvalence_clone.IsDefined(), true, "clone");
      BCL_ExampleCheck
      (
        chemistry::CollectorValence().GetClassIdentifier(),
        GetStaticClassName< chemistry::CollectorValence>()
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // setup input stream
      io::IFStream input_sdf;

      // read in molecule
      const std::string filename_in( AddExampleInputPathToFilename( e_Chemistry, "corina_diazepam.sdf"));
      BCL_ExampleMustOpenInputFile( input_sdf, filename_in);

      chemistry::FragmentComplete diazepam( sdf::FragmentFactory::MakeFragment( input_sdf));

      BCL_ExampleIndirectCheck
      (
        chemistry::CollectorValence().Collect( diazepam).GetSize(), 0, "open valences before removal of hydrogens"
      );

      // Remove all hydrogens to test valence collector
      diazepam.RemoveH();

      BCL_ExampleIndirectCheck
      (
        chemistry::CollectorValence().Collect( diazepam).GetSize(), 13, "open valences after removal of hydrogens"
      );

      io::File::CloseClearFStream( input_sdf);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleChemistryCollectorValence

  const ExampleClass::EnumType ExampleChemistryCollectorValence::s_Instance
  (
    GetExamples().AddEnum( ExampleChemistryCollectorValence())
  );

} // namespace bcl
