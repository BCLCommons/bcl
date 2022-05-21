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
#include "descriptor/bcl_descriptor_cheminfo_properties.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "model/bcl_model_retrieve_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_cheminfo_properties.cpp
  //!
  //! @author kothiwsk
  //! @date Feb 18, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorCheminfoProperties :
    public ExampleInterface
  {
  public:

    ExampleDescriptorCheminfoProperties *Clone() const
    {
      return new ExampleDescriptorCheminfoProperties( *this);
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
      // descriptor::CheminfoProperty is tested implicitly by other descriptor::CheminfoProperty... classes
      // Here, we just verify the basic functionality, and also test the most basic properties, misc property

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor:
      descriptor::CheminfoProperty undef;

      BCL_ExampleAssert( descriptor::CheminfoProperty().IsDefined(), false);

      // construct a property that will retrieve cached values for HAcc
      //descriptor::MoleculeMiscProperty hacc_from_file( "HAcc", 1);
      BCL_ExampleAssert( descriptor::CheminfoProperty( "HbondAcceptor").IsDefined(), true);

      // now make a property that can actually calculate the hydrogen bond acceptors
      descriptor::CheminfoProperty hacc( "HbondAcceptor");
      BCL_ExampleAssert( hacc.IsDefined(), true);
      BCL_ExampleCheck( hacc->GetCachePreference(), descriptor::e_PreferCache);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        descriptor::CheminfoProperty
        (
          descriptor::CheminfoProperty( "HbondAcceptor")
        )->GetSizeOfFeatures(),
        descriptor::CheminfoProperty( "HbondAcceptor")->GetSizeOfFeatures()
      );

      BCL_ExampleCheck
      (
        descriptor::CheminfoProperty( "HbondAcceptor")->GetSizeOfFeatures(),
        size_t( 1)
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleCheck
      (
        TestBCLObjectIOForSymmetry
        (
          descriptor::CheminfoProperty( "HbondAcceptor"),
          descriptor::CheminfoProperty()
        ),
        true
      );

      // set to ignore the terminal's line width; this way this example file won't change depending on the terminal's
      // line width.  This is necessary because currently most of the WriteHelp functions used by various classes do not
      // support taking a io::FixedLineWidthWriter as an argument
      util::GetLogger().SetIgnoreTerminalLineWidth( true);

      io::OFStream write;
      BCL_ExampleMustOpenOutputFile
      (
        write,
        GetExamples().GetExamplePath() + "/../../scripts/machine_learning/descriptors/qsar/QSARDescriptorsHelp.txt"
      );
      // reset the help flag for atom properties so that they too will be output
      descriptor::CheminfoProperty::ResetHaveDisplayedHelp();
      descriptor::CheminfoID::ResetHaveDisplayedHelp();
      util::Implementation< io::RetrieveInterface< chemistry::MoleculeComplete, chemistry::MoleculeEnsemble> >::ResetHaveDisplayedHelp();
      util::Implementation< model::RetrieveInterface>::ResetHaveDisplayedHelp();
      io::FixedLineWidthWriter writer;
      descriptor::CheminfoProperty::WriteInstancesHelp( writer);
      writer << '\n';
      descriptor::CheminfoID::ResetHaveDisplayedHelp();
      descriptor::CheminfoID::WriteInstancesHelp( writer);
      write << writer.String();
      io::File::CloseClearFStream( write);

      // reset set ignore terminal line width
      util::GetLogger().SetIgnoreTerminalLineWidth( false);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorCheminfoProperties

  const ExampleClass::EnumType ExampleDescriptorCheminfoProperties::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorCheminfoProperties())
  );

} // namespace bcl
