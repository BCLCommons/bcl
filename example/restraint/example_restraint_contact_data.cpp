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
#include "restraint/bcl_restraint_contact_data.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "io/bcl_io_file.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "score/bcl_score_restraint_atom_distance.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_contact_data.cpp
  //!
  //! @author weinerbe
  //! @date Mar 19, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintContactData :
    public ExampleInterface
  {
  public:

    ExampleRestraintContactData *Clone() const
    {
      return new ExampleRestraintContactData( *this);
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

      // create object
      restraint::ContactData data_construct;

      // read in restraints
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Biology, "contact.RR"));
      data_construct.ReadRestraints( read);
      io::File::CloseClearFStream( read);

    /////////////////
    // data access //
    /////////////////

      // test GetDataConstruct
      const size_t nr_restraints( 4);
      BCL_ExampleCheck( data_construct.GetAtomDistanceRestraints().GetSize(), nr_restraints);

      // test GetAllFlags
      BCL_ExampleCheck( data_construct.GetLabel().GetNumberArguments() >= 4, true);

      // have a look at the function
      BCL_MessageStd
      (
        "Label of function used for scoring contacts: " + restraint::ContactData::s_ScoreContactDistance->GetString()
        + "\nActual function help: "
      );
      restraint::ContactData::s_ScoreContactRestraintAssignment->WriteHelp( util::GetLogger());

    //////////////////////
    // input and output //
    //////////////////////

      // test WriteRestraints
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, AddExampleOutputPathToFilename( data_construct, "contact_restraints.txt"));
      restraint::HandlerAtomDistanceAssigned().WriteRestraints( write, data_construct.GetAtomDistanceRestraints(), false, true);
      io::File::CloseClearFStream( write);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintContactData

  const ExampleClass::EnumType ExampleRestraintContactData::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintContactData())
  );

} // namespace bcl
