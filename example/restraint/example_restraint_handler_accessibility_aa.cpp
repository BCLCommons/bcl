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
#include "restraint/bcl_restraint_handler_accessibility_aa.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    //! single instance of class Map< restraint::AccessibilityAA::EnvironmentType, double>
    template<> const util::SiPtr< const util::ObjectInterface>
    Map< restraint::AccessibilityAA::EnvironmentEnum, double>::s_Instance
    (
      GetObjectInstances().AddInstance( new Map< restraint::AccessibilityAA::EnvironmentEnum, double>())
    );
  } // namespace storage

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_handler_accessibility_aa.cpp
  //!
  //! @author alexanns
  //! @date
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintHandlerAccessibilityAA :
    public ExampleInterface
  {
  public:

    ExampleRestraintHandlerAccessibilityAA *Clone() const
    {
      return new ExampleRestraintHandlerAccessibilityAA( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "test default constructor");
      // create HandlerAtomDistanceAssigned "handler"
      restraint::HandlerAccessibilityAA handler_def;

      // test constuctor taking member variable parameter
      restraint::HandlerAccessibilityAA handler
      (
        util::ShPtr< assemble::AAExposureInterface>( new assemble::AANeighborVector())
      );

      // create string "restraint_filename" which has path for example restraint file
      const std::string restraint_filename
      (
        AddExampleInputPathToFilename( e_Biology, "restraint_handler_accessibility_aa_restraint_file.txt")
      );

      // create stream to "restraint_filename"
      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, restraint_filename);

      // test ReadRestraints function
      BCL_MessageStd( "test CreateRestraints function");
      // get the restraints from the handler
      const restraint::AccessibilityProfile restraints( handler.ReadRestraints( read));
      io::File::CloseClearFStream( read);

    /////////////////
    // data access //
    /////////////////

      // output restraints
      BCL_MessageStd
      (
        "The restraints in " + restraint_filename + " : \n" + util::Format()( restraints)
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // make sure correct number of restraints was created
      BCL_MessageStd( "make sure created restraints are correct");
      size_t correct_number_of_restraints( 2);
      BCL_ExampleCheck( restraints.GetAccessibilities().GetSize(), correct_number_of_restraints);

      // make sure that the restraints were created correctly

      // check the chain and resi id of the first restraint
      {
        const util::ShPtr< assemble::LocatorAA> restraint( restraints.GetAccessibilities().Begin()->GetAA());
        // make sure it is defined
        BCL_ExampleAssert( restraint.IsDefined(), true);
        BCL_ExampleCheck( restraint->GetLocatorChain().GetChainID(), 'A');
        BCL_ExampleCheck( restraint->GetAAID(), 49);
      }

      // check the chain and resi id of the second restraint
      {
        const util::ShPtr< assemble::LocatorAA> restraint( ( ++restraints.GetAccessibilities().Begin())->GetAA());
        // make sure it is defined
        BCL_ExampleAssert( restraint.IsDefined(), true)
        BCL_ExampleCheck( restraint->GetLocatorChain().GetChainID(), 'B');
        BCL_ExampleCheck( restraint->GetAAID(), 29);
      }

      // check NiEDDA accessibility of first restraint
      {
        const storage::Pair< bool, double> access_type
        (
          restraints.GetAccessibilities().Begin()->GetAccessibilityByEnvironment
          (
            restraint::AccessibilityAA::e_NiEDDA
          )
        );
        BCL_ExampleAssert( access_type.First(), true);
        const double restraint_accessibility( access_type.Second());
        double correct_accessibility( 0.98);
        BCL_ExampleCheck( restraint_accessibility, correct_accessibility);
      }

      // check Oxygen accessibility of first restraint
      {
        const storage::Pair< bool, double> access_type
        (
          restraints.GetAccessibilities().Begin()->GetAccessibilityByEnvironment
          (
            restraint::AccessibilityAA::e_Oxygen
          )
        );
        BCL_ExampleAssert( access_type.First(), true);
        const double restraint_accessibility( access_type.Second());
        double correct_accessibility( 0.43);
        BCL_ExampleCheck( restraint_accessibility, correct_accessibility);
      }

      // check Oxygen accessibility of second restraint
      {
        const storage::Pair< bool, double> access_type
        (
          ( ++restraints.GetAccessibilities().Begin())->GetAccessibilityByEnvironment
          (
            restraint::AccessibilityAA::e_Oxygen
          )
        );
        BCL_ExampleAssert( access_type.First(), true);
        const double restraint_accessibility( access_type.Second());
        double correct_accessibility( 0.12);
        BCL_ExampleCheck( restraint_accessibility, correct_accessibility);
      }

      // check NiEDDA accessibility of second restraint
      {
        const storage::Pair< bool, double> access_type
        (
          ( ++restraints.GetAccessibilities().Begin())->GetAccessibilityByEnvironment
          (
            restraint::AccessibilityAA::e_NiEDDA
          )
        );
        BCL_ExampleAssert( access_type.First(), true);
        const double restraint_accessibility( access_type.Second());
        double correct_accessibility( 0.54);
        BCL_ExampleCheck( restraint_accessibility, correct_accessibility);
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // WriteRestraints
      io::OFStream write;
      const std::string output_file
      (
        AddExampleOutputPathToFilename( handler, "restraint_handler_accessibility_aa_WriteRestraints.access_bcl")
      );

      BCL_ExampleMustOpenOutputFile( write, output_file);
      restraint::HandlerAccessibilityAA::WriteRestraints( write, restraints);
      io::File::CloseClearFStream( write);
      BCL_ExampleCheck( io::File::FilesMatch( output_file, restraint_filename), true);

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDistanceAssigned

  const ExampleClass::EnumType ExampleRestraintHandlerAccessibilityAA::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintHandlerAccessibilityAA())
  );

} // namespace bcl
