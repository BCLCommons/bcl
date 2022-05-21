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
#include "find/bcl_find_locator_coordinates_average.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    template class LocatorCoordinatesKnown< assemble::ProteinModel>;
  } // namespace find
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_find_locator_coordinates_average.cpp
  //! TODO: add an detailed description for this example
  //!
  //! @author alexanns
  //! @date Feb 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFindLocatorCoordinatesAverage :
    public ExampleInterface
  {
  public:

    ExampleFindLocatorCoordinatesAverage *Clone() const
    {
      return new ExampleFindLocatorCoordinatesAverage( *this);
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
      // create protein
      const util::ShPtr< assemble::ProteinModel> ubi
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")).Clone()
      );

      // create ShPtr to LocatorInterface "cb_locator" and initialize with a LocatorAtom object
      const util::ShPtr< assemble::LocatorAtom> cb_locator
      (
        new assemble::LocatorAtom( 'A', 14, biol::GetAtomTypes().CB)
      );

      // create ShPtr to LocatorInterface "ca_locator" and initialize with a LocatorAtom object
      const util::ShPtr< assemble::LocatorAtom> ca_locator
      (
        new assemble::LocatorAtom( 'A', 20, biol::GetAtomTypes().CA)
      );

      // create locator list
      storage::List< util::Implementation< find::LocatorCoordinatesInterface< assemble::ProteinModel> > > locator_list;
      locator_list.PushBack( util::Implementation< find::LocatorCoordinatesInterface< assemble::ProteinModel> >( *cb_locator));
      locator_list.PushBack( util::Implementation< find::LocatorCoordinatesInterface< assemble::ProteinModel> >( *ca_locator));

      // expected coordinate average
      const linal::Vector3D coords_average( 29.138, 25.6125, 7.1485);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const find::LocatorCoordinatesAverage< assemble::ProteinModel> default_constr;
      BCL_ExampleCheck( default_constr.Locate( *ubi).IsDefined(), false);

      // constructor taking parameter
      const find::LocatorCoordinatesAverage< assemble::ProteinModel> param_constr( locator_list);
      BCL_ExampleCheckWithinTolerance( param_constr.Locate( *ubi), coords_average, 0.0001);

      // constructor from atom locator list
      {
        // create atom locator list
        storage::List< util::Implementation< find::LocatorCoordinatesInterface< assemble::ProteinModel> > > locator_list;
        locator_list.PushBack( *cb_locator);
        locator_list.PushBack( *ca_locator);
        const find::LocatorCoordinatesAverage< assemble::ProteinModel> param_constr( locator_list);
        BCL_ExampleCheckWithinTolerance( param_constr.Locate( *ubi), coords_average, 0.0001);
      }

      // clone constructor
      const util::ShPtr< find::LocatorCoordinatesAverage< assemble::ProteinModel> > clone_constr( param_constr.Clone());
      BCL_ExampleCheckWithinTolerance( clone_constr->Locate( *ubi), coords_average, 0.0001);

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd( "test GetStaticClassName GetClassIdentifier");
      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::find::LocatorCoordinatesAverage<bcl::assemble::ProteinModel>");
      BCL_ExampleCheck
      (
        GetStaticClassName< find::LocatorCoordinatesAverage< assemble::ProteinModel> >(), correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< find::LocatorCoordinatesAverage< assemble::ProteinModel> >(),
        clone_constr->GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

      // Locate
      BCL_ExampleCheckWithinTolerance( param_constr.Locate( *ubi), coords_average, 0.0001);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities");
      // write the object
      WriteBCLObject( param_constr);

      // read the object back in
      find::LocatorCoordinatesAverage< assemble::ProteinModel> read;
      ReadBCLObject( read);

      BCL_ExampleCheckWithinTolerance( read.Locate( *ubi), coords_average, 0.0001);
      BCL_MessageDbg( "located average coordinates are " + util::Format()( read.Locate( *ubi)));
      BCL_MessageDbg( "expected average coordinates are " + util::Format()( coords_average));

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFindLocatorCoordinatesAverage

  const ExampleClass::EnumType ExampleFindLocatorCoordinatesAverage::s_Instance
  (
    GetExamples().AddEnum( ExampleFindLocatorCoordinatesAverage())
  );

} // namespace bcl
