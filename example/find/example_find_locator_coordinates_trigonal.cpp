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
#include "find/bcl_find_locator_coordinates_trigonal.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_find_locator_coordinates_trigonal.cpp
  //! @details example for testing the location of coordinates for a point defining trigonal geometry
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by butkiem1 on Nov 6, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFindLocatorCoordinatesTrigonal :
    public ExampleInterface
  {
  public:

    ExampleFindLocatorCoordinatesTrigonal *Clone() const
    {
      return new ExampleFindLocatorCoordinatesTrigonal( *this);
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
      // create string "pdb_filename" which has path for example pdb file
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2CUU_spl_b_0001.pdb"));

      // get the protein model from "pdb_filename"
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

      // create ShPtr to LocatorInterface "ca_locator" and initialize with a LocatorAtom object
      util::ShPtr< assemble::LocatorAtom> ca_locator
      (
        new assemble::LocatorAtom( 'A', 14, biol::GetAtomTypes().CA)
      );

      // create ShPtr to LocatorInterface "c_locator" and initialize with a LocatorAtom object
      util::ShPtr< assemble::LocatorAtom> c_locator
      (
        new assemble::LocatorAtom( 'A', 13, biol::GetAtomTypes().C)
      );

      // create ShPtr to LocatorInterface "n_locator" and initialize with a LocatorAtom object
      util::ShPtr< assemble::LocatorAtom> n_locator
      (
        new assemble::LocatorAtom( 'A', 14, biol::GetAtomTypes().N)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      find::LocatorCoordinatesTrigonal< assemble::ProteinModel> def_constr;

      // test parameter constructor
      // this is checked in the operations section
      find::LocatorCoordinatesTrigonal< assemble::ProteinModel> param_constr
      (
        n_locator,
        ca_locator,
        c_locator,
        0.98
      );

      // create linal::Vector3D "correct_coordinates" initialize with the coordinates that are expected to be located
      linal::Vector3D correct_coordinates
      (
        linal::CoordinatesTrigonal
        (
          linal::Vector3D( 46.698, 15.496, 17.574),
          linal::Vector3D( 47.269, 16.732, 18.110),
          linal::Vector3D( 47.435, 14.573, 16.957),
          0.98
        )
      );

      BCL_MessageDbg( "the correct coordinates are " + util::Format()( correct_coordinates));

      // create const double "x_coord" and initialize with the correct x coordinate
      const double x_coord( correct_coordinates.X());
      // create const double "y_coord" and initialize with the correct y coordinate
      const double y_coord( correct_coordinates.Y());
      // create const double "z_coord" and initialize with the correct z coordinate
      const double z_coord( correct_coordinates.Z());

      // test clone constructor
      util::ShPtr< find::LocatorCoordinatesTrigonal< assemble::ProteinModel> > clone_constr
      (
        param_constr.Clone()
      );

      BCL_ExampleCheck( clone_constr->Locate( protein_model).X(), x_coord);
      BCL_ExampleCheck( clone_constr->Locate( protein_model).Y(), y_coord);
      BCL_ExampleCheck( clone_constr->Locate( protein_model).Z(), z_coord);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( param_constr.GetClassIdentifier(), GetStaticClassName( param_constr));

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test Locate atom function
      BCL_MessageStd( "test Locate function");
      BCL_ExampleCheck( param_constr.Locate( protein_model).X(), x_coord);
      BCL_ExampleCheck( param_constr.Locate( protein_model).Y(), y_coord);
      BCL_ExampleCheck( param_constr.Locate( protein_model).Z(), z_coord);

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( param_constr);

      // read the file back
      find::LocatorCoordinatesTrigonal< assemble::ProteinModel> read_object;
      ReadBCLObject( read_object);
      BCL_ExampleCheck( read_object.Locate( protein_model).X(), x_coord);
      BCL_ExampleCheck( read_object.Locate( protein_model).Y(), y_coord);
      BCL_ExampleCheck( read_object.Locate( protein_model).Z(), z_coord);

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFindLocatorCoordinatesTrigonal

  const ExampleClass::EnumType ExampleFindLocatorCoordinatesTrigonal::s_Instance
  (
    GetExamples().AddEnum( ExampleFindLocatorCoordinatesTrigonal())
  );

} // namespace bcl
