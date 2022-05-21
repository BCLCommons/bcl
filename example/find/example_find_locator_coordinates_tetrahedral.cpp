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
#include "find/bcl_find_locator_coordinates_tetrahedral.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_find_locator_coordinates_tetrahedral.cpp
  //! @details example for class find::LocatorCoordinatesTetrahedral
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by weinerbe on Nov 6, 2010
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFindLocatorCoordinatesTetrahedral :
    public ExampleInterface
  {
  public:

    ExampleFindLocatorCoordinatesTetrahedral *Clone() const
    {
      return new ExampleFindLocatorCoordinatesTetrahedral( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

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
        new assemble::LocatorAtom( 'A', 14, biol::GetAtomTypes().C)
      );

      // create ShPtr to LocatorInterface "cb_locator" and initialize with a LocatorAtom object
      util::ShPtr< assemble::LocatorAtom> cb_locator
      (
        new assemble::LocatorAtom( 'A', 14, biol::GetAtomTypes().CB)
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
      find::LocatorCoordinatesTetrahedral< assemble::ProteinModel> def_constr;

      // test parameter constructor
      // that it was constructed properly is tested in the operations section
      find::LocatorCoordinatesTetrahedral< assemble::ProteinModel> param_constr
      (
        ca_locator,
        c_locator,
        cb_locator,
        n_locator,
        0.98
      );

      // create linal::Vector3D "correct_coordinates" initialize with the coordinates that are expected to be located
      linal::Vector3D correct_coordinates
      (
        linal::CoordinatesTetrahedral
        (
          linal::Vector3D( 30.599,  33.927,   6.546),
          linal::Vector3D( 31.227,  34.987,   7.306),
          linal::Vector3D( 31.487,  32.695,   6.498),
          linal::Vector3D( 30.193,  34.380,   5.135),
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
      util::ShPtr< find::LocatorCoordinatesTetrahedral< assemble::ProteinModel> > clone_constr
      (
        param_constr.Clone()
      );
      BCL_ExampleCheckWithinAbsTolerance( clone_constr->Locate( protein_model).X(), x_coord, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( clone_constr->Locate( protein_model).Y(), y_coord, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( clone_constr->Locate( protein_model).Z(), z_coord, 0.001);

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
      BCL_ExampleCheckWithinAbsTolerance( param_constr.Locate( protein_model).X(), x_coord, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( param_constr.Locate( protein_model).Y(), y_coord, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( param_constr.Locate( protein_model).Z(), z_coord, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( param_constr);

      // read the file back
      find::LocatorCoordinatesTetrahedral< assemble::ProteinModel> read_object;
      ReadBCLObject( read_object);
      BCL_ExampleCheckWithinAbsTolerance( read_object.Locate( protein_model).X(), x_coord, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( read_object.Locate( protein_model).Y(), y_coord, 0.001);
      BCL_ExampleCheckWithinAbsTolerance( read_object.Locate( protein_model).Z(), z_coord, 0.001);

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleFindLocatorCoordinatesTetrahedral

  const ExampleClass::EnumType ExampleFindLocatorCoordinatesTetrahedral::s_Instance
  (
    GetExamples().AddEnum( ExampleFindLocatorCoordinatesTetrahedral())
  );

} // namespace bcl
