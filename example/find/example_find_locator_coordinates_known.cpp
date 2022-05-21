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
#include "find/bcl_find_locator_coordinates_known.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_find_locator_coordinates_known.cpp
  //! @details example for class find::LocatorCoordinatesKnown
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by heinzes1
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleFindLocatorCoordinatesKnown :
    public ExampleInterface
  {
  public:

    ExampleFindLocatorCoordinatesKnown *Clone() const
    {
      return new ExampleFindLocatorCoordinatesKnown( *this);
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
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));

      // get the protein model from "pdb_filename"
      assemble::ProteinModel protein_model( Proteins::GetModel( pdb_filename));

      // create ShPtr to LocatorInterface "atom_locator" and initialize with a LocatorAtom object
      util::ShPtr< assemble::LocatorAtom> atom_locator
      (
        new assemble::LocatorAtom( 'A', 14, biol::GetAtomTypes().CB)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      find::LocatorCoordinatesKnown< assemble::ProteinModel> def_constr;

      // test parameter constructor
      // this is checked in the operations section below
      find::LocatorCoordinatesKnown< assemble::ProteinModel> param_constr( atom_locator);

      // create const double "x_coord" and initialize with the correct x coordinate of the atom
      const double x_coord( 3.259);
      // create const double "y_coord" and initialize with the correct y coordinate of the atom
      const double y_coord( 34.378);
      // create const double "z_coord" and initialize with the correct z coordinate of the atom
      const double z_coord( -26.828);

      // test clone constructor
      util::ShPtr< find::LocatorCoordinatesKnown< assemble::ProteinModel> > clone_constr( param_constr.Clone());
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
      find::LocatorCoordinatesKnown< assemble::ProteinModel> read_object;
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

  }; //end ExampleFindLocatorCoordinatesKnown

  const ExampleClass::EnumType ExampleFindLocatorCoordinatesKnown::s_Instance
  (
    GetExamples().AddEnum( ExampleFindLocatorCoordinatesKnown())
  );

} // namespace bcl
