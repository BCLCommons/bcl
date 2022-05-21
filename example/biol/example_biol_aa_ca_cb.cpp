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
#include "biol/bcl_biol_aa_ca_cb.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_classes.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_ca_cb.cpp
  //!
  //! @author rouvelgh
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAACaCb :
    public ExampleInterface
  {
  public:

    ExampleBiolAACaCb *Clone() const
    {
      return new ExampleBiolAACaCb( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      // example prerequisites
      const linal::Vector3D CA_coordinates( 0.5, 0.2, 0.3);
      const linal::Vector3D CB_coordinates( 1.5, 2.2, 9.3);
      const biol::Atom def_CA( CA_coordinates, biol::AtomType( "CA"), 7, 0.242);
      const biol::Atom def_CB( CB_coordinates, biol::AtomType("CB"), 8, 0.532);
      const util::SiPtrVector< const biol::Atom> siptr_vec_atom( util::SiPtrVector< const biol::Atom>::Create( def_CA, def_CB));
      const biol::AA aa_obj( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 2, 3, ' ')));
      const storage::Set< biol::AtomType> correct_atom_types( biol::GetAtomTypes().CA, biol::GetAtomTypes().CB);
      const linal::Vector3D new_translate( 0.785, -0.325, 0.672);
      const math::TransformationMatrix3D transformation_matrix( 1.0, 0.0, 1.0);
      const math::RotationMatrix3D rotation_matrix( 30, 60, 40);
      const linal::Vector3D expected_coord_ca( 1.285,-0.125, 0.972);
      const linal::Vector3D expected_coord_cb( 2.285, 1.875, 9.972);
      const linal::Vector3D expected_transform_ca( 2.285,-0.125, 1.972);
      const linal::Vector3D expected_transform_cb( 3.285, 1.875, 10.972);
      const linal::Vector3D expected_rotation_ca( -0.0514414, 2.89285, -0.868574);
      const linal::Vector3D expected_rotation_cb( -6.49606, 8.56325, -4.37761);
      const linal::Vector3D expected_center_of_mass( -3.27375, 5.72805, -2.62309);
      const linal::Vector3D ideal_coord_ca( 3.26, 0, 1.661);
      const linal::Vector3D ideal_coord_fsa( 4.06932, 1.03771, 0.849);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      biol::AACaCb def_constr;

      // test constructor from AABase
      biol::AACaCb aa_base_constr( aa_obj);

      // construct from aabase and two atoms
      biol::AACaCb new_constr( aa_obj, def_CA, def_CB);

      // construct from aabase and siptrvec to atoms
      biol::AACaCb aacacb_atoms( aa_obj, siptr_vec_atom);

      //test copy constructor
      biol::AACaCb copy_construct( new_constr);

      // test clone constructor
      util::ShPtr< util::ObjectInterface> ptr( copy_construct.Clone());

      // testing the empty constructor
      util::ShPtr< biol::AACaCb> empty_construct( def_constr.Empty( util::ShPtr< biol::AAData>( new biol::AAData())));

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_Example_Check
      (
        GetStaticClassName( def_constr) == "bcl::biol::AACaCb",
        "unexpected static class name: " + GetStaticClassName( def_constr) + " should be: bcl::biol::AACaCb"
      );

      // class identifier
      BCL_Example_Check
      (
        ptr->GetClassIdentifier() == GetStaticClassName< biol::AACaCb>(),
        "unexpected class identifier class name: " + ptr->GetClassIdentifier() + " should be: " + GetStaticClassName< biol::AACaCb>()
      );

      // GetType
      BCL_MessageStd( "Testing GetType along with class constructors");
      BCL_Example_Check
      (
        !def_constr.GetType().IsDefined() &&
        aa_base_constr.GetType() == biol::GetAATypes().ALA &&
        new_constr.GetType() == biol::GetAATypes().ALA &&
        aacacb_atoms.GetType() == biol::GetAATypes().ALA &&
        copy_construct.GetType() == biol::GetAATypes().ALA &&
        !empty_construct->GetType().IsDefined()
        ,
        "The AAtype of default constructed aacacb should be undefined but is " + util::Format()( def_constr.GetType())
      );

      // check GetNumberOfAtoms
      BCL_MessageStd( "Testing GetNumberOfAtoms");
      BCL_Example_Check
      (
        new_constr.GetNumberOfAtoms() == 2,
        "Number of Atoms are " + util::Format()( new_constr.GetNumberOfAtoms()) + " but should be " + util::Format()( 2)
      );

      // test get types of atoms
      BCL_MessageStd( "Testing GetTypesOfAtoms");
      BCL_Example_Check
      (
        new_constr.GetTypesOfAtoms().InternalData() == correct_atom_types.InternalData(),
        "Types of Atoms are " + util::Format()( new_constr.GetTypesOfAtoms()) + " but should be "
        + util::Format()( correct_atom_types)
      );

      // test get all atoms
      BCL_MessageStd( "Testing GetAtoms");
      BCL_Example_Check
      (
        new_constr.GetAtoms()( 0).GetPointer() == &new_constr.GetCA() &&
        new_constr.GetAtoms()( 1).GetPointer() == &new_constr.GetFirstSidechainAtom() &&
        aacacb_atoms.GetAtoms()( 0).GetPointer() == &aacacb_atoms.GetCA() &&
        aacacb_atoms.GetAtoms()( 1).GetPointer() == &aacacb_atoms.GetFirstSidechainAtom() &&
        aa_base_constr.GetAtoms()( 0).GetPointer() == &aa_base_constr.GetCA()  &&
        aa_base_constr.GetAtoms()( 1).GetPointer() == &aa_base_constr.GetFirstSidechainAtom() &&
        copy_construct.GetAtoms()( 0).GetPointer() == &copy_construct.GetCA()  &&
        copy_construct.GetAtoms()( 1).GetPointer() == &copy_construct.GetFirstSidechainAtom() &&
        def_constr.GetAtoms().IsDefined(),
        "The function returned for new_constr "
        + util::Format()( *new_constr.GetAtoms()( 0).GetPointer())
        + util::Format()( *new_constr.GetAtoms()( 1).GetPointer())
        +" but should have returned "
        + util::Format()( def_CA) + util::Format()( def_CB) + util::Format()( aa_obj)
        +" The function returned for aacacb_atoms"
        + util::Format()( *aacacb_atoms.GetAtoms()( 0).GetPointer())
        + util::Format()( *aacacb_atoms.GetAtoms()( 1).GetPointer())
        +" but should have returned "
        + util::Format()( aa_obj) + util::Format()( siptr_vec_atom)
        +" The function returned for aa_base_constr "
        + util::Format()( *aa_base_constr.GetAtoms()( 0).GetPointer())
        + util::Format()( *aa_base_constr.GetAtoms()( 1).GetPointer())
        +" but should have returned "
        + util::Format()( aa_obj)
        +" The function returned for copy_construct "
        + util::Format()( *copy_construct.GetAtoms()( 0).GetPointer())
        + util::Format()( *copy_construct.GetAtoms()( 1).GetPointer())
        + " but should have returned "
        + util::Format()( def_CA) + util::Format()( def_CB) + util::Format()( aa_obj)
      );

      // test get specific atom type
      BCL_MessageStd( "Testing get GetAtom");
      BCL_Example_Check
      (
        new_constr.GetAtom( biol::GetAtomTypes().CA).GetType() == def_CA.GetType() &&
        new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates() == def_CA.GetCoordinates() &&
        new_constr.GetAtom( biol::GetAtomTypes().CB).GetType() == def_CB.GetType() &&
        new_constr.GetAtom( biol::GetAtomTypes().CB).GetCoordinates() == def_CB.GetCoordinates(),
        "Results from GetAtom are " + util::Format()( new_constr.GetAtom( biol::GetAtomTypes().CA)) +
        util::Format()( new_constr.GetAtom( biol::GetAtomTypes().CB)) + " but should be "
        + util::Format()( def_CA) + util::Format()( def_CB)
      );

      // test get the CA
      BCL_MessageStd( "Testing GetCA");
      BCL_Example_Check
      (
        new_constr.GetCA().GetType() == def_CA.GetType() &&
        new_constr.GetCA().GetCoordinates() == def_CA.GetCoordinates(),
        "The CA atom is " + util::Format()( new_constr.GetCA()) + " but should be "
        + util::Format()( def_CA)
      );

      // test GetFirstSideChainAtom (gets the CB)
      BCL_MessageStd( "Testing GetFirstSideChainAtom");
      BCL_Example_Check
      (
        new_constr.GetFirstSidechainAtom().GetType() == def_CB.GetType() &&
        new_constr.GetFirstSidechainAtom().GetCoordinates() == def_CB.GetCoordinates(),
        "First Sidechain Atom is " + util::Format()( new_constr.GetFirstSidechainAtom()) + " but should be "
        + util::Format()( def_CB)
      );

      // test set atom
      BCL_MessageStd( "Testing SetAtom");
      new_constr.SetAtom( def_CA);
      BCL_Example_Check
      (
        new_constr.GetAtom( biol::GetAtomTypes().CA).GetType() == def_CA.GetType() &&
        new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates() == def_CA.GetCoordinates(),
        "The set atom is" + util::Format()( biol::GetAtomTypes().CA) + " but should be "
        + util::Format()( def_CA)
      );

      // test GetAAClass
      BCL_MessageStd( "Testing GetAAClass");
      BCL_Example_Check
      (
        new_constr.GetAAClass() == def_constr.GetAAClass(),
        "The AAClass is " + util::Format()( new_constr.GetAAClass()) + " but should be "
        + util::Format()( def_constr.GetAAClass())
      );

      // test CalculateOmega
      BCL_MessageStd( "Testing CalculateOmega");
      BCL_Example_Check
      (
        !util::IsDefined( new_constr.CalculateOmega( def_CA, def_CB)),
        "The Omega angle is " + util::Format()( new_constr.CalculateOmega( def_CA, def_CB)) + " but should be "
        + util::Format()( util::GetUndefinedDouble())
      );

      // test CalculatePhi
      BCL_MessageStd( "Testing CalculatePhi");
      BCL_Example_Check
      (
        !util::IsDefined( new_constr.CalculatePhi( def_CA)),
        "The Phi angle is " + util::Format()( new_constr.CalculatePhi( def_CA)) + " but should be "
        + util::Format()( util::GetUndefinedDouble())
      );

      // test CalculatePsi
      BCL_MessageStd( "Testing CalculatePsi");
      BCL_Example_Check
      (
        !util::IsDefined( new_constr.CalculatePsi( def_CA)) &&
        !util::IsDefined( def_constr.CalculatePsi( def_CA)),
        "The Psi angle is " + util::Format()( new_constr.CalculatePsi( def_CA)) + " but should be "
        + util::Format()( util::GetUndefinedDouble())
      );

    ////////////////
    // operations //
    ////////////////

      // check if the SiPtr point to the correct linal::Vector3Ds after using GetAtomCoordinates
      BCL_Example_Check
      (
        new_constr.GetAtomCoordinates()( 0).GetPointer() == &new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates() &&
        new_constr.GetAtomCoordinates()( 1).GetPointer() == &new_constr.GetAtom( biol::GetAtomTypes().CB).GetCoordinates() &&
        *new_constr.GetAtomCoordinates()( 0) == CA_coordinates &&
        *new_constr.GetAtomCoordinates()( 1) == CB_coordinates,
        "GetAtomCoordinates returned" + util::Format()( *new_constr.GetAtomCoordinates()( 0)) + " as the CA coordinates and " +
        util::Format()( *new_constr.GetAtomCoordinates()( 1)) + " as the CB coordinates, but should have returned" +
        util::Format()( def_CA.GetCoordinates())
        + " as the CA coordinates and for the CB coordinates:  " + util::Format()( def_CB.GetCoordinates())
      );

      // test GetAtomCoordinates for specified atom types
      BCL_MessageStd( "Testing GetAtomCoordinates( AtomTypes)");
      BCL_Example_Check
      (
        new_constr.GetAtomCoordinates( correct_atom_types)( 0).GetPointer() ==
        &new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates() &&
        new_constr.GetAtomCoordinates( correct_atom_types)( 1).GetPointer() ==
        &new_constr.GetAtom( biol::GetAtomTypes().CB).GetCoordinates() &&
        *new_constr.GetAtomCoordinates( correct_atom_types)( 0) == CA_coordinates &&
        *new_constr.GetAtomCoordinates( correct_atom_types)( 1) == CB_coordinates,
        "error"
      );

      // test translate
      BCL_MessageStd( "Testing Translate");
      new_constr.Translate( new_translate);
      // check if translate occurred according to translate function
      BCL_Example_Check
      (
         math::EqualWithinTolerance( new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), expected_coord_ca) &&
         math::EqualWithinTolerance( new_constr.GetAtom( biol::GetAtomTypes().CB).GetCoordinates(), expected_coord_cb),
         " Translation: the coordinates resulting from the translate should be " + util::Format()( expected_coord_ca) + " "
         + util::Format()( expected_coord_cb)
         + "but are " + util::Format()( *new_constr.GetAtomCoordinates()( 0)) + util::Format()( *new_constr.GetAtomCoordinates()( 1))
      );

      // test Transform according to transformationMatrix3D
      BCL_MessageStd( "Testing Transform");
      new_constr.Transform( transformation_matrix);
      // check if transformation occurred according to transform function
      BCL_Example_Check
      (
        math::EqualWithinTolerance( new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), expected_transform_ca) &&
        math::EqualWithinTolerance( new_constr.GetAtom( biol::GetAtomTypes().CB).GetCoordinates(), expected_transform_cb),
        " Transformation: the coordinates resulting from the transformation should be " + util::Format()( expected_transform_ca)
        + " for the CA coordinates" + util::Format()( expected_transform_cb) + " for the CB coordinates, but are "
        + util::Format()( *new_constr.GetAtomCoordinates()( 0)) + util::Format()( *new_constr.GetAtomCoordinates()( 1))
      );

      // test Rotation according to transformationMatrix3D
      BCL_MessageStd( "Testing Rotate");
      new_constr.Rotate( rotation_matrix);
      // check if Rotation occurred properly
      BCL_Example_Check
      (
        math::EqualWithinTolerance( new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), expected_rotation_ca) &&
        math::EqualWithinTolerance( new_constr.GetAtom( biol::GetAtomTypes().CB).GetCoordinates(), expected_rotation_cb),
        " Rotation: the coordinates resulting from the rotation should be " + util::Format()( expected_rotation_ca) +
        " for the CA coordinates and " + util::Format()( expected_rotation_cb) + " for the CB coordinates, but are "
        + util::Format()( *new_constr.GetAtomCoordinates()( 0)) + util::Format()( *new_constr.GetAtomCoordinates()( 1))
      );

      // test GetCenter
      BCL_MessageStd( "Testing GetCenter");
      // check if GetCenter returned the proper center of mass
      BCL_Example_Check
      (
        math::EqualWithinTolerance( new_constr.GetCenter(), expected_center_of_mass),
        " Center of Mass: the coordinates resulting from the center of mass calculation should be "
        + util::Format()( expected_center_of_mass) + " but are " + util::Format()( new_constr.GetCenter())
      );

      // test SetToIdealConformation
      BCL_MessageStd( "Testing SetToIdealConformation")
      // invoke the function on the object
      new_constr.SetToIdealConformation( biol::GetSSTypes().HELIX, transformation_matrix);
      // check to see if the AA was properly idealized according to transformation_matrix
      BCL_Example_Check
      (
        math::EqualWithinTolerance( new_constr.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), ideal_coord_ca) &&
        math::EqualWithinTolerance( new_constr.GetFirstSidechainAtom().GetCoordinates(), ideal_coord_fsa),
        " Idealize: the coordinates resulting from the idealizaton should be "
        + util::Format()( ideal_coord_ca)
        + " for the CA coordinates "
        + util::Format()( ideal_coord_fsa)
        + " and for the FirstSideChainAtom coordinates"
        + " but are "
        + util::Format()( *new_constr.GetAtomCoordinates()(0))
        + util::Format()( *new_constr.GetAtomCoordinates()(1))
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for biol::AACaCb");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( def_constr);
      BCL_MessageVrb( "read object");
      biol::AACaCb aacacb_read;
      ReadBCLObject( aacacb_read);

      // compare written and read object
      BCL_Example_Check
      (
        !def_constr.GetCA().GetCoordinates().IsDefined() &&
        !aacacb_read.GetCA().GetCoordinates().IsDefined() &&
        !util::IsDefined( def_constr.GetAtoms()( 0)->GetPdbID()) &&
        !util::IsDefined( aacacb_read.GetAtoms()( 0)->GetPdbID()),
        "read aacacb is different from written: " + util::Format()( def_constr) + " != \n" +
        util::Format()( aacacb_read)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAACaCb

  const ExampleClass::EnumType ExampleBiolAACaCb::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAACaCb())
  );

} // namespace bcl

