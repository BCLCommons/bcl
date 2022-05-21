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
#include "biol/bcl_biol_aa_back_bone.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_back_bone.cpp
  //!
  //! @author rouvelgh
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAABackBone :
    public ExampleInterface
  {
  public:

    ExampleBiolAABackBone *Clone() const
    {
      return new ExampleBiolAABackBone( *this);
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

      // initialize objects required by class functions
      const std::string filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));
      // get the first

      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 4;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 4;
      const assemble::ProteinModel this_model
      (
        Proteins::GetModel( filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );
      // get the first amino acid
      util::SiPtr< const biol::AABase> sp_first_residue( this_model.GetSequences().FirstElement()->GetFirstAA());

      biol::AACaCb aa_ca_cb;

      // these matrices and vectors are used by the translation, rotation, transformation, and center of mass functions
      const linal::Vector3D new_translate( 0.785, -0.325, 0.672);
      const math::TransformationMatrix3D transformation_matrix( 1.0, 0.0, 1.0);
      const math::RotationMatrix3D rotation_matrix( 30, 60, 40);
      const linal::Vector3D expected_coord_n( 28.128, 23.969, 3.355);
      const linal::Vector3D expected_coord_ca( 27.166, 25.036, 3.566);
      const linal::Vector3D expected_transform_n( 29.128, 23.969, 4.355);
      const linal::Vector3D expected_transform_ca( 28.166, 25.036, 4.566);
      const linal::Vector3D expected_rotation_n( 13.1, 15.5435, -32.0735);
      const linal::Vector3D expected_rotation_ca( 12.2353, 14.6433, -32.8153);
      const linal::Vector3D expected_center_of_mass( 12.0517, 15.2011, -33.1466);
      const linal::Vector3D ideal_coord_n( 2.34703, -0.70421, 0.755);
      const linal::Vector3D ideal_coord_ca( 3.26, 0, 1.661);
      const linal::Vector3D ideal_coord_c( 2.50385, 0.744392, 2.747);
      const linal::Vector3D ideal_coord_o( 2.82947, 0.674928, 3.922);
      const linal::Vector3D ideal_coord_fsa( 4.06932, 1.03771, 0.849);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      biol::AABackBone def_constr;

      // test constructor from AABase
      biol::AABackBone aa_base_constr( aa_ca_cb);

      // construct from aabase and siptrvec to atoms
      biol::AABackBone aabackbone_atom
      (
        biol::AA( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().MET, 1, 2))),
        sp_first_residue->GetAtoms()
      );
      storage::Set< biol::AtomType> correct_atom_type
        (
          storage::Set< biol::AtomType>::Create
          (
            biol::GetAtomTypes().N, biol::GetAtomTypes().CA, biol::GetAtomTypes().O, biol::GetAtomTypes().C,
            biol::GetAtomTypes().HA2, biol::GetAtomTypes().CB
          )
        );

      // test copy constructor
      biol::AABackBone copy_construct( aabackbone_atom);

      // test clone constructor
      util::ShPtr< util::ObjectInterface> ptr( copy_construct.Clone());

      // constructor from CaCb amino acids
      biol::AABackBone aa_ca_cb_constr( aa_ca_cb);

      // testing the empty constructor
      util::ShPtr< biol::AABackBone> empty_construct( def_constr.Empty( util::ShPtr< biol::AAData>( new biol::AAData())));

    /////////////////
    // data access //
    /////////////////

      // check the class identifier
      BCL_Example_Check
      (
        GetStaticClassName( def_constr) == "bcl::biol::AABackBone",
        "unexpected static class name: " + GetStaticClassName( def_constr) + " should be: bcl::biol::AABackBone"
      );

      // class identifier
      BCL_Example_Check
      (
        ptr->GetClassIdentifier() == GetStaticClassName< biol::AABackBone>(),
        "unexpected class identifier class name: " + ptr->GetClassIdentifier() + " should be: " + GetStaticClassName< biol::AABackBone>()
      );

      // GetType tests the constructors
      BCL_MessageStd( "Testing GetType along with class constructors");
      BCL_Example_Check
      (
        !def_constr.GetType().IsDefined() &&
        !aa_ca_cb.GetType().IsDefined() &&
        aabackbone_atom.GetType() == biol::GetAATypes().MET &&
        copy_construct.GetType() == biol::GetAATypes().MET &&
        !empty_construct->GetType().IsDefined(),
        "The AAtype of default constructed aabackbone should be undefined but is " + util::Format()( def_constr.GetType())
      );

      // check GetNumberOfAtoms, 5 is the default for aa_back_bone
      BCL_MessageStd( "Testing GetNumberOfAtoms");
      BCL_Example_Check
      (
        aabackbone_atom.GetNumberOfAtoms() == 5,
        "Number of Atoms are " + util::Format()( aabackbone_atom.GetNumberOfAtoms()) + " but should be " + util::Format()( 5)
      );

      // test get types of atoms, correct_atom_type is a set that contains the required AtomTypes
      BCL_MessageStd( "Testing GetTypesOfAtoms");
      BCL_Example_Check
      (
        aabackbone_atom.GetTypesOfAtoms().InternalData() == correct_atom_type.InternalData(),
        "Types of Atoms are " + util::Format()( aabackbone_atom.GetTypesOfAtoms()) + " but should be "
        + util::Format()( correct_atom_type)
      );

      // test get all atoms, asking does the pointer point to the address of the correct atom.  Done for the aabackbone_atom
      // and the copy_constructor, which may not be necessary to test, but is included.
      BCL_MessageStd( "Testing GetAtoms");
      BCL_Example_Check
      (
        aabackbone_atom.GetAtoms()( 0).GetPointer() == &aabackbone_atom.GetAtom( biol::GetAtomTypes().N)   &&
        aabackbone_atom.GetAtoms()( 1).GetPointer() == &aabackbone_atom.GetCA() &&
        aabackbone_atom.GetAtoms()( 2).GetPointer() == &aabackbone_atom.GetAtom( biol::GetAtomTypes().C) &&
        aabackbone_atom.GetAtoms()( 3).GetPointer() == &aabackbone_atom.GetAtom( biol::GetAtomTypes().O) &&
          copy_construct.GetAtoms()( 0).GetPointer() == &copy_construct.GetAtom( biol::GetAtomTypes().N)  &&
          copy_construct.GetAtoms()( 1).GetPointer() == &copy_construct.GetCA() &&
          copy_construct.GetAtoms()( 2).GetPointer() == &copy_construct.GetAtom( biol::GetAtomTypes().C) &&
          copy_construct.GetAtoms()( 3).GetPointer() == &copy_construct.GetAtom( biol::GetAtomTypes().O) &&
            def_constr.GetAtoms().IsDefined(),
         +" The function returned for aabackbone_atom "
         + util::Format()( *aabackbone_atom.GetAtoms()( 0).GetPointer())
         + util::Format()( *aabackbone_atom.GetAtoms()( 1).GetPointer())
         + util::Format()( *aabackbone_atom.GetAtoms()( 2).GetPointer())
         + util::Format()( *aabackbone_atom.GetAtoms()( 3).GetPointer())
         +" but should have returned "
         + util::Format()( biol::GetAtomTypes().N)
         + util::Format()( biol::GetAtomTypes().CA)
         + util::Format()( biol::GetAtomTypes().C)
         + util::Format()( biol::GetAtomTypes().O)
            + " The function returned for copy_construct "
            + util::Format()( *copy_construct.GetAtoms()( 0).GetPointer())
            + util::Format()( *copy_construct.GetAtoms()( 1).GetPointer())
            + util::Format()( *copy_construct.GetAtoms()( 2).GetPointer())
            + util::Format()( *copy_construct.GetAtoms()( 3).GetPointer())
            +" but should have returned "
            + util::Format()( biol::GetAtomTypes().N)
            + util::Format()( biol::GetAtomTypes().CA)
            + util::Format()( biol::GetAtomTypes().C)
            + util::Format()( biol::GetAtomTypes().O)
      );

      // get the specified atom type from the object
      BCL_MessageStd( "Testing get GetAtom");
      BCL_Example_Check
      (
        aabackbone_atom.GetAtom( biol::GetAtomTypes().CA).GetType() == sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetType() &&
        aabackbone_atom.GetAtom( biol::GetAtomTypes().CA).GetCoordinates() == sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetCoordinates() &&
        aabackbone_atom.GetAtom( biol::GetAtomTypes().CB).GetType() == sp_first_residue->GetAtom( biol::GetAtomTypes().CB).GetType() &&
        aabackbone_atom.GetAtom( biol::GetAtomTypes().CB).GetCoordinates() == sp_first_residue->GetAtom( biol::GetAtomTypes().CB).GetCoordinates() &&
        aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetType() == sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetType() &&
        aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetCoordinates() == sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetCoordinates(),
        "Results from GetAtom are " + util::Format()( aabackbone_atom.GetAtom( biol::GetAtomTypes().CA)) +
        util::Format()( aabackbone_atom.GetAtom( biol::GetAtomTypes().CB)) + util::Format()( biol::GetAtomTypes().CB) + " but should be "
        + util::Format()( biol::GetAtomTypes().CA) + util::Format()( biol::GetAtomTypes().CB)
      );

      // test get the CA
      BCL_MessageStd( "Testing GetCA");
      BCL_Example_Check
      (
        aabackbone_atom.GetCA().GetType() == sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetType() &&
        aabackbone_atom.GetCA().GetCoordinates() == sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetCoordinates(),
        "The CA atom is " + util::Format()( aabackbone_atom.GetCA()) + " but should be "
        + util::Format()( sp_first_residue->GetAtom( biol::GetAtomTypes().CA))
      );

      // test GetFirstSideChainAtom (gets the CB)
      BCL_MessageStd( "Testing GetFirstSideChainAtom");
      BCL_Example_Check
      (
        aabackbone_atom.GetFirstSidechainAtom().GetType() == sp_first_residue->GetAtom( biol::GetAtomTypes().CB).GetType() &&
        aabackbone_atom.GetFirstSidechainAtom().GetCoordinates() == sp_first_residue->GetAtom( biol::GetAtomTypes().CB).GetCoordinates(),
        "First Sidechain Atom is " + util::Format()( aabackbone_atom.GetFirstSidechainAtom()) + " but should be "
        + util::Format()( sp_first_residue->GetAtom( biol::GetAtomTypes().CB))
      );

      // Explicitly sets an atom provided by the user
      BCL_MessageStd( "Testing SetAtom");
      aabackbone_atom.SetAtom( sp_first_residue->GetAtom( biol::GetAtomTypes().N));
      BCL_Example_Check
      (
        aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetType() == sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetType() &&
        aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetCoordinates() == sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetCoordinates(),
        "The set atom is" + util::Format()( aabackbone_atom.GetAtom( biol::GetAtomTypes().N)) + " but should be "
        + util::Format()( sp_first_residue->GetAtom( biol::GetAtomTypes().CA))
      );

      // test GetAAClass
      BCL_MessageStd( "Testing GetAAClass");
      BCL_Example_Check
      (
        aabackbone_atom.GetAAClass() == copy_construct.GetAAClass(),
        "The AAClass is " + util::Format()( aabackbone_atom.GetAAClass()) + " but should be "
        + util::Format()( copy_construct.GetAAClass())
      );

      // test CalculateOmega trivially with the copy constructor and the aabackbone_atom objects
      BCL_MessageStd( "Testing CalculateOmega");
      BCL_Example_Check
      (
        aabackbone_atom.CalculateOmega
        (
          sp_first_residue->GetAtom( biol::GetAtomTypes().CA),
          sp_first_residue->GetAtom( biol::GetAtomTypes().C)
        ) ==
          copy_construct.CalculateOmega
          (
            sp_first_residue->GetAtom( biol::GetAtomTypes().CA),
            sp_first_residue->GetAtom( biol::GetAtomTypes().C)
        ),
        "The Omega angle is " + util::Format()
        (
          aabackbone_atom.CalculateOmega
          (
            sp_first_residue->GetAtom( biol::GetAtomTypes().CA),
            sp_first_residue->GetAtom( biol::GetAtomTypes().C)
          )
        )
        + " but should be "
        + util::Format()
        (
          copy_construct.CalculateOmega
          (
            sp_first_residue->GetAtom( biol::GetAtomTypes().CA),
            sp_first_residue->GetAtom( biol::GetAtomTypes().C)
          )
        )
      );

      // test CalculatePhi, trivially with the copy constructor and aabackbone objects
      BCL_MessageStd( "Testing CalculatePhi");
      BCL_Example_Check
      (
        aabackbone_atom.CalculatePhi( sp_first_residue->GetAtom( biol::GetAtomTypes().C)) ==
        copy_construct.CalculatePhi( sp_first_residue->GetAtom( biol::GetAtomTypes().C)),
        "The Phi angle is " + util::Format()
        (
          aabackbone_atom.CalculatePhi( sp_first_residue->GetAtom( biol::GetAtomTypes().C))
        )
        + " but should be " + util::Format()
        (
          copy_construct.CalculatePhi( sp_first_residue->GetAtom( biol::GetAtomTypes().C))
        )
      );

      // test CalculatePsi
      BCL_MessageStd( "Testing CalculatePsi");
      BCL_Example_Check
      (
        aabackbone_atom.CalculatePsi( sp_first_residue->GetAtom( biol::GetAtomTypes().N)) ==
        copy_construct.CalculatePsi( sp_first_residue->GetAtom( biol::GetAtomTypes().N)),
        "The Phi angle is " + util::Format()
        (
          aabackbone_atom.CalculatePsi( sp_first_residue->GetAtom( biol::GetAtomTypes().N))
        )
        + " but should be " + util::Format()
        (
          copy_construct.CalculatePsi( sp_first_residue->GetAtom( biol::GetAtomTypes().N))
        )
      );

    ////////////////
    // operations //
    ////////////////

      // check if the SiPtr point to the correct linal::Vector3Ds after using GetAtomCoordinates
      BCL_Example_Check
      (
        *aabackbone_atom.GetAtomCoordinates()( 0).GetPointer() ==
          sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetCoordinates() &&
        *aabackbone_atom.GetAtomCoordinates()( 1).GetPointer() ==
          sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetCoordinates(),
        "GetAtomCoordinates returned" + util::Format()( *aabackbone_atom.GetAtomCoordinates()( 0).GetPointer()) + " as the N coordinates and " +
        util::Format()( *aabackbone_atom.GetAtomCoordinates()( 1).GetPointer()) + " as the CA coordinates, but should have returned" +
        util::Format()( sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetCoordinates())
        + " as the N coordinates and for the CA coordinates:  " +
        util::Format()( sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetCoordinates())
      );

      // test GetAtomCoordinates for specified atom types
      BCL_MessageStd( "Testing GetAtomCoordinates( AtomTypes)");
      BCL_Example_Check
      (
        *aabackbone_atom.GetAtomCoordinates( correct_atom_type)( 0).GetPointer() ==
          sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetCoordinates() &&
        *aabackbone_atom.GetAtomCoordinates( correct_atom_type)( 1).GetPointer() ==
          sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetCoordinates(),
        "The atom coordinates are" + util::Format()( *aabackbone_atom.GetAtomCoordinates( correct_atom_type)( 0).GetPointer())
        + " for the N coordinates, but should be " +
        util::Format()( sp_first_residue->GetAtom( biol::GetAtomTypes().N).GetCoordinates())
        + "The atom coordinates for the CA are " + util::Format()( *aabackbone_atom.GetAtomCoordinates( correct_atom_type)( 1).GetPointer())
        + " but should be " + util::Format()( sp_first_residue->GetAtom( biol::GetAtomTypes().CA).GetCoordinates())
      );

      // test translate
      BCL_MessageStd( "Testing Translate");
      aabackbone_atom.Translate( new_translate);
      // check if translate occurred according to translate function
      BCL_Example_Check
      (
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetCoordinates(), expected_coord_n) &&
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), expected_coord_ca),
        " Translation: the coordinates resulting from the translate should be " + util::Format()( expected_coord_n) + " "
        + util::Format()( expected_coord_ca)
        + "but are " + util::Format()( *aabackbone_atom.GetAtomCoordinates()( 0)) + util::Format()( *aabackbone_atom.GetAtomCoordinates()( 1))
      );

      // test Transform according to transformationMatrix3D
      BCL_MessageStd( "Testing Transform");
      aabackbone_atom.Transform( transformation_matrix);
      // check if transformation occurred according to transform function
      BCL_Example_Check
      (
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetCoordinates(), expected_transform_n) &&
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), expected_transform_ca),
        " Transformation: the coordinates resulting from the transformation should be " + util::Format()( expected_transform_n)
        + " for the CA coordinates" + util::Format()( expected_transform_ca) + " for the CB coordinates, but are "
        + util::Format()( *aabackbone_atom.GetAtomCoordinates()( 0)) + util::Format()( *aabackbone_atom.GetAtomCoordinates()( 1))
      );

      // test Rotation according to rotation matrix
      BCL_MessageStd( "Testing Rotate");
      aabackbone_atom.Rotate( rotation_matrix);
      // check if rotation occurred according to rotation matrix
      BCL_Example_Check
      (
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetCoordinates(), expected_rotation_n) &&
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), expected_rotation_ca),
        " Rotation: the coordinates resulting from the rotation should be " + util::Format()( expected_rotation_n) +
        " for the N coordinates" + util::Format()( expected_rotation_ca) + " for the CA coordinates, but are "
        + util::Format()( *aabackbone_atom.GetAtomCoordinates()( 0)) + util::Format()( *aabackbone_atom.GetAtomCoordinates()( 1))
      );

      // test GetCenter
      BCL_MessageStd( "Testing GetCenter");
      // check if GetCenter returned the proper center of mass
      BCL_Example_Check
      (
        math::EqualWithinTolerance( aabackbone_atom.GetCenter(), expected_center_of_mass),
        " Center of Mass: the coordinates resulting from the center of mass calculation should be "
        + util::Format()( expected_center_of_mass) + " but are " + util::Format()( aabackbone_atom.GetCenter())
      );

      // test SetToIdealConformation
      BCL_MessageStd( "Testing SetToIdealConformation")
      // invoke the function on the object
      aabackbone_atom.SetToIdealConformation( biol::GetSSTypes().HELIX, transformation_matrix);
      // check to see if the AA was properly idealized according to transformation_matrix
      BCL_Example_Check
      (
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().N).GetCoordinates(), ideal_coord_n) &&
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().CA).GetCoordinates(), ideal_coord_ca) &&
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().C).GetCoordinates(), ideal_coord_c) &&
        math::EqualWithinTolerance( aabackbone_atom.GetAtom( biol::GetAtomTypes().O).GetCoordinates(), ideal_coord_o) &&
        math::EqualWithinTolerance( aabackbone_atom.GetFirstSidechainAtom().GetCoordinates(), ideal_coord_fsa),
        " Idealize: the coordinates resulting from the idealizaton should be "
        + util::Format()( ideal_coord_n)
        + " for the N coordinates"
        + util::Format()( ideal_coord_ca)
        + " for the CA coordinates "
        + util::Format()( ideal_coord_c)
        + " for the C coordinates"
        + util::Format()( ideal_coord_o)
        + " for the O coordinates "
        + util::Format()( ideal_coord_fsa)
        + " for the FirstSideChainAtom coordinates"
        + " but are "
        + util::Format()( *aabackbone_atom.GetAtomCoordinates()(0))
        + util::Format()( *aabackbone_atom.GetAtomCoordinates()(1))
        + util::Format()( *aabackbone_atom.GetAtomCoordinates()(2))
        + util::Format()( *aabackbone_atom.GetAtomCoordinates()(3))
        + util::Format()( *aabackbone_atom.GetAtomCoordinates()(4))
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for biol::AACaCb");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( def_constr);
      BCL_MessageVrb( "read object");
      biol::AABackBone aabackbone_read;
      ReadBCLObject( aabackbone_read);

      // compare written and read object
      BCL_Example_Check
      (
        !def_constr.GetCA().GetCoordinates().IsDefined() &&
        !aabackbone_read.GetCA().GetCoordinates().IsDefined() &&
        !util::IsDefined( def_constr.GetAtoms()( 0)->GetPdbID()) &&
        !util::IsDefined( aabackbone_read.GetAtoms()( 0)->GetPdbID()),
        "read aabackbone_read is different from written: " + util::Format()( def_constr) + " != \n" +
        util::Format()( aabackbone_read)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAABackBone

  const ExampleClass::EnumType ExampleBiolAABackBone::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAABackBone())
  );

} // namespace bcl
