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
#include "biol/bcl_biol_aa_complete.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_back_bone.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_complete.cpp
  //!
  //! @author karakam
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAAComplete :
    public ExampleInterface
  {
  public:

    ExampleBiolAAComplete *Clone() const
    {
      return new ExampleBiolAAComplete( *this);
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

    static const linal::Vector3D CalculateCenterOfMassOfSideChain( const biol::AABase &AMINO_ACID)
    {
      // if the aa passed wasn't an AAComplete, then just return the center of the first side chain atom
      if( AMINO_ACID.GetAAClass() != biol::GetAAClasses().e_AAComplete)
      {
        // calculate the center of the first side chain atom
        return AMINO_ACID.GetFirstSidechainAtom().GetCoordinates();
      }

      // all the side chain atom types of that particular amino acid
      storage::Vector< biol::AtomType> atom_types_sidechain;

      // fill atomtypes with difference of all atom types for that amino acid and the backbone amino acid types
//      std::set_difference
//      (
//        AMINO_ACID.GetTypesOfAtoms().Begin(), AMINO_ACID.GetTypesOfAtoms().End(),
//        biol::GetAtomTypes().GetBackBoneAtomTypes().Begin(), biol::GetAtomTypes().GetBackBoneAtomTypes().End(),
//        atom_types_sidechain.Begin()
//      );

      // obtain all the coordinates for the side chain atoms
      util::SiPtrVector< const linal::Vector3D> side_chain_atom_coordinates
      (
        AMINO_ACID.GetAtomCoordinates( storage::Set< biol::AtomType>( atom_types_sidechain.Begin(), atom_types_sidechain.End()))
      );

      // calculate center of mass of the side chain based on all the side chain coordinates
      return coord::CenterOfMass( side_chain_atom_coordinates);
    }

    int Run() const
    {

    /////////////////
    // preparation //
    /////////////////

      // minimum sse sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 5;

      // initialize pdb filename to be read
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9_shifted.pdb"));
      // initialize read write streams
      io::IFStream read;
      io::OFStream write;

      // create factory and create the model
      pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
      assemble::ProteinModel native_model( factory.ProteinModelFromPDBFilename( pdb_filename, min_sse_sizes));

      // make a copy of the ShPtr to the sequence of chain A
      util::ShPtr< biol::AASequence> sp_sequence( native_model.GetChain( 'A')->GetSequence());

      // create the locator that finds helix 58-79
      const assemble::LocatorSSE helix_58_79_locator( 'A', 58, 79);

      // now locate the SSEs from the protein model
      const util::SiPtr< const assemble::SSE> sp_helix_58_79( helix_58_79_locator.Locate( native_model));

      // make sure the sse was correctly located
      BCL_ExampleIndirectAssert
      (
        sp_helix_58_79.IsDefined(),
        true,
        "Helix 58-79 was not located in 1IE9 pdb " + util::Format()( helix_58_79_locator)
      );

      // create the locator that finds aa 65
      const assemble::LocatorAA aa_65_locator( 'A', 65);

      // now locate the aa from the protein model
      const util::SiPtr< const biol::AABase> sp_aa_65( aa_65_locator.Locate( native_model));

      // create the locator that finds aa 64
      const assemble::LocatorAA aa_64_locator( 'A', 64);

      // now locate the aa from the protein model
      const util::SiPtr< const biol::AABase> sp_aa_64( aa_64_locator.Locate( native_model));

      // make sure the aa was correctly located
      BCL_ExampleIndirectAssert
      (
        sp_aa_65.IsDefined(),
        true,
        "Amino acid 65 was not located in 1IE9 pdb " + util::Format()( aa_65_locator)
      );

      // calculate center of mass of aa 65
      const linal::Vector3D com_aa_65( CalculateCenterOfMassOfSideChain( *sp_aa_65));

      // print out center
      BCL_MessageStd( "center of mass of side chain: " + util::Format()( com_aa_65));

      // create an AABase for use later
      util::ShPtr< biol::AABase> aa_base( new biol::AACaCb( biol::AA( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ARG)))));

      // default AABackBone object
      biol::AABackBone aa_backbone;

      // default CaCb object
      biol::AACaCb aa_ca_cb;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      biol::AAComplete default_aa_complete;

      // test default constructor
      const size_t number_of_atoms( 4);
      util::SiPtrVector< const biol::Atom>::const_iterator this_itr( default_aa_complete.GetAtoms().Begin());
      BCL_MessageStd( "Testing the default constructor");
      BCL_Example_Check
      (
        default_aa_complete.GetAAClass() == biol::GetAAClasses().e_AAComplete &&
        default_aa_complete.GetNumberOfAtoms() == number_of_atoms &&
        ( *this_itr)->GetType() == biol::GetAtomTypes().N &&
        ( *( ++this_itr))->GetType() == biol::GetAtomTypes().CA,
        "The default AAComplete constructor performed unexpectedly!\n" "The number of atoms should be " +
        util::Format()( number_of_atoms) + " but is instead " + util::Format()( default_aa_complete.GetNumberOfAtoms())
      );

      // constructor from AABase (ARG)
      biol::AAComplete aa_comp_from_aa_base( *aa_base);
      const size_t number_of_atoms_2( 5);
      util::SiPtrVector< const biol::Atom>::const_iterator aa_comp_itr( aa_comp_from_aa_base.GetAtoms().Begin());
      BCL_MessageStd( "Testing constructor from an AABase");
      const std::string desc( "Constructor from AABase");
      BCL_ExampleIndirectCheck( aa_comp_from_aa_base.GetAAClass(), biol::GetAAClasses().e_AAComplete, desc);
      BCL_ExampleIndirectCheck( aa_comp_from_aa_base.GetNumberOfAtoms(), 5, desc);
      BCL_ExampleIndirectCheck( ( *aa_comp_itr)->GetType(), biol::GetAtomTypes().N, desc);
      BCL_ExampleIndirectCheck( ( *( ++aa_comp_itr))->GetType(), biol::GetAtomTypes().CA, desc);
      BCL_ExampleIndirectCheck( aa_comp_from_aa_base.GetType(), biol::GetAATypes().ARG, desc);

      // simple ptr vector of trivial atoms
      const util::SiPtrVector< const biol::Atom> atoms
      (
        util::SiPtr< const biol::Atom>
        (
          new biol::Atom( biol::GetAtomTypes().CA)
        ),
        util::SiPtr< const biol::Atom>
        (
          new biol::Atom( biol::GetAtomTypes().N)
        )
      );

      // constructor from AABase and a SiPtrVector to define atoms
      biol::AAComplete aa_comp_from_base_atom_vector( *aa_base, atoms);

      // instantiate iterator
      util::SiPtrVector< const biol::Atom>::const_iterator atom_itr_a( aa_comp_from_base_atom_vector.GetAtoms().Begin());
      util::SiPtrVector< const biol::Atom>::const_iterator atom_itr_b( atom_itr_a + 1);

      // test the constructor for integrity
      BCL_MessageStd( "Testing constructor from an AABase and SiPtrVector to some atoms");
      BCL_Example_Check
      (
        aa_comp_from_base_atom_vector.GetAAClass() == biol::GetAAClasses().e_AAComplete &&
        aa_comp_from_base_atom_vector.GetNumberOfAtoms() == number_of_atoms_2 &&
        ( *atom_itr_a)->GetType() == biol::GetAtomTypes().N &&
        ( *atom_itr_b)->GetType() == biol::GetAtomTypes().CA &&
        aa_comp_from_base_atom_vector.GetType() == biol::GetAATypes().ARG,
        "The constructor from AABase and SiPtrVector to atoms performed unexpectedly!\nThe number of atoms in this object is " +
        util::Format()( aa_comp_from_base_atom_vector.GetNumberOfAtoms()) +
        " but should be " + util::Format()( number_of_atoms_2) +
        " and the atom types are " +
        util::Format()( ( *atom_itr_a)->GetType()) + " and " +
        util::Format()( ( *atom_itr_b)->GetType())
        + " but should be " +
        util::Format()( biol::GetAtomTypes().N) + " and " +
        util::Format()( biol::GetAtomTypes().CA)
      );

      // copy constr copying AABackbone to AAComplete
      biol::AAComplete aa_comp_backbone( aa_backbone);

      // copy constr copying AACaCb to AAComplete
      biol::AAComplete aa_comp_cacb( aa_ca_cb);

      // clone functionality
      util::ShPtr< biol::AAComplete> aa_comp_clone( default_aa_complete.Clone());

      // empty functionality
      util::ShPtr< biol::AAComplete> aa_comp_empty( default_aa_complete.Empty( util::ShPtr< biol::AAData>( new biol::AAData())));

    /////////////////
    // data access //
    /////////////////

      // testing get class identifier
       BCL_Example_Check
       (
         default_aa_complete.GetClassIdentifier() == "bcl::biol::AAComplete",
        "unexpected static class name: " + default_aa_complete.GetClassIdentifier() + " should be: bcl::biol::AAComplete"
       );

      // testing get number of atoms
      const size_t number_of_atoms_4( 2);
      BCL_MessageStd( "Testing GetNumberOfAtoms");
      BCL_Example_Check
      (
        aa_base->GetNumberOfAtoms() == number_of_atoms_4,
        "The number of atoms for this AAComplete as ARG is incorrectly reported as " + util::Format()( aa_base->GetNumberOfAtoms())
        + " but should be " + util::Format()( number_of_atoms_4)
      );

      // storage set for atom types
      const storage::Set< biol::AtomType> atom_set( biol::GetAATypes().ARG->GetAllowedAtomTypes());

      // testing get types of atoms
      BCL_MessageStd( "Testing GetTypesOfAtoms");
      BCL_Example_Check
      (
        aa_comp_from_aa_base.GetTypesOfAtoms().InternalData() == atom_set.InternalData(),
        "This AAComplete object contains atoms of " + util::Format()( aa_comp_from_aa_base.GetTypesOfAtoms())
        + " but should be of type ARG and should have atoms of type: " + util::Format()( atom_set)
      );

      // testing GetAtom, SetAtom && GetCA
      const util::ShPtr< biol::Atom> atom_ca( new biol::Atom( linal::Vector3D( 1.0, 2.0, 3.0), biol::GetAtomTypes().CA));
      aa_comp_from_aa_base.SetAtom( *atom_ca);
      BCL_MessageStd( "Testing SetAtom and GetCA");
      BCL_Example_Check
      (
        aa_comp_from_aa_base.GetAtom( biol::GetAtomTypes().CA).GetType() == atom_ca->GetType() &&
        aa_comp_from_aa_base.GetAtom( biol::GetAtomTypes().CA).GetCoordinates() == atom_ca->GetCoordinates() &&
        aa_comp_from_aa_base.GetCA().GetType() == atom_ca->GetType() &&
        aa_comp_from_aa_base.GetCA().GetCoordinates() == atom_ca->GetCoordinates(),
        "This AAComplete's CA atom is " + util::Format()( aa_comp_from_aa_base.GetAtom( biol::GetAtomTypes().CA))
        + " but should be " + util::Format()( *atom_ca)
      );

      // testing GetFirstSideChainAtom( CB)
      const util::ShPtr< biol::Atom> atom_cb( new biol::Atom( linal::Vector3D( 2.0, 2.7, 2.98), biol::GetAtomTypes().CB));
      aa_comp_from_aa_base.SetAtom( *atom_cb);
      BCL_MessageStd( "Testing GetAtom and SetAtom");
      BCL_Example_Check
      (
        aa_comp_from_aa_base.GetFirstSidechainAtom().GetType() == atom_cb->GetType() &&
        aa_comp_from_aa_base.GetFirstSidechainAtom().GetCoordinates() == atom_cb->GetCoordinates(),
        "This AAComplete's CB atom is " + util::Format()( aa_comp_from_aa_base.GetFirstSidechainAtom())
        + " but should be " + util::Format()( *atom_cb)
      );

      // testing SetAtoms
      BCL_MessageStd( "Testing SetAtoms");
      const util::SiPtr< const biol::Atom> si_ptr_ca( atom_ca);
      const util::SiPtr< const biol::Atom> si_ptr_cb( atom_cb);
      const util::SiPtrVector< const biol::Atom> atoms_to_be_set( si_ptr_ca, si_ptr_cb);
      aa_comp_from_aa_base.SetAtoms( atoms_to_be_set);
      BCL_Example_Check
      (
        aa_comp_from_aa_base.GetFirstSidechainAtom().GetType() == atom_cb->GetType() &&
        aa_comp_from_aa_base.GetFirstSidechainAtom().GetCoordinates() == atom_cb->GetCoordinates() &&
        aa_comp_from_aa_base.GetCA().GetType() == atom_ca->GetType() &&
        aa_comp_from_aa_base.GetCA().GetCoordinates() == atom_ca->GetCoordinates(),
        "This AAComplete's CB atom is " + util::Format()( aa_comp_from_aa_base.GetFirstSidechainAtom())
        + " but should be " + util::Format()( *atom_cb)
      );

      // testing GetAAClass
      BCL_MessageStd( "Testing GetAAClass");
      BCL_Example_Check
      (
        aa_comp_from_aa_base.GetAAClass() == biol::GetAAClasses().e_AAComplete,
        "This object's type is " + util::Format()( aa_comp_from_aa_base.GetAAClass())
        + " but should be " + util::Format()( biol::GetAAClasses().e_AAComplete)
      );

      // testing CalculateOmega
      const double calculated_omega( sp_aa_65->CalculateOmega( sp_aa_64->GetCA(), sp_aa_65->GetAtom( biol::GetAtomTypes().C)));
      const double expected_omega( -2.66141);
      BCL_MessageStd( "Testing CalculateOmega");
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_omega, calculated_omega),
        "This calculated omega angle is " + util::Format()( calculated_omega)
        + " but is expected to be " + util::Format()( expected_omega)
      );

      // testing CalculatePhi
      const double calculated_phi( sp_aa_65->CalculatePhi( sp_aa_64->GetAtom( biol::GetAtomTypes().C)));
      const double expected_phi( -1.06812);
      BCL_MessageStd( "Testing CalculatePhi");
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_phi, calculated_phi),
        "This calculated phi angle is " + util::Format()( calculated_phi)
        + " but is expected to be " + util::Format()( expected_phi)
      );

      // testing CalculatePsi
      const double calculated_psi( sp_aa_64->CalculatePsi( sp_aa_65->GetAtom( biol::GetAtomTypes().N)));
      const double expected_psi( -0.738314);
      BCL_MessageStd( "Testing CalculatePsi");
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_psi, calculated_psi),
        "This calculated psi angle is " + util::Format()( calculated_psi)
        + " but is expected to be " + util::Format()( expected_psi)
      );

    ///////////////
    // operators //
    ///////////////

      // test the = assignment operator
      BCL_MessageStd( "Testing the = assignment operator")
      biol::AAComplete assignment;
      assignment = *sp_aa_64;
      const linal::Vector3D sp_coords( sp_aa_64->GetCA().GetCoordinates());
      const linal::Vector3D assign_coords( assignment.GetCA().GetCoordinates());
      BCL_Example_Check
      (
        sp_coords( 0) == assign_coords( 0) &&
        sp_coords( 1) == assign_coords( 1) &&
        sp_coords( 2) == assign_coords( 2),
        "The assignment operator: the coordinates resulting from this operator are: " + util::Format()( sp_coords)
        + " for the CA coordinates but were expected to be: " + util::Format()( assign_coords)
      );

    ////////////////
    // operations //
    ////////////////

      // testing GetAtomCoordinates
      // set of atoms composing aa_64
      const storage::Set< biol::AtomType> aa_64_atom_set( sp_aa_64->GetType()->GetAllowedAtomTypes());
      util::SiPtrVector< const linal::Vector3D> atom_coordinates( sp_aa_64->GetAtomCoordinates());
      // GetAtomCoordinates for only the atoms in the set
      util::SiPtrVector< const linal::Vector3D> atom_coordinates_from_set( sp_aa_64->GetAtomCoordinates( aa_64_atom_set));
      BCL_MessageStd( "Testing GetAtomCoordinates");
      BCL_Example_Check
      (
        atom_coordinates == atom_coordinates_from_set,
        " The atom coordinates for this AA are " + util::Format()( atom_coordinates) + " but should be "
        + util::Format()( atom_coordinates_from_set)
      );

      // testing Translate
      // translation vector
      const linal::Vector3D translation_vector( 1.755, 3.259, -6.503);
      util::ShPtr< biol::AABase> new_aa_64( new biol::AAComplete( *sp_aa_64));
      linal::Vector3D old_ca_coordinates( new_aa_64->GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
      linal::Vector3D old_n_coordinates( new_aa_64->GetAtom( biol::GetAtomTypes().N).GetCoordinates());
      new_aa_64->Translate( translation_vector);
      const linal::Vector3D new_ca_coordinates( new_aa_64->GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D new_n_coordinates( new_aa_64->GetAtom( biol::GetAtomTypes().N).GetCoordinates());
      BCL_MessageStd( "Testing Translate");
      BCL_Example_Check
      (
        math::EqualWithinTolerance( old_ca_coordinates( 0) + 1.755, new_ca_coordinates( 0)) &&
        math::EqualWithinTolerance( old_ca_coordinates( 1) + 3.259, new_ca_coordinates( 1)) &&
        math::EqualWithinTolerance( old_ca_coordinates( 2) - 6.503, new_ca_coordinates( 2)) &&
        math::EqualWithinTolerance( old_n_coordinates( 0) + 1.755 , new_n_coordinates( 0)) &&
        math::EqualWithinTolerance( old_n_coordinates( 1) + 3.259 , new_n_coordinates( 1)) &&
        math::EqualWithinTolerance( old_n_coordinates( 2) - 6.503 , new_n_coordinates( 2)),
        "Translation function is broken."
      );

      // test Transform
      BCL_MessageStd( "Testing Transform");
      new_aa_64->Transform( math::TransformationMatrix3D( translation_vector));
      const linal::Vector3D transf_ca_coords( new_aa_64->GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D transf_n_coords( new_aa_64->GetAtom( biol::GetAtomTypes().N).GetCoordinates());
      const linal::Vector3D expected_ca_transf_coords( 8.388, 23.506, 16.886);
      const linal::Vector3D expected_n_transf_coords( 8.178, 24.946, 16.847);

      // check if transformation occurred according to transform function
      BCL_Example_Check
      (
        transf_ca_coords == expected_ca_transf_coords &&
        transf_n_coords == expected_n_transf_coords,
        " Transformation: the coordinates resulting from the transformation are: " + util::Format()( transf_ca_coords)
        + " for the CA coordinates but should be " + util::Format()( expected_ca_transf_coords) +
        " and for the N coordiantes " + util::Format()( transf_n_coords) + " but should be "
        + util::Format()( expected_n_transf_coords)
      );

      // test Rotate
      BCL_MessageStd( "Testing Rotate");
      const math::RotationMatrix3D rotation_matrix( translation_vector, 109.5);
      new_aa_64->Rotate( rotation_matrix);
      const linal::Vector3D rot_ca_coords( new_aa_64->GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D rot_n_coords( new_aa_64->GetAtom( biol::GetAtomTypes().N).GetCoordinates());
      const linal::Vector3D expected_ca_rot_coords( -20.8596, -18.1971, -11.9068);
      const linal::Vector3D expected_n_rot_coords( -20.942, -19.0687, -13.0699);

      // check if rotation occurred according to rotation function
      BCL_Example_Check
      (
        math::EqualWithinTolerance( rot_ca_coords, expected_ca_rot_coords) &&
        math::EqualWithinTolerance( rot_n_coords, expected_n_rot_coords),
        " Rotation: the coordinates resulting from the rotation are: " + util::Format()( rot_ca_coords)
        + " for the CA coordinates but should be " + util::Format()( rot_ca_coords) +
        " and for the N coordiantes " + util::Format()( rot_n_coords) + " but should be "
        + util::Format()( expected_n_rot_coords)
      );

      // check if GetCenter occurred according to GetCenter function
      BCL_MessageStd( "Testing GetCenter")
      const linal::Vector3D aa_center( new_aa_64->GetCenter());
      const linal::Vector3D expected_center( -20.1285, -17.7221, -11.9187);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_center, aa_center),
        " The GetCenter function found the center of this AA to be\n" + util::Format()( aa_center)
        + " but was expected to be\n " + util::Format()( expected_center)
      );

      // check if SetToIdealConformation occurred as expected
      BCL_MessageStd( "Testing SetToIdealConformation")
      new_aa_64->SetToIdealConformation( biol::GetSSTypes().HELIX, math::TransformationMatrix3D( translation_vector));
      const linal::Vector3D idl_ca_coords( new_aa_64->GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D idl_n_coords( new_aa_64->GetAtom( biol::GetAtomTypes().N).GetCoordinates());
      const linal::Vector3D expected_ca_idl_coords( 4.015, 3.259, -5.842);
      const linal::Vector3D expected_n_idl_coords( 3.10203, 2.55479, -6.748);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( idl_ca_coords, expected_ca_idl_coords) &&
        math::EqualWithinTolerance( idl_n_coords, expected_n_idl_coords),
        "SetToIdealConformation: the coordinates resulting from the function are: " + util::Format()( idl_ca_coords)
        + " for the CA coordinates but should be " + util::Format()( expected_ca_idl_coords) +
        " and for the N coordiantes " + util::Format()( idl_n_coords) + " but should be "
        + util::Format()( expected_n_idl_coords)
      );

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for biol::AAComplete");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( default_aa_complete);
      BCL_MessageVrb( "read object");
      biol::AAComplete aacomp_read;
      ReadBCLObject( aacomp_read);

      // compare written and read object
      BCL_Example_Check
      (
        !default_aa_complete.GetCA().GetCoordinates().IsDefined() &&
        !aacomp_read.GetCA().GetCoordinates().IsDefined() &&
        !util::IsDefined( default_aa_complete.GetAtoms()( 0)->GetPdbID()) &&
        !util::IsDefined( aacomp_read.GetAtoms()( 0)->GetPdbID()),
        "read AAComplete is different from written: " + util::Format()( default_aa_complete) + " != \n" +
        util::Format()( aacomp_read)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // testing DefaultBackboneAtoms
      // DefaultBackBoneAtoms should be invoked when the constructor for this class is instantiated
      BCL_MessageStd( "Testing DefaultBackBoneAtoms");
      storage::Set< biol::AtomType> atom_types_of_aa;
      util::SiPtrVector< const biol::Atom> default_atoms( default_aa_complete.GetAtoms());
      // insert atom types of default atoms into atom_types_of_aa
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator
         itr( default_atoms.Begin()), itr_end( default_atoms.End());
        itr != itr_end;
        ++itr
      )
      {
        atom_types_of_aa.Insert( ( *itr)->GetType());
      }

      BCL_Example_Check
      (
        biol::GetAtomTypes().GetBackBoneAtomTypes().InternalData() == atom_types_of_aa.InternalData(),
        "The default backbone atom types should be \n" + util::Format()( biol::GetAtomTypes().GetBackBoneAtomTypes()) +
        " but are " + util::Format()( atom_types_of_aa)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAAComplete

  const ExampleClass::EnumType ExampleBiolAAComplete::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAAComplete())
  );

} // namespace bcl
