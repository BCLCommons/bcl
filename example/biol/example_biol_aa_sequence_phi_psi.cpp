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
#include "biol/bcl_biol_aa_sequence_phi_psi.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_biol_aa_sequence_phi_psi.cpp
  //!
  //! @author weinerbe
  //! @date Jan 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleBiolAASequencePhiPsi :
    public ExampleInterface
  {
  public:

    ExampleBiolAASequencePhiPsi *Clone() const
    {
      return new ExampleBiolAASequencePhiPsi( *this);
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

      // create protein model from pdb
      BCL_MessageStd( "building model");
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 9;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 5;
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, ssetype_min_size)
      );

      // get the helix
      const util::SiPtr< const assemble::SSE> sp_helix
      (
        protein_model.GetSSEs( biol::GetSSTypes().HELIX).FirstElement()
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      biol::AASequencePhiPsi def_construct;
      BCL_ExampleIndirectCheck( def_construct.GetAngles().IsEmpty(), true, "Default constructor");

      // test construct from AASequence
      biol::AASequencePhiPsi seq_construct( *sp_helix);

      // test clone
      util::ShPtr< biol::AASequencePhiPsi> clone_construct( seq_construct.Clone());
      BCL_ExampleIndirectCheck( clone_construct->GetCA(), seq_construct.GetCA(), "Clone");

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier
      BCL_ExampleCheck( seq_construct.GetClassIdentifier(), "bcl::biol::AASequencePhiPsi");

      // test GetN
      BCL_ExampleCheck
      (
        seq_construct.GetN(),
        sp_helix->GetData()( sp_helix->GetSize() / 2)->GetAtom( biol::GetAtomTypes().N).GetCoordinates()
      );

      // test GetCa
      BCL_ExampleCheck
      (
        seq_construct.GetCA(),
        sp_helix->GetData()( sp_helix->GetSize() / 2)->GetAtom( biol::GetAtomTypes().CA).GetCoordinates()
      );

      // test GetC
      BCL_ExampleCheck
      (
        seq_construct.GetC(),
        sp_helix->GetData()( sp_helix->GetSize() / 2)->GetAtom( biol::GetAtomTypes().C).GetCoordinates()
      );

      // test GetAngles
      BCL_ExampleIndirectCheck
      (
        seq_construct.GetAngles().GetSize(), 12, "GetAngles"
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

      // write the object
      WriteBCLObject( seq_construct);

      // read the object back in
      biol::AASequencePhiPsi phi_psi_read;
      ReadBCLObject( phi_psi_read);

      // compare them
      BCL_ExampleIndirectCheck
      ( seq_construct.GetN(), phi_psi_read.GetN(), "read and write"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleBiolAASequencePhiPsi

  const ExampleClass::EnumType ExampleBiolAASequencePhiPsi::s_Instance
  (
    GetExamples().AddEnum( ExampleBiolAASequencePhiPsi())
  );

} // namespace bcl
