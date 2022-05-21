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
#include "assemble/bcl_assemble_voxel_grid_aa.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_aa_ca_cb.h"
#include "biol/bcl_biol_aa_complete.h"
#include "biol/bcl_biol_atom.h"
#include "score/bcl_score_protein_model_score_sum.h"
#include "score/bcl_score_protein_model_sse_pairs.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_assemble_voxel_grid_aa.cpp
  //!
  //! @author mendenjl
  //! @date Nov 22, 2016
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleVoxelGridAA :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleAssembleVoxelGridAA
    ExampleAssembleVoxelGridAA *Clone() const
    {
      return new ExampleAssembleVoxelGridAA( *this);
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

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

      BCL_MessageStd( " Starting tests for SlicelistAA. ");

      TestExtractPosition();
      TestGetRelevantNeighbors();
      TestInsertion();
      TestDeletion();
      TestUpdate();

      BCL_MessageStd( " Finished tests for SlicelistAA. ");

      return 0;
    }

    void TestGetRelevantNeighbors() const
    {
      BCL_MessageStd( " Starting relevant neighbor tests. ");

      storage::Vector< biol::AACaCb> aa_vec = GetVectorOfAACaCbs();
      assemble::VoxelGridAA slmaa;
      slmaa.SetObjects( ConvertToSiPtrVector( aa_vec));

      BCL_MessageStd( " Set Objects complete");
      BCL_ExampleCheck( slmaa.GetNeighbors( aa_vec( 0), 3.9).GetSize(), 2);

      BCL_MessageStd( " Finished relevant neighbor tests. ");
    }

    void TestInsertion() const
    {
      BCL_MessageStd( " Starting insertion-tests. ");

      storage::Vector< biol::AACaCb> aa_vec = GetVectorOfAACaCbs();
      assemble::VoxelGridAA slmaa;
      slmaa.SetObjects( ConvertToSiPtrVector( aa_vec));

      const linal::Vector3D CA_coordinates( 0.1, 0.1, 0.1);
      const linal::Vector3D CB_coordinates( 7.5, 7.2, 7.3);
      const biol::Atom def_CA( CA_coordinates, biol::AtomType( "CA"), 7, 0.242);
      const biol::Atom def_CB( CB_coordinates, biol::AtomType( "CB"), 8, 0.532);
      const biol::AA aa_obj( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 11, 11, ' ')));
      biol::AACaCb aa1( aa_obj, def_CA, def_CB);

      slmaa.InsertObject( aa1);

      TestSizes( slmaa, 6);

      storage::Vector< biol::AACaCb> aa_vec2 = GetVectorOfAACaCbs();
      assemble::VoxelGridAA slmaa2;
      slmaa2.SetObjects( ConvertToSiPtrVector( aa_vec2));

      const biol::AA aa_obj2( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 12, 12, ' ')));
      biol::AACaCb aa2( aa_obj2, def_CA, def_CB);

      util::SiPtrVector< biol::AACaCb> sivec_aa;
      sivec_aa.PushBack( aa1);
      sivec_aa.PushBack( aa2);

      slmaa2.InsertObjects( sivec_aa);

      TestSizes( slmaa2, 7);

      BCL_MessageStd( " Finished insertion-tests. ");
    }

    void TestDeletion() const
    {
      BCL_MessageStd( " Starting deletion-tests. ");

      storage::Vector< biol::AACaCb> aa_vec = GetVectorOfAACaCbs();
      assemble::VoxelGridAA slmaa;
      slmaa.SetObjects( ConvertToSiPtrVector( aa_vec));

      slmaa.RemoveObject( aa_vec( 0));

      TestSizes( slmaa, 4);

      slmaa.RemoveObject( aa_vec( 3));

      TestSizes( slmaa, 3);

      storage::Vector< biol::AACaCb> aa_vec2 = GetVectorOfAACaCbs();
      assemble::VoxelGridAA slmaa2;
      slmaa2.SetObjects( ConvertToSiPtrVector( aa_vec2));
      aa_vec2.RemoveElements(0, 1);
      aa_vec2.RemoveElements(1, 1);

      slmaa2.RemoveObjects( util::ConvertToSiPtrVector( aa_vec2));

      TestSizes( slmaa2, 2);

      BCL_MessageStd( " Finished deletion-tests. ");
    }

    void TestUpdate() const
    {
      BCL_MessageStd( " Starting update-tests. ");

      storage::Vector< biol::AACaCb> aa_vec = GetVectorOfAACaCbs();
      assemble::VoxelGridAA slmaa;
      slmaa.SetObjects( ConvertToSiPtrVector( aa_vec));

      aa_vec.PopBack();
      aa_vec.PopBack();

      slmaa.SetObjects( ConvertToSiPtrVector( aa_vec));
      TestSizes( slmaa, 3);

      BCL_MessageStd( " Finished update-tests. ");
    }

    void TestExtractPosition() const
    {
      BCL_MessageStd( " Starting extract-tests. ");

      // Get protein's AAs
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1C1D.pdb"));
      storage::Map< biol::SSType, size_t> ssetype_min_size;
      ssetype_min_size[ biol::GetSSTypes().HELIX] = 4;
      ssetype_min_size[ biol::GetSSTypes().STRAND] = 4;
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AACaCb, ssetype_min_size));
      const util::SiPtrVector< const biol::AACaCb> aas( model.GetAminoAcids());

      assemble::VoxelGridAA slamAACaCb;
      storage::Vector< biol::AACaCb> vecAA = GetVectorOfAACaCbs();

      util::SiPtr< const linal::Vector3D> res = slamAACaCb.ExtractPosition( vecAA( 0));
      BCL_ExampleIndirectCheck( *res, linal::Vector3D( 1.5, 2.2, 6.3), "ALA");

      // only test for the first 10
      size_t counter( 0);
      for( size_t i = 0; i < aas.GetSize(); ++i) {
        ++counter;
        const util::SiPtr< const biol::AACaCb> aab( static_cast< const util::SiPtr< const biol::AACaCb> >( aas( i)));
        const util::SiPtr< const linal::Vector3D> point( aab->GetFirstSidechainAtom().GetCoordinates());
        res = slamAACaCb.ExtractPosition( aab);

        if(res.IsDefined()) {
          BCL_Example_Check
          (
            *res == *point,
            " An incorrect position was extracted."
          );
        }
        if( counter > 8) break;
      }

      BCL_MessageStd( " Finished extract-tests. ");
    }

    void TestSizes( const util::SiPtr< assemble::VoxelGridAA> MANAGER, const size_t EXPECTED_SIZE) const
    {
      BCL_ExampleIndirectCheck( MANAGER->GetNumberItems(), EXPECTED_SIZE, "Deletion");
    }

    storage::Vector< biol::AACaCb> GetVectorOfAACaCbs() const
    {
      const linal::Vector3D CA_coordinates( 0.5, 0.2, 0.3);
      const linal::Vector3D CB_coordinates( 1.5, 2.2, 6.3);
      const biol::Atom def_CA( CA_coordinates, biol::AtomType( "CA"), 7, 0.242);
      const biol::Atom def_CB( CB_coordinates, biol::AtomType( "CB"), 8, 0.532);
      const biol::AA aa_obj( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 1, 3, ' ')));

      const linal::Vector3D CA_coordinates2( 0.8, 1.0, 0.5);
      const linal::Vector3D CB_coordinates2( 1.0, 2.0, 5.5);
      const biol::Atom def_CA2( CA_coordinates2, biol::AtomType( "CA"), 7, 0.242);
      const biol::Atom def_CB2( CB_coordinates2, biol::AtomType( "CB"), 8, 0.532);
      const biol::AA aa_obj2( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 2, 3, ' ')));

      const linal::Vector3D CA_coordinates3( 0.8, 1.0, 0.5);
      const linal::Vector3D CB_coordinates3( 3.0, 1.0, 4.5);
      const biol::Atom def_CA3( CA_coordinates3, biol::AtomType( "CA"), 7, 0.242);
      const biol::Atom def_CB3( CB_coordinates3, biol::AtomType( "CB"), 8, 0.532);
      const biol::AA aa_obj3( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ALA, 3, 3, ' ')));

      const linal::Vector3D CA_coordinates4( 0.8, 1.0, 0.5);
      const linal::Vector3D CB_coordinates4( 7.0, 0.5, 2.0);
      const biol::Atom def_CA4( CA_coordinates4, biol::AtomType( "CA"), 7, 0.242);
      const biol::Atom def_CB4( CB_coordinates4, biol::AtomType( "CB"), 8, 0.532);
      const biol::AA aa_obj4( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ILE, 4, 3, ' ')));

      const linal::Vector3D CA_coordinates5( 0.8, 1.0, 0.5);
      const linal::Vector3D CB_coordinates5( 2.0, 4.0, 1.0);
      const biol::Atom def_CA5( CA_coordinates5, biol::AtomType( "CA"), 7, 0.242);
      const biol::Atom def_CB5( CB_coordinates5, biol::AtomType( "HA2"), 8, 0.532);
      const biol::AA aa_obj5( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().GLY, 5, 3, ' ')));

      biol::AACaCb aa1( aa_obj, def_CA, def_CB);
      biol::AACaCb aa2( aa_obj2, def_CA2, def_CB2);
      biol::AACaCb aa3( aa_obj3, def_CA3, def_CB3);
      biol::AACaCb aa4( aa_obj4, def_CA4, def_CB4);
      biol::AACaCb aa5( aa_obj5, def_CA5, def_CB5);

      storage::Vector< biol::AACaCb> v_aacacb;
      v_aacacb.PushBack( aa1);
      v_aacacb.PushBack( aa2);
      v_aacacb.PushBack( aa3);
      v_aacacb.PushBack( aa4);
      v_aacacb.PushBack( aa5);

      return v_aacacb;
    }

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleOctree

  //! single instance of this class
  const ExampleClass::EnumType ExampleAssembleVoxelGridAA::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleVoxelGridAA())
  );

} // namespace bcl
