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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "fold/bcl_fold_protein_geometry.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_matrix3x3.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    const util::SiPtr< const util::ObjectInterface> ProteinGeometry::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinGeometry())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinGeometry::ProteinGeometry()
    {
    }

    //! @brief clone function
    //! @return pointer to a new ProteinGeometry
    ProteinGeometry *ProteinGeometry::Clone() const
    {
      return new ProteinGeometry( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &ProteinGeometry::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns the local coordinate system of the given amino acid
    //! @param AMINO_ACID amino aid to return the local coordinate system for
    //! @return local coordinate system of the given amino acid
    linal::Matrix3x3< double> ProteinGeometry::GetLocalCoordinateSystem( const biol::AABase &AMINO_ACID)
    {
      // get the coordinates to compute the local coordinate system
      const linal::Vector3D &ca_coord( AMINO_ACID.GetAtom( biol::GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D &cb_coord( AMINO_ACID.GetAtom( biol::GetAtomTypes().CB).GetCoordinates());
      const linal::Vector3D &n_coord( AMINO_ACID.GetAtom( biol::GetAtomTypes().N).GetCoordinates());

      // compute the orthonormal base of the local coordinate system
      const linal::Vector3D y_axis( ( n_coord - ca_coord).Normalize());
      const linal::Vector3D ca_cb( cb_coord - ca_coord);
      const linal::Vector3D ca_cb_y( ( ca_cb * y_axis) * y_axis);
      const linal::Vector3D x_axis( ( ca_cb - ca_cb_y).Normalize());
      const linal::Vector3D z_axis( linal::CrossProduct( x_axis, y_axis).Normalize());

      // create the matrix defining the coordinate system
      linal::Matrix3x3< double> base( util::GetUndefinedDouble());
      for( size_t i( 0); i < 3; ++i)
      {
        base( 0, i) = x_axis( i);
      }
      for( size_t i( 0); i < 3; ++i)
      {
        base( 1, i) = y_axis( i);
      }
      for( size_t i( 0); i < 3; ++i)
      {
        base( 2, i) = z_axis( i);
      }

      return base;
    }

    //! @brief fits the given loop to the given template and returns a new protein model containing the loop
    //! @param PROTEIN_MODEL protein model to containing the loop
    //! @param LOOP loop to be the given template
    //! @param TEMPLATE template used to fit the loop to
    //! @return new protein model containing the fitted loop
    util::ShPtr< assemble::ProteinModel> ProteinGeometry::FitToTemplate
    (
      const assemble::ProteinModel &PROTEIN_MODEL, const LoopParameters &LOOP, const LoopParameters &TEMPLATE
    )
    {
      // get the amino acids which need to be fitted
      const assemble::Chain &chain( *PROTEIN_MODEL.GetChain( LOOP.GetChainID()));
      const util::SiPtrVector< const biol::AABase> res_tmp( chain.GetAminoAcids());
      const util::SiPtrVector< const biol::AABase> residues
      (
        res_tmp[ LOOP.GetAnchors()( 0)], // first residue in the loop to be fitted
        res_tmp[ LOOP.GetAnchors()( 1)]  // last residue in the loop to be fitted
      );

      // get an initial conformation of the loop
      const biol::AASequence &sequence( *chain.GetSequence());
      const size_t length( LOOP.GetAnchors()( 1) - LOOP.GetAnchors()( 0) - 1);
      util::ShPtr< assemble::SSE> sp_loop
      (
        new assemble::SSE( sequence.SubSequence( LOOP.GetAnchors()( 0), length), biol::GetSSTypes().COIL)
      );
      sp_loop->SetType( biol::GetSSTypes().STRAND);
      sp_loop->SetToIdealConformationInPlace();
      sp_loop->SetType( biol::GetSSTypes().COIL);

      // fit the sequence to the template
      const storage::Vector< double> &angles( TEMPLATE.GetAngles());
      int seq_id( LOOP.GetAnchors()( 0) + 1);
      for( auto it( angles.Begin() + 1), it_end( angles.End() - 1); it != it_end; it += 2, ++seq_id)
      {
        storage::VectorND< 2, double> phi_psi( *it, *( it + 1));
        biol::AASequenceFlexibility::SetPhiPsi
        (
          *sp_loop, seq_id, phi_psi, biol::AASequenceFlexibility::SequenceDirection::e_CTerminal
        );
      }

      // attach the sequence to the anchor residue
      const int n_term_id( sp_loop->GetFirstMember()->GetSeqID());
      const biol::AABase &n_term_loop( *sp_loop->GetFirstMember());
      util::ShPtrVector< biol::AABase> n_term_res( 1, *chain.GetSequence()->FindAABySeqID( n_term_id - 1));
      biol::AASequence n_term_seq( n_term_res, n_term_res( 0)->GetChainID());
      const math::TransformationMatrix3D transform
      (
        biol::AASequenceFactory::TransformationAppend( n_term_seq, n_term_loop, angles( 1))
      );
      sp_loop->Transform( transform);

      util::ShPtr< assemble::ProteinModel> sp_new_model( PROTEIN_MODEL.HardCopy());
      for( auto res_it( sp_loop->Begin()), res_it_end( sp_loop->End()); res_it != res_it_end; ++res_it)
      {
        const biol::AABase &new_residue( **res_it);
        const int seq_id( new_residue.GetSeqID());
        biol::AABase &residue( **sp_new_model->GetChain( sp_loop->GetChainID())->GetSequence()->FindAABySeqID( seq_id));
        const util::SiPtrVector< const biol::Atom> new_atoms( new_residue.GetAtoms());
        residue.SetAtoms( new_atoms);
      }
      sp_new_model->ConnectSSEToChainData();

      return sp_new_model;
    }

    //! @brief combines the two given sequence at the given merging point
    //! @detail the given sequences are merged at the given merging points, which denote the residue index.
    //! merging points 5 and 3 results in the fifth residue of the n-terminal sequence being connected top the
    //! third residue of the c-terminal sequence
    //! @param SEQ_N n-terminal sequence to be merged
    //! @param SEQ_C c-terminal sequence to be merged
    //! @param MERGE_N merging point of the n-terminal sequence
    //! @param MERGE_C merging point of the c-terminal sequence
    util::ShPtr< biol::AASequence> ProteinGeometry::CombineSequences
    (
      const biol::AASequence &SEQ_N,
      const biol::AASequence &SEQ_C,
      int MERGE_N,
      int MERGE_C
     )
    {
      // truncate the sequences to be consistent with the merging points
      const biol::AASequence seq_n( SEQ_N.SubSequence( 0, MERGE_N + 1));
      biol::AASequence seq_c( SEQ_C.SubSequence( MERGE_C, SEQ_C.GetSize() - MERGE_C));

      // make chain and sequence ids of c-terminal sequence consistent with n-terminal sequence
      int n_new_id( ( **seq_n.GetMembers().Last()).GetSeqID());
      const char n_chain_id( SEQ_N.GetChainID());
      for( auto aa_it( seq_c.Begin()); aa_it != seq_c.End(); ++aa_it)
      {
        util::ShPtr< biol::AAData> sp_res_data( new biol::AAData( *( **aa_it).GetData()));
        sp_res_data->SetSeqID( ++n_new_id);
        sp_res_data->SetChainID( n_chain_id);
        ( **aa_it).SetData( sp_res_data);
      }

      // determine the phi- and psi-angles of the merging points in the n- and c-terminal sequences
      double phi( 0.0);
      if( MERGE_C != 0)
      {
        // compute the phi-angle for the residue located at the merging point
        const biol::AABase &merge_res_prev( *SEQ_C.GetMembers()( MERGE_C - 1));
        const biol::AABase &merge_res( *SEQ_C.GetMembers()( MERGE_C));
        phi = merge_res.CalculatePhi( merge_res_prev.GetAtom( biol::GetAtomTypes().C));
      }
      else
      {
        // cannot compute the phi-angle for the first residue in the sequence
        BCL_MessageVrb( "Cannot compute phi-angle of c-terminal sequence. Using random angle instead.");
        phi = random::GetGlobalRandom().Double( math::Range< double>( 0, 2 * math::g_Pi));
      }
      double psi( 0.0);
      if( MERGE_N < ( int) SEQ_N.GetSize() - 1)
      {
        // compute the phi-angle for the residue located at the merging point
        const biol::AABase &merge_res_suc( *SEQ_N.GetMembers()( MERGE_N + 1));
        const biol::AABase &merge_res( *SEQ_N.GetMembers()( MERGE_N));
        psi = merge_res.CalculatePsi( merge_res_suc.GetAtom( biol::GetAtomTypes().N));
      }
      else
      {
        // cannot compute the phi-angle for the first residue in the sequence
        BCL_MessageVrb( "Cannot compute psi-angle of n-terminal sequence. Using random angle instead.");
        psi = random::GetGlobalRandom().Double( math::Range< double>( 0, 2 * math::g_Pi));
      }

      // append the c-terminal sequence to the n-terminal sequence
      util::ShPtr< biol::AASequence> sp_new_sequence( new biol::AASequence( seq_n));
      biol::AASequenceFactory::AppendSequence( *sp_new_sequence, seq_c, phi, psi);

      return sp_new_sequence;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads members from a given input stream
    //! @param ISTREAM input stream to read members from
    //! @return input stream which members were read from
    std::istream &ProteinGeometry::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief writes members into a given output stream
    //! @param OSTREAM output stream to write members into
    //! @param INDENT number of indentations
    //! @return output stream into which members were written
    std::ostream &ProteinGeometry::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
