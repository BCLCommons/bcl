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
#include "biol/bcl_biol_aa_sequence_factory.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_aa_back_bone_completer.h"
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "biol/bcl_biol_aa_sequence_phi_psi.h"
#include "fold/bcl_fold_mutate_aa_set_phi.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically
#include <iterator>

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AASequenceFactory::s_Instance
    (
      GetObjectInstances().AddInstance( new AASequenceFactory())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequenceFactory::AASequenceFactory()
    {
    }

    //! @brief Clone function
    //! @return pointer to new AASequenceFactory
    AASequenceFactory *AASequenceFactory::Clone() const
    {
      return new AASequenceFactory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequenceFactory::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief appends a residue to the C-terminal end of a sequence
    //! @param SEQUENCE AASequence to be extended
    //! @param AMINO_ACID to be appended
    //! @param PHI phi angle for added residue
    //! @param PSI psi angle for last residue in original sequence
    void AASequenceFactory::AppendAA
    (
      AASequence &SEQUENCE,
      const AABase &AMINO_ACID,
      const double PHI,
      const double PSI
    )
    {
      AppendSequence
      (
        SEQUENCE,
        AASequence
        (
          util::ShPtrVector< AABase>( 1, util::ShPtr< AABase>( AMINO_ACID.Clone())),
          SEQUENCE.GetChainID(),
          SEQUENCE.GetFastaHeader()
        ),
        PHI,
        PSI
      );
    }

    //! @brief appends a sequence to the C-terminal end of a sequence
    //! @param N_TERMINAL_SEQUENCE AASequence to be extended
    //! @param C_TERMINAL_SEQUENCE AASequence to be appended
    //! @param PHI phi angle for first residue in C-terminal sequence
    //! @param PSI psi angle for last residue in N-terminal sequence
    void AASequenceFactory::AppendSequence
    (
      AASequence &N_TERMINAL_SEQUENCE,
      const AASequence &C_TERMINAL_SEQUENCE,
      const double PHI,
      const double PSI
    )
    {
      // change psi angle
      {
        const int seq_id( N_TERMINAL_SEQUENCE.GetLastAA()->GetSeqID());
        const storage::VectorND< 2, double> new_phi_psi( 0.0, PSI);
        storage::VectorND< 2, double> new_phi_psi_change( AASequenceFlexibility::CalculatePhiPsiChange( N_TERMINAL_SEQUENCE, seq_id, new_phi_psi));
        new_phi_psi_change.First() = 0;
        AASequenceFlexibility::ChangePhiPsi( N_TERMINAL_SEQUENCE, seq_id, new_phi_psi_change, AASequenceFlexibility::e_CTerminal);
      }

      // calculate necessary transformation, apply and prepend
      const math::TransformationMatrix3D trans_c_term( TransformationAppend( N_TERMINAL_SEQUENCE, *C_TERMINAL_SEQUENCE.GetFirstAA(), PHI));
      AASequence c_term_seq( C_TERMINAL_SEQUENCE);
      c_term_seq.Transform( trans_c_term);
      N_TERMINAL_SEQUENCE.AppendSequence( c_term_seq);
    }

    //! @brief transformation to append a c term sequence onto the given n-term sequence
    //! @param N_TERMINAL_SEQUENCE AASequence fixed in space
    //! @param C_TERMINAL_AA amino acid to be appended
    //! @param PHI angle for first residue in cterm
    math::TransformationMatrix3D AASequenceFactory::TransformationAppend
    (
      const AASequence &N_TERMINAL_SEQUENCE,
      const AABase &C_TERMINAL_AA,
      const double PHI
    )
    {
      // make sure that the cterminal sequence has its coordinates defined
      BCL_Assert
      (
        N_TERMINAL_SEQUENCE.HasDefinedCoordinates(),
        "N-terminal aa sequence does not have all defined coordinates " + util::Format()( N_TERMINAL_SEQUENCE)
      );

      // make sure the n-terminal sequence has its coordinates defined
      BCL_Assert
      (
        C_TERMINAL_AA.HasDefinedCoordinates(),
        "C-terminal aa does not have all defined coordinates " + util::Format()( C_TERMINAL_AA)
      );

      // get the coordinates of the c atom of the last residue in the n-terminal sequence
      const linal::Vector3D &anchor_aa_coord_c
      (
        N_TERMINAL_SEQUENCE.GetLastAA()->GetAtom( GetAtomTypes().C).GetCoordinates()
      );

      // overall transformation
      math::TransformationMatrix3D overall_trans;

      // create ShPtr to a new c-terminal sequence cloned from "C_TERMINAL_SEQUENCE"
      util::ShPtr< AABase> new_c_terminal_aa( C_TERMINAL_AA.Clone());

    ///////////////////
    // superimpose N //
    ///////////////////

      {
        // calculate the coordinates of where the N of first residue in the cterminal sequence should go according to
        // "PSI"
        const Atom desired_nitrogen
        (
          AABackBoneCompleter::GenerateN( *N_TERMINAL_SEQUENCE.GetLastAA())
        );

        // get the current position of the nitrogen of the first residue of the c-terminal sequence
        const linal::Vector3D &current_coord_n
        (
          new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates()
        );

        // translate c-terminal sequence to correct position (at desired N position)
        const linal::Vector3D c_super_trans( desired_nitrogen.GetCoordinates() - current_coord_n);
        new_c_terminal_aa->Translate( c_super_trans);

        // accumulate overall transformation
        overall_trans( c_super_trans);
      }

    ///////////////////
    // superimpose C //
    ///////////////////

      {
        // desired n
        const Atom desired_carbon( AABackBoneCompleter::GenerateC( *new_c_terminal_aa, PHI));

        // cross product
        const linal::Vector3D rotation_axis
        (
          linal::CrossProduct
          (
            new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates(),
            desired_carbon.GetCoordinates(),
            anchor_aa_coord_c
          )
        );

        // rotation angle
        const double angle
        (
          linal::ProjAngle
          (
            new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates(),
            desired_carbon.GetCoordinates(),
            anchor_aa_coord_c
          )
        );

        // transformation matrix
        math::TransformationMatrix3D transformation( -new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates());

        // apply rotation to superimpose n onto n
        transformation( math::RotationMatrix3D( rotation_axis, -angle));

        // move back
        transformation( new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates());

        // apply to amino acid and overall transformation
        new_c_terminal_aa->Transform( transformation);
        overall_trans( transformation);
      }

    ///////////////
    // set omega //
    ///////////////

      {
        // calculate the rotation necessary in order to properly set omega to 180 degrees
        const double omega_rotation
        (
          AABackBoneCompleter::s_OmegaCACNCA - // optimal omega
          new_c_terminal_aa->CalculateOmega( N_TERMINAL_SEQUENCE.GetLastAA()->GetCA(), N_TERMINAL_SEQUENCE.GetLastAA()->GetAtom( GetAtomTypes().C))
        );

        // create transformation matrix and add the necessary transformations to set the omega angle
        math::TransformationMatrix3D omega_transform
        (
          -new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates()
        );

        // add the rotation to "c_n_ca_transform" to set the omega angle
        omega_transform
        (
          math::RotationMatrix3D
          (
            new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates() - anchor_aa_coord_c,
            -omega_rotation
          )
        );

        // add the translation to move the c-terminal sequence back to its original position but rotated
        omega_transform( new_c_terminal_aa->GetAtom( GetAtomTypes().N).GetCoordinates());

        // now transform the c-terminal sequence so that omega is properly set
        new_c_terminal_aa->Transform( omega_transform);

        // accumulate overall transformation
        overall_trans( omega_transform);
      }

    /////////////
    // set phi //
    /////////////

      {
        // transform the c-terminal sequence so that the phi of the first residue is set correctly to "PHI"
        const math::TransformationMatrix3D phi_trans
        (
          fold::MutateAASetPhi
          (
            N_TERMINAL_SEQUENCE.GetLastAA()->GetAtom( GetAtomTypes().C), PHI
          ).GetTransformationMatrix( *new_c_terminal_aa)
        );
//        new_c_terminal_aa->Transform
//        (
//          phi_trans
//        );

        // accumulate overall transformation
        overall_trans( phi_trans);
      }

      // return the overall transformation
      return overall_trans;
    }

    //! @brief prepends a residue to the N-terminal end of a sequence
    //! @param AMINO_ACID to be prepended
    //! @param SEQUENCE AASequence to be extended
    //! @param PHI phi angle for first residue in original sequence
    //! @param PSI psi angle for added residue
    void AASequenceFactory::PrependAA
    (
      const AABase &AMINO_ACID,
      AASequence &SEQUENCE,
      const double PHI,
      const double PSI
    )
    {
      PrependSequence
      (
        AASequence
        (
          util::ShPtrVector< AABase>( 1, util::ShPtr< AABase>( AMINO_ACID.Clone())),
          SEQUENCE.GetChainID(),
          SEQUENCE.GetFastaHeader()
        ),
        SEQUENCE,
        PHI,
        PSI
      );
    }

    //! @brief prepends a sequence to the N-terminal end of a sequence
    //! @param N_TERMINAL_SEQUENCE AASequence to be prepended
    //! @param C_TERMINAL_SEQUENCE AASequence to be extended
    //! @param PHI phi angle for first residue in C-terminal sequence
    //! @param PSI psi angle for last residue in N-terminal sequence
    void AASequenceFactory::PrependSequence
    (
      const AASequence &N_TERMINAL_SEQUENCE,
      AASequence &C_TERMINAL_SEQUENCE,
      const double PHI,
      const double PSI
    )
    {
      // create ShPtr to a new n-terminal sequence cloned from "N_TERMINAL_SEQUENCE"
      AASequence new_n_terminal_sequence( N_TERMINAL_SEQUENCE);

      // change psi angle
      {
        const int seq_id( new_n_terminal_sequence.GetLastAA()->GetSeqID());
        const storage::VectorND< 2, double> new_phi_psi( 0.0, PSI);
        storage::VectorND< 2, double> new_phi_psi_change( AASequenceFlexibility::CalculatePhiPsiChange( new_n_terminal_sequence, seq_id, new_phi_psi));
        new_phi_psi_change.First() = 0;
        AASequenceFlexibility::ChangePhiPsi( new_n_terminal_sequence, seq_id, new_phi_psi_change, AASequenceFlexibility::e_NTerminal);
      }

      // transformation matrix
      const math::TransformationMatrix3D trans_n_term( TransformationPrepend( *new_n_terminal_sequence.GetLastAA(), C_TERMINAL_SEQUENCE, PHI));
      new_n_terminal_sequence.Transform( trans_n_term);
      C_TERMINAL_SEQUENCE.PrependSequence( new_n_terminal_sequence);
    }

    //! @brief transformation to append a c term sequence onto the given n-term sequence
    //! @param N_TERMINAL_AA amino acid to be appended
    //! @param C_TERMINAL_SEQUENCE AASequence fixed in space
    //! @param PHI angle for first residue in cterm
    math::TransformationMatrix3D AASequenceFactory::TransformationPrepend
    (
      const AABase &N_TERMINAL_AA,
      const AASequence &C_TERMINAL_SEQUENCE,
      const double PHI
    )
    {
      // make sure that the n-terminal residue has its coordinates defined
      BCL_Assert
      (
        N_TERMINAL_AA.HasDefinedCoordinates(),
        "N-terminal aa does not have all defined coordinates " + util::Format()( N_TERMINAL_AA)
      );

      // make sure that the c-terminal sequence has its coordinates defined
      BCL_Assert
      (
        C_TERMINAL_SEQUENCE.HasDefinedCoordinates(),
        "C-terminal aa sequence does not have all defined coordinates " + util::Format()( C_TERMINAL_SEQUENCE)
      );

      // get the coordinates of the n atom of the first residue in the c-terminal sequence
      const linal::Vector3D &anchor_aa_coord_n
      (
        C_TERMINAL_SEQUENCE.GetFirstAA()->GetAtom( GetAtomTypes().N).GetCoordinates()
      );

      // create ShPtr to a new n-terminal aa
      util::ShPtr< AABase> new_n_term_aa( N_TERMINAL_AA.Clone());

      // overall transformation
      math::TransformationMatrix3D overall_trans;

    ///////////////////
    // superimpose C //
    ///////////////////

      {
        // calculate the coordinates of where the C of last residue in the nterminal sequence should go according to
        // "PHI"
        const Atom desired_carbon
        (
          AABackBoneCompleter::GenerateC( *C_TERMINAL_SEQUENCE.GetFirstAA(), PHI)
        );

        // get the current position of the carbon of the last residue of the n-terminal sequence
        const linal::Vector3D &current_coord_c
        (
          new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates()
        );

        // translate n-terminal sequence to correct position (at desired C position)
        const linal::Vector3D super_c_trans( desired_carbon.GetCoordinates() - current_coord_c);
        new_n_term_aa->Translate( super_c_trans);

        // accumulate overall transformation
        overall_trans( super_c_trans);
      }

    ///////////////////
    // superimpose N //
    ///////////////////

      {
        // desired n
        const Atom desired_nitrogen( AABackBoneCompleter::GenerateN( *new_n_term_aa));

        // cross product
        const linal::Vector3D rotation_axis
        (
          linal::CrossProduct
          (
            new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates(),
            desired_nitrogen.GetCoordinates(),
            anchor_aa_coord_n
          )
        );

        // rotation angle
        const double angle
        (
          linal::ProjAngle
          (
            new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates(),
            desired_nitrogen.GetCoordinates(),
            anchor_aa_coord_n
          )
        );

        // transformation matrix
        math::TransformationMatrix3D transformation( -new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates());

        // apply rotation to superimpose n onto n
        transformation( math::RotationMatrix3D( rotation_axis, -angle));

        // move back
        transformation( new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates());

        // apply to amino acid and overall transformation
        new_n_term_aa->Transform( transformation);
        overall_trans( transformation);
      }

    ///////////////
    // set omega //
    ///////////////

      {
        // calculate the rotation necessary in order to properly set omega to 180 degrees
        const double omega_rotation
        (
          AABackBoneCompleter::s_OmegaCACNCA - // optimal omega
          C_TERMINAL_SEQUENCE.GetFirstAA()->CalculateOmega( new_n_term_aa->GetCA(), new_n_term_aa->GetAtom( GetAtomTypes().C))
        );

        // create transformation matrix and add the necessary transformations to set the omega angle
        math::TransformationMatrix3D omega_transform
        (
          -anchor_aa_coord_n
        );

        // creat rotation matrix around peptide bond  and add to omega tranformation
        omega_transform
        (
          math::RotationMatrix3D
          (
            anchor_aa_coord_n - new_n_term_aa->GetAtom( GetAtomTypes().C).GetCoordinates(),
            omega_rotation //< rotate in the correct direction
          )
        );

        // add the translation to move the n-terminal sequence back to its original position but rotated
        omega_transform( anchor_aa_coord_n);

        // now transform the n-terminal sequence so that omega is properly set
        new_n_term_aa->Transform( omega_transform);

        // accumulate overall transformation
        overall_trans( omega_transform);
      }

    ////////////////////////
    // phi transformation //
    ////////////////////////

      {
        // transform the n-terminal aa so that the phi is correct relative to the c terminal sequence
        const math::TransformationMatrix3D phi_trans
        (
          fold::MutateAASetPhi
          (
            new_n_term_aa->GetAtom( GetAtomTypes().C), PHI
          ).GetTransformationMatrix( *C_TERMINAL_SEQUENCE.GetFirstAA())
        );
//        new_n_term_aa->Transform
//        (
//          phi_trans
//        );

        // accumulate overall transformation
        overall_trans( phi_trans);
      }

      // return the overall transformation applied
      return overall_trans;
    }

    //! @brief fits the passed sequence to the phi/psi angles
    //! @param SEQUENCE AASequence to be used
    //! @param PHI_PSI phi and psi values to be applied
    //! @param SS_TYPE sstype to be applied if the sequence is longer than the given phi/psi information
    void AASequenceFactory::FitSequence
    (
      AASequence &SEQUENCE,
      const AASequencePhiPsi &PHI_PSI,
      const SSType &SS_TYPE
    )
    {
      // if the sequence is empty or no phi/psi's are given
      if( SEQUENCE.GetSize() == 0 || PHI_PSI.GetAngles().IsEmpty())
      {
        // return with no changes
        return;
      }

      // idealize the sequence according to the given type
      IdealizeSequence( SEQUENCE, SS_TYPE);

      // get the middle index
      const size_t middle_seq_index( SEQUENCE.GetSize() / 2);

      // move the sequence so that the middle residues are superimposed
      SEQUENCE.Transform
      (
        quality::RMSD::SuperimposeCoordinates
        (
          util::SiPtrVector< const linal::Vector3D>::Create
          (
            PHI_PSI.GetN(),
            PHI_PSI.GetCA(),
            PHI_PSI.GetC()
          ),
          util::SiPtrVector< const linal::Vector3D>::Create
          (
            SEQUENCE.GetData()( middle_seq_index)->GetAtom( GetAtomTypes().N).GetCoordinates(),
            SEQUENCE.GetData()( middle_seq_index)->GetAtom( GetAtomTypes().CA).GetCoordinates(),
            SEQUENCE.GetData()( middle_seq_index)->GetAtom( GetAtomTypes().C).GetCoordinates()
          )
        )
      );

      // determine the corresponding indeces to be used
      const size_t middle_angle_index( PHI_PSI.GetAngles().GetSize() / 2);
      size_t first_index( 0);
      size_t last_index( PHI_PSI.GetAngles().GetSize() - 1);
      int first_seq_id( SEQUENCE.GetFirstAA()->GetSeqID());

      // if the sequence has more residues than phi/psi angles given
      if( SEQUENCE.GetSize() > PHI_PSI.GetAngles().GetSize())
      {
        // move up the first seq id to match with the beginning of the angle vector
        first_seq_id += ( SEQUENCE.GetSize() - PHI_PSI.GetAngles().GetSize()) / 2;
      }

      // if the angle vector is larger than the sequence
      if( PHI_PSI.GetAngles().GetSize() > SEQUENCE.GetSize())
      {
        // adjust the indeces to match the vector size
        first_index += ( PHI_PSI.GetAngles().GetSize() - SEQUENCE.GetSize()) / 2;
        last_index = first_index + SEQUENCE.GetSize() - 1;

        // move the first seq id to compensate
        first_seq_id -= first_index;
      }

      // iterate backwards from the residue
      for( int i( middle_angle_index - 1); i >= int( first_index); --i)
      {
        // set the previous phi/psi's
        AASequenceFlexibility::SetPhiPsi
        (
          SEQUENCE, first_seq_id + i, PHI_PSI.GetAngles()( i), AASequenceFlexibility::e_NTerminal
        );
      }

      // set the middle residue
      AASequenceFlexibility::SetPhiPsi
      (
        SEQUENCE,
        first_seq_id + middle_angle_index,
        PHI_PSI.GetAngles()( middle_angle_index),
        AASequenceFlexibility::e_Bidirectional
      );

      // iterate forwards from the residue
      for( int i( middle_angle_index + 1); i <= int( last_index); ++i)
      {
        // set the next phi/psi's
        AASequenceFlexibility::SetPhiPsi
        (
          SEQUENCE, first_seq_id + i, PHI_PSI.GetAngles()( i), AASequenceFlexibility::e_CTerminal
        );
      }
    }

    //! @brief Idealize the coordinates for the given AASequence according to given SSType
    //! @param AA_SEQUENCE AASequence to be idealized
    //! @param SS_TYPE SSType to be used for idealization
    void AASequenceFactory::IdealizeSequence( AASequence &AA_SEQUENCE, const SSType &SS_TYPE)
    {
      // find the transformation matrix that transforms one residue to the following residue in this SSType,
      // including the translation related to z axis( half the height of this SSE)
      math::TransformationMatrix3D transform
      (
        linal::Vector3D( 0.0, 0.0, -1.0 * AA_SEQUENCE.GetSize() * 0.5 * SS_TYPE->GetRiseInZPerResidue())
      );

      // iterate over aa's
      for
      (
        AASequence::iterator
          aa_itr( AA_SEQUENCE.GetData().Begin()), aa_itr_end( AA_SEQUENCE.GetData().End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // first set the residues to ideal transformed conformation
        ( *aa_itr)->SetToIdealConformation( SS_TYPE, transform);

        // transform the transformation matrix
        transform( SS_TYPE->GetTransformationMatrixForResidues());
      }
    }

    //! @brief read fasta from given ISTREAM and initialize the sequence with instances of AA_CLASS
    //! @param ISTREAM input stream
    //! @param AA_CLASS type of amino acid to be used
    //! @param CHAIN_ID chain id to be used for the sequence
    //! @return newly constructed AASequence
    AASequence AASequenceFactory::BuildSequenceFromFASTA
    (
      std::istream &ISTREAM,
      const AAClass &AA_CLASS,
      const char CHAIN_ID
    )
    {
      // local variables
      char c;

      // read header
      // whitespace is skipped automatically (unless you add "<< noskipws")
      std::string fasta_header( ">");
      ISTREAM >> c;

      // if ">" is found as first character, assume there is a header
      if( c == AASequence::s_FastaHeaderChar)
      {
        // read the complete first line and save in m_FastaHeader
        std::string line;
        std::getline( ISTREAM, line);
        fasta_header += line;
      }
      // if there is no ">"
      else
      {
        // put back the read character and create a new header
        ISTREAM.putback( c);
        fasta_header += AASequence::GetDefaultFastaHeader();
      }

      // initialize iterators on the stream to read the sequence
      std::istream_iterator< char> stream_itr( ISTREAM);
      std::istream_iterator< char> stream_itr_end;
      size_t id = 1;

      // initialize vector of aabase
      util::ShPtrVector< AABase> aa_data;

      // while there are chars left in STREAM
      while( stream_itr != stream_itr_end)
      {
        // if current_oneletter_code is "*", we are at the end of the sequence, then break
        if( *stream_itr == '*')
        {
          break;
        }

        // if current_oneletter_code is ">", we found the next sequence, so putback the char and end this sequence
        if( *stream_itr == AASequence::s_FastaHeaderChar)
        {
          ISTREAM.unget(); // unget the ">" character that the next sequence can be read correctly
          break;
        }

        // get aatype from current_oneletter_code
        const AAType current_aatype( GetAATypes().AATypeFromOneLetterCode( toupper( *stream_itr)));

        // if current aa type is undefined
        if
        (
          current_aatype == GetAATypes().e_Undefined &&
          ( *stream_itr) != GetAATypes().e_Undefined->GetOneLetterCode()
        )
        {
          // give user a message on every unkown aatype we found
          BCL_MessageCrt
          (
            std::string( "Unknown one letter code: ") + std::string( 1, *stream_itr)
            + " (ASCII: " + util::Format()( int( *stream_itr)) + ") supplied!"
          );
          break;
        }

        util::ShPtr< AAData> sp_aa_data( new AAData( current_aatype, int( id), int( id), AAData::s_DefaultPdbICode, CHAIN_ID));

        // initialize ShPtr to an amino acid of type AA_CLASS corresponding AAData
        util::ShPtr< AABase> sp_aa( ( *AA_CLASS)->Empty( sp_aa_data));

        // insert this aminoacid into the sequence
        aa_data.PushBack( sp_aa);

        // move further in the stream and in the sequence ids
        ++stream_itr;
        ++id;
      }

      // end
      return AASequence( aa_data, CHAIN_ID, fasta_header);
    }

    //! @brief read fasta from provided SEQUENCE_STRING and initialize the sequences with amino acids of type AA_CLASS
    //! @param SEQUENCE_STRING sequence string that has the one letter code sequence
    //! @param AA_CLASS type of amino acid to be used
    //! @param CHAIN_ID chain id to be used for the sequence
    //! @return AA sequence built from the string
    AASequence AASequenceFactory::BuildSequenceFromFASTAString
    (
      const std::string &SEQUENCE_STRING,
      const AAClass &AA_CLASS,
      const char CHAIN_ID
    )
    {
      // read the string into a stream
      std::stringstream ss;
      ss << SEQUENCE_STRING;

      // call BuildSequenceFromFASTA function with this stream
      return BuildSequenceFromFASTA( ss, AA_CLASS, CHAIN_ID);
    }

    //! @brief read sequence data, fasta and pdb files, from the given filenames and returns their sequences
    //! @param FILENAMES a vector of sequence filenames
    //! @return a ShPtrVector of sequences
    util::ShPtrVector< AASequence>
    AASequenceFactory::ReadSequenceData( const storage::Vector< std::string> &FILENAMES)
    {
      util::ShPtrVector< AASequence> sequences; // instantiate ShPtrVector to collect all sequences
      io::IFStream read; // IFStream for reading files
      pdb::Factory factory; // pdb::Factory to reading pdbs

      for
      (
        storage::Vector< std::string>::const_iterator itr( FILENAMES.Begin()), itr_end( FILENAMES.End());
        itr != itr_end;
        ++itr
      )
      {
        io::File::MustOpenIFStream( read, *itr);
        std::string extension( io::File::GetFullExtension( io::File::RemovePath( *itr)));

        // read sequences from different file formats
        // searching for the format extension in all extensions allow also compressed files
        if( extension.find( "fasta") != std::string::npos)
        {
          while( !read.eof()) // read all sequences in this file
          {
            // new AASequence to read next sequence
            util::ShPtr< AASequence> sequence( new AASequence( AASequenceFactory::BuildSequenceFromFASTA( read)));
            sequences.PushBack( sequence); // save sequence
          }
        }
        else if( extension.find( pdb::GetDefaultFileExtension()) != std::string::npos)
        {
          pdb::Handler pdb( read);
          util::ShPtrVector< AASequence> pdb_sequences( factory.AASequencesFromPDB( pdb));
          sequences.InsertElements( sequences.End(), pdb_sequences);
        }

        io::File::CloseClearFStream( read);
      }

      return sequences;
    }

    //! @brief calculates the transformation needed for superimposing SEQUENCE_A to SEQUENCE_B
    //! @param AA_SEQUENCE_A AASequence which will be superimposed to AA_SEQUENCE_B
    //! @param AA_SEQUENCE_B AASequence to which superimposition will be calculated against
    //! @return transformation matrix that gives the best superimposition
    math::TransformationMatrix3D AASequenceFactory::CalculateSuperimposition
    (
      const AASequence &AA_SEQUENCE_A,
      const AASequence &AA_SEQUENCE_B
    )
    {
      // if this is the same AASequence then return identity matrix
      if( &AA_SEQUENCE_A == &AA_SEQUENCE_B)
      {
        return math::TransformationMatrix3D();
      }

      // calculate transformation matrix for super imposition
      math::TransformationMatrix3D trans_matrix
      (
        quality::RMSD::SuperimposeCoordinates
        (
          AA_SEQUENCE_B.GetAtomCoordinates( GetAtomTypes().GetBackBoneAtomTypes()),
          AA_SEQUENCE_A.GetAtomCoordinates( GetAtomTypes().GetBackBoneAtomTypes())
        )
      );

      // problem with superimposition? try to superimpose all coordinates
      if( trans_matrix == math::TransformationMatrix3D())
      {
        trans_matrix = quality::RMSD::SuperimposeCoordinates
        (
          AA_SEQUENCE_B.GetAtomCoordinates(),
          AA_SEQUENCE_A.GetAtomCoordinates()
        );
      }

      // end
      return trans_matrix;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequenceFactory::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AASequenceFactory::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
