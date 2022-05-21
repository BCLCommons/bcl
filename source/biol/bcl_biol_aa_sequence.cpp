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
#include "biol/bcl_biol_aa_sequence.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_compare.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_atom.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AASequence::s_Instance
    (
      GetObjectInstances().AddInstance( new AASequence())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASequence::AASequence() :
      m_Data(),
      m_ChainID( s_DefaultChainID),
      m_FastaHeader( GetDefaultFastaHeader())
    {
      Initialize( GetAAClasses().e_AA, 0);
    }

    //! @brief construct from AACLASS, LENGTH, CHAIN_ID and FASTA_HEADER
    //! @param AACLASS AAClass type for m_Data
    //! @param LENGTH number of residues in this sequence
    //! @param CHAIN_ID chain id of sequence
    //! @param FASTA_HEADER fasta header
    AASequence::AASequence
    (
      const AAClass &AACLASS,
      const size_t LENGTH,
      const char CHAIN_ID,
      const std::string &FASTA_HEADER
    ) :
      m_Data(),
      m_ChainID( CHAIN_ID)
    {
      SetFastaHeader( FASTA_HEADER);
      Initialize( AACLASS, LENGTH);
    }

    // TODO: can not ensure the people pass DATA from a sequence with AAClass of same type with this->m_AAClass
    //! @brief construct from ShPtrVector of AABases DATA, CHAIN_ID and FASTA_HEADER
    //! @param DATA ShPtrVector of AABases
    //! @param CHAIN_ID chain id of sequence
    //! @param FASTA_HEADER fasta header
    AASequence::AASequence
    (
      const util::ShPtrVector< AABase> &DATA,
      const char CHAIN_ID,
      const std::string &FASTA_HEADER
    ) :
      m_Data( DATA.HardCopy()),
      m_ChainID( CHAIN_ID),
      m_FastaHeader()
    {
      BCL_Assert( IsContinuous( m_Data), "given data are not continuous amino acids");
      BCL_Assert
      (
        m_Data.IsEmpty() || m_Data.FirstElement()->GetChainID() == m_ChainID,
        "first amino acid has different chain id than passed chainid! first amino acid chain id |" +
        m_Data.FirstElement()->GetIdentification() + "| The given chain id is  |" + m_ChainID + "|"
      );
      SetFastaHeader( FASTA_HEADER);
    }

    //! @brief copy constructor from another AA_SEQUENCE
    //! @param AA_SEQUENCE AASequence to be copied
    AASequence::AASequence( const AASequence &AA_SEQUENCE) :
      m_Data( AA_SEQUENCE.m_Data.HardCopy()),
      m_ChainID( AA_SEQUENCE.m_ChainID),
      m_FastaHeader( AA_SEQUENCE.m_FastaHeader)
    {
    }

    //! @brief virtual copy constructor
    AASequence *AASequence::Clone() const
    {
      return new AASequence( *this);
    }

    //! @brief hard copy all amino acids and their data
    //! @return aa sequence with independent hard copied AADatas
    AASequence *AASequence::HardCopy() const
    {
      // copy for this sequence
      AASequence *new_sequence( Clone());

      // iterate over all amino acids to make hard copy of the Data
      for
      (
        AASequence::iterator aa_itr( new_sequence->Begin()), aa_itr_end( new_sequence->End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // hard copy the data for that aa
        ( *aa_itr)->SetData( ( *aa_itr)->GetData().HardCopy());
      }

      // end
      return new_sequence;
    }

    //! @brief construct a new SequenceInterface from an AASequence interface implementation, sequence id and members
    //! @param ID sequence identifier
    //! @param MEMBERS sequence members
    //! @param CHAIN_ID chain id used to build the sequence
    //! @return pointer to a SequenceInterface
    align::SequenceInterface< AABase> *AASequence::Construct
    (
      const std::string &ID,
      const std::string &MEMBERS,
      const char CHAIN_ID
    ) const
    {
      AASequence *p_aa_sequence
      (
        AASequenceFactory::BuildSequenceFromFASTAString( MEMBERS, GetAAClasses().e_AA, CHAIN_ID).Clone()
      );
      p_aa_sequence->SetFastaHeader( ">" + ID);
      return p_aa_sequence;
    }

    //! @brief destructor
    AASequence::~AASequence()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AASequence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get default fasta header
    //! @return default fasta header
    const std::string &AASequence::GetDefaultFastaHeader()
    {
      // initialize static const default fasta header
      static const std::string s_default_fasta_header( "BCL_AASequence");

      // return
      return s_default_fasta_header;
    }

    //! @brief get identification of this AASequence
    //! @return string with identification for this AASequence
    std::string AASequence::GetSequenceIdentification() const
    {
      // initialize identification composed of chain ID, first seq id, the sequence one-letter codes, last seq-id
      // example: "chain: 'A' 23 VLALLV 28"
      std::string identification( "chain: " + util::Format()( m_ChainID));

      // add identifiers for the first and last residues this sse is spanning
      if( GetSize() > 0)
      {
        identification += ' ' + util::Format()( GetFirstAA()->GetSeqID())
          + ' ' + Sequence() + ' ' + util::Format()( GetLastAA()->GetSeqID());
      }

      // end
      return identification;
    }

    //! @brief add a AABase to the sequence with type based on ONE_LETTER_CODE
    //! @param ONE_LETTER_CODE char encoding the type of the AA
    void AASequence::AddMember( const char &ONE_LETTER_CODE)
    {
      // if m_Data is empty, set a default type of AAClass, otherwise get the type from the first AA
      const AAClass new_member_type( m_Data.IsEmpty() ? GetAAClasses().e_AA : GetFirstAA()->GetAAClass());

      // create ShPtr to AAData and set aatype, seq_id, chain_id
      util::ShPtr< AAData> new_member_data
      (
        new AAData
        (
          GetAATypes().AATypeFromOneLetterCode( toupper( ONE_LETTER_CODE)),
          m_Data.IsEmpty() ? AAData::s_DefaultSeqID : m_Data.LastElement()->GetSeqID() + 1,
          AAData::s_DefaultPdbID,
          AAData::s_DefaultPdbICode,
          m_ChainID
        )
      );

      // create ShPtr to AABase
      util::ShPtr< AABase> new_member( ( *new_member_type)->Empty( new_member_data));

      // push back in m_Data
      m_Data.PushBack( new_member);
    }

    //! @brief initializes the sequence with a new type of AAClass and fills with LENGTH number of residues
    //! deletes all existing aminoacids and builds up a new sequence with empty AABase interfaces
    //! @param AA_CLASS AAClass that identifies the amino acid type to be used
    //! @param LENGTH number of residues to be inserted
    void AASequence::Initialize( const AAClass &AA_CLASS, const size_t LENGTH)
    {
      // reset m_Data
      m_Data.Reset();

      // allocate memory for LENGTH number of residues
      m_Data.AllocateMemory( LENGTH);

      // iterate until LENGTH
      for( size_t i( 0); i < LENGTH; ++i)
      {
        util::ShPtr< AAData> sp_aa_data
        (
          new AAData( GetAATypes().e_Undefined, i + 1, i + 1, AAData::s_DefaultPdbICode, m_ChainID)
        );

        // create a ShPtr of AABase type with instance of AA_CLASS behind it
        util::ShPtr< AABase> sp_aa( ( *AA_CLASS)->Empty( sp_aa_data));

        // insert this new amino acid ShPtr into m_Data
        m_Data.PushBack( sp_aa);
      }
    }

    //! @brief get all atoms of all residues in this sequence
    //! @return SiPtrVector of Atoms of all residues in this sequence
    util::SiPtrVector< const Atom> AASequence::GetAtoms() const
    {
      // initialize atoms SiPtrVector to be returned
      util::SiPtrVector< const Atom> atoms;

      // allocate memory that is enough for total number of atoms to be inserted
      if( !m_Data.IsEmpty())
      {
        atoms.AllocateMemory( GetFirstAA()->GetNumberOfAtoms() * GetSize());
      }

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        //get atoms from the amino acid
        atoms.Append( ( *aa_itr)->GetAtoms());
      }

      // return
      return atoms;
    }

    //! @brief get all atoms for the specified atom types for all residues in the sequence
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return atoms of types specified in ATOM_TYPES  for all residues in the sequence
    util::SiPtrVector< const Atom> AASequence::GetAtoms
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // initialize atoms SiPtrVector to be returned
      util::SiPtrVector< const Atom> atoms;

      // allocate memory that is enough for storing all atoms of specified ATOM_TYPES
      atoms.AllocateMemory( ATOM_TYPES.GetSize() * GetSize());

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get atoms from the amino acid for specified ATOM_TYPES
        atoms.Append( ( *aa_itr)->GetAtoms( ATOM_TYPES));
      }

      // return
      return atoms;
    }

    //! @brief get all atom coordinates for all residues in the sequence
    //! @return all atom coordinates for all residues in the sequence
    util::SiPtrVector< const linal::Vector3D> AASequence::GetAtomCoordinates() const
    {
      // initialize coordinates SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // allocate memory that is enough for storing all atom coordinates
      if( !m_Data.IsEmpty())
      {
        coordinates.AllocateMemory( GetFirstAA()->GetNumberOfAtoms() * GetSize());
      }

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
         aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get coordinates from the amino acid
        coordinates.Append( ( *aa_itr)->GetAtomCoordinates());
      }

      // return
      return coordinates;
    }

    //! @brief get all atom coordinates for the specified atom types for all residues in the sequence
    //! @param ATOM_TYPES AtomTypes of interest
    //! @return coordinates of atoms of types specified in ATOM_TYPES  for all residues in the sequence
    util::SiPtrVector< const linal::Vector3D> AASequence::GetAtomCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES
    ) const
    {
      // initialize coordinates SiPtrVector to be returned
      util::SiPtrVector< const linal::Vector3D> coordinates;

      // allocate memory that is enough for storing all atom coordinates of specified ATOM_TYPES
      coordinates.AllocateMemory( ATOM_TYPES.GetSize() * GetSize());

      // iterate over all amino acids
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // get coordinates from the amino acid for specified ATOM_TYPES
        coordinates.Append( ( *aa_itr)->GetAtomCoordinates( ATOM_TYPES));
      }

      // return
      return coordinates;
    }

    //! @brief find and return the iterator to the amino acid with the given sequence id
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @return the iterator to the amino acid with the given sequence id
    AASequence::iterator AASequence::FindAABySeqID( const int SEQ_ID)
    {
      // find and return by seq_id
      return std::find_if( Begin(), End(), AACompareBySeqID( SEQ_ID));
    }

    //! @brief find and return the iterator to the amino acid with the given sequence id
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @return the iterator to the amino acid with the given sequence id
    AASequence::const_iterator AASequence::FindAABySeqID( const int SEQ_ID) const
    {
      // find and return by seq_id
      return std::find_if( Begin(), End(), AACompareBySeqID( SEQ_ID));
    }

    //! @brief calculates the phi and psi for the amino acid with the given seq id
    //! @param SEQ_ID sequence id of the amino acid of interest
    //! @return pair of phi and psi values for the specified amino acid
    storage::VectorND< 2, double> AASequence::CalculatePhiPsi( const int SEQ_ID) const
    {
      // static undefined vector
      static const storage::VectorND< 2, double> s_undefined_vector( util::GetUndefined< double>());

      // find the amino acid
      const AASequence::const_iterator aa_itr( FindAABySeqID( SEQ_ID));

      // if not found return undefined
      if( aa_itr == m_Data.End())
      {
        return s_undefined_vector;
      }

      // construct the pair
      storage::VectorND< 2, double> phi_psi_pair( s_undefined_vector);

      // if the first residue
      if( aa_itr != m_Data.Begin())
      {
        phi_psi_pair.First() = ( *aa_itr)->CalculatePhi( ( *( aa_itr - 1))->GetAtom( GetAtomTypes().C));
      }
      else
      {
        phi_psi_pair.First() = ( *aa_itr)->Phi();
      }

      // if not the last residue
      if( aa_itr != m_Data.Last())
      {
        phi_psi_pair.Second() = ( *aa_itr)->CalculatePsi( ( *( aa_itr + 1))->GetAtom( GetAtomTypes().N));
      }
      else
      {
        phi_psi_pair.Second() = ( *aa_itr)->Psi();
      }

      // return
      return phi_psi_pair;
    }

    //! @brief set chain id of the sequence and all residues in the sequence to given CHAIN_ID
    //! @param CHAIN_ID new chain id of the sequence
    void AASequence::SetChainID( const char CHAIN_ID)
    {
      // first update the chain id of the sequence to CHAIN_ID
      m_ChainID = CHAIN_ID;

      // iterate over all residues in the sequence
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // create a new amino acid with new aa data with the new CHAIN_ID
        util::ShPtr< AABase> sp_aa( ( *aa_itr)->Clone());
        util::ShPtr< AAData> sp_aa_data
        (
          new AAData( sp_aa->GetType(), sp_aa->GetSeqID(), sp_aa->GetPdbID(), sp_aa->GetPdbICode(), CHAIN_ID)
        );

        sp_aa->SetData( sp_aa_data);
        if( ( *aa_itr)->GetBlastProfilePtr().IsDefined())
        {
          sp_aa->SetBlastProfile( ( *aa_itr)->GetBlastProfile());
        }
        sp_aa->SetSSPredictions( ( *aa_itr)->GetSSPredictions());
        *aa_itr = sp_aa;
      }
    }

    //! @brief returns the amino acids with the given type
    //! @param AA_TYPE type of the amino acid to return
    //! @return shared pointers to the amino acids of the the given type
    const util::ShPtrList< AABase> AASequence::GetData( const AAType &AA_TYPE) const
    {
      // iterate over the sequence to find the amino acids of the given type
      util::ShPtrList< AABase> aa_list;
      typedef util::ShPtrVector< AABase>::const_iterator aa_it;
      for( aa_it it( Begin()), it_end( End()); it != it_end; ++it)
      {
        if( ( *it)->GetType() == AA_TYPE)
        {
          aa_list.Append( *it);
        }
      }

      return aa_list;
    }

    //! @brief change the fasta header of the sequence to given FASTA_HEADER
    //! @param FASTA_HEADER new fasta header of the sequence
    void AASequence::SetFastaHeader( const std::string &FASTA_HEADER)
    {
      m_FastaHeader = util::TrimString( FASTA_HEADER);

      if( !m_FastaHeader.empty() && m_FastaHeader[ 0] == s_FastaHeaderChar)
      {
        m_FastaHeader = m_FastaHeader.substr( 1, m_FastaHeader.size() - 1);
      }
      m_FastaHeader = util::TrimString( m_FastaHeader);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief translate the object along a given TRANSLATION vector
    //! @param TRANSLATION translation vector to be applied
    void AASequence::Translate( const linal::Vector3D &TRANSLATION)
    {
      //iterate over all amino acids
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // translate the amino acid coordinates
        ( *aa_itr)->Translate( TRANSLATION);
      }
    }

    //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
    //! @param TRANSFORMATIONMATRIX3D TransformationMatrix3D to be applied
    void AASequence::Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D)
    {
      //iterate over all amino acids
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // transform the amino acid coordinates
        ( *aa_itr)->Transform( TRANSFORMATIONMATRIX3D);
      }
    }

    //! @brief rotate the object by a given ROTATION_MATRIX_3D
    //! @param ROTATIONMATRIX3D RotationMatrix3D to be applied
    void AASequence::Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D)
    {
      //iterate over all amino acids
      for
      (
        AASequence::iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // rotate the amino acid coordinates
        ( *aa_itr)->Rotate( ROTATIONMATRIX3D);
      }
    }

    //! return util::ShPtrVector of AASequences chopped to pieces of size SIZE with one aa gap
    util::ShPtrVector< AASequence> AASequence::ChopSequence( const size_t SIZE) const
    {
      BCL_Assert( util::IsDefined( SIZE), "undefined size supplied");

      //calculate the numer of pieces
      size_t number_of_sequences = ( GetSize() + 1) / ( SIZE + 1);

      //instantiate util::ShPtrVector to pieces
      util::ShPtrVector< AASequence> chopped_sequence;

      //nothing to be chopped
      if( number_of_sequences < 2)
      {
        chopped_sequence.PushBack( util::ShPtr< AASequence>( new AASequence( *this)));
        return chopped_sequence;
      }

      // pushback last sequence which does not provide a gap
      chopped_sequence.PushBack
                       (
                         util::ShPtr< AASequence>
                         (
                           new AASequence
                           (
                             SubSequence
                             (
                               0, //position of subsequence
                               size_t( math::Absolute( double( 1) * double( GetSize()) / double( number_of_sequences))) //length of subsequence
                             )
                           )
                         )
                       );

      //collect chopped sequences and make an one aagap between sequences
      for( size_t i( 1); i < number_of_sequences - 1; ++i)
      {
        chopped_sequence.PushBack
        (
          util::ShPtr< AASequence>
          (
            new AASequence
            (
              SubSequence
              (
                size_t( math::Absolute( double( i    ) * double( GetSize()) / double( number_of_sequences))) + 1,  //position of subsequence
                size_t( math::Absolute( double( i + 1) * double( GetSize()) / double( number_of_sequences)))       //length of subsequence
                - size_t( math::Absolute( double( i    ) * double( GetSize()) / double( number_of_sequences))) - 1 //length of subsequence
              )
            )
          )
        );
      }

      // pushback last sequence which does not provide a gap
      chopped_sequence.PushBack
      (
        util::ShPtr< AASequence>
        (
          new AASequence
          (
            SubSequence
            (
              size_t( math::Absolute( double( number_of_sequences - 1) * double( GetSize()) / double( number_of_sequences))) + 1, //position of subsequence
              GetSize() - size_t( math::Absolute( double( number_of_sequences - 1) * double( GetSize()) / double( number_of_sequences))) - 1 //length of subsequence
            )
          )
        )
      );

      //return chopped sequences
      return chopped_sequence;
    }

    //! @brief clips number of amino acidss from the beginning and the end of the sequence
    //! @param NUMBER_AMINO_ACIDS number of amino acids to be clipped from each end
    //! @return this sequence after clipping NUMBER_AMINO_ACIDS from each end
    AASequence &AASequence::ClipEnds( const size_t NUMBER_AMINO_ACIDS)
    {
      // check that after clipping we are left with at least 1 amino acid in the sequence
      if( m_Data.GetSize() <= 2 * NUMBER_AMINO_ACIDS + 1)
      {
        BCL_MessageCrt
        (
          "sequence too short: " + util::Format()( m_Data.GetSize()) +
          " <= " + util::Format()( 2 * NUMBER_AMINO_ACIDS + 1)
        );

        // end
        return *this;
      };
      // create a subsequence accordingly
      operator =( SubSequence( NUMBER_AMINO_ACIDS, m_Data.GetSize() - 2 * NUMBER_AMINO_ACIDS));

      // return
      return *this;
    }

    //! @brief prepends new amino acid to the beginning of the sequence
    //! @param AMINO_ACID amino acid to be prepended to the sequence
    void AASequence::PushFront( const AABase &AMINO_ACID)
    {
      if( m_Data.IsEmpty())
      {
        m_ChainID = AMINO_ACID.GetChainID();
      }
      else
      {
        // check that amino acid precedes sequence
        BCL_Assert
        (
          AMINO_ACID.DoesPrecede( *GetFirstAA()),
          "amino acid does not preced this: " + AMINO_ACID.GetIdentification() + " .. " + GetSequenceIdentification()
        );
      }

      // clone the AMINO_ACID and pushback it into the m_Data
      m_Data.InsertElement( m_Data.Begin(), util::ShPtr< AABase>( AMINO_ACID.Clone()));
    }

    //! @brief appends util::ShPtr to amino acid to the end of the sequence
    //! @param SP_AMINO_ACID ShPtr to amino acid to be appended to the end of the sequence
    void AASequence::PushBack( const util::ShPtr< AABase> &SP_AMINO_ACID)
    {
      if( m_Data.IsEmpty())
      {
        m_ChainID = SP_AMINO_ACID->GetChainID();
      }
      else
      {
        // check that amino acids follows this sequence
        BCL_Assert
        (
          GetLastAA()->DoesPrecede( *SP_AMINO_ACID),
          "amino acid does not follow this: " + GetSequenceIdentification() + " .. " + SP_AMINO_ACID->GetIdentification()
        );
      }

      // push back SP_AMINO_ACID
      m_Data.PushBack( SP_AMINO_ACID);
    }

    //! @brief check if a second sequence follows this by seqid and chain id (this precedes given sequence)
    //! @param SEQUENCE AASequence that should follow this
    //! @return true if neither seq is empty and this last aa precedes arguments first aa
    bool AASequence::DoesPrecede( const AASequence &SEQUENCE) const
    {
      return
           !m_Data.IsEmpty() && !SEQUENCE.m_Data.IsEmpty() // at least one aa in each data
        && m_Data.LastElement()->DoesPrecede( *SEQUENCE.m_Data.FirstElement()); // checks seqid and chainid
    }

    //! @brief appends SEQUENCE to the end of this sequence
    //! @param SEQUENCE AASequence to be appended to this sequence
    void AASequence::AppendSequence( const AASequence &SEQUENCE)
    {
      if( SEQUENCE.m_Data.IsEmpty())
      {
        return;
      }

      if( m_Data.IsEmpty())
      {
        m_ChainID = SEQUENCE.GetChainID();
        m_Data = SEQUENCE.m_Data.HardCopy();
      }
      else
      {
        // assure that chain ids do match
        BCL_Assert
        (
          DoesPrecede( SEQUENCE),
          "given sequence does not follow this: " + GetSequenceIdentification() +
          " .. " + SEQUENCE.GetSequenceIdentification()
        );
        // insert the aminoacids of the SEQUENCE
        m_Data.InsertElements( m_Data.GetSize(), SEQUENCE.m_Data.HardCopy());
      }
    }

    //! @brief prepends SEQUENCE to the end of this sequence
    //! @param SEQUENCE AASequence to be prepended to this sequence
    void AASequence::PrependSequence( const AASequence &SEQUENCE)
    {
      if( SEQUENCE.m_Data.IsEmpty())
      {
        return;
      }

      if( m_Data.IsEmpty())
      {
        m_ChainID = SEQUENCE.GetChainID();
        m_Data = SEQUENCE.m_Data.HardCopy();
      }
      else
      {
        // assure that chain ids do match
        BCL_Assert
        (
          SEQUENCE.DoesPrecede( *this),
          "given sequence does precede this: " + SEQUENCE.GetSequenceIdentification() +
          " .. " + GetSequenceIdentification()
        );

        // insert the amino acids of the SEQUENCE
        m_Data.InsertElements( 0, SEQUENCE.m_Data.HardCopy());
      }
    }

    //! @brief return sequence as string
    //! @return sequence as a string
    std::string AASequence::Sequence() const
    {
      // instantiate string for one letter code sequence
      std::string sequence;

      // iterate over all amino acids in the sequence
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()),
          aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // concatenate one letter code for this amino acid
        sequence += ( *aa_itr)->GetType()->GetOneLetterCode();
      }

      // end
      return sequence;
    }

    //! @brief let the aadata of the given sequence point to the aadata of this sequence - determined by pdbID
    //! @param AA_SEQUENCE sequence of amino acids which aas will be connected to the aadata of this sequence
    //! @param SET_ATOMS set to true if atoms should also be connected by this function
    void AASequence::ConnectAADataByPdbID
    (
      AASequence &AA_SEQUENCE,
      const bool &SET_ATOMS
    ) const
    {
      // if given sequence is empty, return
      if( AA_SEQUENCE.m_Data.IsEmpty())
      {
        // set the sequence's chain id to this chain id
        AA_SEQUENCE.m_ChainID = m_ChainID;

        return;
      }

      // instantiate iterator to the begin and end of this sequence and given AA_SEQUENCE
      AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
      AASequence::iterator arg_aa_itr( AA_SEQUENCE.Begin()), arg_aa_itr_end( AA_SEQUENCE.End());
      BCL_Assert( !m_Data.IsEmpty(), "m_Data.IsEmpty()");
      BCL_Assert( !AA_SEQUENCE.m_Data.IsEmpty(), "AA_SEQUENCE.IsEmpty()");

      // find first matching amino acid
      while
      (
        // iterate as long PdbID or ICode are not equal
        ( *aa_itr)->GetPdbID() != ( *arg_aa_itr)->GetPdbID() ||
        ( *aa_itr)->GetPdbICode() != ( *arg_aa_itr)->GetPdbICode()
      )
      {
        // move further in this sequence until the match is found
        ++aa_itr;

        // if the end has been reached
        if( aa_itr == aa_itr_end)
        {
          // issue warning and return
          BCL_MessageCrt( "could not find any matching starting amino acid");
          return;
        }
      }

      // iterate over the residues in the portion of the full sequence that is equal to the given sequence
      while( aa_itr != aa_itr_end && arg_aa_itr != arg_aa_itr_end)
      {
        // update the data of the residue in AASEQUENCE to corresponding residues data in this sequence
        ( *arg_aa_itr)->SetData( ( *aa_itr)->GetData());
        if( SET_ATOMS)
        {
          ( **arg_aa_itr).SetAtoms( ( **aa_itr).GetAtoms());
        }

        // move further in both sequences
        ++aa_itr;
        ++arg_aa_itr;

        // if end has been reached at either of the sequence
        if( aa_itr == aa_itr_end && arg_aa_itr != arg_aa_itr_end)
        {
          // exit
          BCL_Exit( "there are not enough amino acids in this sequence to be connected to the aas in the argument", -1);
        }
      }

      // set the sequence's chain id to this chain id
      AA_SEQUENCE.m_ChainID = m_ChainID;

      // end
      return;
    }

    //! @brief returns whether all coordinates for the given atom types are defined
    //! @param ATOM_TYPES Atom Types of interest
    bool AASequence::HasDefinedCoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES// = biol::GetAtomTypes().GetBackBoneAtomTypes()
    ) const
    {
      // iterate over the residues
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // if there are undefined coords
        if( !( *aa_itr)->HasDefinedCoordinates())
        {
          return false;
        }
      }

      // this point is reached only if all coordinates were defined, therefore return true
      return true;
    }

    //! @brief counts the number of residues that have defined coordinates for all the given atom types
    //! @param ATOM_TYPES Atom Types of interest
    size_t AASequence::CountDefinedAACoordinates
    (
      const storage::Set< AtomType> &ATOM_TYPES // = GetAtomTypes().GetBackBoneAtomTypes()
    ) const
    {
      size_t number_defined( 0);
      // iterate over the residues
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end; ++aa_itr
      )
      {
        // if there are undefined coords
        if( ( *aa_itr)->HasDefinedCoordinates( ATOM_TYPES))
        {
          ++number_defined;
        }
      }

      // this point is reached only if all coordinates were defined, therefore return true
      return number_defined;
    }

  ///////////////
  // operators //
  ///////////////

    //! operator = for assigning this sequence to AA_SEQUENCE
    //! @param AA_SEQUENCE AASequence to be copied
    //! @return this sequence after being assigned to AA_SEQUENCE
    AASequence &AASequence::operator =( const AASequence &AA_SEQUENCE)
    {
      // assign data members
      m_Data = AA_SEQUENCE.m_Data.HardCopy();
      m_ChainID = AA_SEQUENCE.m_ChainID;
      m_FastaHeader = AA_SEQUENCE.m_FastaHeader;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write fasta to given OSTREAM in blocks of BLOCK_SIZE
    //! @param OSTREAM output stream
    //! @param BLOCK_SIZE number of amino acids to be outputted in each block
    //! @return std::ostream which was written to
    std::ostream &AASequence::WriteFasta( std::ostream &OSTREAM, const size_t BLOCK_SIZE) const
    {
      // output the fasta header
      OSTREAM << s_FastaHeaderChar << m_FastaHeader << '\n';

      // initialize counter
      size_t counter( 0);

      // iterate over the amino acids in the sequence
      for
      (
        AASequence::const_iterator aa_itr( m_Data.Begin()), aa_itr_end( m_Data.End());
        aa_itr != aa_itr_end;
        ++aa_itr
      )
      {
        // output the one letter code of this amino acid
        ( *aa_itr)->WriteFasta( OSTREAM);

        // if BLOCK_SIZE has been reached
        if( !( ( counter + 1) % BLOCK_SIZE))
        {
          // output endline
          OSTREAM << '\n';
        }

        // increment counter
        ++counter;
      }

      // line break at end of output
      OSTREAM << '\n';

      // end
      return OSTREAM;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AASequence::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Data, ISTREAM);
      io::Serialize::Read( m_ChainID, ISTREAM);
      io::Serialize::Read( m_FastaHeader, ISTREAM);
      SetFastaHeader( m_FastaHeader);

      // check that chain id is correct and the amino acids are continuous
      BCL_Assert
      (
           m_Data.IsEmpty()
        || ( IsContinuous( m_Data) && ( m_Data.FirstElement()->GetChainID() == m_ChainID)),
        "read sequence is not continuous"
      );

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AASequence::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Data       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChainID    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FastaHeader, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if a given vector of amino acids is continuous
    //! @param DATA amino acids
    //! @return true if two consecutive amino acids follow each other by chain id and seq id
    bool AASequence::IsContinuous( const util::ShPtrVector< AABase> &DATA)
    {
      if( DATA.IsEmpty())
      {
        return true;
      }

      // iterate over amino acids
      const_iterator itr( DATA.Begin()), itr_follow( itr), itr_end( DATA.End());
      ++itr_follow;
      for( ; itr_follow != itr_end; ++itr, ++itr_follow)
      {
        if( !( *itr)->DoesPrecede( **itr_follow))
        {
          return false;
        }
      }

      return true;
    }

    //! checks whether two sequences have overlap
    bool DoOverlap( const AASequence &SEQUENCE_A, const AASequence &SEQUENCE_B)
    {
      //check for identical chain IDs
      if( SEQUENCE_A.GetChainID() != SEQUENCE_B.GetChainID())
      {
        BCL_MessageDbg( "different chains");
        return false;
      }

      // empty sequences do not overlaps
      if( SEQUENCE_A.GetSize() == 0 || SEQUENCE_B.GetSize() == 0)
      {
        return false;
      }

      //check if the range between beginning and ending seqid of each sequence are not overlapping
      return SEQUENCE_A.GetLastAA()->GetSeqID() >= SEQUENCE_B.GetFirstAA()->GetSeqID()
             && SEQUENCE_A.GetFirstAA()->GetSeqID() <= SEQUENCE_B.GetLastAA()->GetSeqID();
    }

    //! @brief calculates the distance in sequence between SEQUENCE_A and SEQUENCE_B
    //! @param SEQUENCE_A first sequence
    //! @param SEQUENCE_B second sequence
    //! @return the distance between SEQUENCE_A and SEQUENCE_B, Undefined if they overlap or are in different chains
    size_t CalculateSequenceDistance( const AASequence &SEQUENCE_A, const AASequence &SEQUENCE_B)
    {
      // if different chains or sequences overlap return undefined size_t
      if( SEQUENCE_A.GetChainID() != SEQUENCE_B.GetChainID() || DoOverlap( SEQUENCE_A, SEQUENCE_B))
      {
        return util::GetUndefinedSize_t();
      }
      // if SEQUENCE_A comes first
      if( SEQUENCE_A.GetLastAA()->GetSeqID() < SEQUENCE_B.GetFirstAA()->GetSeqID())
      {
        // return the distance from last residue of SEQUENCE_A to first residue of SEQUENCE_B
        return SEQUENCE_B.GetFirstAA()->GetSeqID() - SEQUENCE_A.GetLastAA()->GetSeqID() - 1;
      }
      // else if SEQUENCE_B comes first
      else
      {
        // return the distance from last residue of SEQUENCE_A to first residue of SEQUENCE_B
        return SEQUENCE_A.GetFirstAA()->GetSeqID() - SEQUENCE_B.GetLastAA()->GetSeqID() - 1;
      }
    }

    //! @brief gives the distance between the c atom of the n terminal aa and the n atom of the c terminal aa
    //! @param N_TERMINAL_AA the aa that provides the c in the peptide bond
    //! @param C_TERMINAL_AA the aa that provides the n in the peptide bond
    //! @return distance between the c atom of the n terminal aa and the n atom of the c terminal aa
    double GetPeptideBondLength( const AABase &N_TERMINAL_AA, const AABase &C_TERMINAL_AA)
    {
      if
      (
        !N_TERMINAL_AA.GetAtom( GetAtomTypes().C).GetCoordinates().IsDefined() ||
        !C_TERMINAL_AA.GetAtom( GetAtomTypes().N).GetCoordinates().IsDefined()
      )
      {
        return util::GetUndefinedDouble();
      }

      return linal::Distance
        (
          N_TERMINAL_AA.GetAtom( GetAtomTypes().C).GetCoordinates(),
          C_TERMINAL_AA.GetAtom( GetAtomTypes().N).GetCoordinates()
        );
    }

    //! @brief generates subsequences of specified radius around every residue and returns them
    //! @param SEQUENCE AASequence of interest
    //! @param WINDOW_RADIUS radius of the window
    //! @param UNDEFINED_AA undefined amino acid to be used for windows of border amino acids
    //! @return map of window sequences
    storage::Map< util::SiPtr< const AABase>, util::SiPtrVector< const AABase> > CreateWindowsFromAminoAcids
    (
      const AASequence &SEQUENCE,
      const size_t WINDOW_RADIUS,
      const AABase &UNDEFINED_AA
    )
    {
      // initialize window map
      storage::Map< util::SiPtr< const AABase>, util::SiPtrVector< const AABase> > window_map;

      // initialize SiPtr to undefined AA
      util::SiPtr< const AABase> undefined_aa( UNDEFINED_AA);

      // initialize window length
      const size_t window_length( 2 * WINDOW_RADIUS + 1);

      // make sure sequence is long enough for at least one window
      if( SEQUENCE.GetSize() <= window_length)
      {
        BCL_MessageStd
        (
          "The given sequence's length is smaller than the window radius " + util::Format()( SEQUENCE.GetSize()) +
            " < " + util::Format()( WINDOW_RADIUS)
         );
        return window_map;
      }

      // initialize amino acid counter
      size_t aa_count( 0);

      // now iterate over residues
      for
      (
        AASequence::const_iterator aa_itr( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr != aa_itr_end;
        ++aa_itr, ++aa_count
      )
      {
        // initialize the SiPtrVector and store the reference for this amino acid
        util::SiPtrVector< const AABase> &current_window( window_map[ **aa_itr]);

        // if there are not enough residues on the left of this residue fill with undefined residues
        if( aa_count < WINDOW_RADIUS)
        {
          // calculate how many undefined residues are different
          const size_t number_undefined( WINDOW_RADIUS - aa_count);

          // insert number_undefined counts of given UNDEFINED_AA
          current_window.Append
          (
            storage::Vector< util::SiPtr< const AABase> >( number_undefined, undefined_aa)
          );

          // append the rest of the the sequence
          current_window.Append
          (
            SEQUENCE.GetData().SubShPtrVector( 0, ( window_length - number_undefined))
          );
        }
        // if there are not enough residues on the right side of this residue fill with undefined residues
        else if( SEQUENCE.GetSize() <= WINDOW_RADIUS + aa_count)
        {
          // calculate how many undefined residues are different
          const size_t number_undefined( aa_count + WINDOW_RADIUS + 1 - SEQUENCE.GetSize());

          // append the rest of the the sequence
          current_window.Append
          (
            SEQUENCE.GetData().SubShPtrVector( aa_count - WINDOW_RADIUS, ( window_length - number_undefined))
          );

          // insert number_undefined counts of given UNDEFINED_AA
          current_window.Append
          (
            storage::Vector< util::SiPtr< const AABase> >( number_undefined, undefined_aa)
          );
        }
        // else it is safe to make a simple subsequence
        else
        {
          // append the subsequence
          current_window.Append
          (
            SEQUENCE.GetData().SubShPtrVector( aa_count - WINDOW_RADIUS, window_length)
          );
        }
      }

      // end
      return window_map;
    }

  } // namespace biol
} // namespace bcl
