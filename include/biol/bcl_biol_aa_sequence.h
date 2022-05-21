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

#ifndef BCL_BIOL_AA_SEQUENCE_H_
#define BCL_BIOL_AA_SEQUENCE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_base.h"
#include "align/bcl_align_sequence_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AASequence
    //! @brief This is a class for storing an amino acid sequence of a selected AAClass type
    //! @details This class allows representation of an peptide sequence. It has a internal ShPtrVector that can store
    //! any kind of AABase derived amino acid class, in addition to knowing its fasta header and the chain ID.
    //! In addition to data storage, it provides many convenience functions. It should be noted that the copy structor
    //! of this class actually hard-copies the the amino acid storage.
    //!
    //! @see @link example_biol_aa_sequence.cpp @endlink
    //! @author meilerj, staritrd, woetzen, karakam
    //! @date 21.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASequence :
      public virtual coord::MovableInterface,
      public align::SequenceInterface< AABase>
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      util::ShPtrVector< AABase> m_Data;        //!< amino acid data
      char                       m_ChainID;     //!< chain identifier
      std::string                m_FastaHeader; //!< fasta header line

    public:

      //! @brief typedef to reverse iterate over non-const members of a sequence
      typedef util::ShPtrVector< AABase>::reverse_iterator       reverse_iterator;
      //! @brief typedef to reverse iterate over const members of a sequence
      typedef util::ShPtrVector< AABase>::const_reverse_iterator const_reverse_iterator;

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! default chain id for chains
      static const char s_DefaultChainID = 'A';

      //! starting character for fasta header
      static const char s_FastaHeaderChar = '>';

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AASequence();

      //! @brief construct from AACLASS, LENGTH, CHAIN_ID and FASTA_HEADER
      //! @param AACLASS AAClass type for m_Data
      //! @param LENGTH number of residues in this sequence
      //! @param CHAIN_ID chain id of sequence
      //! @param FASTA_HEADER fasta header
      AASequence
      (
        const AAClass &AACLASS,
        const size_t LENGTH,
        const char CHAIN_ID,
        const std::string &FASTA_HEADER = GetDefaultFastaHeader()
      );

      // TODO: can not ensure the people pass DATA from a sequence with AAClass of same type with this->m_AAClass
      //! @brief construct from ShPtrVector of AABases DATA, CHAIN_ID and FASTA_HEADER
      //! @param DATA ShPtrVector of AABases
      //! @param CHAIN_ID chain id of sequence
      //! @param FASTA_HEADER fasta header
      AASequence
      (
        const util::ShPtrVector< AABase> &DATA,
        const char CHAIN_ID = s_DefaultChainID,
        const std::string &FASTA_HEADER = GetDefaultFastaHeader()
      );

      //! @brief copy constructor from another AA_SEQUENCE
      //! @param AA_SEQUENCE AASequence to be copied
      AASequence( const AASequence &AA_SEQUENCE);

      //! @brief virtual copy constructor
      virtual AASequence *Clone() const;

      //! @brief virtual hard copy all amino acids and their data
      //! @return aa sequence with independent hard copied AADatas
      virtual AASequence *HardCopy() const;

      //! @brief construct a new SequenceInterface from an AASequence interface implementation, sequence id and members
      //! @param ID sequence identifier
      //! @param MEMBERS sequence members
      //! @param CHAIN_ID chain id used to build the sequence
      //! @return pointer to a SequenceInterface
      virtual align::SequenceInterface< AABase> *Construct
      (
        const std::string &ID,
        const std::string &MEMBERS,
        const char CHAIN_ID = 0
      ) const;

      //! @brief destructor
      virtual ~AASequence();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      virtual const std::string &GetClassIdentifier() const;

      //! @brief Get default fasta header
      //! @return default fasta header
      static const std::string &GetDefaultFastaHeader();

      //! @brief get identification of this AASequence
      //! @return string with identification for this AASequence
      std::string GetSequenceIdentification() const;

      //! @brief return reference to const Sequence information
      //! @return reference to const Sequence information
      virtual const util::ShPtrVector< AABase> &GetData() const
      {
        return m_Data;
      }

      //! @brief return reference to non-const Sequence information
      //! @return reference to non-const Sequence information
      virtual util::ShPtrVector< AABase> &GetData()
      {
        return m_Data;
      }

      //! @brief returns the amino acids with the given type
      //! @param AA_TYPE type of the amino acid to return
      //! @return shared pointers to the amino acids of the the given type
      const util::ShPtrList< AABase> GetData( const AAType &AA_TYPE) const;

    protected:

      //! @brief Return non-const reference to all AA
      //! @return non-const reference to all AA
      virtual util::ShPtrVector< AABase> &GetAminoAcids()
      {
        return m_Data;
      }

    public:

      //! @brief set sequence information to the supplied AA_SEQUENCE_DATA
      //! @param AA_SEQUENCE_DATA ShPtrVector of AABases to be copied
      virtual void SetData( const util::ShPtrVector< AABase> &AA_SEQUENCE_DATA)
      {
        m_Data = AA_SEQUENCE_DATA;
      }

      //! @brief returns the number of amino acids contained in sequence
      //! @return the number of amino acids contained in sequence
      virtual size_t GetSize() const
      {
        return m_Data.GetSize();
      }

      //! @brief returns a vector of members
      //! @return a vector of SiPtrs of type t_Member contained in the sequence
      virtual util::SiPtrVector< const AABase> GetMembers() const
      {
        return util::SiPtrVector< const AABase>( m_Data);
      }

      //! @brief returns a vector of members
      //! @return a vector of SiPtrs of type t_Member contained in the sequence
      virtual util::SiPtrVector< AABase> GetMembers()
      {
        return util::SiPtrVector< AABase>( m_Data);
      }

      //! @brief returns a SiPtr to the first element
      //! @return SiPtr to the first AABase
      virtual util::SiPtr< const AABase> GetFirstMember() const
      {
        return util::SiPtr< const AABase>( m_Data.FirstElement());
      }

      //! @brief returns a SiPtr to the last element
      //! @return SiPtr to the last AABase
      virtual util::SiPtr< const AABase> GetLastMember() const
      {
        return util::SiPtr< const AABase>( m_Data.LastElement());
      }

      //! @brief add a AABase to the sequence with type based on ONE_LETTER_CODE
      //! @param ONE_LETTER_CODE char encoding the type of the AA
      virtual void AddMember( const char &ONE_LETTER_CODE);

      //! @brief initializes the sequence with a new type of AAClass and fills with LENGTH number of residues
      //! deletes all existing aminoacids and builds up a new sequence with empty AABase interfaces
      //! @param AA_CLASS AAClass that identifies the amino acid type to be used
      //! @param LENGTH number of residues to be inserted
      void Initialize( const AAClass &AA_CLASS, const size_t LENGTH);

      //! @brief get all atoms of all residues in this sequence
      //! @return SiPtrVector of Atoms of all residues in this sequence
      util::SiPtrVector< const Atom> GetAtoms() const;

      //! @brief get all atoms for the specified atom types for all residues in the sequence
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return atoms of types specified in ATOM_TYPES  for all residues in the sequence
      util::SiPtrVector< const Atom> GetAtoms
      (
        const storage::Set< AtomType> &ATOM_TYPES
      ) const;

      //! @brief get all atom coordinates for all residues in the sequence
      //! @return all atom coordinates for all residues in the sequence
      util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates() const;

      //! @brief get all atom coordinates for the specified atom types for all residues in the sequence
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return coordinates of atoms of types specified in ATOM_TYPES  for all residues in the sequence
      util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates
      (
        const storage::Set< AtomType> &ATOM_TYPES
      ) const;

      //! @brief return iterator to Begin of Sequence
      //! @return iterator to Begin of Sequence
      iterator Begin()
      {
        return m_Data.Begin();
      }

      //! @brief return const_iterator to Begin of Sequence
      //! @return const_iterator to Begin of Sequence
      const_iterator Begin() const
      {
        return m_Data.Begin();
      }

      //! @brief return iterator to End of Sequence
      //! @return iterator to End of Sequence
      iterator End()
      {
        return m_Data.End();
      }

      //! @brief return const_iterator to End of Sequence
      //! @return const_iterator to End of Sequence
      const_iterator End() const
      {
        return m_Data.End();
      }

      //! @brief return iterator to reverse begin of Sequence
      //! @return iterator to reverse begin of Sequence
      reverse_iterator ReverseBegin()
      {
        return m_Data.ReverseBegin();
      }

      //! @brief return const_iterator to reverse begin of Sequence
      //! @return const_iterator to reverse begin of Sequence
      const_reverse_iterator ReverseBegin() const
      {
        return m_Data.ReverseBegin();
      }

      //! @brief return iterator to reverse end of Sequence
      //! @return iterator to reverse end of Sequence
      reverse_iterator ReverseEnd()
      {
        return m_Data.ReverseEnd();
      }

      //! @brief return const_iterator to reverse end of Sequence
      //! @return const_iterator to reverse end of Sequence
      const_reverse_iterator ReverseEnd() const
      {
        return m_Data.ReverseEnd();
      }

      //! @brief set the amino acid residing at index INDEX, to given AMINO_ACID
      //! @param INDEX index of the residue of interest in the sequence
      //! @param AMINO_ACID amino acid to be copied
      void SetAA( const size_t INDEX, const AABase &AMINO_ACID)
      {
        *m_Data( INDEX) = AMINO_ACID;
      }

      //! @brief return pointer to residue at index INDEX
      //! @param INDEX index of the residue of interest in the sequence
      const util::ShPtr< AABase> &GetAA( const size_t INDEX) const
      {
        return m_Data( INDEX);
      }

      //! @brief find and return the iterator to the amino acid with the given sequence id
      //! @param SEQ_ID sequence id of the amino acid of interest
      //! @return the iterator to the amino acid with the given sequence id
      AASequence::iterator FindAABySeqID( const int SEQ_ID);

      //! @brief find and return the iterator to the amino acid with the given sequence id
      //! @param SEQ_ID sequence id of the amino acid of interest
      //! @return the iterator to the amino acid with the given sequence id
      AASequence::const_iterator FindAABySeqID( const int SEQ_ID) const;

      //! @brief calculates the phi and psi for the amino acid with the given seq id
      //! @param SEQ_ID sequence id of the amino acid of interest
      //! @return pair of phi and psi values for the specified amino acid
      storage::VectorND< 2, double> CalculatePhiPsi( const int SEQ_ID) const;

      //! @brief return pointer to first amino acid in the sequence
      //! @return pointer to first amino acid in the sequence
      const util::ShPtr< AABase> &GetFirstAA() const
      {
        return m_Data.FirstElement();
      }

      //! @brief return pointer to last amino acid in the sequence
      //! @return pointer to last amino acid in the sequence
      const util::ShPtr< AABase> &GetLastAA() const
      {
        return m_Data.LastElement();
      }

      //! @brief get chain id of the sequence
      //! @return chain id of the sequence
      const char &GetChainID() const
      {
        return m_ChainID;
      }

      //! @brief set chain id of the sequence and all residues in the sequence to given CHAIN_ID
      //! @param CHAIN_ID new chain id of the sequence
      void SetChainID( const char CHAIN_ID);

      //! @brief return sequence identifier
      //! @return sequence identifier of the sequence
      std::string GetSequenceId() const
      {
        return m_FastaHeader;
      }

      //! @brief return fasta header of the sequence
      //! @return fasta header of the sequence
      const std::string &GetFastaHeader() const
      {
        return m_FastaHeader;
      }

      //! @brief change the fasta header of the sequence to given FASTA_HEADER
      //! @param FASTA_HEADER new fasta header of the sequence
      void SetFastaHeader( const std::string &FASTA_HEADER);

    ////////////////
    // operations //
    ////////////////

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION translation vector to be applied
      virtual void Translate( const linal::Vector3D &TRANSLATION);

      //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
      //! @param TRANSFORMATIONMATRIX3D TransformationMatrix3D to be applied
      virtual void Transform( const math::TransformationMatrix3D &TRANSFORMATIONMATRIX3D);

      //! @brief rotate the object by a given ROTATION_MATRIX_3D
      //! @param ROTATIONMATRIX3D RotationMatrix3D to be applied
      virtual void Rotate( const math::RotationMatrix3D &ROTATIONMATRIX3D);

      //! @brief returns the geometric center of the object
      //! @return geometric center of the object
      virtual linal::Vector3D GetCenter() const
      {
        return coord::CenterOfMass( GetAtomCoordinates( storage::Set< AtomType>( GetAtomTypes().CA)), true);
      }

      //! @brief return util::ShPtrVector of AASequences formed after chopping in lengths of SIZE with one gap in between
      //! @param SIZE length of each resultant aasequence
      //! @return ShPtrVector of AASequence formed after chopping
      util::ShPtrVector< AASequence> ChopSequence( const size_t SIZE) const;

      //! @brief return subsequence starting from position INDEX, and has length LENGTH
      //! @param INDEX index of the first residue of the subsequence
      //! @param LENGTH length of the subsequence
      //! @return subsequence starting from position INDEX, and has length LENGTH
      AASequence SubSequence( const size_t &INDEX, const size_t LENGTH) const
      {
        // update the fasta header to indicate first and last indices of the subsequence
        std::ostringstream header
        (
          m_FastaHeader + "[" + util::Format()( INDEX + 1) + ".." + util::Format()( INDEX + LENGTH) + "]"
        );

        // construct a sub ShPtrVector of m_Data starting at INDEX and of length LENGTH and return it
        return AASequence( m_Data.SubShPtrVector( INDEX, LENGTH), m_ChainID, header.str());
      }

      //! @brief clips number of amino acids from the beginning and the end of the sequence
      //! @param NUMBER_AMINO_ACIDS number of amino acids to be clipped from each end
      //! @return this sequence after clipping NUMBER_AMINO_ACIDS from each end
      AASequence &ClipEnds( const size_t NUMBER_AMINO_ACIDS);

      //! @brief appends new amino acid to the end of the sequence
      //! @param AMINO_ACID amino acid to be appended to the sequence
      void PushBack( const AABase &AMINO_ACID)
      {
        PushBack( util::ShPtr< AABase>( AMINO_ACID.Clone()));
      }

      //! @brief prepends new amino acid to the beginning of the sequence
      //! @param AMINO_ACID amino acid to be prepended to the sequence
      void PushFront( const AABase &AMINO_ACID);

      //! @brief appends util::ShPtr to amino acid to the end of the sequence
      //! @param SP_AMINO_ACID ShPtr to amino acid to be appended to the end of the sequence
      void PushBack( const util::ShPtr< AABase> &SP_AMINO_ACID);

      //! @brief check if a second sequence follows this by seqid and chain id (this precedes given sequence)
      //! @param SEQUENCE AASequence that should follow this
      //! @return true if neither seq is empty and this last aa precedes arguments first aa
      bool DoesPrecede( const AASequence &SEQUENCE) const;

      //! @brief appends SEQUENCE to the end of this sequence
      //! @param SEQUENCE AASequence to be appended to this sequence
      void AppendSequence( const AASequence &SEQUENCE);

      //! @brief prepends SEQUENCE to the end of this sequence
      //! @param SEQUENCE AASequence to be prepended to this sequence
      void PrependSequence( const AASequence &SEQUENCE);

      //! @brief removes last amino acid
      void PopBack()
      {
        m_Data.PopBack();
      }

      //! @brief reset complete sequence
      void Reset()
      {
        m_Data.Reset();
      }

      //! @brief return sequence as string
      //! @return sequence as a string
      std::string Sequence() const;

      //! @brief let the aadata of the given sequence point to the aadata of this sequence - determined by pdbID
      //! @param AA_SEQUENCE sequence of amino acids which aas will be connected to the aadata of this sequence
      //! @param SET_ATOMS set to true if atoms should also be connected by this function
      void ConnectAADataByPdbID( AASequence &AA_SEQUENCE, const bool &SET_ATOMS = true) const;

      //! @brief returns whether all coordinates for the given atom types are defined
      //! @param ATOM_TYPES Atom Types of interest
      bool HasDefinedCoordinates
      (
        const storage::Set< AtomType> &ATOM_TYPES = GetAtomTypes().GetBackBoneAtomTypes()
      ) const;

      //! @brief counts the number of residues that have defined coordinates for all the given atom types
      //! @param ATOM_TYPES Atom Types of interest
      size_t CountDefinedAACoordinates
      (
        const storage::Set< AtomType> &ATOM_TYPES = GetAtomTypes().GetBackBoneAtomTypes()
      ) const;

    ///////////////
    // operators //
    ///////////////

      //! operator = for assigning this sequence to AA_SEQUENCE
      //! @param AA_SEQUENCE AASequence to be copied
      //! @return this sequence after being assigned to AA_SEQUENCE
      virtual AASequence &operator =( const AASequence &AA_SEQUENCE);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write fasta to given OSTREAM in blocks of BLOCK_SIZE
      //! @param OSTREAM output stream
      //! @param BLOCK_SIZE number of amino acids to be outputted in each block
      //! @return std::ostream which was written to
      std::ostream &WriteFasta( std::ostream &OSTREAM, const size_t BLOCK_SIZE = 50) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief check if a given vector of amino acids is continuous
      //! @param DATA amino acids
      //! @return true if two consecutive amino acids follow each other by chain id and seq id
      static bool IsContinuous( const util::ShPtrVector< AABase> &DATA);

    }; // class AASequence

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief operator == checking if SEQUENCE_A is equal to SEQUENCE_B
    //! @param SEQUENCE_A first AASequence
    //! @param SEQUENCE_B second AASequence
    //! @return whether SEQUENCE_A and SEQUENCE_B have same sequence, chain id and fasta header
    inline bool operator ==( const AASequence &SEQUENCE_A, const AASequence &SEQUENCE_B)
    {
      return
      (
           SEQUENCE_A.Sequence()       == SEQUENCE_B.Sequence()
        && SEQUENCE_A.GetChainID()     == SEQUENCE_B.GetChainID()
        && SEQUENCE_A.GetFastaHeader() == SEQUENCE_B.GetFastaHeader()
      );
    }

    //! @brief DoOverlap checks if two amino acid sequences overlap
    //! chain id needs to be the same and the first or the last amino acid of sequence b must be within the first and
    //! last sequence ids of sequence a
    //! @param SEQUENCE_A first sequence which will be checked against SEQUENCE_B
    //! @param SEQUENCE_B second sequence which will be checked against SEQUENCE_A
    //! @return returns a bool true if they overlap, otherwise false
    BCL_API bool DoOverlap( const AASequence &SEQUENCE_A, const AASequence &SEQUENCE_B);

    //! @brief calculates the distance in sequence between SEQUENCE_A and SEQUENCE_B
    //! @param SEQUENCE_A first sequence
    //! @param SEQUENCE_B second sequence
    //! @return the distance between SEQUENCE_A and SEQUENCE_B, Undefined if they overlap or are in different chains
    BCL_API size_t CalculateSequenceDistance( const AASequence &SEQUENCE_A, const AASequence &SEQUENCE_B);

    //! @brief gives the distance between the c atom of the n terminal aa and the n atom of the c terminal aa
    //! @param N_TERMINAL_AA the aa that provides the c in the peptide bond
    //! @param C_TERMINAL_AA the aa that provides the n in the peptide bond
    //! @return distance between the c atom of the n terminal aa and the n atom of the c terminal aa
    BCL_API double GetPeptideBondLength( const AABase &N_TERMINAL_AA, const AABase &C_TERMINAL_AA);

    //! @brief generates subsequences of specified radius around every residue and returns them
    //! @param SEQUENCE AASequence of interest
    //! @param WINDOW_RADIUS radius of the window
    //! @param UNDEFINED_AA undefined amino acid to be used for windows of border amino acids
    //! @return map of window sequences
    BCL_API storage::Map< util::SiPtr< const AABase>, util::SiPtrVector< const AABase> > CreateWindowsFromAminoAcids
    (
      const AASequence &SEQUENCE,
      const size_t WINDOW_RADIUS,
      const AABase &UNDEFINED_AA
    );

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_AA_SEQUENCE_H_
