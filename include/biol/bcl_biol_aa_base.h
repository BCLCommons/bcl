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

#ifndef BCL_BIOL_AA_BASE_H_
#define BCL_BIOL_AA_BASE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_data.h"
#include "bcl_biol_aa_types.h"
#include "coord/bcl_coord_movable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AABase
    //! @brief Base class to be used for deriving different AAClass types from: AA, AACaCb, AABackBone, etc..
    //! @details AABase provides a unique interface for all AA classes preventing the usage of templates in upper classes
    //! such as AASequence. Each AABase derived class should be also added to AAClasses enumerator.
    //!
    //! @see @link example_biol_aa_base.cpp @endlink
    //! @author woetzen, karakam
    //! @date 11/25/2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AABase :
      public coord::MovableInterface
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< AAData> m_Data; //!< amino acid data

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AABase();

      //! @brief construct AABase from util::ShPtr to AAData
      //! @brief SP_AA_DATA ShPTr to AAData to be copied
      AABase( const util::ShPtr< AAData> &SP_AA_DATA);

      //! @brief copy constructor - makes just a soft copy of m_Data
      //! @param AA_BASE_RHS AABase to be copied
      AABase( const AABase &AA_BASE_RHS);

      //! @brief virtual empty constructor with AAData
      //! @param SP_AA_DATA AAData object with information for the new AA
      //! this function is designed to be used in cases where AAClass is used for production multiple AA's
      //! which do not share one common AAData
      virtual AABase *Empty( const util::ShPtr< AAData> &SP_AA_DATA) const = 0;

      //! @brief destructor
      virtual ~AABase();

    /////////////////
    // data access //
    /////////////////

      //! @brief get the pointer to AAData
      //! @return SharedPointer to AAData
      const util::ShPtr< AAData> &GetData() const
      {
        return m_Data;
      }

      //! @brief set the pointer to AAData
      //! @param SP_AA_DATA ShPtr to AAData to be copied
      void SetData( const util::ShPtr< AAData> &SP_AA_DATA)
      {
        m_Data = SP_AA_DATA;
      }

      //! @brief return AAType
      //! @return AAType
      const AAType &GetType() const
      {
        return m_Data->m_Type;
      }

      //! @brief return SeqID
      //! @return SeqID
      const int &GetSeqID() const
      {
        return m_Data->m_SeqID;
      }

      //! @brief return PdbID
      //! @return PdbID
      int GetPdbID() const
      {
        return m_Data->m_PdbID;
      }

      //! @brief return pdb insertion code
      //! @return pdb insertion code
      char GetPdbICode() const
      {
        return m_Data->m_PdbICode;
      }

      //! @brief return ChainID
      //! @return chain id
      char GetChainID() const
      {
        return m_Data->m_ChainID;
      }

      //! @brief get number of atoms
      //! @return number of atoms
      virtual size_t GetNumberOfAtoms() const = 0;

      //! @brief get types of atoms
      //! @return set of AtomTypes
      virtual const storage::Set< AtomType> &GetTypesOfAtoms() const = 0;

      //! @brief get all atoms
      //! @return SiPtrVector of Atoms
      virtual const util::SiPtrVector< const Atom> &GetAtoms() const = 0;

      //! @brief get all atoms of specified types
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return SiPtrVector of Atoms with specified type
      virtual util::SiPtrVector< const Atom> GetAtoms( const storage::Set< AtomType> &ATOM_TYPES) const = 0;

      //! @brief get all atoms of the specified element type
      //! @param ELEMENT_TYPE element types of interest
      //! @return all atoms of the specified element type
      virtual util::SiPtrVector< const Atom> GetAtoms( const chemistry::ElementType &ELEMENT_TYPE) const;

      //! @brief get the specified atom
      //! @brief ATOM_TYPE AtomType of interest
      //! @return atom with type ATOM_TYPE
      virtual const Atom &GetAtom( const AtomType &ATOM_TYPE) const = 0;

      //! @brief get CA atom
      //! @return CA atom
      virtual const Atom &GetCA() const = 0;

      //! @brief get CB atom
      //! @return CB atom
      virtual const Atom &GetFirstSidechainAtom() const = 0;

      //! @brief set the aa type
      //! @return AAType
      void SetType( const AAType &TYPE)
      {
        m_Data = m_Data.HardCopy();
        m_Data->m_Type = TYPE;
      }

      //! @brief set the specified atom
      //! @brief ATOM Atom to be set
      virtual void SetAtom( const Atom &ATOM) = 0;

      //! @brief set all atoms
      //! @param ATOMS SiPtrVector of Atoms to be set
      virtual void SetAtoms( const util::SiPtrVector< const Atom> &ATOMS) = 0;

      //! @brief get AAClass
      //! @return AAClass
      virtual const AAClass &GetAAClass() const = 0;

      //! @brief calculate Omega backbone angle
      //! @param PREVIOUS_CA previous CA atom
      //! @param PREVIOUS_C previous carbon atom
      //! @return omega angle
      virtual double CalculateOmega( const Atom &PREVIOUS_CA, const Atom &PREVIOUS_C) const = 0;

      //! @brief calculate phi backbone angle if hydrogen is part of this aa
      //! @return phi angle
      virtual double Phi() const = 0;

      //! @brief calculate phi backbone angle
      //! @param PREVIOUS_C previous carbon atom
      //! @return phi angle
      virtual double CalculatePhi( const Atom &PREVIOUS_C) const = 0;

      //! @brief calculate psi backbone angle if oxygen is part of this aa
      //! @return psi angle
      virtual double Psi() const = 0;

      //! @brief calculate psi backbone angle
      //! @param FOLLOWING_N following nitrogen atom
      //! @return psi angle
      virtual double CalculatePsi( const Atom &FOLLOWING_N) const = 0;

      //! @brief calculate phi and psi
      //! @param PREVIOUS_C previous carbon atom
      //! @param FOLLOWING_N followin nitrogen atom
      //! @return phi and psi angles
      const storage::Pair< double, double> CalculatePhiPsi( const Atom &PREVIOUS_C, const Atom &FOLLOWING_N) const;

      //! @brief create a locator for that aa
      //! @param USE_PDB_ID
      //! @return locator aa to locate an amino acid in a protein model
      assemble::LocatorAA Locator( const bool USE_PDB_ID) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get all atom coordinates
      //! @return all atom coordinates
      virtual util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates() const = 0;

      //! @brief get all atom coordinates for the specified atom types
      //! @param ATOM_TYPES AtomTypes of interest
      //! @return coordinates of atoms of types specified in ATOM_TYPES
      virtual util::SiPtrVector< const linal::Vector3D> GetAtomCoordinates
      (
        const storage::Set< AtomType> &ATOM_TYPES
      ) const = 0;

      //! @brief calculate the center of mass for all side chain atoms
      //! @return center of mass (not weighted by mass) of all side chain atoms - undefined if there are no atoms
      linal::Vector3D CalculateCenterOfMassOfSideChain() const;

      //! @brief translate the object along a given TRANSLATION vector
      //! @param TRANSLATION translation vector to be applied
      virtual void Translate( const linal::Vector3D &TRANSLATION) = 0;

      //! @brief transform the object by a given TRANSFORMATION_MATRIX_3D
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      virtual void Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D) = 0;

      //! @brief rotate the object by a given TRANSFORMATION_MATRIX_3D
      //! @param TRANSFORMATION_MATRIX_3D RotationMatrix3D to be applied
      virtual void Rotate( const math::RotationMatrix3D &TRANSFORMATION_MATRIX_3D) = 0;

      //! @brief returns the geometric center of the object
      //! @return geometric center of the object
      virtual linal::Vector3D GetCenter() const = 0;

      //! @brief set the aminoacid to the ideal conformation according to the SS_TYPE with given TRANSFORMATION_MATRIX_3D
      //! @param SS_TYPE SSType this AABase derived class is in
      //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
      virtual void SetToIdealConformation
      (
        const SSType &SS_TYPE,
        const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D
      ) = 0;

      //! @brief returns whether all coordinates for the given atom types are defined
      //! @param ATOM_TYPES Atom Types of interest
      bool HasDefinedCoordinates
      (
        const storage::Set< AtomType> &ATOM_TYPES = GetAtomTypes().GetBackBoneAtomTypes()
      ) const;

      //! @brief check if this aa precedes a given amino acid
      //! @param AMINO_ACID second amino acid, that could precede this
      //! @return true if both amino acids have same chain id and if (other seq id) - (this seq_id) == 1
      bool DoesPrecede( const AABase &AMINO_ACID) const;

      //! @brief gives dihedral angles of the side chain - starting closest to backbone and moving out along side chain
      //! @return vector holding the dihedral angles - first element is dihedral angle closest to backbone
      Rotamer CalculateSideChainDihedralAngles() const;

    ////////////////////////////////////////
    // data access - SSPrediction related //
    ////////////////////////////////////////

      //! @brief return map of secondary structure predictions stored for various method
      //! @return map of secondary structure predictions stored for various method
      const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &GetSSPredictions() const;

      //! @brief return secondary structure predictions stored for given SS_METHOD
      //! @brief METHOD sspred::Method of interest
      //! @return secondary structure predictions stored for SS_METHOD
      util::SiPtr< const sspred::MethodInterface> GetSSPrediction
      (
        const sspred::Method &SS_METHOD
      ) const;

      //! @brief provide sspredictions from different methods in PREDICTIONS_MAP to overwrite existing ones
      //! @param PREDICTIONS_MAP Map of predictions with different methods
      void SetSSPredictions
      (
        const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &PREDICTIONS_MAP
      );

      //! @brief Sets the predictions for the given SS_METHOD
      //! @param SS_METHOD SSMethod of interest
      //! @param PREDICTION MethodInterface derived class that stores the predictions
      //! For thread safety, PREDICTION should not be a shared pointer because otherwise classes would be free to use
      //! static global shared pointers, which would then require either 1. thread-safe shared pointers or 2. muteces
      //! around every AABase copy constructor and destructor.
      void SetSSPrediction
      (
        const sspred::Method &SS_METHOD,
        const sspred::MethodInterface &PREDICTION
      );

      //! @brief remove structure-based secondary-structure and TM methods from the map for this AA
      void RemoveStructureBasedSSTMInfo();

      //! @brief get the exposure prediction
      //! @return exposure prediction
      const double &GetExposurePrediction() const;

      //! @brief set the exposure prediction
      //! @param EXPOSURE exposure to set
      void SetExposurePrediction( const double &EXPOSURE);

    /////////////////////////////////
    // data access - Blast related //
    /////////////////////////////////

      //! @brief returns BlastProfile Object
      //! @param ISTREAM input stream to read the blast profile from
      void ReadBlastProfile( std::istream &ISTREAM);

      //! @brief returns const BlastProfile Object
      //! @return const BlastProfile Object
      const BlastProfile &GetBlastProfile() const;

      //! @brief returns const BlastProfile Object
      //! @return const BlastProfile Object
      const util::ShPtr< BlastProfile> &GetBlastProfilePtr() const;

      //! @brief set BlastProfile and Probabilities accordingly to given BLAST_PROFILE object
      //! @param BLAST_PROFILE blast profile object to be set to
      void SetBlastProfile( const BlastProfile &BLAST_PROFILE);

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator for AABase class
      //! @param AA_BASE_RHS AABase to be assigned from
      //! @return this AABase after being updated to AA_BASE_RHS
      virtual AABase &operator =( const AABase &AA_BASE_RHS);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read fasta from given ISTREAM with given SEQ_ID
      //! @param ISTREAM input stream
      //! @param SEQ_ID seq id for this AABase
      //! @return std::istream from which was read
      std::istream &ReadFasta( std::istream &ISTREAM, const size_t SEQ_ID = 1);

      //! @brief write fasta to provided OSTREAM
      //! @param OSTREAM output stream
      //! @return std::ostream to which was written
      std::ostream &WriteFasta( std::ostream &OSTREAM) const;

      //! @brief get Identification of this amino acid
      //! @return string with identification
      std::string GetIdentification() const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief check if two amino acids are in range to be peptide bonded
      //! @param AA_LEFT left in sequence
      //! @param AA_RIGHT right in sequence
      //! @param TEST_OMEGA check the angle also
      //! @return true if amino acids are close enough for a peptide bond
      static bool AreAminoAcidsPeptideBonded
      (
        const AABase &AA_LEFT,
        const AABase &AA_RIGHT,
        const bool TEST_OMEGA
      );

      //! @brief calculate peptide bond length and angle
      //! @param AA_LEFT left in sequence
      //! @param AA_RIGHT right in sequence
      //! @return vector of peptide bond length (C->N) [0, +inf] and peptide bond angle (ca->c->n->ca) [-pi, pi]
      //!         if any coord is undefined, the result will be undefined
      static storage::VectorND< 2, double> PeptideBondLengthAndAngle
      (
        const AABase &AA_LEFT,
        const AABase &AA_RIGHT
      );

      //! @brief construct single letter code and seqid
      //! @param ONE_LETTER_CODE aa one letter code
      //! @param SEQ_ID sequence id
      //! @return pointer to new amino acid
      static AABase *Construct( const char ONE_LETTER_CODE, const int SEQ_ID);

    }; // class AABase

  ///////////////
  // operators //
  ///////////////

    //! @brief compare to amino acids by identity
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return whether AMINO_ACID_A and AMINO_ACID_B have the same identity
    inline bool operator ==( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B)
    {
      return AMINO_ACID_A.GetType() == AMINO_ACID_B.GetType();
    }

    //! @brief Calculate the CB distance between two aminoacids by calling AtomDistance function
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return CB distance between AMINO_ACID_A and AMINO_ACID_B
    BCL_API double FirstSidechainAtomDistance( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B);

    //! @brief Compute the nearest atom separation between two AAs (nearest of atom distances - VdW radii of atoms involved)
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @param IGNORE_BB_BB ignore backbone-backbone interaction
    //! @return the nearest atom separation between two AAs (nearest of atom distances - VdW radii of atoms involved)
    BCL_API double NearestAtomVdWSphereSeparation( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B, const bool &IGNORE_BB_BB);

    //! @brief Calculate the sequence separation between two amino acids
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return sequence separation between AMINO_ACID_A and AMINO_ACID_B, will be undefined if they are from the same
    //!          chain or if the sequence ids are the same
    BCL_API size_t SequenceSeparation( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B);

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_AA_BASE_H_

