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
#include "biol/bcl_biol_aa_base.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "align/bcl_align.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_classes.h"
#include "biol/bcl_biol_rotamer.h"
#include "chemistry/bcl_chemistry_atom_types.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AABase::AABase() :
      m_Data( util::ShPtr< AAData>( new AAData()))
    {
    }

    //! @brief construct AABase from util::ShPtr to AAData
    //! @brief SP_AA_DATA ShPTr to AAData to be copied
    AABase::AABase( const util::ShPtr< AAData> &SP_AA_DATA) :
      m_Data( SP_AA_DATA)
    {
    }

    //! @brief copy constructor - makes just a soft copy of m_Data
    //! @param AA_BASE_RHS AABase to be copied
    AABase::AABase( const AABase &AA_BASE_RHS) :
      m_Data( AA_BASE_RHS.m_Data)
    {
    }

    //! destructor
    AABase::~AABase()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get all atoms of the specified element type
    //! @param ELEMENT_TYPE element types of interest
    //! @return all atoms of the specified element type
    util::SiPtrVector< const Atom> AABase::GetAtoms( const chemistry::ElementType &ELEMENT_TYPE) const
    {
      // get all the atoms
      const util::SiPtrVector< const Atom> all_atoms( GetAtoms());
      util::SiPtrVector< const Atom> found_atoms;

      // iterate over the atoms
      for
      (
        util::SiPtrVector< const Atom>::const_iterator atom_itr( all_atoms.Begin()), atom_itr_end( all_atoms.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // if the atom is the right element type
        if( ( *atom_itr)->GetType()->GetElementType() == ELEMENT_TYPE)
        {
          // add it to the vector
          found_atoms.PushBack( *atom_itr);
        }
      }

      // end
      return found_atoms;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the center of mass for all side chain atoms
    //! @return center of mass (not weighted by mass) of all side chain atoms - undefined if there are no atoms
    linal::Vector3D AABase::CalculateCenterOfMassOfSideChain() const
    {
       // if the aa passed wasn't an AAComplete, then just return the center of the first side chain atom
      if( GetAAClass() != GetAAClasses().e_AAComplete)
      {
        // center of the first side chain atom is center of mass for the sidechain
        return GetFirstSidechainAtom().GetCoordinates();
      }

      // all the side chain atom types of that particular amino acid
      std::vector< AtomType> atom_types_sidechain;

      // fill atom types with difference of all atom types for that amino acid and the backbone amino acid types
      std::set_difference
      (
        GetTypesOfAtoms().Begin(), GetTypesOfAtoms().End(),
        GetAtomTypes().GetBackBoneAtomTypes().Begin(), GetAtomTypes().GetBackBoneAtomTypes().End(),
        std::back_insert_iterator< std::vector< AtomType> >( atom_types_sidechain)
      );

      // obtain all the coordinates for the side chain atoms
      util::SiPtrVector< const linal::Vector3D> side_chain_atom_coordinates
      (
        GetAtomCoordinates( storage::Set< AtomType>( atom_types_sidechain.begin(), atom_types_sidechain.end()))
      );

      // calculate center of mass of the side chain based on all the side chain coordinates
      return coord::CenterOfMass( side_chain_atom_coordinates);
    }

    //! @brief calculate phi and psi
    //! @param PREVIOUS_C previous carbon atom
    //! @param FOLLOWING_N followin nitrogen atom
    //! @return phi and psi angles
    const storage::Pair< double, double> AABase::CalculatePhiPsi( const Atom &PREVIOUS_C, const Atom &FOLLOWING_N) const
    {
      return storage::Pair< double, double>( CalculatePhi( PREVIOUS_C), CalculatePsi( FOLLOWING_N));
    }

    //! @brief create a locator for that aa
    //! @param USE_PDB_ID
    //! @return locator aa to locate an amino acid in a protein model
    assemble::LocatorAA AABase::Locator( const bool USE_PDB_ID) const
    {
      return assemble::LocatorAA( m_Data->GetChainID(), USE_PDB_ID ? m_Data->GetPdbID() : m_Data->GetSeqID(), USE_PDB_ID);
    }

  ////////////////////////////////////////
  // data access - SSPrediction related //
  ////////////////////////////////////////

    //! @brief return map of secondary structure predictions stored for various method
    //! @return map of secondary structure predictions stored for various method
    const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &AABase::GetSSPredictions() const
    {
      return m_Data->m_SSPredictions;
    }

    //! @brief secondary structure predictions stored for given SS_METHOD
    //! @brief METHOD sspred::Method of interest
    //! @return secondary structure predictions stored for SS_METHOD
    util::SiPtr< const sspred::MethodInterface> AABase::GetSSPrediction
    (
      const sspred::Method &SS_METHOD
    ) const
    {
      return m_Data->GetSSPrediction( SS_METHOD);
    }

    //! @brief provide sspredictions from different methods in PREDICTIONS_MAP to overwrite existing ones
    //! @param PREDICTIONS_MAP Map of predictions with different methods
    void AABase::SetSSPredictions
    (
      const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &PREDICTIONS_MAP
    )
    {
      m_Data->m_SSPredictions = PREDICTIONS_MAP;
    }

    //! @brief Sets the predictions for the given SS_METHOD
    //! @param SS_METHOD SSMethod of interest
    //! @param PREDICTION MethodInterface derived class that stores the predictions
    void AABase::SetSSPrediction
    (
      const sspred::Method &SS_METHOD,
      const sspred::MethodInterface &PREDICTION
    )
    {
      m_Data->m_SSPredictions[ SS_METHOD] = util::ShPtr< sspred::MethodInterface>( PREDICTION.Clone());
    }

    //! @brief remove structure-based secondary-structure and TM methods from the map for this AA
    void AABase::RemoveStructureBasedSSTMInfo()
    {
      for
      (
        storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> >::iterator
          itr( m_Data->m_SSPredictions.Begin()), itr_end( m_Data->m_SSPredictions.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr->first)->GetIsDeterminedFromSturcture())
        {
          itr->second.Reset();
        }
      }
    }

    //! @brief get the exposure prediction
    //! @return exposure prediction
    const double &AABase::GetExposurePrediction() const
    {
      return m_Data->m_ExposurePrediction;
    }

    //! @brief set the exposure prediction
    //! @param EXPOSURE exposure to set
    void AABase::SetExposurePrediction( const double &EXPOSURE)
    {
      m_Data->m_ExposurePrediction = EXPOSURE;
    }

  /////////////////////////////////
  // data access - Blast related //
  /////////////////////////////////

    //! @brief returns BlastProfile Object
    //! @param ISTREAM input stream to read the blast profile from
    void AABase::ReadBlastProfile( std::istream &ISTREAM)
    {
      if( !m_Data->m_BlastProfile.IsDefined())
      {
        m_Data->m_BlastProfile = util::ShPtr< BlastProfile>( new BlastProfile);
      }
      m_Data->m_BlastProfile->ReadProfile( ISTREAM);
    }

    //! @brief returns const BlastProfile Object
    //! @return const BlastProfile Object
    const BlastProfile &AABase::GetBlastProfile() const
    {
      static const BlastProfile s_undefined_profile;
      return m_Data->m_BlastProfile.IsDefined() ? *m_Data->m_BlastProfile : s_undefined_profile;
    }

    //! @brief returns const BlastProfile Object
    //! @return const BlastProfile Object
    const util::ShPtr< BlastProfile> &AABase::GetBlastProfilePtr() const
    {
      return m_Data->m_BlastProfile;
    }

    //! @brief set BlastProfile and Probabilities accordingly to given BLAST_PROFILE object
    //! @param BLAST_PROFILE blast profile object to be set to
    void AABase::SetBlastProfile( const BlastProfile &BLAST_PROFILE)
    {
      m_Data->m_BlastProfile = util::ShPtr< BlastProfile>( BLAST_PROFILE.Clone());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns whether all coordinates for the given atom types are defined
    //! @param ATOM_TYPES Atom Types of interest
    bool AABase::HasDefinedCoordinates( const storage::Set< AtomType> &ATOM_TYPES) const
    {
      // get all the atom coordinates
      const util::SiPtrVector< const linal::Vector3D> coordinates( GetAtomCoordinates( ATOM_TYPES));

      // iterate over all coordinates
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coord_itr( coordinates.Begin()), coord_itr_end( coordinates.End());
        coord_itr != coord_itr_end; ++coord_itr
      )
      {
        // if any one of them is undefined
        if( !( *coord_itr)->IsDefined())
        {
          // then return failure
          return false;
        }
      }

      // this point is reached only if all coordinates were defined, therefore return true
      return true;
    }

    //! @brief check if this aa precedes a given amino acid
    //! @param AMINO_ACID second amino acid, that could precede this
    //! @return true if both amino acids have same chain id and if (other seq id) - (this seq_id) == 1
    bool AABase::DoesPrecede( const AABase &AMINO_ACID) const
    {
      return AMINO_ACID.GetChainID() == GetChainID() && ( AMINO_ACID.GetSeqID() - GetSeqID()) == 1;
    }

    //! @brief gives dihedral angles of the side chain - starting closest to backbone and moving out along side chain
    //! @return vector holding the dihedral angles - first element is dihedral angle closest to backbone
    Rotamer AABase::CalculateSideChainDihedralAngles() const
    {
      // get the atom types that are involved in the side chain dihedral angles for this residue
      const storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> > atom_types
      (
        m_Data->GetType()->GetSideChainDihedralAngleAtomTypes()
      );

      // will hold the dihedral angles
      Rotamer dihedral_angles;

      // true if there are not enough atom types to calculate a dihedral angle
      if( atom_types.IsEmpty())
      {
        BCL_MessageCrt
        (
          "cannot calculate dihedral angles of residue " + GetIdentification() +
          " because the atom types specified for calculating dihedral angles are " + util::Format()( atom_types)
        );
        // return empty vector
        return dihedral_angles;
      }

      // iterate through the atom types
      for
      (
        storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> >::const_iterator
          chi_itr( atom_types.Begin()), chi_itr_end( atom_types.End());
        chi_itr != chi_itr_end;
        ++chi_itr
      )
      {
        // get the atoms
        const linal::Vector3D &coords_a( GetAtom( chi_itr->second( 0)).GetCoordinates());
        const linal::Vector3D &coords_b( GetAtom( chi_itr->second( 1)).GetCoordinates());
        const linal::Vector3D &coords_c( GetAtom( chi_itr->second( 2)).GetCoordinates());
        const linal::Vector3D &coords_d( GetAtom( chi_itr->second( 3)).GetCoordinates());

        // true if any of the coordinates are not defined
        if( !coords_a.IsDefined() || !coords_b.IsDefined() || !coords_c.IsDefined() || !coords_d.IsDefined())
        {
          // add the dihedral to vector
          dihedral_angles.Insert( ChiAngle( chi_itr->first, util::GetUndefinedDouble(), math::Angle::e_Radian));

          // go to to next dihedral
          continue;
        }

        // get the dihedral angle
        const double dihedral( linal::Dihedral( coords_a, coords_b, coords_c, coords_d));

        // add the dihedral to vector
        dihedral_angles.Insert( ChiAngle( chi_itr->first, dihedral, math::Angle::e_Radian));
      }

      return dihedral_angles;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator for AABase class
    //! @param AA_BASE_RHS AABase to be assigned from
    //! @return this AABase after being updated to AA_BASE_RHS
    AABase &AABase::operator =( const AABase &AA_BASE_RHS)
    {
      // update data
      m_Data = AA_BASE_RHS.m_Data;

      // return
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read fasta from given ISTREAM with given SEQ_ID
    //! @param ISTREAM input stream
    //! @param SEQ_ID seq id for this AABase
    //! @return std::istream from which was read
    std::istream &AABase::ReadFasta( std::istream &ISTREAM, const size_t SEQ_ID)
    {
      // initialize one_letter_code
      char one_letter_code;

      // skip newline and or empty characters
      do
      {
        // if the end of stream is reached
        if( !( ISTREAM >> one_letter_code))
        {
          // end
          return ISTREAM;
        }
      }
      while( one_letter_code == '\n' || one_letter_code == ' ');

      // change this aa to the AAType that come from the one letter code and also update the seqids.
      m_Data = util::ShPtr< AAData>
      (
        new AAData
        (
          GetAATypes().AATypeFromOneLetterCode( one_letter_code),
          SEQ_ID,
          SEQ_ID,
          m_Data->m_PdbICode,
          m_Data->m_ChainID
        )
      );

      //end
      return ISTREAM;
    }

    //! @brief write fasta to provided OSTREAM
    //! @param OSTREAM output stream
    //! @return std::ostream to which was written
    std::ostream &AABase::WriteFasta( std::ostream &OSTREAM) const
    {
      // write one letter code to OSTREAM
      OSTREAM << m_Data->GetType()->GetOneLetterCode();

      // return
      return OSTREAM;
    }

    //! @brief get Identification of this amino acid
    //! @return string with identification
    std::string AABase::GetIdentification() const
    {
      // initialize identification with sequence id, one letter and three letter code and one state ss prediction

      std::string identification
      (
        util::Format().W( 5)( GetSeqID()) + " " +
        GetType()->GetOneLetterCode() + " " +
        GetType()->GetThreeLetterCode() + " "
      );

      // check if the PDB SSType was set
      if( GetSSPrediction( sspred::GetMethods().e_PDB).IsDefined())
      {
        identification += GetSSPrediction( sspred::GetMethods().e_PDB)->GetOneStateSSPrediction()->GetOneLetterCode();
      }
      // if jufo was provided
      else if( GetSSPrediction( sspred::GetMethods().e_JUFO).IsDefined())
      {
        identification += GetSSPrediction( sspred::GetMethods().e_JUFO)->GetOneStateSSPrediction()->GetOneLetterCode();
      }
      // else print
      else
      {
        identification += "U";
      }

      // end
      return identification;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AABase::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Data, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write Identification to given OSTREAM
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return std::ostream to which was written
    std::ostream &AABase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Data, OSTREAM, INDENT);

      //end
      return OSTREAM;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief Calculate the distance between two aminoacids by calling AtomDistance function for the first side chain atoms
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return CB distance between AMINO_ACID_A and AMINO_ACID_B
    double FirstSidechainAtomDistance( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B)
    {
      return Distance( AMINO_ACID_A.GetFirstSidechainAtom(), AMINO_ACID_B.GetFirstSidechainAtom());
    }

    //! @brief Calculate the sequence separation between two amino acids
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return sequence separation between AMINO_ACID_A and AMINO_ACID_B, will be undefined if they are from the same
    //!          chain or if the sequence ids are the same
    size_t SequenceSeparation( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B)
    {
      // if the chain ids are different return undefined
      if( AMINO_ACID_A.GetChainID() != AMINO_ACID_B.GetChainID())
      {
        return util::GetUndefined< size_t>();
      }
      // if the sequence ids are the same return undefined
      if( AMINO_ACID_A.GetSeqID() == AMINO_ACID_B.GetSeqID())
      {
        return util::GetUndefined< size_t>();
      }
      // if AMINO_ACID_A comes first
      if( AMINO_ACID_A.GetSeqID() < AMINO_ACID_B.GetSeqID())
      {
        return AMINO_ACID_B.GetSeqID() - AMINO_ACID_A.GetSeqID() - 1;
      }
      // else
      return AMINO_ACID_A.GetSeqID() - AMINO_ACID_B.GetSeqID() - 1;
    }

    //! @brief check if two amino acids are in range to be peptide bonded
    //! @param AA_LEFT left in sequence
    //! @param AA_RIGHT right in sequence
    //! @param TEST_OMEGA check the angle also
    //! @return true if amino acids are close enough for a peptide bond
    bool AABase::AreAminoAcidsPeptideBonded
    (
      const AABase &AA_LEFT,
      const AABase &AA_RIGHT,
      const bool TEST_OMEGA
    )
    {
      static const double s_peptide_bond_length( GetAtomTypes().C->GetBondLength( GetAtomTypes().N));
      static const double s_peptide_bond_tolerance( 0.02 * s_peptide_bond_length);

      const storage::VectorND< 2, double> bond_length_angle( PeptideBondLengthAndAngle( AA_LEFT, AA_RIGHT));

      // check the distance
      if( !util::IsDefined( bond_length_angle.First()) || !math::EqualWithinAbsoluteTolerance( s_peptide_bond_length, bond_length_angle.First(), s_peptide_bond_tolerance))
      {
//        if( bond_length_angle.First() < 1.5)
//        {
//          BCL_Message
//          (
//            util::Message::e_Verbose,
//            "barely missed peptide bond: " + util::Format()( bond_length_angle.First()) + " "
//            + AA_LEFT.GetIdentification() + " " + AA_RIGHT.GetIdentification()
//          );
//        }
        return false;
      }

      // distance matched, need to test angle also?
      if( !TEST_OMEGA)
      {
        return true;
      }

      // tolerance of 10 degrees
      static const double s_angle_tolerance( math::g_Pi / 18.0);

      // check angle
      if
      (
           !util::IsDefined( bond_length_angle.Second())                                             // will be undefined if one of the coordinates is not defined
        ||
           (
                 !math::EqualWithinAbsoluteTolerance(    -math::g_Pi, bond_length_angle.Second(), s_angle_tolerance) // check for cis peptide bond
              && !math::EqualWithinAbsoluteTolerance(            0.0, bond_length_angle.Second(), s_angle_tolerance) // check for cis peptide bond
              && !math::EqualWithinAbsoluteTolerance(     math::g_Pi, bond_length_angle.Second(), s_angle_tolerance) // check for trans peptide bond
           )
      )
      {
//        BCL_Message
//        (
//          util::Message::e_Verbose,
//          "missed peptide bond angle: " + util::Format()( bond_length_angle.Second()) + " " +
//          AA_LEFT.GetIdentification() + " " + AA_RIGHT.GetIdentification()
//        );
        return false;
      }

      // meets all conditions
      return true;
    }

    //! @brief calculate peptide bond length and angle
    //! @param AA_LEFT left in sequence
    //! @param AA_RIGHT right in sequence
    //! @return vector of peptide bond length (C->N) and peptide bond angle ca->c->n->ca)
    storage::VectorND< 2, double> AABase::PeptideBondLengthAndAngle
    (
      const AABase &AA_LEFT,
      const AABase &AA_RIGHT
    )
    {
      const linal::Vector3D &coord_ca_left( AA_LEFT.GetAtom( GetAtomTypes().CA).GetCoordinates());
      const linal::Vector3D &coord_c( AA_LEFT.GetAtom( GetAtomTypes().C).GetCoordinates());
      const linal::Vector3D &coord_n( AA_RIGHT.GetAtom( GetAtomTypes().N).GetCoordinates());
      const linal::Vector3D &coord_ca_right( AA_RIGHT.GetAtom( GetAtomTypes().CA).GetCoordinates());

      return
        storage::VectorND< 2, double>
        (
          linal::Distance( coord_c, coord_n), // bond
          linal::Dihedral( coord_ca_left, coord_c, coord_n, coord_ca_right) // angle
        );
    }

    //! @brief Compute the nearest atom separation between two AAs (nearest of atom distances - VdW radii of atoms involved)
    //! @param AMINO_ACID_A first amino acid
    //! @param AMINO_ACID_B second amino acid
    //! @return the nearest atom separation between two AAs (nearest of atom distances - VdW radii of atoms involved)
    double NearestAtomVdWSphereSeparation( const AABase &AMINO_ACID_A, const AABase &AMINO_ACID_B, const bool &IGNORE_BB_BB)
    {
      double nearest_so_far( std::numeric_limits< double>::max());
      storage::Vector< double> vdw_radii_b;
      vdw_radii_b.AllocateMemory( AMINO_ACID_B.GetAtoms().GetSize());
      const AAType type_a( AMINO_ACID_A.GetType()), type_b( AMINO_ACID_B.GetType());
      for( auto itr_b( AMINO_ACID_B.GetAtoms().Begin()), itr_b_end( AMINO_ACID_B.GetAtoms().End()); itr_b != itr_b_end; ++itr_b)
      {
        vdw_radii_b.PushBack
        (
          type_b->GetChemistryAtomType( ( *itr_b)->GetType())->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD)
        );
      }
      for( auto itr_a( AMINO_ACID_A.GetAtoms().Begin()), itr_a_end( AMINO_ACID_A.GetAtoms().End()); itr_a != itr_a_end; ++itr_a)
      {
        if( !( *itr_a)->GetCoordinates().IsDefined())
        {
          continue;
        }
        if( IGNORE_BB_BB && ( *itr_a)->GetType()->IsBackBone() && ( *itr_a)->GetType() != GetAtomTypes().CA)
        {
          continue;
        }
        const double vdw_radius_a
        (
          type_a->GetChemistryAtomType( ( *itr_a)->GetType())->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD)
        );
        auto itr_vdw_radii_b( vdw_radii_b.Begin());
        for( auto itr_b( AMINO_ACID_B.GetAtoms().Begin()), itr_b_end( AMINO_ACID_B.GetAtoms().End()); itr_b != itr_b_end; ++itr_b, ++itr_vdw_radii_b)
        {
          if( !( *itr_b)->GetCoordinates().IsDefined())
          {
            continue;
          }
          if( IGNORE_BB_BB && ( *itr_b)->GetType()->IsBackBone() && ( *itr_b)->GetType() != GetAtomTypes().CA)
          {
            continue;
          }
          nearest_so_far
            = std::min
              (
                nearest_so_far,
                std::max
                (
                  0.0,
                  linal::Distance( ( *itr_b)->GetCoordinates(), ( *itr_a)->GetCoordinates())
                  - *itr_vdw_radii_b - vdw_radius_a
                )
              );
        }
      }
      return nearest_so_far;
    }

    //! @brief construct single letter code and seqid
    //! @param ONE_LETTER_CODE aa one letter code
    //! @param SEQ_ID sequence id
    //! @return pointer to new amino acid
    AABase *AABase::Construct( const char ONE_LETTER_CODE, const int SEQ_ID)
    {
      // create ShPtr to AAData and set aatype, seq_id, chain_id
      util::ShPtr< AAData> new_member_data
      (
        new AAData
        (
          GetAATypes().AATypeFromOneLetterCode( toupper( ONE_LETTER_CODE)),
          SEQ_ID,
          AAData::s_DefaultPdbID,
          AAData::s_DefaultPdbICode,
          'A'
        )
      );

      // create AA
      return new AA( new_member_data);
    }

  } // namespace biol

  namespace align
  {
  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief specialization of the template function t_Member objects of type AABase
    //! @param MEMBER the t_Member object to extract the character id from
    //! @return the character id
    template<>
    char GetCharId< biol::AABase>( const biol::AABase &MEMBER)
    {
      return MEMBER.GetType()->GetOneLetterCode();
    }

    //! @brief specialization of the template function t_Member objects of type AABase
    //! @param MEMBER the t_Member object to extract the complete identifier from
    //! @return the complete identifier
    template<>
    std::string GetCompleteId< biol::AABase>( const biol::AABase &MEMBER)
    {
      return MEMBER.GetIdentification();
    }

  } // namespace align
} // namespace bcl
