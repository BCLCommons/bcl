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

#ifndef BCL_BIOL_DSSP_H_
#define BCL_BIOL_DSSP_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_mutate_interface.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically
#include <deque>

namespace bcl
{
  namespace biol
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DSSP
    //! @brief use DSSP to assign secondary structure to a protein model
    //! @details using DSSP (W. KABSCH AND C.SANDER, BIOPOLYMERS 22 (1983) 2577-2637) hydrogen bonds between amino acids
    //!          are calculated, beta bridges are determined and according to that, secondary structure is assigned to
    //!          the amino acids
    //!
    //! @see @link example_biol_dssp.cpp @endlink
    //! @author woetzen
    //! @date Aug 22, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DSSP :
      public math::MutateInterface< assemble::ProteinModel>
    {

    private:

      //! @
      enum BridgeType
      {
        e_BTNoBridge, e_BTParallel, e_BTAntiParallel
      };

      enum SecondaryStructureType
      {
        e_Loop,   //' '
        e_Alphahelix, // H
        e_Betabridge, // B
        e_Strand,   // E
        e_Helix_3,  // G
        e_Helix_5,  // I
        e_Turn,   // T
        e_Bend    // S
      };

      enum HelixFlagType
      {
        e_HelixNone, e_HelixStart, e_HelixEnd, e_HelixStartAndEnd, e_HelixMiddle
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class Bridge
      //! @brief TODO
      //! @author woetzen
      //! @date Aug 22, 2012
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct Bridge
      {
        BridgeType m_Type;
        size_t     m_Sheet;
        size_t     m_Ladder;
        std::set< Bridge *> m_Link;
        std::deque< util::SiPtr< const AABase> > m_I;
        std::deque< util::SiPtr< const AABase> > m_J;
        char      m_ChainI;
        char      m_ChainJ;

        //! @brief compare bridge by chainid and seqid of residues in bridge
        //! @param BRIDGE bridge to comparet his to
        //! @return true if this bridge comes before the argument BRIDGE in sequence
        bool operator<( const Bridge &BRIDGE) const;

      }; // struct Bridge

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class BridgePartner
      //! @brief TODO
      //! @author woetzen
      //! @date Aug 22, 2012
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class BridgePartner :
        public util::ObjectInterface
      {
      public:
        util::SiPtr< const AABase> m_Partner;
        size_t                     m_Ladder;
        bool                       m_Parallel;

        //! @brief constructor from members
        BridgePartner();

        //! @brief constructor from members
        BridgePartner( const util::SiPtr< const AABase> &PARTNER, const size_t LADDER, const bool PARALLEL);

        //! @brief Clone function
        //! @return pointer to new BridgePartner
        BridgePartner *Clone() const
        {
          return new BridgePartner( *this);
        }

        //! @brief returns class name
        //! @return the class name as const ref std::string
        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName( *this);
        }

      //////////////////////
      // input and output //
      //////////////////////

      protected:

        //! @brief read from std::istream
        //! @param ISTREAM input stream
        //! @return istream which was read from
        std::istream &Read( std::istream &ISTREAM)
        {
          return ISTREAM;
        }

        //! @brief write to std::ostream
        //! @param OSTREAM outputstream to write to
        //! @param INDENT number of indentations
        //! @return outputstream which was written to
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
        {
          return OSTREAM;
        }

      }; // class BridgePartner

    //////////
    // data //
    //////////

      //! model that was complete with required N-H atoms
      mutable util::ShPtr< assemble::ProteinModel> m_Model;

      //! @brief max hydrogen bond energy
      double m_MaxHydrogenBondEnergy;

      typedef storage::Map< util::SiPtr< const AABase>, storage::VectorND< 2, storage::Pair< util::SiPtr< const AABase>, double> > > HydrogenBondContainerType;
      //! for each amino acid in the model, the best two hbond donors are stored
      mutable HydrogenBondContainerType m_HydrogenBondDonor;
      //! for each amino acid in the model, the best two hbond acceptors are stored
      mutable HydrogenBondContainerType m_HydrogenBondAcceptor;

      typedef storage::Map< util::SiPtr< const AABase>, storage::VectorND< 2, BridgePartner> > BridgePartnerContainerType;
      //! store the bridge partner for each amino acid in the model
      mutable BridgePartnerContainerType m_BridgePartner;

      //! store the sheet number each amino acid belongs to
      mutable storage::Map< util::SiPtr< const AABase>, size_t> m_Sheet;

      typedef std::map< util::SiPtr< const AABase>, SecondaryStructureType> SSContainerType;
      //! store the secondary structure assignment of each amino acid
      mutable SSContainerType m_SecondaryStructure;

      typedef std::map< util::SiPtr< const AABase>, HelixFlagType> HelixFlagContainerType;
      //! store the helix flag for each stride and amino acid
      mutable HelixFlagContainerType m_HelixFlag[ 3];

      typedef storage::Set< util::SiPtr< const AABase> > BendContainerType;
      //! store if the amino acid is bend
      mutable BendContainerType m_Bend;

      typedef std::map< util::SiPtr< const AABase>, std::pair< double, char> > AlphaContainerType;
      //! store alpha and chirlaity for each amino acid
      mutable AlphaContainerType m_Alpha;

      typedef std::map< util::SiPtr< const AABase>, double> AngleContainerType;
      //! store kappa angle for each amino acid
      mutable AngleContainerType m_Kappa;
      //! store phi angle for each amino acid
      mutable AngleContainerType m_Phi;
      //! store psi angle for each amino acid
      mutable AngleContainerType m_Psi;
      //! store tco for each amino acid
      mutable AngleContainerType m_TCO;

      // statistics
      static const size_t s_HistogramSize = 30;
      mutable size_t m_NrOfHydrogenBondsInParallelBridges;
      mutable size_t m_NrOfHydrogenBondsInAntiparallelBridges;
      mutable size_t m_ParallelBridgesPerLadderHistogram[     s_HistogramSize];
      mutable size_t m_AntiparallelBridgesPerLadderHistogram[ s_HistogramSize];
      mutable size_t m_LaddersPerSheetHistogram[              s_HistogramSize];

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! @brief default max HBond energy for an HBond to be considered
      //! @return default max energy
      static const double &GetDefaultHBondMaxEnergy();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from max hydrogen bond energy
      //! @param MAX_HYDROGEN_BOND_ENERGY maximal hydrogen bond energy to be considered an H-bond
      DSSP( const double MAX_HYDROGEN_BOND_ENERGY = GetDefaultHBondMaxEnergy());

      //! @brief Clone function
      //! @return pointer to new DSSP
      DSSP *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this mutate
      //! @return the scheme for this mutate
      const std::string &GetScheme() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief get various statistics
      void GetStatistics
      (
        size_t &NR_OF_RESIDUES,
        size_t &NR_OF_CHAINS,
        size_t &NR_OF_SS_BRIDGES,
        size_t &NR_OF_INTRA_CHAIN_SS_BRIDGES,
        size_t &NR_OF_HYDROGEN_BONDS,
        size_t NR_OF_HYDROGEN_BONDS_PER_DISTANCE[ s_HistogramSize],
        size_t NR_RES_PER_ALPHA_HELIX_HISTOGRAM[ s_HistogramSize]
      ) const;

      //! @brief reset all members
      void Reset() const;

      //! @brief Set PDB SS Predictions using the DSSP algorithm
      //! @param PROTEIN_MODEL the protein model with sses (can be just one single loop)
      void SetPDBSSPred( assemble::ProteinModel &PROTEIN_MODEL) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that uses the coordinates of the secondary structure elements within them to generate a new
      //!        new secondary structure assignment using DSSP
      //! @param PROTEIN_MODEL the protein model with sses (can be just one single loop)
      //! @return a mutate result with a protein model that has new secondary structure elements based on dssp
      math::MutateResult< assemble::ProteinModel> operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to file in standard DSSP format
      //! @param OSTREAM stream to write to
      //! @return the ostream written to
      std::ostream &WriteToFile( std::ostream &OSTREAM) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      //! @brief write information for a single residue
      //! @param SP_AMINO_ACID SiPtr to amino acid residue
      //! @param OSTREAM stream to write to
      //! @return the stream written to
      std::ostream &WriteResidue( const util::SiPtr< const AABase> &SP_AMINO_ACID, std::ostream &OSTREAM) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculate all pairwise hydrogen bonding terms
      void CalculateHydrogenBondEnergies() const;

      //! @brief using all hydrogen bonds, determine all beta sheets
      void CalculateBetaSheets() const;

      //! @brief using all hydrogen bonds, determine all helix types
      void CalculateAlphaHelices() const;

      //! @brief test and calculate hydrogen bond energy between donor and acceptor, store in member
      //! @param AA_DONOR (hydrogen, nitrogen)
      //! @param AA_ACCEPTOR (oxygen, carbon)
      //! @return the energy calculated - if relevant interaction was found
      double CalculateHydrogenBondEnergy( const AABase &AA_DONOR, const AABase &AA_ACCEPTOR) const;

      //! @brief calculate alpha angle, store in member
      //! @param PREV previous amino acid
      //! @param CENTER center amino acid
      //! @param NEXT next amino acid
      //! @param NEXT_NEXT next amino acid after NEXT
      std::pair< double, char> CalculateAlpha( const AABase &PREV, const AABase &CENTER, const AABase &NEXT, const AABase &NEXT_NEXT) const;

      //! @brief test if three amino acids in one stretch - I is interacting via a beta bridge with a second stretch - J
      //! @param I_PREV previous amino acid in i
      //! @param I_CURRENT center amino acid in i
      //! @param I_NEXT next amino acid in i
      //! @param J_PREV previous amino acid in j
      //! @param J_CURRENT center amino acid in j
      //! @param J_NEXT next amino acid in j
      //! @return parallel or anti parallel or no bridge was found
      BridgeType TestBridge
      (
        const AABase &I_PREV, const AABase &I_CURR, const AABase &I_NEXT,
        const AABase &J_PREV, const AABase &J_CURR, const AABase &J_NEXT
      ) const;

      //! @brief test if donor and acceptor are connected through a hydrogen bond
      //! @param AA_DONOR donor amino acid
      //! @param AA_ACCEPTOR acceptor amino acid
      //! @return true if acceptor is forming a hydrogen bond to donor
      bool TestBond( const AABase &AA_DONOR, const AABase &AA_ACCEPTOR) const;

      //! @brief test if any of the residues in bridge a is identical to any of the residues in bridge b in case they
      //!        need to be joined
      //! @param BRIDGE_A bridge a potentially containing any AA of bridge b
      //! @param BRIDGE_B bridge b potentially containing any AA of bridge a
      //! @return true, if bridge and b are overlapping
      static bool Linked( const Bridge &BRIDGE_A, const Bridge &BRIDGE_B);

      //! @brief calculate kappa angle over 5 residues
      //! @param PREV_PREV center - 2 amino acid
      //! @param CENTER center amino acid
      //! @param NEXT_NEXT center + 2 amino acid
      //! @return the kappa angle
      static double Kappa( const AABase &PREV_PREV, const AABase &CENTER, const AABase &NEXT_NEXT);

      //! @brief test if given amino acid is a start residue for given helix stride
      //! @param AMINO_ACID the amino acid in question
      //! @param STRIDE the stride of the helix
      //! @return true if amino acid was flagged start compatible for the stride
      bool IsHelixStart( const AABase &AMINO_ACID, const size_t STRIDE) const;

    }; // class DSSP

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_DSSP_H_
