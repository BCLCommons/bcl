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

#ifndef BCL_BIOL_AA_DATA_H_
#define BCL_BIOL_AA_DATA_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_types.h"
#include "bcl_biol_blast_profile.h"
#include "sspred/bcl_sspred_method_interface.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAData
    //! @brief This is a class for storing amino acid data
    //! @details Amino acid data class for storage of sequence data only. It stores
    //! <ul>
    //!   <li> amino acid type
    //!   <li> Sequence id ( bcl internal, starts with 1 for the first amino acid in a sequence
    //!   <li> PDB sequence id, can even be negative
    //!   <li> PDB insertion code
    //!   <li> chain id
    //!   <li> blast profile
    //!   <li> secondary structure predictions
    //! </ul>
    //! It only gives read access - once the aadata object is initialized, you can only get the data - you might use the
    //! Read function.
    //!
    //! It is a member of the AA, which is the only class that has direct access to the members and may change them. The
    //! AA class is the only object that uses this class. The advantage is, that when you copy an AA, the AAData will
    //! not be copied - just the pointer to the AAData object.
    //!
    //! @see @link example_biol_aa_data.cpp @endlink
    //! @author meilerj, woetzen
    //! @date 10.10.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAData :
      public util::ObjectInterface
    {

    /////////////
    // friends //
    /////////////

      friend class AA;
      friend class AABase;

    private:

    //////////
    // data //
    //////////

      AAType                     m_Type;          //!< AminoAcid Type
      int                        m_SeqID;         //!< Sequence ID
      int                        m_PdbID;         //!< Pdb-File ID
      char                       m_PdbICode;      //!< Insertion code for pdb-residues
      char                       m_ChainID;       //!< Chain id of the aminoacid
      util::ShPtr< BlastProfile> m_BlastProfile;  //!< Position based scoring matrix

      //! Secondary Structure Predictions
      storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > m_SSPredictions;

      double m_ExposurePrediction; //!< Exposure (rSASA) prediction

    public:

    //////////
    // data //
    //////////

      static const int s_DefaultSeqID = 1; //!< default seq id
      static const int s_DefaultPdbID = 1; //!< default pdb id
      static const char s_DefaultPdbICode = ' '; //!< default pdb insertion id
      static const char s_DefaultChainID = 'A'; //!< default chain id

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct AAData from optional AAType and other information
      //! @param AA_TYPE amino acid type
      //! @param SEQ_ID sequence id
      //! @param PDB_ID pdb-file ID
      //! @param PDB_I_CODE insertion code for pdb-residues
      //! @param CHAIN_ID chain id of the amino acid
      AAData
      (
        const AAType &AA_TYPE = GetAATypes().e_Undefined,
        const int SEQ_ID = s_DefaultSeqID,
        const int PDB_ID = s_DefaultPdbID,
        const char PDB_I_CODE = s_DefaultPdbICode,
        const char CHAIN_ID = s_DefaultChainID
      );

      //! @brief copy constructor
      AAData *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return Type
      const AAType &GetType() const;

      //! @brief return SeqID
      //! @return SeqID
      const int &GetSeqID() const;

      //! @brief sets the sequence id to the given value
      //@ @param ID new sequence id
      void SetSeqID( int ID);

      //! @brief return PdbID
      //! @return PdbID
      const int &GetPdbID() const;

      //! @brief return PdbICode
      //! @return PdbICode
      const char &GetPdbICode() const;

      //! @brief return chain id
      //! @return chain id
      char GetChainID() const;

      //! @brief sets the chain id to the given value
      //! @param ID new chain id
      void SetChainID( char ID);

      //! @brief return blast profile
      //! @return blast profile
      const util::ShPtr< BlastProfile> &GetBlastProfile() const;

      //! @brief return SSPredictions set
      //! @return SSPredictions set
      const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &GetSSPredictions() const;

      //! @brief return SSPredictions for given SS_METHOD
      //! @param SS_METHOD method of interest
      //! @return SSPredictions for given SS_METHOD
      util::SiPtr< const sspred::MethodInterface> GetSSPrediction( const sspred::Method &SS_METHOD) const;

      //! @brief get the exposure prediction
      //! @return exposure prediction
      const double &GetExposurePrediction() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief smaller than, that compares the position in the sequence
      //! @param RHS_AA_DATA
      //! @return true if given AAData comes after this aa data by chain id, or for the same chain by seqid
      bool operator <( const AAData &RHS_AA_DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class AAData

  } // namespace biol
} // namespace bcl

#endif //BCL_BIOL_AA_DATA_H_
