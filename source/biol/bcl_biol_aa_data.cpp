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
#include "biol/bcl_biol_aa_data.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAData::s_Instance
    (
      GetObjectInstances().AddInstance( new AAData())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct AAData from optional AAType and other information
    //! @param AA_TYPE amino acid type
    //! @param SEQ_ID sequence id
    //! @param PDB_ID pdb-file ID
    //! @param PDB_I_CODE insertion code for pdb-residues
    //! @param CHAIN_ID chain id of the amino acid
    AAData::AAData
    (
      const AAType &AA_TYPE,
      const int SEQ_ID,
      const int PDB_ID,
      const char PDB_I_CODE,
      const char CHAIN_ID
    ) :
      m_Type( AA_TYPE),
      m_SeqID( SEQ_ID),
      m_PdbID( PDB_ID),
      m_PdbICode( PDB_I_CODE),
      m_ChainID( CHAIN_ID),
      m_BlastProfile(),
      m_SSPredictions(),
      m_ExposurePrediction( util::GetUndefined< double>())
    {
    }

    //! @brief copy constructor
    AAData *AAData::Clone() const
    {
      return new AAData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return Type
    const AAType &AAData::GetType() const
    {
      return m_Type;
    }

    //! @brief return SeqID
    //! @return SeqID
    const int &AAData::GetSeqID() const
    {
      return m_SeqID;
    }

    //! @brief sets the sequence id to the given value
    //@ @param ID new sequence id
    void AAData::SetSeqID( int ID)
    {
      m_SeqID = ID;
    }

    //! @brief return PdbID
    //! @return PdbID
    const int &AAData::GetPdbID() const
    {
      return m_PdbID;
    }

    //! @brief return PdbICode
    //! @return PdbICode
    const char &AAData::GetPdbICode() const
    {
      return m_PdbICode;
    }

    //! @brief return chain id
    //! @return chain id
    char AAData::GetChainID() const
    {
      return m_ChainID;
    }

    //! @brief sets the chain id to the given value
    //! @param ID new chain id
    void AAData::SetChainID( char ID)
    {
      m_ChainID = ID;
    }

    //! @brief return blast profile
    //! @return blast profile
    const util::ShPtr< BlastProfile> &AAData::GetBlastProfile() const
    {
      return m_BlastProfile;
    }

    //! @brief return SSPredictions set
    //! @return SSPredictions set
    const storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> > &AAData::GetSSPredictions() const
    {
      return m_SSPredictions;
    }

    //! @brief return SSPredictions for given SS_METHOD
    //! @param SS_METHOD method of interest
    //! @return SSPredictions for given SS_METHOD
    util::SiPtr< const sspred::MethodInterface> AAData::GetSSPrediction( const sspred::Method &SS_METHOD) const
    {
      // search for the predictions for SS_METHOD in this map and store the iterator
      storage::Map< sspred::Method, util::ShPtr< sspred::MethodInterface> >::const_iterator method_itr
      (
        m_SSPredictions.Find( SS_METHOD)
      );

      // if equal to end
      if( method_itr == m_SSPredictions.End())
      {
        // return undefined
        return util::SiPtr< const sspred::MethodInterface>();
      }

      // if found return itr dereferenced
      return method_itr->second;
    }

    //! @brief get the exposure prediction
    //! @return exposure prediction
    const double &AAData::GetExposurePrediction() const
    {
      return m_ExposurePrediction;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief smaller than, that compares the position in the sequence
    //! @param RHS_AA_DATA
    //! @return true if given AAData comes after this aa data by chain id, or for the same chain by seqid
    bool AAData::operator <( const AAData &RHS_AA_DATA) const
    {
      // compare chain id
      if( m_ChainID < RHS_AA_DATA.m_ChainID)
      {
        return true;
      }
      else if( m_ChainID > RHS_AA_DATA.m_ChainID)
      {
        return false;
      }

      // compare sequence id
      if( m_SeqID < RHS_AA_DATA.m_SeqID)
      {
        return true;
      }

      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AAData::Read( std::istream &ISTREAM)
    {
      // read data
      io::Serialize::Read( m_Type, ISTREAM);
      io::Serialize::Read( m_SeqID, ISTREAM);
      io::Serialize::Read( m_PdbID, ISTREAM);
      io::Serialize::Read( m_PdbICode, ISTREAM);
      io::Serialize::Read( m_ChainID, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT indentation
    //! @return ostream which was written to
    std::ostream &AAData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write data
      io::Serialize::Write( m_Type, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SeqID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PdbID, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PdbICode, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ChainID, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace biol
} // namespace bcl
