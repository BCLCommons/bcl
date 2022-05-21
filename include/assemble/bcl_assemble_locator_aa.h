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

#ifndef BCL_ASSEMBLE_LOCATOR_AA_H_
#define BCL_ASSEMBLE_LOCATOR_AA_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_locator_chain.h"
#include "biol/bcl_biol_aa_types.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorAA
    //! @brief is used for locating a specified amino acid from a protein model
    //! @details This class is used for locating a specified amino acid from a given protein model. It uses first the
    //! member chain locator to find the corresponding chain and then iterates over the amino acids of the SSEs
    //! to find the amino acid with the specified amino acid seq id.
    //!
    //! @see @link example_assemble_locator_aa.cpp @endlink
    //! @author alexanns
    //! @date 01/16/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorAA :
      public find::LocatorInterface< util::SiPtr< const biol::AABase>, ProteinModel>
    {

    private:

    //////////
    // data //
    //////////

      LocatorChain m_LocatorChain; //!< Chain Locator class
      int          m_AA_ID;        //!< the SeqID of the amino acid
      bool         m_UsePDBID;     //!< use the pdb id instead of the seq id to locate the amino acid
      biol::AAType m_AAType;       //!< Amino acid type; used for verification

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      LocatorAA() :
        m_LocatorChain( 'A'),
        m_AA_ID( 1),
        m_UsePDBID( false)
      {
      }

      //! constructor from ChainID and amino acid ID
      //! @param CHAIN_ID char which indicates the chain
      //! @param AA_ID int which indicates the SeqID of the amino acid
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorAA( const char CHAIN_ID, const int AA_ID, const bool USE_PDB_ID = false) :
        m_LocatorChain( CHAIN_ID),
        m_AA_ID( AA_ID),
        m_UsePDBID( USE_PDB_ID)
      {
      }

      //! constructor from ChainID, amino acid ID, and aa type
      //! @param CHAIN_ID char which indicates the chain
      //! @param AA_ID int which indicates the SeqID of the amino acid
      //! @param TYPE the expected aa type
      //! @param USE_PDB_ID true if the pdd id should be used to locate the residue instead of the seqid
      LocatorAA( const char CHAIN_ID, const int AA_ID, const biol::AAType &TYPE, const bool USE_PDB_ID = false) :
        m_LocatorChain( CHAIN_ID),
        m_AA_ID( AA_ID),
        m_UsePDBID( USE_PDB_ID),
        m_AAType( TYPE)
      {
      }

      //! clone constructor
      LocatorAA *Clone() const
      {
        return new LocatorAA( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives formatted string description
      //! @return gives formatted string description
      std::string GetIdentification() const;

      //! @brief reads formatted string describing the locator
      //! @return formatted string describing the locator
      std::istream &ReadIdentification( std::istream &ISTREAM);

      //! @brief returns a bool indicating if the pdb ids will be used instead of the seq ids
      //! @return boolean true if pdb ids should be used - false otherwise
      bool GetUsePDBID() const;

      //! @brief returns aa type, if known
      //! @return aa type, if known
      const biol::AAType &GetAAType() const;

      //! @brief GetAAID gives the seq id of the amino acid to be located
      //! @return returns "m_AA_ID" int
      int const &GetAAID() const;

      //! @brief SetAAID changes the amino acid id
      //! @param AA_ID int which indicates the new amino acid id
      void SetAAID( const int AA_ID);

      //! @brief returns the chain locator
      //! @return returns the const reference to the chain locator
      const LocatorChain &GetLocatorChain() const;

      //! @brief SetLocatorChain the chain locator to a new one
      //! @param LOCATOR_CHAIN the new locator chain which "m_LocatorChain" will be set to
      void SetLocatorChain( const LocatorChain &LOCATOR_CHAIN);

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief GetAA finds the amino acid denoted by the LocatorAA
      //! @param SSE SSE which the LocatorAA refers to
      //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
      util::SiPtr< const biol::AABase> Locate( const SSE &SSE) const;

      //! @brief GetAA finds the amino acid denoted by the LocatorAA
      //! @param SSE SSE which the LocatorAA refers to
      //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
      util::SiPtr< biol::AABase> Locate( SSE &SSE) const;

      //! @brief Locate finds the amino acid denoted by the LocatorAA in a list of residues
      //! @param RESIDUE_LIST SSE the list of residues which will be searched for the residue of interest
      //! @return returns SiPtr to the amino acid denoted by the LocatorAA
      util::SiPtr< const biol::AABase> Locate( const util::SiPtrVector< const biol::AABase> &RESIDUE_LIST) const;

      //! @brief Locate finds the amino acid denoted by the LocatorAA in a list of residues
      //! @param RESIDUE_LIST SSE the list of residues which will be searched for the residue of interest
      //! @return returns SiPtr to the amino acid denoted by the LocatorAA
      util::SiPtr< biol::AABase> Locate( util::SiPtrVector< biol::AABase> &RESIDUE_LIST) const;

      //! @brief GetSSE finds the SSE which contains the amino acid denoted by the LocatorAA
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
      //! @return returns SharedPointer to the SSE containing the amino acid denoted by the LocatorAA
      util::SiPtr< const SSE> LocateSSE( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief GetAA finds the amino acid denoted by the LocatorAA
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
      //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
      util::SiPtr< const biol::AABase> Locate( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief GetSSE finds the SSE which contains the amino acid denoted by the LocatorAA
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
      //! @return returns SharedPointer to the SSE containing the amino acid denoted by the LocatorAA
      util::SiPtr< SSE> LocateSSE( ProteinModel &PROTEIN_MODEL) const;

      //! @brief GetAA finds the amino acid denoted by the LocatorAA
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
      //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
      util::SiPtr< biol::AABase> Locate( ProteinModel &PROTEIN_MODEL) const;

      //! @brief GetAA finds the amino acid denoted by the LocatorAA
      //! @param PROTEIN_MODEL ProteinModel which the LocatorAA refers to
      //! @return returns SharedPointer to the amino acid denoted by the LocatorAA
      util::SiPtr< const biol::AABase> Locate( const DomainInterface &PROTEIN_MODEL) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LocatorAA

    //! @brief less than operator for comparing two locator aa objects by seq and chain id
    //! @param LOCATOR_A the first locator
    //! @param LOCATOR_B the second locator
    //! @return boolean true if LOCATOR_A is less than LOCATOR_B - false otherwise
    BCL_API
    bool operator <( const LocatorAA &LOCATOR_A, const LocatorAA &LOCATOR_B);

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_AA_H_
