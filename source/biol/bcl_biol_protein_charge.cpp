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
#include "biol/bcl_biol_protein_charge.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_protein_params.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    const AATypeData::PropertyType ProteinCharge::s_FirstpKAAProperty   = AATypeData::e_pK_EMBOSS;
    const AATypeData::PropertyType ProteinCharge::s_LastpKAAProperty    = AATypeData::e_pK_ProMoST;

    const double ProteinCharge::s_pK_NTerm[ s_LastpKAAProperty - s_FirstpKAAProperty + 1] =
    {
      8.6, 8.0, 9.6, 8.2, 8.0, 11.2, 8.2, 9.69, 7.7, 0.0, 0.0
    };

    const double ProteinCharge::s_pK_CTerm[ s_LastpKAAProperty - s_FirstpKAAProperty + 1] =
    {
      3.6, 3.1, 2.4, 3.2, 3.1, 4.2, 3.65, 2.34, 3.3, 0.0, 0.0
    };

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinCharge::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinCharge())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinCharge::ProteinCharge() :
      m_AACount(),
      m_NTermAA(),
      m_CTermAA(),
      m_pKProperty( s_LastpKAAProperty)
    {
    }

    //! @brief constructor from amino acid sequence
    //! @param SEQUENCE the amino acid seuqence of the protein of interest
    ProteinCharge::ProteinCharge( const AASequence &SEQUENCE) :
      m_AACount( ProteinParams::CountAAs( SEQUENCE)),
      m_pKProperty( s_LastpKAAProperty)
    {
      if( SEQUENCE.GetSize() < 1)
      {
        return;
      }
      m_NTermAA = SEQUENCE.GetFirstAA()->GetType();
      m_CTermAA = SEQUENCE.GetLastAA()->GetType();
    }

    //! @brief construct from amino acid count
    //! @param AA_COUNT the number of amino acids
    //! @param N_TERM_AA n terminal amino acid
    //! @param C_TERM_AA c terminal amino acid
    ProteinCharge::ProteinCharge
    (
      const storage::Map< AAType, size_t> &AA_COUNT,
      const AAType &N_TERM_AA,
      const AAType &C_TERM_AA
    ) :
      m_AACount( AA_COUNT),
      m_NTermAA( N_TERM_AA),
      m_CTermAA( C_TERM_AA),
      m_pKProperty( s_LastpKAAProperty)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinCharge
    ProteinCharge *ProteinCharge::Clone() const
    {
      return new ProteinCharge( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinCharge::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinCharge::GetScheme() const
    {
      return this->GetClassIdentifier();
    }

    //! @brief set the pk aa property to be used
    //! @brief AA_PROPERTY the aa property to be used; needs to be within the allowed pk property range
    void ProteinCharge::SetPKProperty( const AATypeData::PropertyType &PROPERTY)
    {
      if( ( PROPERTY >= s_FirstpKAAProperty && PROPERTY <= s_LastpKAAProperty))
      {
        m_pKProperty = PROPERTY;
      }
      else
      {
        BCL_MessageStd
        (
          "setting the pk property to the default property: " +
          AATypeData::PropertyTypeEnum( s_LastpKAAProperty).GetString()
        );
        m_pKProperty = s_LastpKAAProperty;
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief for a given pH value clauclate the charge
    //! @param PH the pH of the solution
    //! @return the net charge of the protein within this pH
    double ProteinCharge::operator()( const double &PH) const
    {
      double net_charge( 0.0);

      for
      (
        storage::Map< AAType, size_t>::const_iterator itr( m_AACount.Begin()), itr_end( m_AACount.End());
        itr != itr_end;
        ++itr
      )
      {
        net_charge += double( itr->second) * HendersonHasselbachEquation( PH, itr->first->GetAAProperty( m_pKProperty), itr->first->GetAAProperty( AATypeData::e_Charge));
      }

      net_charge = CorrectTerminalCharge( net_charge, PH);

      return net_charge;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinCharge::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_AACount, ISTREAM);
      io::Serialize::Read( m_NTermAA, ISTREAM);
      io::Serialize::Read( m_CTermAA, ISTREAM);
      AATypeData::PropertyTypeEnum property;
      io::Serialize::Read( property, ISTREAM);
      SetPKProperty( property);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinCharge::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_AACount, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NTermAA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CTermAA, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( AATypeData::PropertyTypeEnum( m_pKProperty), OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief correct the given charge by the terminal residue charges for that ph
    //! @param PH the desired pH
    //! @param CHARGE ignoring the termini
    //! @return the correct charge by the terminal amino acids, that usually have different pK value
    double ProteinCharge::CorrectTerminalCharge( const double &CHARGE, const double &PH) const
    {
      double corrected_charge( CHARGE);

      // correct termini for ProMoST
      if( m_pKProperty == AATypeData::e_pK_ProMoST)
      {
//        corrected_charge -= HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( m_pKProperty), m_NTermAA->GetAAProperty( AATypeData::e_Charge));
//        corrected_charge -= HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( m_pKProperty), m_CTermAA->GetAAProperty( AATypeData::e_Charge));
        corrected_charge += HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( AATypeData::e_pK_ProMoST_NTerm),  1.0);
        corrected_charge += HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( AATypeData::e_pK_ProMoST_CTerm), -1.0);
      }

      // correct termini for Bjellqvist
      else if( m_pKProperty == AATypeData::e_pK_Bjellqvist)
      {
//        corrected_charge -= HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( m_pKProperty), m_NTermAA->GetAAProperty( AATypeData::e_Charge));
//        corrected_charge -= HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( m_pKProperty), m_CTermAA->GetAAProperty( AATypeData::e_Charge));
        corrected_charge += HendersonHasselbachEquation( PH, m_NTermAA->GetAAProperty( AATypeData::e_pK_Bjellqvist_NTerm),  1.0);
        corrected_charge += HendersonHasselbachEquation( PH, m_CTermAA->GetAAProperty( AATypeData::e_pK_Bjellqvist_CTerm), -1.0);
      }

      // for all other pK vlaue sources, C and N termini are the same for every amino acid
      else
      {
        corrected_charge += HendersonHasselbachEquation( PH, s_pK_CTerm[ m_pKProperty - s_FirstpKAAProperty],  1.0);
        corrected_charge += HendersonHasselbachEquation( PH, s_pK_NTerm[ m_pKProperty - s_FirstpKAAProperty], -1.0);
      }

      return corrected_charge;
    }

    //! @brief Henderson-Hasselbach-Equation
    //! @param PH the ph
    //! @param PK the pK
    //! @param AA_CHARGE the charge of the amino acid (positive or negative)
    //! @return the net charge - if pk is 0 charge returned is 0
    double ProteinCharge::HendersonHasselbachEquation( const double &PH, const double &PK, const double &AA_CHARGE)
    {
      if( PK == 0.0 || AA_CHARGE == 0.0)
      {
        return 0;
      }

      if( AA_CHARGE < 0.0)
      {
        return -1.0 / ( 1.0 + math::Pow( 10.0, PK - PH));
      }

      return 1.0 / ( 1.0 + math::Pow( 10.0, PH - PK));
    }

  } // namespace biol
} // namespace bcl
