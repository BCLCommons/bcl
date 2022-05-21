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

#ifndef BCL_BIOL_PROTEIN_CHARGE_H_
#define BCL_BIOL_PROTEIN_CHARGE_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_types.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinCharge
    //! @brief calculate the net charge for a given pH
    //! @details this function is constructed from an amino acid sequence, and with the count of the positive and
    //!          negative amino acids, using the pK values and the Henderson-Hasselbach-Equation, the charge is calculated
    //!          this function can be used with a root finder to determine the pH where the net charge is 0 which equals
    //!          to the pI of the protein - the pH value where the protein would not migrate in a pH gradient gel (2D gel)
    //! @see @link example_biol_protein_charge.cpp @endlink
    //! @author woetzen
    //! @date Jun 5, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinCharge :
      public math::FunctionInterfaceSerializable< double, double>
    {

    private:

    //////////
    // data //
    //////////

      storage::Map< AAType, size_t> m_AACount;    //!< count for each aa type in a sequence
      AAType                  m_NTermAA;    //!< n terminal amino acid since their charges are treated differently
      AAType                  m_CTermAA;    //!< C terminal amino acid since their charges are treated differently
      AATypeData::PropertyType      m_pKProperty; //!< pK scale/values to use

      //! N term charge - amino acid type independent for some of the scales - the scales, where it is amino acid
      //! dependent, values are set to 0
      static const double s_pK_NTerm[];
      //! C term charge - amino acid type independent for some of the scales - the scales, where it is amino acid
      //! dependent, values are set to 0
      static const double s_pK_CTerm[];

    public:

      //! first AATypeData::PropertyType enum which is a pK value
      static const AATypeData::PropertyType s_FirstpKAAProperty;
      //! last AATypeData::PropertyType enum which is a pK value
      static const AATypeData::PropertyType s_LastpKAAProperty;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinCharge();

      //! @brief constructor from amino acid sequence
      //! @param SEQUENCE the amino acid seuqence of the protein of interest
      ProteinCharge( const AASequence &SEQUENCE);

      //! @brief construct from amino acid count
      //! @param AA_COUNT the number of amino acids
      //! @param N_TERM_AA n terminal amino acid
      //! @param C_TERM_AA c terminal amino acid
      ProteinCharge
      (
        const storage::Map< AAType, size_t> &AA_COUNT,
        const AAType &N_TERM_AA,
        const AAType &C_TERM_AA
      );

      //! @brief Clone function
      //! @return pointer to new ProteinCharge
      ProteinCharge *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief set the pk aa property to be used
      //! @brief AA_PROPERTY the aa property to be used; needs to be within the allowed pk property range
      void SetPKProperty( const AATypeData::PropertyType &PROPERTY);

    ///////////////
    // operators //
    ///////////////

      //! @brief for a given pH value clauclate the charge
      //! @param PH the pH of the solution
      //! @return the net charge of the protein within this pH
      double operator()( const double &PH) const;

    //////////////////////
    // input and output //
    //////////////////////

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief correct the given charge by the terminal residue charges for that ph
      //! @param PH the desired pH
      //! @param CHARGE ignoring the termini
      //! @return the correct charge by the terminal amino acids, that usually have different pK value
      double CorrectTerminalCharge( const double &CHARGE, const double &PH) const;

      //! @brief Henderson-Hasselbach-Equation
      //! @param PH the ph
      //! @param PK the pK
      //! @param AA_CHARGE the charge of the amino acid (positive or negative)
      //! @return the net charge - if pk is 0 charge returned is 0
      static double HendersonHasselbachEquation( const double &PH, const double &PK, const double &AA_CHARGE);

    }; // class ProteinCharge

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_PROTEIN_CHARGE_H_
