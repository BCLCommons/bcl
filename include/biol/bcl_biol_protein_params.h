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

#ifndef BCL_BIOL_PROTEIN_PARAMS_H_
#define BCL_BIOL_PROTEIN_PARAMS_H_

// include the namespace header
#include "bcl_biol.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_biol_aa_type_data.h"
#include "util/bcl_util_si_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinParams
    //! @brief give parameters for a given protein sequence like ExPASy ProtParam program
    //! @details counts all amino acids, report molecular weight, estimated extinction coefficient and other helpful
    //!          params
    //!
    //! @see @link example_biol_protein_params.cpp @endlink
    //! @author woetzen
    //! @date May 21, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinParams :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinParams();

      //! @brief Clone function
      //! @return pointer to new ProteinParams
      ProteinParams *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate protein parameters
      //! @param SEQUENCE amino acid sequence
      //! @return table containing all information
      storage::Table< double> operator()( const AASequence &SEQUENCE) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculate the molecular weight from the number of amino acids
      //! @param AA_COUNT a count for each amino acid type
      //! @return molecular weight in Dalton [Da]
      static double CalcualteMolecularWeight( const storage::Map< AAType, size_t> &AA_COUNT);

      //! @brief calculate the pI value (the pH where the protein is not charged)
      //! @param AA_COUNT a count for each amino acid type
      //! @param N_TERM_TYPE type of n terminal amino acid
      //! @param C_TERM_TYPE type of c terminal amino acid
      //! @param PK_PROPERTY the pk value scale to use
      //! @return the PI value - the pH at which the protein would have no net-charge
      static double CalculatePI
      (
        const storage::Map< AAType, size_t> &AA_COUNT,
        const AAType &N_TERM_AATYPE,
        const AAType &C_TERM_AATYPE,
        const AATypeData::PropertyType &PK_PROPERTY
      );

      //! @brief calculate the extinction coefficient from the number of amino acids in units of  M-1 cm-1, at 280 nm measured in water
      //! @param AA_COUNT a count for each amino acid type
      //! @param RESULT result table
      static void CalcualteExtinctionCoefficient
      (
        const storage::Map< AAType, size_t> &AA_COUNT, storage::Table< double> &RESULT
      );

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

    public:

      //! @brief count the number of each amino acid in the given sequence
      //! @param SEQUENCE amino acid sequence
      //! @return the number of each amino acid type
      static storage::Map< AAType, size_t> CountAAs( const AASequence &SEQUENCE);

    private:

      //! @brief extinction coefficient map at 280 nm in water
      //! @return map of amino acid extinction coefficients
      static storage::Map< AAType, double> ExtinctionCoefficientMap();

      //! @brief add the aa count to the result table
      //! @param AA_COUNT map of aatype and their counts
      //! @param RESULT result table
      static void AddAACountToResultTable
      (
        const storage::Map< AAType, size_t> &AA_COUNT,
        storage::Table< double> &RESULT
      );

    }; // class ProteinParams

  } // namespace biol
} // namespace bcl

#endif // BCL_BIOL_PROTEIN_PARAMS_H_
