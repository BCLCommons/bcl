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

#ifndef BCL_SCORE_AA_ASSIGNMENT_BLOSUM_H_
#define BCL_SCORE_AA_ASSIGNMENT_BLOSUM_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "function/bcl_function_binary_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAAssignmentBLOSUM
    //! @brief This is a Function derived class for scoring BLOSUM Table for pairs of AA
    //!
    //! @see @link example_score_aa_assignment_blosum.cpp @endlink
    //! @author meilerj, woetzen
    //! @date 21.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAAssignmentBLOSUM :
      public function::BinaryInterface< const biol::AABase, const biol::AABase, double>
    {

    public:

    ///////////
    // types //
    ///////////

      //! types of Blosum tables
      enum TableType
      {
        e_BLOSUM_90,
        e_BLOSUM_80,
        e_BLOSUM_62,
        e_BLOSUM_45,
        s_NumberTableTypes
      };

      //! @brief TableType as string
      //! @param TABLE_TYPE the TableType
      //! @return the string for the TableType
      static const std::string &GetTableDescriptor( const TableType &TABLE_TYPE);

      //! @brief TableTypeEnum is used for I/O of TableType
      typedef util::WrapperEnum< TableType, &GetTableDescriptor, s_NumberTableTypes> TableTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! Blosum90, Blosum80, Blosum62, Blosum45 substitution matrix for global alignments
      //! taken from David Mount. "Bioinformatics" and
      //! http://eta.embl-heidelberg.de:8000/misc/mat/
      //! log of probability of replacing aa i with aa j divided by freqency of aa i
      static const double s_BLOSUMTable[][ biol::AATypes::s_NumberStandardAATypes + 4][ biol::AATypes::s_NumberStandardAATypes + 4];

      TableTypeEnum m_BLOSUMTable; //!< TableType to use

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a BLOSUM table
      //! @param BLOSUM_TABLE TableType to be used
      AAAssignmentBLOSUM( const TableType &BLOSUM_TABLE = e_BLOSUM_62);

      //! @brief virtual copy constructor
      //! @return pointer to a new AAAssignmentBLOSUM copied from this one
      AAAssignmentBLOSUM *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the TableType
      //! @return returns m_BLOSUMTable
      const TableType GetTableType() const
      {
        return m_BLOSUMTable;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief return the mutation probability for two AAs given a table
      //! @param BLOSUM_TABLE TableType to be used
      //! @param AA_TYPE_A first AAType of interest
      //! @param AA_TYPE_B second AAType of interest
      //! @return the mutation probability for two AAs given a table
      static double Probability
      (
        const TableType BLOSUM_TABLE,
        const biol::AAType &AA_TYPE_A,
        const biol::AAType &AA_TYPE_B
      );

      //! @brief returns the requested BLOSUM table as a matrix
      //! @param BLOSUM_TABLE requested BLOSUM table name
      //! @return the requested BLOSUM table as a matrix
      static linal::Matrix< double> GetBLOSUMMatrix( const TableType &BLOSUM_TABLE);

      //! @brief returns the requested BLOSUM table as a matrix
      //! @param BLOSUM_TABLE requested BLOSUM table name
      //! @return the requested BLOSUM table as a matrix
      static linal::VectorConstReference< double> GetBLOSUMRow
      (
        const TableType &BLOSUM_TABLE,
        const biol::AAType &AATYPE
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the score between two assigned members
      //! @param MEMBER_A amino acid A that is compared
      //! @param MEMBER_B amino acid A that is compared
      //! @return value from the m_BLOSUMTable for this combination of amino acids
      double operator()( const biol::AABase &MEMBER_A, const biol::AABase &MEMBER_B) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

    }; // class AAAssignmentBLOSUM

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_AA_ASSIGNMENT_BLOSUM_H_
