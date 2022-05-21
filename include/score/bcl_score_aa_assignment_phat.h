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

#ifndef BCL_SCORE_AA_ASSIGNMENT_PHAT_H_
#define BCL_SCORE_AA_ASSIGNMENT_PHAT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "function/bcl_function_binary_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAAssignmentPHAT
    //! @brief This is a Function derived class for scoring PHAT Table for pairs of AA
    //!
    //! @see @link example_score_aa_assignment_phat.cpp @endlink
    //! @author dongen
    //! @date 08.24.2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAAssignmentPHAT :
      public function::BinaryInterface< const biol::AABase, const biol::AABase, double>
    {

    public:

    ///////////
    // enums //
    ///////////

      //! types of Phat tables
      enum TableType
      {
        e_PHAT_85,
        e_PHAT_80,
        e_PHAT_75,
        e_PHAT_70,
        s_NumberTableTypes
      };

      //! @brief TableType as string
      //! @param PHAT_TABLE the message level
      //! @return the TableType as string
      static const std::string &GetTableTypeString( const TableType &PHAT_TABLE);

      //! @brief TableTypeEnum is used for I/O of TableType
      typedef util::WrapperEnum< TableType, &GetTableTypeString, s_NumberTableTypes> TableTypeEnum;

    private:

    //////////
    // data //
    //////////

      //! PHAT T70_B66, PHAT T75_B73, PHAT T80_B78, PHAT T85_B82 substitution matrix for transmembrane regions
      //! taken from Ng, Henikoff, Henikoff 2000. "Bioinformatics" and
      //! http://blocks.fhcrc.org/sift/PHAT.html
      static const double s_PHATTable[][ biol::AATypes::s_NumberStandardAATypes + 4][ biol::AATypes::s_NumberStandardAATypes + 4];

      TableTypeEnum m_PHATTable; //!< which TableType to use

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a PHAT table
      //! @param PHAT_TABLE TableType to be used
      AAAssignmentPHAT( const TableType &PHAT_TABLE = e_PHAT_70);

      //! @brief virtual copy constructor
      //! @return pointer to a new AAAssignmentPHAT copied from this one
      AAAssignmentPHAT *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetTableType gives "m_PHATTable"
      //! @return return "m_PHATTable"
      const TableType &GetTableType() const
      {
        return m_PHATTable;
      }

      //! @brief returns the requested PHAT table as a matrix
      //! @param PHAT_TABLE requested Phat table name
      //! @return the requested PHAT table as a matrix
      static linal::Matrix< double> GetPHATMatrix( const TableType &PHAT_TABLE);

    ////////////////
    // operations //
    ////////////////

      //! @brief return the mutation probability for two AAs given a table
      //! @param PHAT_TABLE TableType to be used
      //! @param AA_TYPE_A first AAType of interest
      //! @param AA_TYPE_B second AAType of interest
      //! @return the mutation probability for two AAs given a table
      static double Probability
      (
        const TableType &PHAT_TABLE,
        const biol::AAType &AA_TYPE_A,
        const biol::AAType &AA_TYPE_B
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the score between two assigned members
      //! @param MEMBER_A amino acid A that is compared
      //! @param MEMBER_B amino acid A that is compared
      //! @return value from the m_PHATTable for this combination of amino acids
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

    }; // class AAAssignmentPHAT

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_AA_ASSIGNMENT_PHAT_H_
