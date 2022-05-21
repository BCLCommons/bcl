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

#ifndef BCL_SCORE_AA_ASSIGNMENT_PAM_H_
#define BCL_SCORE_AA_ASSIGNMENT_PAM_H_

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
    //! @class AAAssignmentPAM
    //! @brief This is a Function derived class for scoring PAM Table for pairs of AA
    //!
    //! @see @link example_score_aa_assignment_pam.cpp @endlink
    //! @author meilerj, woetzen
    //! @date 21.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAAssignmentPAM :
      public function::BinaryInterface< const biol::AABase, const biol::AABase, double>
    {

    public:

    ///////////
    // enums //
    ///////////

      //! types of PAM tables
      enum TableType
      {
        e_PAM_100,
        e_PAM_120,
        e_PAM_160,
        e_PAM_250,
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

      //! Pam100, Pam120, Pam160, Pam250 substitution matrix for global alignments
      //! taken from David Mount. "Bioinformatics" and
      //! http://eta.embl-heidelberg.de:8000/misc/mat/
      //! log of probability of replacing aa i with aa j divided by frequency of aa i
      static const double s_PAMTable[][ biol::AATypes::s_NumberStandardAATypes + 4][ biol::AATypes::s_NumberStandardAATypes + 4];

      TableTypeEnum m_PAMTable; //!< which TableType to use

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a PAM_TABLE
      //! @param PAM_TABLE Pam table to be used
      AAAssignmentPAM( const TableType PAM_TABLE = e_PAM_250);

      //! @brief virtual copy constructor
      //! @return pointer to a new AAAssignmentPAM copied from this one
      AAAssignmentPAM *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief return the mutation probability for two AAs given a table
      //! @param PAM_TABLE Pam table to be used
      //! @param AA_TYPE_A first AAType of interest
      //! @param AA_TYPE_B second AAType of interest
      //! @return the mutation probability for two AAs given a table
      static double Probability( const TableType PAM_TABLE, const biol::AAType &AA_TYPE_A, const biol::AAType &AA_TYPE_B)
      {
        return pow( 10, s_PAMTable[ PAM_TABLE][ AA_TYPE_A][ AA_TYPE_B]) / double( 20);
      }

      //! @brief returns the requested PAM table as a matrix
      //! @param PAM_TABLE requested PAM table name
      //! @return the requested PAM table as a matrix
      static linal::Matrix< double> GetPAMMatrix( const TableType PAM_TABLE);

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that calculates the score between two assigned members
      //! @param MEMBER_A amino acid A that is compared
      //! @param MEMBER_B amino acid A that is compared
      //! @return value from the m_PAMTable for this combination of amino acids
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

    }; // class AAAssignmentPAM

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_AA_ASSIGNMENT_PAM_H_
