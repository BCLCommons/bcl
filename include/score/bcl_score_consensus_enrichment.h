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

#ifndef BCL_SCORE_CONSENSUS_ENRICHMENT_H_
#define BCL_SCORE_CONSENSUS_ENRICHMENT_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConsensusEnrichment
    //! @brief TODO: add a brief comment
    //! @details TODO: add an general comment to this class
    //!
    //! @see @link example_score_consensus_enrichment.cpp @endlink
    //! @author woetzen
    //! @date Jul 29, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConsensusEnrichment :
      public math::FunctionInterfaceSerializable< linal::Vector< double>, double>
    {

    public:
    //////////
    // type //
    //////////

      enum ROCProperty { e_Enrichment, e_ROCIntegeral};

    private:

    //////////
    // data //
    //////////

      storage::Table< double> m_ScoresRmsdTable;
      //! table column to use for classification true/false of enrichment
      size_t                  m_IndexEnrichmentCol;
      //! binary predicate function to classify true/false
      util::ShPtr< util::FunctionInterface< double, bool> > m_UnaryPredicate;

      //! actual roc property to be returned by operator
      ROCProperty m_Property;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConsensusEnrichment()
      {
      }

      //! @brief constructor form table and enrichment column name
      ConsensusEnrichment
      (
        const storage::Table< double> &SCORES_RMSD_TABLE,
        const std::string &ENRICHMENT_COL,
        const ROCProperty &ROC_PROPERTY,
        const util::ShPtr< math::FunctionInterfaceSerializable< double, bool> > &SP_UNARY_PREDICATE
      ) :
        m_ScoresRmsdTable( SCORES_RMSD_TABLE),
        m_IndexEnrichmentCol( SCORES_RMSD_TABLE.GetHeader()[ ENRICHMENT_COL]),
        m_UnaryPredicate( SP_UNARY_PREDICATE),
        m_Property( ROC_PROPERTY)
      {
      }

      //! @brief Clone function
      //! @return pointer to new ConsensusEnrichment
      ConsensusEnrichment *Clone() const
      {
        return new ConsensusEnrichment( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief operator
      //! @param WEIGHTS vector< double> containing weights to be used to calculate sum over cols of table and classify results
      //! @return enrichment
      double operator()( const linal::Vector< double> &WEIGHTS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class ConsensusEnrichment

  } // namespace score
} // namespace bcl

#endif // BCL_SCORE_CONSENSUS_ENRICHMENT_H_ 
