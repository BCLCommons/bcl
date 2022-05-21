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
#include "score/bcl_score_consensus_enrichment.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_roc_curve.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConsensusEnrichment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief operator
    //! @param WEIGHTS vector< double> containing weights to be used to calculate enrichment
    //! @return enrichment
    double ConsensusEnrichment::operator()( const linal::Vector< double> &WEIGHTS) const
    {
      storage::List< storage::Pair< double, bool> > classified_results;

      for
      (
        storage::Table< double>::const_iterator itr( m_ScoresRmsdTable.Begin()), itr_end( m_ScoresRmsdTable.End());
        itr != itr_end;
        ++itr
      )
      {
        const linal::Vector< double> scores( itr->Second().GetData());
        classified_results.PushBack
        (
          storage::Pair< double, bool>
          (
            scores * WEIGHTS,
            m_UnaryPredicate->operator()( itr->Second()( m_IndexEnrichmentCol))
          )
        );
      }

      classified_results.Sort
      (
        storage::PairBinaryPredicateFirst< double, bool>
        (
          util::BinaryFunctionSTLWrapper< std::less< double> >()
        )
      );

      const math::ROCCurve roc_curve( classified_results);
      switch( m_Property)
      {
        case e_Enrichment:
        {
          const double enrichment
          (
            roc_curve.ContingencyMatrixFraction
            (
              double( roc_curve.GetNumberActualPositives()) /
              double( roc_curve.GetNumberResults())
            ).GetEnrichment()
          );
          return util::IsDefined( enrichment) ? enrichment : 0.0;
        }
        case e_ROCIntegeral:
        {
          return roc_curve.Integral();
        }
        default:
        {
          return util::GetUndefinedDouble();
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConsensusEnrichment::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ScoresRmsdTable   , ISTREAM);
      io::Serialize::Read( m_IndexEnrichmentCol, ISTREAM);
      io::Serialize::Read( m_UnaryPredicate    , ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ConsensusEnrichment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ScoresRmsdTable   , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IndexEnrichmentCol, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UnaryPredicate    , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
