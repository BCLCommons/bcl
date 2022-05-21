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
#include "score/bcl_score_read_histograms.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_environment_types.h"
#include "math/bcl_math_histogram.h"
#include "score/bcl_score_energy_distribution.h"
#include "util/bcl_util_si_ptr_vector.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  ////////////////
  // operations //
  ////////////////

    //! @brief read histograms for membrane dependent score from stream
    //! also removes additional empty bins from the end of the histograms - but keeping them all of the same length
    //! @param ISTREAM stream containing the histograms for each membrane region 20 amino acid histograms
    //! @return vector for all membrane regions containing vectors of histograms containing all amino acids
    storage::Vector< storage::Vector< math::Histogram> >
    ReadHistograms::ReadMembraneDependentEnvironmentHistograms
    (
      std::istream &ISTREAM
    )
    {
      // vector for all membrane regions containing vectors of histograms containing all amino acids
      storage::Vector< storage::Vector< math::Histogram> >
        membrane_aa_histograms
        (
          biol::GetEnvironmentTypes().GetEnumCount(),
          storage::Vector< math::Histogram>
          (
            biol::GetAATypes().GetEnumCount()
          )
        );

      // initialize vector of histograms
      util::SiPtrVector< math::Histogram> all_histograms;

      // iterate over environment types
      while( ISTREAM.good() && !ISTREAM.eof())
      {
        // read environment type and make sure it is correct
        std::string tmp_env;
        ISTREAM >> tmp_env;
        if( !ISTREAM.good())
        {
          break;
        }
        ISTREAM >> tmp_env;
        biol::EnvironmentType current_env( util::Strip( tmp_env, "\""));

        // iterate over AATypes
        for( biol::AATypes::const_iterator
               aa_itr( biol::GetAATypes().Begin()),
               aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
             aa_itr != aa_itr_end;
             ++aa_itr)
        {
          // get the one letter code and make sure it is correct
          std::string tmp;
          ISTREAM >> tmp;
          biol::AAType aatype( biol::GetAATypes().AATypeFromOneLetterCode( tmp[ 0]));
          BCL_Assert
          (
            *aa_itr == aatype,
            "unexpected aatype read from file! " + util::Format()( *aa_itr) + " != " +   util::Format()( aatype)
          );

          // read the corresponding histogram
          ISTREAM >> membrane_aa_histograms( current_env)( aatype);

          // pushback the histogram
          all_histograms.PushBack( util::SiPtr< math::Histogram>( membrane_aa_histograms( current_env)( aatype)));
        }
      }

      // remove all empty bins - but keep the length of all histograms the same
      EnergyDistribution::RemoveAdditionalEmptyBinsExceptOne( all_histograms);

      // end
      return membrane_aa_histograms;
    }

    //! read histograms from stream, containing 20 histograms for each aminoacid
    storage::Vector< math::Histogram>
    ReadHistograms::ReadEnvironmentHistograms
    (
      std::istream &ISTREAM
    )
    {
      // histograms for all aminoacids
      storage::Vector< math::Histogram> histograms( biol::AATypes::s_NumberStandardAATypes);

      for( biol::AATypes::const_iterator
             aa_itr( biol::GetAATypes().Begin()),
             aa_itr_end( biol::GetAATypes().GetEnumIteratorFromIndex( biol::AATypes::s_NumberStandardAATypes));
           aa_itr != aa_itr_end;
           ++aa_itr)
      {
        math::Histogram current_aa_env_histogram;
        biol::AAType current_aatype;

        // read one letter code and make sure it is correct amino acid
        std::string tmp;
        ISTREAM >> tmp;
        current_aatype = biol::GetAATypes().AATypeFromOneLetterCode( tmp[ 0]);
        BCL_Assert
        (
          *aa_itr == current_aatype,
          "unexpected aatype read from file! " + util::Format()( *aa_itr) + " != " +   util::Format()( current_aatype)
        );

        // read the corresponding histogram for this AAType
        ISTREAM >> histograms( current_aatype);
      }

      util::SiPtrVector< math::Histogram> all_histograms( util::ConvertToSiPtrVector( histograms));
      // remove all empty bins - but keep the length of all histograms the same
      EnergyDistribution::RemoveAdditionalEmptyBinsExceptOne( all_histograms);

      // end
      return histograms;
    }

  } // namespace score
} // namespace bcl
