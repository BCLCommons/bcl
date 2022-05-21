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
#include "score/bcl_score_contact_order.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_average_sd.h"
#include "score/bcl_score_energy_distribution.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @param NORMALIZATION_TYPE type of contact order normalization used
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &ContactOrder::GetDefaultHistogramFilename( const contact::Order::NormalizationType &NORMALIZATION_TYPE)
    {
      // static string
      static const std::string s_default_histogram_filename[ contact::Order::s_NumberNormalizationType + 1] =
      {
        "", // nothing for absolute
        "contact_order_chain_relative_sses.histogram",
        "contact_order_chain_relative_sequence.histogram",
        "contact_order_chain_relative_sses_sqr.histogram",
        "contact_order_chain_relative_sequence_sqr.histogram",
        ""
      };

      // end
      return s_default_histogram_filename[ NORMALIZATION_TYPE];
    }

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &ContactOrder::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "co_score");

      // end
      return s_default_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ContactOrder::ContactOrder() :
      m_Normalize( false),
      m_Scheme(),
      m_HistogramFileName(),
      m_EnergyFunction(),
      m_ContactOrder( contact::Order::e_Absolute, "co", false)
    {
    }

    //! @brief default constructor
    //! @brief NORMALIZATION_TYPE normalization for contact order
    //! @param NORMALIZE flag to enable normalization
    //! @param SCHEME scheme to be used
    //! @param CACHE whether to cache neighbor list generation or not
    ContactOrder::ContactOrder
    (
      const contact::Order::NormalizationType NORMALIZATION_TYPE,
      const bool NORMALIZE,
      const std::string &SCHEME,
      const bool CACHE
    ) :
      m_Normalize( NORMALIZE),
      m_Scheme( SCHEME),
      m_HistogramFileName( GetDefaultHistogramFilename( NORMALIZATION_TYPE)),
      m_EnergyFunction(),
      m_ContactOrder( NORMALIZATION_TYPE, "co", CACHE)
    {
      // read energy function
      ReadEnergyFunction();
    }

    //! @brief virtual copy constructor
    ContactOrder *ContactOrder::Clone() const
    {
      return new ContactOrder( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ContactOrder::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ContactOrder::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief returns filename of the histogram being used
    //! @return filename of the histogram being used
    const std::string &ContactOrder::GetHistogramFilename() const
    {
      return m_HistogramFileName;
    }

    //! @brief returns the energy function
    //! @return energy function
    const util::ShPtr< math::CubicSplineDamped> &ContactOrder::GetEnergyFunction() const
    {
      return m_EnergyFunction;
    }

    //! @brief get a more readable score scheme
    //! @return a more readable score scheme
    const std::string &ContactOrder::GetReadableScheme() const
    {
      static const std::string s_readable_scheme( "Contact order");
      return s_readable_scheme;
    }

    //! @brief get score type
    //! @return score type
    ProteinModel::Type ContactOrder::GetType() const
    {
      return ProteinModel::e_Structure;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of radius of gyration for the given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return the score of radius of gyration for the given ProteinModel
    double ContactOrder::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // if there is a single chain in the model
      if( PROTEIN_MODEL.GetChains().GetSize() == 1)
      {
        // then use the operator to calculate the contact order
        double rel_co( m_ContactOrder( PROTEIN_MODEL));

        // score the co
        double score( m_EnergyFunction->operator ()( rel_co));

        // normalize
        if( !m_Normalize)
        {
          score *= PROTEIN_MODEL.GetNumberAAs();
        }

        // end
        return score;
      }

      // if more than one chain, we need to iterate one by one over the chains and score them individually
      // initialize score
      double score( 0);

      // iterate over chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator itr( PROTEIN_MODEL.GetChains().Begin()),
        itr_end( PROTEIN_MODEL.GetChains().End());
        itr != itr_end;
        ++itr
      )
      {
        // calculate the contact order
        double rel_co( m_ContactOrder.ContactOrder( **itr));

        // score for chain
        double current_score( m_EnergyFunction->operator ()( rel_co));

        // normalize
        if( !m_Normalize)
        {
          current_score *= ( *itr)->GetNumberAAs();
        }

        // sum up the score
        score += current_score;
      }

      // end
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ContactOrder::GetSerializer() const
    {
      io::Serializer serial( m_ContactOrder.GetSerializer());
      serial.SetClassDescription
      (
        "Measure of locality of protein contacts, see https://en.wikipedia.org/wiki/Contact_order "
        "Higher values indicate more non-local structure"
      );
      serial.AddInitializer
      (
        "normalize",
        "True - score will be independent of the size of the protein. "
        "false - the score is multiplied by the number of AAs",
        io::Serialization::GetAgent( &m_Normalize),
        "false"
      );
      return serial;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ContactOrder::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( !m_ContactOrder.ReadInitializerSuccessHook( LABEL, ERR_STREAM))
      {
        return false;
      }
      m_HistogramFileName = GetDefaultHistogramFilename( m_ContactOrder.GetNormalization());
      ReadEnergyFunction();
      return true;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief generate histograms from a table of relative contact orders
    //! @param TABLE each col is a relative contact order
    //! @return Map of col names and histograms
    storage::Map< std::string, math::Histogram> ContactOrder::HistogramsFromColumns( const storage::Table< double> &TABLE)
    {
      storage::Vector< math::RunningAverageSD< double> > mean_sd( TABLE.GetHeader().GetSize());

      // iterate over rows
      for( storage::Table< double>::const_iterator itr( TABLE.Begin()), itr_end( TABLE.End()); itr != itr_end; ++itr)
      {
        // iterate over cols
        for( size_t i( 0), i_max( mean_sd.GetSize()); i < i_max; ++i)
        {
          if( util::IsDefined( itr->Second()( i)))
          {
            mean_sd( i) += itr->Second()( i);
          }
        }
      }

      storage::Vector< math::Histogram> histograms;

      const size_t nr_bins( 30);
      // iterate over mean sd
      for( size_t i( 0), i_max( mean_sd.GetSize()); i < i_max; ++i)
      {
        // binsize and start
        const double binsize( mean_sd( i).GetStandardDeviation() / ( nr_bins / 6));
        const double start( 0.0);

        histograms.PushBack( math::Histogram( start, binsize, nr_bins));
      }

      // iterate over rows
      for( storage::Table< double>::const_iterator itr( TABLE.Begin()), itr_end( TABLE.End()); itr != itr_end; ++itr)
      {
        // iterate over cols
        for( size_t i( 0), i_max( mean_sd.GetSize()); i < i_max; ++i)
        {
          histograms( i).PushBack( itr->Second()( i));
        }
      }

      // fill histograms into map
      storage::Map< std::string, math::Histogram> histogram_map;

      // iterate over histograms
      for( size_t i( 0), i_max( histograms.GetSize()); i < i_max; ++i)
      {
        histogram_map[ TABLE.GetHeader()( i)] = histograms( i);
      }

      // end
      return histogram_map;
    }

    //! @brief read energy function for scoring radius of gyration
    void ContactOrder::ReadEnergyFunction()
    {
      // initialize read
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( m_HistogramFileName));

      // read histogram and reset the stream
      math::Histogram rel_co_histogram;
      read >> rel_co_histogram;
      io::File::CloseClearFStream( read);

      // remove unnecessary bins
      BCL_MessageDbg( util::Format()( rel_co_histogram));
      rel_co_histogram.RemoveBinsAfterIndex( rel_co_histogram.GetIndexOfLastInformationContainingBin());
      BCL_MessageDbg( util::Format()( rel_co_histogram));
      rel_co_histogram.RemoveBinsBeforeIndex( rel_co_histogram.GetIndexOfFirstInformationContainingBin());
      BCL_MessageDbg( util::Format()( rel_co_histogram));
      if( rel_co_histogram.GetHistogram().IsEmpty())
      {
        return;
      }
      rel_co_histogram.SetCount( 0, rel_co_histogram.GetHistogram()( 1) / 2);

      // create a cubic spline from the histogram and store it in the data member m_EnergyFunction
      m_EnergyFunction =
        util::ShPtr< math::CubicSplineDamped>
        (
          new math::CubicSplineDamped
          (
            EnergyDistribution::GeneratePotentialFromHistogram
            (
              rel_co_histogram, 1.0, math::e_FirstDer, storage::Pair< double, double>( 0.0, 0.0)
            )
          )
        );
    }

  } // namespace score
} // namespace bcl
