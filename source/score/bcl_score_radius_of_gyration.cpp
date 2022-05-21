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
#include "score/bcl_score_radius_of_gyration.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"
#include "score/bcl_score_energy_distribution.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> RadiusOfGyration::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new RadiusOfGyration())
    );

  //////////
  // data //
  //////////

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return the default filename for the default phi_psi_angle statistics
    const std::string &RadiusOfGyration::GetDefaultHistogramFilename()
    {
      // static string
      static const std::string s_default_histogram_filename( "radius_of_gyration_chains.histogram");

      // end
      return s_default_histogram_filename;
    }

    //! @brief returns default file where the statistics and in consequence the energy potentials are read from
    //! @return default file where the statistics and in consequence the energy potentials are read from
    const std::string &RadiusOfGyration::GetDefaultHistogramFilenameMembrane()
    {
      // static string
      static const std::string s_default_histogram_filename( "radius_of_gyration_model_membrane.histogram");

      // end
      return s_default_histogram_filename;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    //! @param NORMALIZE flag to enable normalization
    //! @param SCHEME scheme to be used
    //! @param SOLUBLE_FILENAME filename for soluble histogram
    //! @param MEMBRANE_FILENAME filename for membrane histogram
    RadiusOfGyration::RadiusOfGyration
    (
      const bool NORMALIZE,
      const bool RAW,
      const std::string &SCHEME,
      const std::string &SOLUBLE_FILENAME,
      const std::string &MEMBRANE_FILENAME
    ) :
      m_Normalize( NORMALIZE),
      m_Raw( RAW),
      m_Scheme( SCHEME.empty() ? std::string( RAW ? "rgyr_raw" : "rgyr") + std::string( NORMALIZE ? "_norm" : "") : SCHEME),
      m_HistogramFileNameSoluble( SOLUBLE_FILENAME),
      m_HistogramFileNameMembrane( MEMBRANE_FILENAME),
      m_EnergyFunctionSoluble(),
      m_EnergyFunctionMembrane()
    {
      std::ostringstream oss;
      ReadInitializerSuccessHook( util::ObjectDataLabel(), oss);
    }

    //! @brief virtual copy constructor
    RadiusOfGyration *RadiusOfGyration::Clone() const
    {
      return new RadiusOfGyration( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RadiusOfGyration::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &RadiusOfGyration::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief returns filename of the histogram being used for soluble proteins
    //! @return filename of the histogram being used for soluble proteins
    const std::string &RadiusOfGyration::GetHistogramFilenameSoluble() const
    {
      return m_HistogramFileNameSoluble;
    }

    //! @brief returns filename of the histogram being used for membrane proteins
    //! @return filename of the histogram being used for membrane proteins
    const std::string &RadiusOfGyration::GetHistogramFilenameMembrane() const
    {
      return m_HistogramFileNameMembrane;
    }

    //! @brief returns the energy function
    //! @return energy function
    const math::CubicSplineDamped &RadiusOfGyration::GetEnergyFunctionSoluble() const
    {
      return *m_EnergyFunctionSoluble;
    }

    //! @brief returns the membrane energy function
    //! @return membrane energy function
    const math::CubicSplineDamped &RadiusOfGyration::GetEnergyFunctionMembrane() const
    {
      return *m_EnergyFunctionMembrane;
    }

    //! @brief get a more readable score scheme
    //! @return a more readable score scheme
    const std::string &RadiusOfGyration::GetReadableScheme() const
    {
      static const std::string s_readable_scheme( "Radius of gyration");
      return s_readable_scheme;
    }

    //! @brief get score type
    //! @return score type
    ProteinModel::Type RadiusOfGyration::GetType() const
    {
      return ProteinModel::e_Structure;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RadiusOfGyration::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Calculates the radius of gyration score of a protein model."
      );
      serializer.AddInitializer
      (
        "raw",
        "if set, no energy function (from histograms) will be used. Instead, the raw radius of gyration will be returned. ",
        io::Serialization::GetAgent( &m_Raw),
        "False"
      );
      serializer.AddInitializer
      (
        "normalize",
        "normalize by sequence length",
        io::Serialization::GetAgent( &m_Normalize),
        "False"
      );
      serializer.AddInitializer
      (
        "histogram file name soluble",
        "file name of the histogram for the radius of gyration of soluble proteins.",
        io::Serialization::GetAgent( &m_HistogramFileNameSoluble),
        GetDefaultHistogramFilename()
      );
      serializer.AddInitializer
      (
        "histogram file name membrane",
        "file name of the histogram for the radius of gyration of membrane proteins.",
        io::Serialization::GetAgent( &m_HistogramFileNameMembrane),
        GetDefaultHistogramFilenameMembrane()
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates the radius of gyration for CB atoms for ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return the radius of gyration for CB atoms for ProteinModel
    double RadiusOfGyration::SquareRadiusOfGyration( const assemble::ProteinModel &PROTEIN_MODEL)
    {
      // calculate the radius of gyration
      return coord::SquareRadiusOfGyration
      (
        PROTEIN_MODEL.GetAtomCoordinates( biol::GetAtomTypes().GetFirstSidechainAtomTypes())
      );
    }

    //! @brief calculates the radius of gyration for specified atom type for ProteinModel with a collapsing normal
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @param ATOM_TYPES Atom types that will be considered when doing the calculations
    //! @return the radius of gyration for specified atom type for ProteinModel with a collapsing normal
    double RadiusOfGyration::SquareRadiusOfGyrationCollapsed
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      const storage::Set< biol::AtomType> &ATOM_TYPES
    )
    {
      // try to get the membrane
      const util::ShPtr< biol::Membrane> &sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // calculate the radius of gyration
      return SquareRadiusOfGyrationCollapsed
      (
        PROTEIN_MODEL.GetAtomCoordinates( ATOM_TYPES),
        sp_membrane
      );
    }

    //! @brief compute the square radius of gyration for a set of coordinates and a collapsing normal
    //! @param COORDINATES coordinates to be collapsed
    //! @param MEMBRANE membrane object to collapse
    //! @return square radius of gyration collapsed
    double RadiusOfGyration::SquareRadiusOfGyrationCollapsed
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::ShPtr< biol::Membrane> &MEMBRANE
    )
    {
      // if membrane not defined
      if( !MEMBRANE.IsDefined() || MEMBRANE->GetThickness( biol::GetEnvironmentTypes().e_MembraneCore) == 0.0)
      {
        // return normal radius of gyration
        return coord::SquareRadiusOfGyration( COORDINATES);
      }

      //collapsed coordinates
      storage::Vector< linal::Vector3D> collapsed_coordinates;
      collapsed_coordinates.AllocateMemory( COORDINATES.GetSize());

      // collapse distance
      const double collapse_distance( MEMBRANE->GetThickness( biol::GetEnvironmentTypes().e_MembraneCore));
      const linal::Vector3D normal( MEMBRANE->GetNormal());

      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          coordinates_itr( COORDINATES.Begin()), coordinates_itr_end( COORDINATES.End());
        coordinates_itr != coordinates_itr_end;
        ++coordinates_itr
      )
      {
        // skip undefined coordinates
        if( !( *coordinates_itr)->IsDefined())
        {
          continue;
        }

        linal::Vector3D current_coord( **coordinates_itr);
        double distance( linal::ScalarProduct( normal, current_coord - MEMBRANE->GetCenter()));

        // if distance is above collapse zone
        if( distance > collapse_distance)
        {
          // move toward plane by collapse distance
          current_coord -= collapse_distance * normal;
        }
        // if distance is below collapse zone
        else if( distance < -collapse_distance)
        {
          // move toward plane by collapse distance
          current_coord += collapse_distance * normal;
        }
        // if distance is within collapse zone
        else
        {
          // move onto plane
          current_coord -= distance * normal;
        }

        // insert
        collapsed_coordinates.PushBack( current_coord);
      }

      //end
      return coord::SquareRadiusOfGyration( util::ConvertToSiPtrVector( collapsed_coordinates));
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    //! @return return code indicating success or failure
    bool RadiusOfGyration::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        m_EnergyFunctionSoluble = ReadEnergyFunction( m_HistogramFileNameSoluble);
        m_EnergyFunctionMembrane = ReadEnergyFunction( m_HistogramFileNameMembrane);
      }
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the score of radius of gyration for the given ProteinModel
    //! @param PROTEIN_MODEL ProteinModel of interest
    //! @return the score of radius of gyration for the given ProteinModel
    double RadiusOfGyration::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      BCL_MessageDbg( "score radius of gyration");
      if( PROTEIN_MODEL.GetChains().IsEmpty())
      {
        return double( 0);
      }

      if( PROTEIN_MODEL.GetNumberSSE( biol::GetSSTypes().COIL) > size_t( 0))
      {
        assemble::ProteinModel copy( PROTEIN_MODEL);
        auto coils( copy.GetSSEs( biol::GetSSTypes().COIL));
        for( auto itr( coils.Begin()), itr_end( coils.End()); itr != itr_end; ++itr)
        {
          copy.Remove( **itr);
        }
        copy.SetProteinModelData( PROTEIN_MODEL.GetProteinModelData());
        return operator()( copy);
      }
      // collect all first side chain atom coordinates
      const util::SiPtrVector< const linal::Vector3D> first_sc_atom_coordinates
      (
        PROTEIN_MODEL.GetAtomCoordinates( biol::GetAtomTypes().GetFirstSidechainAtomTypes())
      );

      // get the total number of residues between
      size_t number_res( 0);
      for
      (
        auto chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        if( ( *chain_itr)->GetNumberSSEs())
        {
          auto sses( ( *chain_itr)->GetSSEs());
          number_res += sses.LastElement()->GetFirstAA()->GetSeqID() - sses.FirstElement()->GetFirstAA()->GetSeqID() + 1;
        }
      }

      // if no atom coordinates
      if( first_sc_atom_coordinates.IsEmpty())
      {
        return double( 0);
      }

      // initialize values to calculate
      double square_radius_of_gyration;
      double ratio_square_rgyr_residues;
      double score;

      // try to get the membrane
      const util::ShPtr< biol::Membrane> &sp_membrane
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Membrane)
      );

      // if this is a membrane protein
      if( sp_membrane.IsDefined() && sp_membrane->GetThickness( biol::GetEnvironmentTypes().e_MembraneCore) != 0.0)
      {
        square_radius_of_gyration = SquareRadiusOfGyrationCollapsed( first_sc_atom_coordinates, sp_membrane);

        // calculate the ratio of square radius of gyration to number of coordinates (residues) used
        ratio_square_rgyr_residues = square_radius_of_gyration / double( number_res);

        // return the value from the spline
        score = util::IsDefined( ratio_square_rgyr_residues)
                ? ( m_Raw ? ratio_square_rgyr_residues - 1.3 : m_EnergyFunctionMembrane->operator ()( ratio_square_rgyr_residues))
                : 0;
      }
      // a soluble protein
      else
      {
        square_radius_of_gyration = coord::SquareRadiusOfGyration( first_sc_atom_coordinates);

        // calculate the ratio of square radius of gyration to number of coordinates (residues) used
        ratio_square_rgyr_residues = square_radius_of_gyration / double( number_res);

        // return the value from the spline
        score = util::IsDefined( ratio_square_rgyr_residues)
                ? ( m_Raw ? ratio_square_rgyr_residues - 1.3 : m_EnergyFunctionSoluble->operator ()( ratio_square_rgyr_residues))
                : 0;
      }

      // if normalization flag is not given
      if( !m_Normalize)
      {
        score *= double( first_sc_atom_coordinates.GetSize());
      }

      // end
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief create histograms from a table generated from statistics application
    //! @param TABLE containing rows: rgyr_sqr, nr_coordinates, subunits
    //! @return map containg histograms for proteins with a single subunit and one over all with the filenames as key
    storage::Map< std::string, math::Histogram> RadiusOfGyration::HistogramsFromTable( const storage::Table< double> &TABLE)
    {
      // calculate mean and sd of rgyr2 and number residues
      math::RunningAverageSD< double> single_mean_sd;
      math::RunningAverageSD< double> multi_mean_sd;
      for( storage::Table< double>::const_iterator itr( TABLE.Begin()), itr_end( TABLE.End()); itr != itr_end; ++itr)
      {
        const double ratio( itr->Second()[ "rgyr_sqr"] / itr->Second()[ "nr_aa"]);
        if( util::IsDefined( ratio) && ratio > 0.0)
        {
          if( itr->Second()[ "subunits"] == 1.0)
          {
            single_mean_sd += ratio;
          }
          multi_mean_sd += ratio;
        }
      }

      // binsize and start
      const size_t nr_bins( 30);
      const double binsize_single( 0.2);
      const double start_single( 0.0);
      const double binsize_multi( 0.2);
      const double start_multi( 0.0);

      math::Histogram radius_of_gryation_chains_histogram( start_single, binsize_single, nr_bins);
      math::Histogram radius_of_gryation_model_histogram( start_multi, binsize_multi, nr_bins);

      for( storage::Table< double>::const_iterator itr( TABLE.Begin()), itr_end( TABLE.End()); itr != itr_end; ++itr)
      {
        const double ratio( itr->Second()[ "rgyr_sqr"] / itr->Second()[ "nr_aa"]);
        if( util::IsDefined( ratio) && ratio > 0.0)
        {
          if( itr->Second()[ "subunits"] == 1.0)
          {
            radius_of_gryation_chains_histogram.PushBack( ratio);
          }
          radius_of_gryation_model_histogram.PushBack( ratio);
        }
      }

      storage::Map< std::string, math::Histogram> histograms;

      // write rgyr histogram for chains
      histograms[ "radius_of_gyration_chains"] = radius_of_gryation_chains_histogram;
      histograms[ "radius_of_gyration_model"] = radius_of_gryation_model_histogram;

      // end
      return histograms;
    }

    //! @brief read individual energy functions for scoring radius of gyration
    //! @param FILENAME filename to be read in
    //! @return read in energy function
    util::ShPtr< math::CubicSplineDamped> RadiusOfGyration::ReadEnergyFunction( const std::string &FILENAME)
    {
      // initialize read
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( FILENAME));

      // read histogram and reset the stream
      math::Histogram rgyr_histogram;
      read >> rgyr_histogram;
      io::File::CloseClearFStream( read);

      // remove unnecessary bins
      // iterate over the all bins and find the second positive slope
      bool downwards( false);
      rgyr_histogram.RemoveBinsBeforeIndex( rgyr_histogram.GetIndexOfFirstInformationContainingBin());
      size_t i( 0), i_max( rgyr_histogram.GetNumberOfBins() - 1);
      for( ; i < i_max; ++i)
      {
        // still first slope?
        if( !downwards)
        {
          downwards = rgyr_histogram.GetHistogram()( i) > rgyr_histogram.GetHistogram()( i + 1);
        }
        // first slope overcome - ist there another positive slope
        else if( rgyr_histogram.GetHistogram()( i) < rgyr_histogram.GetHistogram()( i + 1))
        {
          rgyr_histogram.RemoveBinsAfterIndex( i);
          break;
        }
      }

      // was there a second positive slope?
      if( i >= i_max)
      {
        rgyr_histogram.RemoveBinsAfterIndex( rgyr_histogram.GetIndexOfLastInformationContainingBin());
      }

      // create a cubic spline from the histogram and return
      return util::ShPtr< math::CubicSplineDamped>
      (
        new math::CubicSplineDamped
        (
          EnergyDistribution::GeneratePotentialFromHistogram
          (
            rgyr_histogram, 1.0, math::e_FirstDer, storage::Pair< double, double>( 0.0, 0.0)
          )
        )
      );
    }

  } // namespace score
} // namespace bcl
