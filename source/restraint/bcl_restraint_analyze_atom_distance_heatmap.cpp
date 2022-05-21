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
#include "restraint/bcl_restraint_analyze_atom_distance_heatmap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"
#include "restraint/bcl_restraint_epr_distance_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeAtomDistanceHeatmap::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeAtomDistanceHeatmap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeAtomDistanceHeatmap::AnalyzeAtomDistanceHeatmap() :
      m_OutFilePostFix( ".AtomDistanceHeatmap"),
      m_ProteinModelDataType( EPRDistanceData::GetDefaultHandler()),
      m_HistogramMinimum( 10),
      m_HistogramBinSize( 2),
      m_HistogramNumberOfBins( 25),
      m_ShowAtomDistance( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeAtomDistanceHeatmap
    AnalyzeAtomDistanceHeatmap *AnalyzeAtomDistanceHeatmap::Clone() const
    {
      return new AnalyzeAtomDistanceHeatmap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeAtomDistanceHeatmap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeAtomDistanceHeatmap::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeAtomDistanceHeatmap::GetAlias() const
    {
      static const std::string s_Name( "AtomDistanceHeatmap");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeAtomDistanceHeatmap::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // string to hold analysis
      std::string analysis;

      // get the protein model data from one of the models of the ensemble and make sure it could be cast
      util::ShPtrVector< AtomDistance> atom_distances( m_ProteinModelDataType->ReadRestraintsFromFile());

      // get the histograms of distances
      const storage::Vector< math::Histogram> score_histograms( GetDistanceHistograms( atom_distances, ENSEMBLE));

      // get the names of the restraints that will be used as tics
      const storage::Vector< std::string> restraint_tics( GetRestraintNameTics( atom_distances));

      // get the gnuplot script object
      math::GnuplotHeatmap heatmap( GetHeatMap( score_histograms, restraint_tics));

      // true if the atom distance restraints should be displayed
      if( m_ShowAtomDistance)
      {
        const storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > boxes
        (
          GetDistanceBoxes( atom_distances)
        );
        const double x_range( score_histograms.GetSize());
        const double x_min( 0);
        const double y_range( score_histograms.FirstElement().GetBoundaries().Second() - score_histograms.FirstElement().GetBoundaries().First());
        const double y_min( score_histograms.FirstElement().GetBoundaries().First());
        heatmap.SetBoxes( boxes, x_range, y_range, x_min, y_min);
      }

      // write heat map to string stream
      std::stringstream stream;
      heatmap.WriteScript( stream);

      // return the string
      return stream.str();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeAtomDistanceHeatmap::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Makes a heat map showing the frequency with which a distance occurs for each atom distance restraint."
        "Can optionally also show in the heat map the restraint distance with upper and lower bounds."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".AtomDistanceHeatmap"
      );

      parameters.AddInitializer
      (
        "restraint",
        "the type of atom distance related restraint needed for analysis",
        io::Serialization::GetAgent( &m_ProteinModelDataType),
        EPRDistanceData::GetDefaultHandler().GetString()
      );

      parameters.AddInitializer
      (
        "histogram_minimum",
        "the minimal value representing the left boundary of the score histogram",
        io::Serialization::GetAgent( &m_HistogramMinimum),
        "-1.0"
      );

      parameters.AddInitializer
      (
        "histogram_binsize",
        "the width of one bin of the score histograms",
        io::Serialization::GetAgent( &m_HistogramBinSize),
        "0.1"
      );

      parameters.AddInitializer
      (
        "histogram_num_bins",
        "the number of bins in the score histograms",
        io::Serialization::GetAgent( &m_HistogramNumberOfBins),
        "10"
      );

      parameters.AddInitializer
      (
        "show_atom_distance",
        "one if the atom distance restraint should be indicated on the heatmap with upper and lower bounds - 0 otherwise",
        io::Serialization::GetAgent( &m_ShowAtomDistance),
        "1"
      );

      return parameters;
    }

    //! @brief creates histograms of how frequently a distance is observed in models for each restraint
    //! @param DATA the restraints that will be used to get distances
    //! @param ENSEMBLE the ensemble of models that will be used to get distances
    //! @return vector of histograms - one histogram for each restraint
    storage::Vector< math::Histogram> AnalyzeAtomDistanceHeatmap::GetDistanceHistograms
    (
      const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // hold the frequency of distance observed for each atom distance restraint
      storage::Vector< math::Histogram> distance_histograms
      (
        DATA.GetSize(),
        math::Histogram( m_HistogramMinimum, m_HistogramBinSize, m_HistogramNumberOfBins)
      );

      // iterator for the vector of distance histograms
      storage::Vector< math::Histogram>::iterator histogram_itr
      (
        distance_histograms.Begin()), histogram_itr_end( distance_histograms.End()
      );

      // iterate through the atom distance restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
        restraint_itr != restraint_itr_end && histogram_itr != histogram_itr_end;
        ++restraint_itr, ++histogram_itr
      )
      {
        const DataPairwise &data_pair( ( *restraint_itr)->GetData());

        // get the distance from the ensemble fro the current data pair
        const storage::Vector< double> distances( ENSEMBLE.GetDistances( data_pair));

        // add the distances to the current histogram
        histogram_itr->CalculateHistogram( distances);

        histogram_itr->Normalize();
      }

      return distance_histograms;
    }

    //! @brief creates the tics that will be used in the heatmap - converts each restraint into a string to use as tic
    //! @param DATA the restraints that will be used to create tics
    //! @return vector of strings which are the tics representing each restraint
    storage::Vector< std::string> AnalyzeAtomDistanceHeatmap::GetRestraintNameTics
    (
      const util::ShPtrVector< AtomDistance> &DATA
    )
    {
      storage::Vector< std::string> tics;

      // iterate through the data to get the tics from the restraints
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator data_itr( DATA.Begin()), data_itr_end( DATA.End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_a( *( *data_itr)->GetData().First());
        const assemble::LocatorAtomCoordinatesInterface &atom_locator_b( *( *data_itr)->GetData().Second());
        tics.PushBack( assemble::LocatorAtomCoordinatesInterface::GetNameFromPair( atom_locator_a, atom_locator_b));
      }

      return tics;
    }

    //! @brief creates gnuplot heat map object from the distance histograms and the restraint names as tics
    //! @param HISTOGRAMS the histograms that will be used to make the heat map
    //! @param TICS the names of the restraints
    //! @return heat map object which represents the distribution of distances for each restraint
    math::GnuplotHeatmap AnalyzeAtomDistanceHeatmap::GetHeatMap
    (
      const storage::Vector< math::Histogram> &HISTOGRAMS, const storage::Vector< std::string> &TICS
    )
    {
      const util::SiPtrVector< const math::Histogram> histograms
      (
        util::ConvertToConstSiPtrVector( HISTOGRAMS)
      );

      math::GnuplotHeatmap heatmap;

      heatmap.SetFromHistograms( histograms, true, true);
      heatmap.SetTitleAndLabel( "Title", "restraint", "distance", "Frequency (fraction of models)");
      heatmap.SetTicsX( TICS, true, 1);
      heatmap.SetRotationXTics( 90);
      heatmap.SetFilename( "AnalyzeAtomDistanceHeatmap.gnuplot");
      heatmap.SetFont( "/usr/share/fonts/dejavu-lgc/DejaVuLGCSansMono.ttf", 10);

      return heatmap;
    }

    //! @brief gives the coordinates of boxes representing the atom distance restraint data
    //!        each atom distance restraint gives two boxes, one for the upper portion of the restraint, and
    //!        one for the lower portion of the restraint. The actual distance is given by the junction of the two
    //!        boxes
    //! @param DATA the restraints that will be used to create the boxes
    //! @return vector of [(x,y),(x,y)] coordinates - the lower left corner and upper right corner of the box,
    //!         respectively
    storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > >
    AnalyzeAtomDistanceHeatmap::GetDistanceBoxes( const util::ShPtrVector< AtomDistance> &DATA)
    {
      storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > box_coords;

      // iterate through the atom distance restraint data
      for
      (
        util::ShPtrVector< AtomDistance>::const_iterator restraint_itr( DATA.Begin()), restraint_itr_end( DATA.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        const double distance( ( *restraint_itr)->GetDistance()->GetDistance());
        const double upper_bound( ( *restraint_itr)->GetUpperBound());
        const double lower_bound( ( *restraint_itr)->GetLowerBound());
        const double lower_left_x( restraint_itr - DATA.Begin());
        const double upper_right_x( restraint_itr - DATA.Begin() + 1);

        // lower bound box
        {
          const double lower_left_y( lower_bound);
          const double upper_right_y( distance);

          storage::VectorND< 2, storage::VectorND< 2, double> > corner_coords;
          corner_coords.First() = storage::VectorND< 2, double>( lower_left_x, lower_left_y);
          corner_coords.Second() = storage::VectorND< 2, double>( upper_right_x, upper_right_y);
          box_coords.PushBack( corner_coords);
        }

        // upper bound box
        {
          const double lower_left_y( distance);
          const double upper_right_y( upper_bound);

          storage::VectorND< 2, storage::VectorND< 2, double> > corner_coords;
          corner_coords.First() = storage::VectorND< 2, double>( lower_left_x, lower_left_y);
          corner_coords.Second() = storage::VectorND< 2, double>( upper_right_x, upper_right_y);
          box_coords.PushBack( corner_coords);
        }
      }

      return box_coords;
    }

  } // namespace restraint
} // namespace bcl
